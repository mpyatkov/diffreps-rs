//! Step 3 of diffReps: slide a window along each chromosome, test it, and
//! merge significant windows into regions.
//!
//! Port of the per-chromosome loop in diffReps.pl (lines ~394-480), the
//! region logic of `diffReps::ChromaModSite` and the `WindowCoverage` class.
//!
//! For window `w` (0-based), the window covers `[w*step + 1, w*step + window]`
//! (1-based, inclusive), i.e. bins `w .. w + window/step` of size `step`.
//!
//! ```text
//! for each window w:
//!     if all treatment counts > cutoff, or all control counts > cutoff:
//!         coverage += bp not counted yet                  (WindowCoverage)
//!         test the window (NB)
//!         if p < --pval:
//!             if an open region ends within --gap bp and has the same direction:
//!                 extend it to this window; keep the best (lowest p) window
//!             else:
//!                 close the open region (re-test it if it has >1 window)
//!                 open a new region = this window
//! close the last region
//! ```

use crate::background::Background;
use crate::bed::{BinCounts, Sample};
use crate::inputs::Chrom;
use crate::nbtest::{Event, Stat, calc_stat};

/// Settings for the scan.
#[derive(Debug, Clone, Copy)]
pub struct ScanParams {
    pub window: u64,
    pub step: u64,
    pub gap: u64,
    pub nsd: f64,
    pub pval: f64,
}

/// A differential region, ready for the report.
#[derive(Debug, Clone)]
pub struct Region {
    pub chrom: usize,
    /// 1-based inclusive coordinates.
    pub start: u64,
    pub end: u64,
    /// Normalized counts (raw / norm factor) of every sample, treatment first.
    pub norm_counts: Vec<f64>,
    /// Statistics of the whole region.
    pub stat: Stat,
    /// Best (lowest p) window: its start and statistics.
    pub win_sta: u64,
    pub win_stat: Stat,
}

/// Regions and window coverage of one chromosome.
#[derive(Debug, Clone, Default)]
pub struct ChromResult {
    pub regions: Vec<Region>,
    /// Genome coverage (bp) of windows that passed the count cutoff.
    pub window_cov: u64,
}

/// Region being extended. Mirrors the fields of `ChromaModSite`.
struct Open {
    start: u64,
    end: u64,
    /// Direction of the first window; windows in the other direction start a new region.
    dirn: Event,
    /// `Some` while the region is a single window (its stats are the window's);
    /// `None` after it has been extended (Perl: `pval` undefined -> `needStat`).
    stat: Option<(Stat, Vec<f64>)>,
    win_sta: u64,
    win_stat: Stat,
}

/// Counts of all samples on one chromosome, in `step`-size bins.
struct ChromCounts<'a> {
    bins: Vec<BinCounts>,
    n_tr: usize,
    norm: &'a [f64],
    step: u64,
}

impl ChromCounts<'_> {
    /// Raw counts of every sample over bins `first..=last`.
    fn raw(&self, first: usize, last: usize) -> Vec<f64> {
        self.bins.iter().map(|b| b.range(first, last) as f64).collect()
    }

    /// `calc_stat` on raw counts; also returns the normalized counts.
    fn stat(&self, raw: &[f64]) -> (Stat, Vec<f64>) {
        let (tr, co) = raw.split_at(self.n_tr);
        let (nt, nc) = self.norm.split_at(self.n_tr);
        let s = calc_stat(tr, co, nt, nc);
        let normalized = raw.iter().zip(self.norm).map(|(c, f)| c / f).collect();
        (s, normalized)
    }

    /// Statistics of the region `[start, end]` (`get_region_count`, SRBed.pm:240).
    fn region_stat(&self, start: u64, end: u64) -> (Stat, Vec<f64>) {
        let first = ((start - 1) / self.step) as usize;
        let last = (end / self.step) as usize - 1;
        self.stat(&self.raw(first, last))
    }
}

/// Scan one chromosome. `norm` and `bg` hold one value per sample, treatment first.
pub fn scan_chrom(
    samples: &[Sample],
    n_tr: usize,
    ci: usize,
    chrom: &Chrom,
    p: &ScanParams,
    bg: &Background,
    norm: &[f64],
) -> ChromResult {
    let bpw = (p.window / p.step) as usize;
    let counts = ChromCounts {
        bins: samples.iter().map(|s| BinCounts::new(&s.pos2[ci], p.step, chrom.len)).collect(),
        n_tr,
        norm,
        step: p.step,
    };
    let Some(n_windows) = counts.bins[0].n_windows(bpw) else {
        eprintln!("Create count vector for treatment at: {} failed! Skip.", chrom.name);
        return ChromResult::default();
    };
    let cutoff: Vec<f64> = (0..samples.len()).map(|j| bg.cutoff(j, p.nsd)).collect();
    // isAboveCutoff: every replicate of a group must have a nonzero count above its cutoff.
    let group_passes = |w: usize, range: std::ops::Range<usize>| {
        range.into_iter().all(|j| {
            let c = counts.bins[j].window(w, bpw);
            c != 0 && c as f64 > cutoff[j]
        })
    };

    let mut result = ChromResult::default();
    // WindowCoverage: bp of passing windows, counting overlaps once.
    let mut moving_cov = p.window;
    let mut open: Option<Open> = None;

    for w in 0..n_windows {
        if group_passes(w, 0..n_tr) || group_passes(w, n_tr..samples.len()) {
            result.window_cov += moving_cov;
            moving_cov = 0;
            let raw: Vec<f64> = counts.bins.iter().map(|b| b.window(w, bpw) as f64).collect();
            let (ws, wnorm) = counts.stat(&raw);
            if ws.pval < p.pval {
                let win_sta = w as u64 * p.step + 1;
                let win_end = w as u64 * p.step + p.window;
                // isWithinGap: winSta - end - 1 <= gap, same direction.
                let within = open.as_ref().is_some_and(|r| {
                    win_sta as i64 - r.end as i64 - 1 <= p.gap as i64 && ws.event == r.dirn
                });
                if within {
                    let r = open.as_mut().unwrap();
                    r.end = win_end;
                    r.stat = None;
                    if ws.pval < r.win_stat.pval {
                        r.win_sta = win_sta;
                        r.win_stat = ws;
                    }
                } else {
                    if let Some(r) = open.take() {
                        result.regions.push(close(r, ci, &counts));
                    }
                    open = Some(Open {
                        start: win_sta,
                        end: win_end,
                        dirn: ws.event,
                        stat: Some((ws, wnorm)),
                        win_sta,
                        win_stat: ws,
                    });
                }
            }
        }
        // WindowCoverage::move: add one step, but never more than a window.
        if moving_cov < p.window {
            moving_cov += p.step;
        }
    }
    if let Some(r) = open.take() {
        result.regions.push(close(r, ci, &counts));
    }
    result
}

/// Finish a region: a multi-window region is re-tested on its total counts.
fn close(r: Open, ci: usize, counts: &ChromCounts) -> Region {
    let (stat, norm_counts) = match r.stat {
        Some(s) => s,
        None => counts.region_stat(r.start, r.end),
    };
    Region { chrom: ci, start: r.start, end: r.end, norm_counts, stat, win_sta: r.win_sta, win_stat: r.win_stat }
}
