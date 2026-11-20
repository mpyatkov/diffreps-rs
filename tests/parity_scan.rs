//! M5 parity: the full scan (cutoff, window tests, region merging, coverage)
//! against all four reference runs. The background (mean/deviation) is taken
//! from each report's header, which removes diffReps' random sampling.

mod common;

use common::*;
use diffreps::background::{Background, estimate_norm};
use diffreps::inputs;
use diffreps::report::fmt_f2;
use diffreps::scan::{ScanParams, scan_chrom};
use rayon::prelude::*;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * b.abs()
}

fn check(run: &str) {
    let r = Report::load(run);
    let (samples, n_tr) = samples();
    let chroms = chroms();
    let bg = Background {
        mean: [r.floats("Treatment count mean"), r.floats("Control count mean")].concat(),
        dev: [r.floats("Treatment count deviation"), r.floats("Control count deviation")].concat(),
    };
    let window = r.number("Window size") as u64;
    let norm = if run.starts_with("RIPPM") {
        let (t, c) = inputs::read_norm(&test_dir().join(run).join("norm.txt"), *n_tr, samples.len() - n_tr).unwrap();
        [t, c].concat()
    } else {
        estimate_norm(samples, &chroms, window, &bg, Some(2.0)).unwrap()
    };
    let p = ScanParams { window, step: r.number("Step size") as u64, gap: 0, nsd: 2.0, pval: 1e-4 };
    let results: Vec<_> = chroms
        .par_iter()
        .enumerate()
        .map(|(ci, c)| scan_chrom(samples, *n_tr, ci, c, &p, &bg, &norm))
        .collect();

    let cov: u64 = results.iter().map(|r| r.window_cov).sum();
    let regions: Vec<_> = results.into_iter().flat_map(|r| r.regions).collect();
    eprintln!(
        "{run}: coverage ours {cov} diffReps {}; regions ours {} diffReps {}",
        r.number("Genome coverage for windows that passed cutoff"),
        regions.len(),
        r.rows.len()
    );
    assert_eq!(cov as f64, r.number("Genome coverage for windows that passed cutoff"));
    assert_eq!(regions.len(), r.rows.len());
    let mut mismatches = 0;
    for (g, row) in regions.iter().zip(&r.rows) {
        let ours_counts: Vec<String> = g.norm_counts.iter().map(|v| format!("{v:.2}")).collect();
        let theirs_counts: Vec<String> = row.tr_cnt.iter().chain(&row.co_cnt).map(|v| format!("{v:.2}")).collect();
        let ok = chroms[g.chrom].name == row.chrom
            && g.start == row.start
            && g.end == row.end
            && ours_counts == theirs_counts
            && g.stat.event.as_str() == row.event
            && fmt_f2(g.stat.log2fc) == row.log2fc
            && close(g.stat.pval, row.pval)
            && g.win_sta == row.win_sta
            && g.win_sta + window - 1 == row.win_end
            && fmt_f2(g.win_stat.log2fc) == row.win_fc
            && close(g.win_stat.pval, row.win_p);
        if !ok {
            mismatches += 1;
            if mismatches <= 5 {
                eprintln!("  diffReps: {}\n  ours:     {} {} {} {:?} {} {:e} win {} {:e}", row.line,
                    chroms[g.chrom].name, g.start, g.end, ours_counts, g.stat.event.as_str(), g.stat.pval,
                    g.win_sta, g.win_stat.pval);
            }
        }
    }
    assert_eq!(mismatches, 0, "{run}: {mismatches} regions differ");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn scan_rippm_1000() {
    check("RIPPM_1000");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn scan_rippm_200() {
    check("RIPPM_200");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn scan_diffreps_1000() {
    check("DIFFREPS_1000");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn scan_diffreps_200() {
    check("DIFFREPS_200");
}
