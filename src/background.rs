//! Step 2 of diffReps: per-sample background level and normalization factors.
//!
//! Port of `mean_std_norm` (diffReps.pl:594). Everything here works on bins
//! of size **window** (not step): diffReps creates its SRBed objects with
//! `bin_size = window` and only switches to `step` for the scan.
//!
//! # Background: mean and deviation of bin counts
//! * robust (default): draw a random subset of about 100,000 bins (the same
//!   bins for every sample); mean = right-trimmed mean (drop the top `alpha`
//!   fraction), deviation = MAD around that mean, scaled by 1.4826.
//! * classic (`--std`): mean and population SD over all bins.
//!
//! A window is later tested only if every replicate of one group has
//! `count > mean + nsd * deviation`.
//!
//! # Normalization factors (when `--norm` is not given)
//! DESeq-style median of ratios: for every bin where all samples pass
//! `count >= mean + nsd * deviation` (only with `--nrpass 1`, the default)
//! and all counts are > 0, take `count_j / geometric_mean(bin)`; the factor of
//! sample `j` is the median of those ratios. Counts are later divided by it.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::bed::{BinCounts, Sample};
use crate::inputs::Chrom;
use crate::stats::{geomean, mad, mean_right_trimmed, median};

#[derive(Debug, Clone, PartialEq)]
pub struct Background {
    /// One value per sample, treatment samples first.
    pub mean: Vec<f64>,
    pub dev: Vec<f64>,
}

impl Background {
    /// Count cutoff of sample `j`: `mean + nsd * dev`.
    pub fn cutoff(&self, j: usize, nsd: f64) -> f64 {
        self.mean[j] + nsd * self.dev[j]
    }
}

/// How the background is estimated.
#[derive(Debug, Clone, Copy)]
pub enum Method {
    /// Trimmed mean + MAD on a random sample of bins.
    Robust { alpha: f64, sample_size: u64, seed: u64 },
    /// Mean and SD over all bins.
    Classic,
}

/// Bins of size `window` for every sample on one chromosome.
fn chrom_bins(samples: &[Sample], ci: usize, chrom: &Chrom, window: u64) -> Vec<BinCounts> {
    samples.iter().map(|s| BinCounts::new(&s.pos2[ci], window, chrom.len)).collect()
}

pub fn estimate_background(samples: &[Sample], chroms: &[Chrom], window: u64, method: Method) -> Background {
    let n = samples.len();
    let total_bins: u64 = chroms.iter().map(|c| c.len.div_ceil(window)).sum();
    match method {
        Method::Classic => {
            // diffReps sums count and count^2 over all bins, then
            // mean = sum / nbin, sd = sqrt(sumsq / nbin - mean^2).
            let mut sum = vec![0f64; n];
            let mut sumsq = vec![0f64; n];
            for (ci, chrom) in chroms.iter().enumerate() {
                for (j, b) in chrom_bins(samples, ci, chrom, window).iter().enumerate() {
                    for i in 0..b.n_bins() {
                        let c = b.bin(i) as f64;
                        sum[j] += c;
                        sumsq[j] += c * c;
                    }
                }
            }
            let nb = total_bins as f64;
            let mean: Vec<f64> = sum.iter().map(|s| s / nb).collect();
            let dev = sumsq.iter().zip(&mean).map(|(sq, m)| (sq / nb - m * m).sqrt()).collect();
            Background { mean, dev }
        }
        Method::Robust { alpha, sample_size, seed } => {
            // Each bin is kept with probability sample_size / total_bins,
            // visiting chromosomes in diffReps order. sample_size = 0 keeps all.
            let rate = if sample_size == 0 { 1.0 } else { sample_size as f64 / total_bins as f64 };
            let mut rng = SplitMix64::new(seed);
            let mut picked: Vec<Vec<f64>> = vec![Vec::new(); n];
            for (ci, chrom) in chroms.iter().enumerate() {
                let bins = chrom_bins(samples, ci, chrom, window);
                for i in 0..bins[0].n_bins() {
                    if sample_size == 0 || rng.next_f64() < rate {
                        for (j, b) in bins.iter().enumerate() {
                            picked[j].push(b.bin(i) as f64);
                        }
                    }
                }
            }
            let mut mean = Vec::with_capacity(n);
            let mut dev = Vec::with_capacity(n);
            for v in &picked {
                let m = mean_right_trimmed(v, alpha);
                mean.push(m);
                dev.push(mad(v, m));
            }
            Background { mean, dev }
        }
    }
}

/// Median-of-ratios normalization factors. `nsd = None` means `--nrpass 0`
/// (use every bin with all counts > 0).
pub fn estimate_norm(
    samples: &[Sample],
    chroms: &[Chrom],
    window: u64,
    bg: &Background,
    nsd: Option<f64>,
) -> Result<Vec<f64>> {
    let n = samples.len();
    let cutoffs: Option<Vec<f64>> = nsd.map(|k| (0..n).map(|j| bg.cutoff(j, k)).collect());
    let mut ratios: Vec<Vec<f64>> = vec![Vec::new(); n];
    let mut counts = vec![0f64; n];
    for (ci, chrom) in chroms.iter().enumerate() {
        let bins = chrom_bins(samples, ci, chrom, window);
        'bin: for i in 0..bins[0].n_bins() {
            for (j, b) in bins.iter().enumerate() {
                let c = b.bin(i) as f64;
                // Perl: `last if $c < mean + nsd*std` (note: `<`, the scan uses `<=`).
                if let Some(cut) = &cutoffs
                    && c < cut[j]
                {
                    continue 'bin;
                }
                counts[j] = c;
            }
            let gm = geomean(&counts);
            if gm > 0.0 {
                for j in 0..n {
                    ratios[j].push(counts[j] / gm);
                }
            }
        }
    }
    if ratios[0].is_empty() {
        bail!("no bins passed the cutoff in all samples; cannot estimate normalization factors (try --nrpass 0 or --norm)");
    }
    Ok(ratios.iter().map(|r| median(r)).collect())
}

/// Read "# Treatment count mean" etc. from a diffReps report header
/// (developer option `--bg-from-report`).
pub fn background_from_report(path: &Path, n_tr: usize, n_co: usize, window: u64) -> Result<Background> {
    let text = fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
    let get = |key: &str| -> Result<Vec<f64>> {
        let prefix = format!("# {key}: ");
        let line = text
            .lines()
            .find(|l| l.starts_with(&prefix))
            .with_context(|| format!("'{key}' not found in {}", path.display()))?;
        line[prefix.len()..]
            .split('\t')
            .map(|v| v.trim().parse::<f64>().with_context(|| format!("bad number in '{line}'")))
            .collect()
    };
    let bin = get("Read count statistics based on bin size")?;
    if bin.first() != Some(&(window as f64)) {
        bail!("report statistics use bin size {:?}, but --window is {window}", bin);
    }
    let (tm, cm) = (get("Treatment count mean")?, get("Control count mean")?);
    let (td, cd) = (get("Treatment count deviation")?, get("Control count deviation")?);
    if tm.len() != n_tr || td.len() != n_tr || cm.len() != n_co || cd.len() != n_co {
        bail!("sample counts in {} do not match the command line", path.display());
    }
    Ok(Background { mean: [tm, cm].concat(), dev: [td, cd].concat() })
}

/// SplitMix64: a tiny, well-mixed 64-bit generator. Plenty for picking bins.
pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        SplitMix64(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1) with 53 random bits.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bed::parse_bed;

    fn toy() -> (Vec<Sample>, Vec<Chrom>) {
        let chroms = vec![Chrom { name: "chr1".into(), len: 400 }];
        // Two samples, 4 bins of 100 bp. Sample A: 2,4,0,6 reads; sample B: 1,2,0,3.
        let mk = |counts: &[usize]| {
            let mut bed = String::new();
            for (b, &c) in counts.iter().enumerate() {
                for _ in 0..c {
                    bed += &format!("chr1\t0\t{}\t.\t0\t.\n", b * 100 + 50);
                }
            }
            parse_bed(bed.as_bytes(), &chroms, 0).unwrap()
        };
        (vec![mk(&[2, 4, 0, 6]), mk(&[1, 2, 0, 3])], chroms)
    }

    #[test]
    fn classic_background() {
        let (s, c) = toy();
        let bg = estimate_background(&s, &c, 100, Method::Classic);
        assert_eq!(bg.mean, [3.0, 1.5]);
        // population SD of 2,4,0,6 = sqrt(56/4 - 9) = sqrt(5)
        assert!((bg.dev[0] - 5f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn robust_background_all_bins() {
        let (s, c) = toy();
        let m = Method::Robust { alpha: 0.0, sample_size: 0, seed: 1 };
        let bg = estimate_background(&s, &c, 100, m);
        assert_eq!(bg.mean, [3.0, 1.5]);
        // |2-3|,|4-3|,|0-3|,|6-3| = 1,1,3,3 -> median 2 -> 2*1.4826
        assert!((bg.dev[0] - 2.0 * 1.4826).abs() < 1e-12);
    }

    #[test]
    fn norm_is_median_of_ratios() {
        let (s, c) = toy();
        let bg = Background { mean: vec![0.0, 0.0], dev: vec![0.0, 0.0] };
        // Bins with all counts > 0: (2,1), (4,2), (6,3); ratio A/gm = sqrt(2) each.
        let f = estimate_norm(&s, &c, 100, &bg, Some(2.0)).unwrap();
        assert!((f[0] - 2f64.sqrt()).abs() < 1e-12 && (f[1] - 0.5f64.sqrt()).abs() < 1e-12);
        // A cutoff that removes every bin is an error.
        let high = Background { mean: vec![100.0, 0.0], dev: vec![0.0, 0.0] };
        assert!(estimate_norm(&s, &c, 100, &high, Some(2.0)).is_err());
    }
}
