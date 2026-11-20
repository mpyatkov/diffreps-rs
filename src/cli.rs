//! Command-line interface.
//!
//! The option names mirror `diffReps.pl` so existing pipelines keep working.
//! Perl's `Getopt::Long` accepts any unambiguous prefix of a long option
//! (the pipelines use `--frag` for `--fragsize`); clap's `infer_long_args`
//! gives the same behaviour.
//!
//! Annotation (`--noanno`, `--gname`) and hotspots (`--nohs`) are not ported.
//! Background samples (`--btr`, `--bco`, `--bkg`) and the non-NB tests
//! (`--meth tt|gt|cs`) are not ported either.

use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Parser;

/// Raw arguments as typed by the user. Turned into a validated [`Config`]
/// by [`Args::resolve`].
#[derive(Parser, Debug)]
#[command(
    name = "diffreps",
    version,
    about = "Detect differential sites from ChIP-seq/ATAC-seq with biological replicates \
             (Rust port of diffReps 1.55.6, negative binomial test only).",
    infer_long_args = true
)]
pub struct Args {
    /// Treatment group BED6 files (one per replicate).
    #[arg(long, num_args = 1.., required = true)]
    pub treatment: Vec<PathBuf>,

    /// Control group BED6 files (one per replicate).
    #[arg(long, num_args = 1.., required = true)]
    pub control: Vec<PathBuf>,

    /// Output report file.
    #[arg(long)]
    pub report: PathBuf,

    /// Chromosome length file: two whitespace-separated columns, name and length.
    #[arg(long)]
    pub chrlen: PathBuf,

    /// Scanning mode, used only when --window and --step are not given:
    /// (p)eak = 1000/100, (n)ucleosome = 200/20, (b)lock = 10000/1000.
    #[arg(long, default_value = "peak")]
    pub mode: String,

    /// Window size in bp (default: from --mode).
    #[arg(long)]
    pub window: Option<u64>,

    /// Window moving step in bp (default: window / 10).
    #[arg(long)]
    pub step: Option<u64>,

    /// Gap in bp allowed between two consecutive significant windows.
    #[arg(long, default_value_t = 0)]
    pub gap: u64,

    /// Cutoff for low counts, in deviations above the mean:
    /// (b)road = 2, (s)harp = 20, or a number.
    #[arg(long, default_value = "broad")]
    pub nsd: String,

    /// Use the classic mean/standard deviation instead of the robust
    /// right-trimmed mean and median absolute deviation.
    #[arg(long)]
    pub std: bool,

    /// Right-trim fraction for the robust mean, in [0, 0.5).
    #[arg(long, default_value_t = 0.05)]
    pub alpha: f64,

    /// P-value cutoff for significant windows.
    #[arg(long, default_value_t = 1e-4)]
    pub pval: f64,

    /// Normalization constants file (lines `treatment f1 f2 ...` and
    /// `control f1 f2 ...`). Counts are DIVIDED by these factors.
    /// Default: estimate them from the data.
    #[arg(long)]
    pub norm: Option<PathBuf>,

    /// When estimating norm factors, use only bins that pass the --nsd cutoff (1) or all bins (0).
    #[arg(long, default_value_t = 1)]
    pub nrpass: u8,

    /// Fragment size used to shift read positions (half of it is added to
    /// `start` on the + strand, subtracted from `end` otherwise).
    #[arg(long, default_value_t = 100)]
    pub fragsize: u64,

    /// Number of threads.
    #[arg(long, default_value_t = 1)]
    pub nproc: usize,

    /// Statistical test. Only `nb` (negative binomial) is implemented.
    #[arg(long, default_value = "nb")]
    pub meth: String,

    /// Seed for the random bin sample used by the robust mean/deviation estimate.
    /// diffReps does not seed this, so its runs are not reproducible; we are.
    #[arg(long, default_value_t = 20131204)]
    pub seed: u64,

    /// Expected number of bins sampled for the robust estimate (diffReps uses 100000).
    /// 0 = use every bin (deterministic, no sampling).
    #[arg(long, default_value_t = 100_000)]
    pub sample_size: u64,

    /// Developer option: take the treatment/control count mean and deviation
    /// from the header of an existing diffReps report instead of estimating them.
    /// Removes the random sampling, so results can be compared with that report.
    #[arg(long, value_name = "REPORT")]
    pub bg_from_report: Option<PathBuf>,
}

/// Validated settings used by the rest of the program.
#[derive(Debug, Clone)]
pub struct Config {
    pub treatment: Vec<PathBuf>,
    pub control: Vec<PathBuf>,
    pub report: PathBuf,
    pub chrlen: PathBuf,
    pub window: u64,
    pub step: u64,
    pub gap: u64,
    pub nsd: f64,
    pub robust: bool,
    pub alpha: f64,
    pub pval: f64,
    pub norm: Option<PathBuf>,
    pub nrpass: bool,
    pub fragsize: u64,
    pub nproc: usize,
    pub seed: u64,
    pub sample_size: u64,
    pub bg_from_report: Option<PathBuf>,
    /// The command line, recorded in the report header.
    pub full_command: String,
}

impl Args {
    pub fn resolve(self, full_command: String) -> Result<Config> {
        if self.meth != "nb" {
            bail!("only --meth nb (negative binomial) is implemented, got '{}'", self.meth);
        }
        // diffReps.pl: "To use Negative Binomial, you must have at least two
        // replicates per condition." The NB test needs a sample variance.
        if self.treatment.len() < 2 || self.control.len() < 2 {
            bail!("the negative binomial test needs at least two replicates per group");
        }
        if !(0.0..0.5).contains(&self.alpha) {
            bail!("--alpha must be in [0, 0.5), got {}", self.alpha);
        }
        if self.nproc == 0 {
            bail!("--nproc must be at least 1");
        }
        let (window, step) = choose_winstep(&self.mode, self.window, self.step)?;
        let nsd = parse_nsd(&self.nsd)?;
        Ok(Config {
            treatment: self.treatment,
            control: self.control,
            report: self.report,
            chrlen: self.chrlen,
            window,
            step,
            gap: self.gap,
            nsd,
            robust: !self.std,
            alpha: self.alpha,
            pval: self.pval,
            norm: self.norm,
            nrpass: self.nrpass != 0,
            fragsize: self.fragsize,
            nproc: self.nproc,
            seed: self.seed,
            sample_size: self.sample_size,
            bg_from_report: self.bg_from_report,
            full_command,
        })
    }
}

/// Port of `choose_winstep` from diffReps.pl.
///
/// * neither given: the preset for `--mode`;
/// * only window given: step = window / 10 if window is a multiple of 10,
///   otherwise step = window;
/// * both given: window must be a multiple of step, otherwise step = window.
///
/// (diffReps treats "only step given" like "only window given" with
/// window = -1, which is nonsense; we reject it instead.)
pub fn choose_winstep(mode: &str, window: Option<u64>, step: Option<u64>) -> Result<(u64, u64)> {
    match (window, step) {
        (None, None) => Ok(match mode {
            "peak" | "p" => (1000, 100),
            "nucleosome" | "n" => (200, 20),
            "block" | "b" => (10000, 1000),
            other => {
                eprintln!("Unspecified scanning mode '{other}'. Use peak instead!");
                (1000, 100)
            }
        }),
        (Some(0), _) | (_, Some(0)) => bail!("--window and --step must be positive"),
        (Some(w), None) => {
            if w % 10 == 0 {
                Ok((w, w / 10))
            } else {
                eprintln!("Window size is not a multiple of 10. Set step size to be equal to window size.");
                Ok((w, w))
            }
        }
        (Some(w), Some(s)) => {
            if w % s != 0 {
                eprintln!("Window size is not a multiple of step size. Set step size to be equal to window size.");
                Ok((w, w))
            } else {
                Ok((w, s))
            }
        }
        (None, Some(_)) => bail!("--step requires --window"),
    }
}

/// `--nsd`: `broad`/`b` = 2, `sharp`/`s` = 20, or a number.
/// diffReps accepts only whole numbers and silently falls back to 2 otherwise;
/// we accept any non-negative number and reject the rest.
pub fn parse_nsd(s: &str) -> Result<f64> {
    match s {
        "broad" | "b" => Ok(2.0),
        "sharp" | "s" => Ok(20.0),
        other => {
            let v: f64 = other
                .parse()
                .with_context(|| format!("--nsd must be broad, sharp or a number, got '{other}'"))?;
            if !(v >= 0.0 && v.is_finite()) {
                bail!("--nsd must be a non-negative number, got '{other}'");
            }
            Ok(v)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn winstep_matches_perl() {
        assert_eq!(choose_winstep("peak", None, None).unwrap(), (1000, 100));
        assert_eq!(choose_winstep("n", None, None).unwrap(), (200, 20));
        assert_eq!(choose_winstep("block", None, None).unwrap(), (10000, 1000));
        assert_eq!(choose_winstep("peak", Some(1000), None).unwrap(), (1000, 100));
        assert_eq!(choose_winstep("peak", Some(200), None).unwrap(), (200, 20));
        assert_eq!(choose_winstep("peak", Some(205), None).unwrap(), (205, 205));
        assert_eq!(choose_winstep("peak", Some(1000), Some(250)).unwrap(), (1000, 250));
        assert_eq!(choose_winstep("peak", Some(1000), Some(300)).unwrap(), (1000, 1000));
        assert!(choose_winstep("peak", None, Some(100)).is_err());
    }

    #[test]
    fn nsd_presets() {
        assert_eq!(parse_nsd("broad").unwrap(), 2.0);
        assert_eq!(parse_nsd("s").unwrap(), 20.0);
        assert_eq!(parse_nsd("3").unwrap(), 3.0);
        assert_eq!(parse_nsd("2.5").unwrap(), 2.5);
        assert!(parse_nsd("wide").is_err());
        assert!(parse_nsd("-1").is_err());
    }

    #[test]
    fn prefix_options_like_getopt_long() {
        // The user's pipelines pass `--frag 0` for `--fragsize 0`.
        let a = Args::try_parse_from([
            "diffreps", "--treatment", "a.bed", "b.bed", "--control", "c.bed", "d.bed",
            "--report", "out", "--chrlen", "sizes", "--window", "1000", "--frag", "0",
        ])
        .unwrap();
        assert_eq!(a.fragsize, 0);
        assert_eq!(a.treatment.len(), 2);
    }
}
