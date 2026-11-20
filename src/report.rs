//! The report file.
//!
//! The layout copies diffReps: `#` header lines written with Perl's `spit()`
//! (`# description: v1<TAB>v2...`), then a 19-column table.

use std::io::Write;

use anyhow::Result;

use crate::background::Background;
use crate::cli::Config;

pub const DIFFREPS_VERSION: &str = "1.55.6";

/// Perl `spit($h, $desc, @values)`.
pub fn spit<W: Write>(out: &mut W, desc: &str, values: &[String]) -> Result<()> {
    if values.is_empty() {
        writeln!(out, "{desc}")?;
    } else {
        writeln!(out, "{desc}: {}", values.join("\t"))?;
    }
    Ok(())
}

fn one<T: ToString>(v: T) -> Vec<String> {
    vec![v.to_string()]
}

/// The "Parameter Settings" block written before any computation starts.
pub fn write_parameters<W: Write>(out: &mut W, cfg: &Config) -> Result<()> {
    let names = |v: &[std::path::PathBuf]| v.iter().map(|p| p.display().to_string()).collect::<Vec<_>>();
    spit(
        out,
        "# diffReps version",
        &one(format!("{DIFFREPS_VERSION} (diffreps-rs {})", env!("CARGO_PKG_VERSION"))),
    )?;
    spit(out, "## Parameter Settings ##", &[])?;
    spit(out, "# Treatment files", &names(&cfg.treatment))?;
    spit(out, "# Control files", &names(&cfg.control))?;
    spit(out, "# Treatment background", &one("None Specified"))?;
    spit(out, "# Control background", &one("None Specified"))?;
    spit(out, "# Output file", &one(cfg.report.display()))?;
    spit(out, "# Statistical Test", &one("Negative Binomial"))?;
    spit(out, "# Use robust estimation", &one(if cfg.robust { "yes" } else { "no" }))?;
    spit(out, "# Number of deviation cutoff", &one(cfg.nsd))?;
    spit(out, "# Alpha in right-trimmed mean", &one(cfg.alpha))?;
    spit(out, "# P-value cutoff", &one(cfg.pval))?;
    spit(out, "# Window size", &one(cfg.window))?;
    spit(out, "# Step size", &one(cfg.step))?;
    spit(out, "# Gap size", &one(cfg.gap))?;
    spit(out, "# Fragment size", &one(cfg.fragsize))?;
    spit(out, "# Number of processors used", &one(cfg.nproc))?;
    spit(out, "# Full command used", &one(&cfg.full_command))?;
    Ok(())
}

/// Perl `fprecision(2, @values)`: `sprintf "%.2f"` each value.
pub fn f2(values: &[f64]) -> Vec<String> {
    values.iter().map(|v| fmt_f2(*v)).collect()
}

/// `sprintf "%.2f"` as Perl prints it, including `Inf`/`-Inf`.
pub fn fmt_f2(v: f64) -> String {
    if v.is_infinite() {
        if v > 0.0 { "Inf".into() } else { "-Inf".into() }
    } else {
        format!("{v:.2}")
    }
}

/// The "Important values estimated from data" block.
pub fn write_estimates<W: Write>(
    out: &mut W,
    window: u64,
    n_tr: usize,
    bg: &Background,
    norm: &[f64],
) -> Result<()> {
    spit(out, "## Important values estimated from data ##", &[])?;
    spit(out, "# Read count statistics based on bin size", &one(window))?;
    spit(out, "# Treatment count mean", &f2(&bg.mean[..n_tr]))?;
    spit(out, "# Control count mean", &f2(&bg.mean[n_tr..]))?;
    spit(out, "# Treatment count deviation", &f2(&bg.dev[..n_tr]))?;
    spit(out, "# Control count deviation", &f2(&bg.dev[n_tr..]))?;
    spit(out, "# Treatment normalization constant", &f2(&norm[..n_tr]))?;
    spit(out, "# Control normalization constant", &f2(&norm[n_tr..]))?;
    Ok(())
}
