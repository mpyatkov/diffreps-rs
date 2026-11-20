//! The report file.
//!
//! The layout copies diffReps: `#` header lines written with Perl's `spit()`
//! (`# description: v1<TAB>v2...`), then a 19-column table.

use std::io::Write;

use anyhow::Result;

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
