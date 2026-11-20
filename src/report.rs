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

/// Perl's default number-to-string conversion (`%.15g`), used for p-values.
pub fn fmt_g15(v: f64) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    if v == 0.0 {
        return if v.is_sign_negative() { "-0".into() } else { "0".into() };
    }
    // Exponent X of the value rounded to 15 significant digits.
    let e_form = format!("{v:.14e}");
    let (mantissa, exp) = e_form.split_once('e').unwrap();
    let x: i32 = exp.parse().unwrap();
    if !(-4..15).contains(&x) {
        let m = trim_fraction(mantissa);
        let sign = if x < 0 { '-' } else { '+' };
        format!("{m}e{sign}{:02}", x.abs())
    } else {
        trim_fraction(&format!("{v:.*}", (14 - x) as usize)).to_string()
    }
}

fn trim_fraction(s: &str) -> &str {
    if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.') } else { s }
}

/// The table header line.
pub const TABLE_HEADER: &str = "Chrom\tStart\tEnd\tLength\tTreatment.cnt\tControl.cnt\tTreatment.avg\tControl.avg\t\
Treatment.enr\tControl.enr\tEvent\tlog2FC\tpval\tpadj\twinSta\twinEnd\twinFC\twinP\twinQ";

/// One table row (`print_reg`, ChromaModSite.pm:114).
pub fn write_region<W: Write>(
    out: &mut W,
    chrom: &str,
    r: &crate::scan::Region,
    n_tr: usize,
    window: u64,
    padj: f64,
    win_q: f64,
) -> Result<()> {
    let (tr, co) = r.norm_counts.split_at(n_tr);
    writeln!(
        out,
        "{chrom}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNA\tNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        r.start,
        r.end,
        r.end - r.start + 1,
        f2(tr).join(";"),
        f2(co).join(";"),
        fmt_f2(r.stat.tr_mean),
        fmt_f2(r.stat.co_mean),
        r.stat.event.as_str(),
        fmt_f2(r.stat.log2fc),
        fmt_g15(r.stat.pval),
        fmt_g15(padj),
        r.win_sta,
        r.win_sta + window - 1,
        fmt_f2(r.win_stat.log2fc),
        fmt_g15(r.win_stat.pval),
        fmt_g15(win_q),
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g15_like_perl() {
        // Strings taken from diffReps reports.
        for s in ["0.000203170734106721", "2.90676945224698e-09", "0.00767328870850612",
                  "3.16430352616811e-06", "1.84249089429946e-05", "0.0001"] {
            assert_eq!(fmt_g15(s.parse().unwrap()), s);
        }
        assert_eq!(fmt_g15(1.0), "1");
        assert_eq!(fmt_g15(0.5), "0.5");
        assert_eq!(fmt_g15(123456789012345678.0), "1.23456789012346e+17");
        assert_eq!(fmt_g15(1e-300), "1e-300");
    }

    #[test]
    fn f2_like_perl_sprintf() {
        assert_eq!(fmt_f2(f64::INFINITY), "Inf");
        assert_eq!(fmt_f2(f64::NEG_INFINITY), "-Inf");
        assert_eq!(fmt_f2(-0.38), "-0.38");
        assert_eq!(fmt_f2(0.125), "0.12"); // round-half-even on the binary value, like glibc
    }
}
