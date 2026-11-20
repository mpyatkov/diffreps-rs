//! diffreps: command-line entry point. The pipeline follows diffReps.pl:
//! read BED files -> background and normalization -> sliding-window scan ->
//! Benjamini-Hochberg -> report.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;
use rayon::prelude::*;

use diffreps::background::{self, Method};
use diffreps::cli::Args;
use diffreps::padj::padj_bh;
use diffreps::scan::{ScanParams, scan_chrom};
use diffreps::{bed, inputs, report};

fn main() -> Result<()> {
    let full_command = std::env::args().collect::<Vec<_>>().join(" ");
    let cfg = Args::parse().resolve(full_command)?;
    rayon::ThreadPoolBuilder::new().num_threads(cfg.nproc).build_global()?;

    let chroms = inputs::read_chrlen(&cfg.chrlen)?;
    let (n_tr, n_co) = (cfg.treatment.len(), cfg.control.len());
    let user_norm = match &cfg.norm {
        Some(p) => Some(inputs::read_norm(p, n_tr, n_co)?),
        None => None,
    };

    let file = File::create(&cfg.report)
        .with_context(|| format!("cannot create report file {}", cfg.report.display()))?;
    let mut out = BufWriter::new(file);
    report::write_parameters(&mut out, &cfg)?;
    out.flush()?;

    // Step 1: read all BED files (one task per file).
    let t = Instant::now();
    let files: Vec<_> = cfg.treatment.iter().chain(&cfg.control).collect();
    let samples = files
        .par_iter()
        .map(|p| bed::read_sample(p, &chroms, cfg.fragsize))
        .collect::<Result<Vec<_>>>()?;
    for s in &samples {
        if s.unknown_chrom + s.malformed + s.negative > 0 {
            eprintln!(
                "{}: skipped {} reads on chromosomes not in --chrlen, {} malformed lines, {} negative positions",
                s.path.display(),
                s.unknown_chrom,
                s.malformed,
                s.negative
            );
        }
    }
    report::spit(&mut out, &format!("# Split BED files done in {} seconds.", t.elapsed().as_secs()), &[])?;
    out.flush()?;

    // Step 2: background level (mean, deviation) and normalization factors.
    let t = Instant::now();
    let bg = match &cfg.bg_from_report {
        Some(p) => background::background_from_report(p, n_tr, n_co, cfg.window)?,
        None => {
            let method = if cfg.robust {
                Method::Robust { alpha: cfg.alpha, sample_size: cfg.sample_size, seed: cfg.seed }
            } else {
                Method::Classic
            };
            background::estimate_background(&samples, &chroms, cfg.window, method)
        }
    };
    let norm: Vec<f64> = match user_norm {
        Some((tr, co)) => [tr, co].concat(),
        None => {
            let nsd = cfg.nrpass.then_some(cfg.nsd);
            background::estimate_norm(&samples, &chroms, cfg.window, &bg, nsd)?
        }
    };
    // As in diffReps, this time includes the normalization factors.
    report::spit(
        &mut out,
        &format!("# Estimate mean and deviation for read counts done in {} seconds", t.elapsed().as_secs()),
        &[],
    )?;
    report::write_estimates(&mut out, cfg.window, n_tr, &bg, &norm)?;
    out.flush()?;

    // Step 3: scan every chromosome (one task per chromosome).
    let t = Instant::now();
    let params = ScanParams { window: cfg.window, step: cfg.step, gap: cfg.gap, nsd: cfg.nsd, pval: cfg.pval };
    let per_chrom: Vec<_> = chroms
        .par_iter()
        .enumerate()
        .map(|(ci, c)| scan_chrom(&samples, n_tr, ci, c, &params, &bg, &norm))
        .collect();
    let coverage: u64 = per_chrom.iter().map(|r| r.window_cov).sum();
    let regions: Vec<_> = per_chrom.into_iter().flat_map(|r| r.regions).collect();

    // Step 4: Benjamini-Hochberg with N = covered bp / window.
    let n = (coverage / cfg.window) as f64;
    let pvals: Vec<f64> = regions.iter().map(|r| r.stat.pval).collect();
    let win_pvals: Vec<f64> = regions.iter().map(|r| r.win_stat.pval).collect();
    let padj = padj_bh(&pvals, n);
    let win_q = padj_bh(&win_pvals, n);

    report::spit(&mut out, "## FDR control information ##", &[])?;
    report::spit(&mut out, "# Genome coverage for windows that passed cutoff", &[coverage.to_string()])?;
    report::spit(&mut out, "# N used in BH adjustment", &[(coverage / cfg.window).to_string()])?;
    report::spit(
        &mut out,
        &format!("# Detecting differential sites done in {} seconds", t.elapsed().as_secs()),
        &[],
    )?;
    writeln!(out, "{}", report::TABLE_HEADER)?;
    for (i, r) in regions.iter().enumerate() {
        report::write_region(&mut out, &chroms[r.chrom].name, r, n_tr, cfg.window, padj[i], win_q[i])?;
    }
    out.flush()?;
    eprintln!("{} differential regions written to {}", regions.len(), cfg.report.display());
    Ok(())
}
