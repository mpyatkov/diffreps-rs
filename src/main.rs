
use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};
use clap::Parser;

use diffreps::cli::Args;
use diffreps::background::{self, Method};
use diffreps::{bed, inputs, report};
use rayon::prelude::*;
use std::time::Instant;

fn main() -> Result<()> {
    let full_command = std::env::args().collect::<Vec<_>>().join(" ");
    let cfg = Args::parse().resolve(full_command)?;

    let chroms = inputs::read_chrlen(&cfg.chrlen)?;
    let user_norm = match &cfg.norm {
        Some(p) => Some(inputs::read_norm(p, cfg.treatment.len(), cfg.control.len())?),
        None => None,
    };

    let file = File::create(&cfg.report)
        .with_context(|| format!("cannot create report file {}", cfg.report.display()))?;
    let mut out = BufWriter::new(file);
    report::write_parameters(&mut out, &cfg)?;
    out.flush()?;

    // Step 1: read all BED files (one rayon task per file).
    rayon::ThreadPoolBuilder::new().num_threads(cfg.nproc).build_global()?;
    let t = Instant::now();
    let files: Vec<_> = cfg.treatment.iter().chain(&cfg.control).collect();
    let samples = files
        .par_iter()
        .map(|p| bed::read_sample(p, &chroms, cfg.fragsize))
        .collect::<anyhow::Result<Vec<_>>>()?;
    for s in &samples {
        eprintln!(
            "{}: {} reads ({} unknown chromosome, {} malformed, {} negative position)",
            s.path.display(), s.nread, s.unknown_chrom, s.malformed, s.negative
        );
    }
    report::spit(&mut out, &format!("# Split BED files done in {} seconds.", t.elapsed().as_secs()), &[])?;
    out.flush()?;

    // Step 2: background level (mean, deviation) and normalization factors.
    let t = Instant::now();
    let (n_tr, n_co) = (cfg.treatment.len(), cfg.control.len());
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
    report::spit(
        &mut out,
        &format!("# Estimate mean and deviation for read counts done in {} seconds", t.elapsed().as_secs()),
        &[],
    )?;
    let norm: Vec<f64> = match user_norm {
        Some((tr, co)) => [tr, co].concat(),
        None => {
            let nsd = cfg.nrpass.then_some(cfg.nsd);
            background::estimate_norm(&samples, &chroms, cfg.window, &bg, nsd)?
        }
    };
    report::write_estimates(&mut out, cfg.window, n_tr, &bg, &norm)?;
    out.flush()?;
    Ok(())
}
