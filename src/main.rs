mod cli;
mod inputs;
mod report;

use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};
use clap::Parser;

use cli::Args;

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

    eprintln!(
        "{} chromosomes, window {} / step {}, user norm: {}",
        chroms.len(),
        cfg.window,
        cfg.step,
        user_norm.is_some()
    );
    Ok(())
}
