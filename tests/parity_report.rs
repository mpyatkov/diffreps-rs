//! M5/M6 parity: the full pipeline after background estimation (cutoff,
//! window tests, region merging, coverage, Benjamini-Hochberg, row format)
//! against all four reference runs. The background (mean/deviation) is taken
//! from each report's header, which removes diffReps' random sampling.
//!
//! Every column must be string-identical, except the four p-value columns,
//! which must agree within a relative 1e-13: about 0.3% of diffReps' p-values
//! differ from ours in the last printed digits (see docs/porting-notes.md, M6).

mod common;

use common::*;
use diffreps::background::{Background, estimate_norm};
use diffreps::inputs;
use diffreps::padj::padj_bh;
use diffreps::report::write_region;
use diffreps::scan::{ScanParams, scan_chrom};
use rayon::prelude::*;

const P_COLS: [usize; 4] = [12, 13, 17, 18]; // pval, padj, winP, winQ

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
    let n = (cov / window) as f64;
    let padj = padj_bh(&regions.iter().map(|g| g.stat.pval).collect::<Vec<_>>(), n);
    let winq = padj_bh(&regions.iter().map(|g| g.win_stat.pval).collect::<Vec<_>>(), n);
    let (mut mismatches, mut exact_p, mut worst) = (0, 0, 0f64);
    for (k, (g, row)) in regions.iter().zip(&r.rows).enumerate() {
        let mut buf = Vec::new();
        write_region(&mut buf, &chroms[g.chrom].name, g, *n_tr, window, padj[k], winq[k]).unwrap();
        let ours = String::from_utf8(buf).unwrap();
        let ours: Vec<&str> = ours.trim_end().split('\t').collect();
        let theirs: Vec<&str> = row.line.split('\t').collect();
        let mut ok = ours.len() == theirs.len();
        for (c, (a, b)) in ours.iter().zip(&theirs).enumerate() {
            if P_COLS.contains(&c) {
                let (x, y): (f64, f64) = (a.parse().unwrap(), b.parse().unwrap());
                let rel = (x - y).abs() / y;
                worst = worst.max(rel);
                if a == b { exact_p += 1 }
                ok &= rel <= 1e-13;
            } else {
                ok &= a == b;
            }
        }
        if !ok {
            mismatches += 1;
            if mismatches <= 5 {
                eprintln!("  diffReps: {}\n  ours:     {}", row.line, ours.join("\t"));
            }
        }
    }
    eprintln!(
        "{run}: {} of {} p-value fields identical as printed, worst relative difference {worst:.1e}",
        exact_p,
        4 * regions.len()
    );
    assert_eq!(mismatches, 0, "{run}: {mismatches} regions differ");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn report_rippm_1000() {
    check("RIPPM_1000");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn report_rippm_200() {
    check("RIPPM_200");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn report_diffreps_1000() {
    check("DIFFREPS_1000");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn report_diffreps_200() {
    check("DIFFREPS_200");
}
