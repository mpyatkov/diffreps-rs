//! M2 parity: raw read counts of every reported region.
//!
//! diffReps prints normalized counts (`raw / norm`, 2 decimals). For the
//! RIPPM runs the norm factors are known exactly (norm.txt), so
//! `round(cnt * norm)` recovers the raw count, which must equal the sum of
//! our step-size bins over the region.

mod common;

use common::*;
use diffreps::bed::BinCounts;
use diffreps::inputs;

fn check_run(run: &str) {
    let report = Report::load(run);
    let (samples, n_tr) = samples();
    let chroms = chroms();
    let step = report.number("Step size") as u64;
    let (norm_tr, norm_co) =
        inputs::read_norm(&test_dir().join(run).join("norm.txt"), *n_tr, samples.len() - n_tr).unwrap();
    let norms: Vec<f64> = norm_tr.iter().chain(&norm_co).copied().collect();

    let mut checked = 0;
    for (ci, chrom) in chroms.iter().enumerate() {
        let bins: Vec<BinCounts> =
            samples.iter().map(|s| BinCounts::new(&s.pos2[ci], step, chrom.len)).collect();
        for row in report.rows.iter().filter(|r| r.chrom == chrom.name) {
            // Region [start, end] (1-based, inclusive) covers bins (start-1)/step ..= end/step - 1.
            let first = ((row.start - 1) / step) as usize;
            let last = (row.end / step) as usize - 1;
            let printed: Vec<f64> = row.tr_cnt.iter().chain(&row.co_cnt).copied().collect();
            for (k, b) in bins.iter().enumerate() {
                let ours = b.range(first, last) as f64;
                let theirs = (printed[k] * norms[k]).round();
                assert_eq!(ours, theirs, "{run} sample {k}: {}", row.line);
                // And the printed value is our count / norm, rounded to 2 decimals.
                assert_eq!(format!("{:.2}", ours / norms[k]), format!("{:.2}", printed[k]));
            }
            checked += 1;
        }
    }
    assert_eq!(checked, report.rows.len());
    eprintln!("{run}: {checked} regions match");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn region_counts_rippm_1000() {
    check_run("RIPPM_1000");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn region_counts_rippm_200() {
    check_run("RIPPM_200");
}
