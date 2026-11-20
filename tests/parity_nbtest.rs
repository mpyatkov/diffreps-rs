//! M4 parity: recompute every region p-value printed by diffReps.
//!
//! Each row gives the normalized counts of the region; for the RIPPM runs the
//! norm factors are exact, so raw counts are `round(cnt * norm)` and the row's
//! `pval` (and `log2FC`) can be recomputed with `calc_stat`. No BED files needed.

mod common;

use common::*;
use diffreps::inputs;
use diffreps::nbtest::calc_stat;
use diffreps::report::fmt_f2;

fn check(run: &str) {
    let r = Report::load(run);
    let (nt, nc) = inputs::read_norm(&test_dir().join(run).join("norm.txt"), 5, 3).unwrap();
    let mut rels: Vec<f64> = Vec::new();
    let mut worst_line = String::new();
    for row in &r.rows {
        let rt: Vec<f64> = row.tr_cnt.iter().zip(&nt).map(|(c, f)| (c * f).round()).collect();
        let rc: Vec<f64> = row.co_cnt.iter().zip(&nc).map(|(c, f)| (c * f).round()).collect();
        let s = calc_stat(&rt, &rc, &nt, &nc);
        assert_eq!(s.event.as_str(), row.event, "{}", row.line);
        assert_eq!(fmt_f2(s.log2fc), row.log2fc, "{}", row.line);
        let rel = (s.pval - row.pval).abs() / row.pval;
        if rel > rels.iter().copied().fold(0.0, f64::max) {
            worst_line = format!("{} -> ours {:e}", row.line, s.pval);
        }
        rels.push(rel);
    }
    rels.sort_by(f64::total_cmp);
    let q = |f: f64| rels[((rels.len() - 1) as f64 * f) as usize];
    let exact = rels.iter().filter(|&&x| x < 1e-12).count();
    eprintln!(
        "{run}: {} rows, {} within 1e-12, median rel diff {:e}, 99% {:e}, max {:e}\n  worst: {worst_line}",
        rels.len(), exact, q(0.5), q(0.99), q(1.0)
    );
}

#[test]
#[ignore = "needs the TEST reports"]
fn region_pvals_rippm_1000() {
    check("RIPPM_1000");
}

#[test]
#[ignore = "needs the TEST reports"]
fn region_pvals_rippm_200() {
    check("RIPPM_200");
}
