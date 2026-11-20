//! Multiple testing correction: port of `padjBH` (Common.pm:47).
//!
//! Benjamini-Hochberg with a caller-supplied number of tests `n`. diffReps
//! uses `n = floor(window coverage / window size)`, i.e. the number of
//! non-overlapping windows that passed the count cutoff, rather than the
//! number of regions or of (overlapping) windows tested.
//!
//! The algorithm, literally:
//! 1. rank p-values ascending; `q_r = min(1, p_r * n / r)`;
//! 2. walk the q values in ascending order; whenever a q belongs to a higher
//!    raw rank than seen so far, copy it down to all lower raw ranks not yet
//!    fixed. This is the usual step-up `q_r = min_{s >= r} q_s` ...
//!
//! ... except for one quirk we keep for parity: if the *smallest* p-value
//! also has the smallest q, step 2 skips it (its rank 0 is not "higher"
//! than the initial 0) and the next copy overwrites it with a larger q.
//! Standard BH would leave it alone. Only the top-ranked region is affected.

/// `padjBH(\@p, $n)`. Ties are ordered by input position (Perl's order is
/// hash-key order, which is random, but tied p-values end up with equal q
/// either way).
pub fn padj_bh(p: &[f64], n: f64) -> Vec<f64> {
    let len = p.len();
    let mut q = p.to_vec();
    let mut sorted_r: Vec<usize> = (0..len).collect();
    sorted_r.sort_by(|&a, &b| p[a].total_cmp(&p[b]));
    let mut rev_r = vec![0usize; len];
    for r in 1..=len {
        let i = sorted_r[r - 1];
        let f = q[i] * n / r as f64;
        q[i] = if f > 1.0 { 1.0 } else { f };
        rev_r[i] = r - 1;
    }
    let mut sorted_rr: Vec<usize> = (0..len).collect();
    sorted_rr.sort_by(|&a, &b| q[a].total_cmp(&q[b]));
    let mut sta_r = 0;
    let mut raw_r = 0;
    for &k in &sorted_rr {
        if rev_r[k] > raw_r {
            raw_r = rev_r[k];
            let v = q[k];
            for &j in &sorted_r[sta_r..raw_r] {
                q[j] = v;
            }
            sta_r = raw_r + 1;
        }
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Plain BH (R's p.adjust(method = "BH", n = n)).
    fn bh(p: &[f64], n: f64) -> Vec<f64> {
        let mut idx: Vec<usize> = (0..p.len()).collect();
        idx.sort_by(|&a, &b| p[b].total_cmp(&p[a]));
        let mut q = vec![0.0; p.len()];
        let mut run = f64::INFINITY;
        for (k, &i) in idx.iter().enumerate() {
            let rank = (p.len() - k) as f64;
            run = run.min((p[i] * n / rank).min(1.0));
            q[i] = run;
        }
        q
    }

    #[test]
    fn equals_bh_when_top_p_is_not_the_smallest_q() {
        let p = [0.01, 0.04, 0.03, 0.2, 0.011];
        assert_eq!(padj_bh(&p, 10.0), bh(&p, 10.0));
        // R: p.adjust(c(0.01, 0.04, 0.03, 0.2, 0.011), "BH", n = 10)
        //    = 0.055 0.100 0.100 0.400 0.055
        let want = [0.055, 0.1, 0.1, 0.4, 0.055];
        for (a, b) in padj_bh(&p, 10.0).iter().zip(want) {
            assert!((a - b).abs() < 1e-15);
        }
    }

    #[test]
    fn diffreps_quirk_for_the_top_p_value() {
        // Raw ranks: 1e-9 -> 4e-9 (the smallest q), 0.2 -> 0.4, 0.25 -> 0.333,
        // 0.3 -> 0.3. Standard BH keeps 4e-9 for the top p-value; diffReps
        // overwrites it with 0.3 when it copies 0.3 down to the lower ranks.
        let p = [1e-9, 0.3, 0.2, 0.25];
        let ours = padj_bh(&p, 4.0);
        let std = bh(&p, 4.0);
        assert_eq!(std[0], 4e-9);
        assert!(ours[0] > std[0], "{ours:?}");
        assert_eq!(&ours[1..], &std[1..]);
    }
}
