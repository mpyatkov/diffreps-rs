//! Small statistics helpers, ported from MyBioinfo/Common.pm with the same
//! definitions (these details change results, so they are spelled out).

/// Arithmetic mean (`mean`, Common.pm:371). Panics on empty input.
pub fn mean(v: &[f64]) -> f64 {
    assert!(!v.is_empty(), "mean of empty slice");
    v.iter().sum::<f64>() / v.len() as f64
}

/// Sample variance with n-1 denominator (`var`, Common.pm:430).
pub fn var(v: &[f64]) -> f64 {
    let m = mean(v);
    v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (v.len() as f64 - 1.0)
}

/// Median; for an even count, the mean of the two middle values
/// (`median`, Common.pm:440).
pub fn median(v: &[f64]) -> f64 {
    assert!(!v.is_empty(), "median of empty slice");
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    let n = s.len();
    if n % 2 == 1 { s[n / 2] } else { (s[n / 2] + s[n / 2 - 1]) / 2.0 }
}

/// Right-trimmed mean (`mean_r($v, 0, $alpha)`, Common.pm:380): sort, drop
/// the largest `floor(alpha * n)` values, average the rest.
pub fn mean_right_trimmed(v: &[f64], alpha: f64) -> f64 {
    if v.is_empty() {
        eprintln!("Empty array encountered! Return zero.");
        return 0.0;
    }
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    let keep = s.len() - (alpha * s.len() as f64).floor() as usize;
    s[..keep].iter().sum::<f64>() / keep as f64
}

/// Median absolute deviation around `center`, times 1.4826
/// (`mad`, Common.pm:448). Returns 0 for empty input.
pub fn mad(v: &[f64], center: f64) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let dev: Vec<f64> = v.iter().map(|x| (x - center).abs()).collect();
    1.4826 * median(&dev)
}

/// Geometric mean `exp(mean(log x))`; 0 if any value is 0
/// (`geomean`, Common.pm:417).
pub fn geomean(v: &[f64]) -> f64 {
    let mut s = 0.0;
    for &x in v {
        if x == 0.0 {
            return 0.0;
        }
        s += x.ln();
    }
    (s / v.len() as f64).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers() {
        assert_eq!(mean(&[1.0, 2.0, 6.0]), 3.0);
        assert_eq!(var(&[1.0, 2.0, 6.0]), 7.0);
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 2.0, 3.0]), 2.5);
        // 20 values 1..=20, alpha 0.05 -> drop 1 value (20) -> mean of 1..=19 = 10
        let v: Vec<f64> = (1..=20).map(f64::from).collect();
        assert_eq!(mean_right_trimmed(&v, 0.05), 10.0);
        // alpha 0.05 of 19 values drops floor(0.95) = 0 values
        assert_eq!(mean_right_trimmed(&v[..19], 0.05), 10.0);
        assert_eq!(mad(&[1.0, 2.0, 3.0, 4.0, 100.0], 3.0), 1.4826);
        assert!((geomean(&[1.0, 4.0]) - 2.0).abs() < 1e-15);
        assert_eq!(geomean(&[0.0, 4.0]), 0.0);
    }
}
