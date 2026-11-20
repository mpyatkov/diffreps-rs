//! The exact negative binomial test (port of MyBioinfo/NBTest.pm) and the
//! per-window/per-region statistics (`calc_stat`, ChromaModSite.pm:518).
//!
//! # The test
//! It is the exact test of DESeq (Anders & Huber 2010), as adapted by diffReps.
//! Under H0 both groups share one normalized mean `μ`. Group sums of raw
//! counts are modelled as `K_A ~ NB(mean m_A, var v_A)` and `K_B ~ NB(m_B, v_B)`.
//! Conditioning on the total `kS = k_A + k_B`, the p-value is the probability
//! of all splits `(a, kS - a)` that are at most as likely as the observed one:
//!
//! ```text
//! p = Σ_{a: P(a)P(kS-a) <= P(kA)P(kB)} P(a)P(kS-a)  /  Σ_a P(a)P(kS-a)
//! ```
//!
//! The sum is walked from both ends towards the expected split
//! `kexp = kS * m_A / (m_A + m_B)` (**not rounded**, so probabilities are
//! evaluated at fractional counts), and in flat tails the step grows by 1.5x,
//! approximating skipped terms by a rectangle. We reproduce this walk exactly,
//! including diffReps' `+1e-8` added to every probability, and compute each
//! probability the way diffReps does, as a difference of two CDF values from
//! DCDFLIB. Whether the step grows depends on the rounding noise of that
//! difference, so an "exact" pmf would give different p-values.
//!
//! # Mean and variance of the group sums (`raw_sum_mean`, `raw_sum_var`)
//! With normalized counts `y = raw / f`, pooled mean `μ = (mean_A + mean_B)/2`:
//! ```text
//! m_A = μ · Σ f_i
//! v_A = max(var_A(y) - μ·mean(1/f), 1e-8·μ) · Σ f_i²  +  m_A
//! ```
//! (the "biological" part of the variance, floored, plus the Poisson part).

use std::cell::Cell;

use crate::cdflib;
use crate::stats::{mean, var};

/// diffReps' `EPSILON`.
pub const EPS: f64 = 1e-8;

/// NB parameters from mean and variance (`nb_r_p`): `r = m²/(v-m)`, `p = r/(r+m)`.
fn nb_r_p(m: f64, v: f64, eps: f64) -> (f64, f64) {
    if m == 0.0 {
        return (eps, eps);
    }
    let v = if v < m { m + eps } else { v };
    let r = m * m / (v - m);
    (r, r / (r + m))
}

/// One negative binomial distribution, as diffReps uses it: `(r, p)` with
/// `pnbinom(x) = I_p(r, x + 1)` from DCDFLIB (see [`crate::cdflib`]).
///
/// The walk asks for `F(x)` and `F(x-1)` at neighbouring points, so the last
/// few CDF values are cached (exact `f64` keys; the results are unchanged).
struct Nb {
    r: f64,
    p: f64,
    cache: [Cell<(u64, f64)>; 4],
    next: Cell<usize>,
}

impl Nb {
    fn new(r: f64, p: f64) -> Self {
        // NaN bits never equal a real key, so the cache starts empty.
        let empty = Cell::new((f64::NAN.to_bits(), 0.0));
        Nb { r, p, cache: [empty.clone(), empty.clone(), empty.clone(), empty], next: Cell::new(0) }
    }

    /// `Math::CDF::pnbinom(x, r, p)`; `undef` (invalid arguments) counts as 0 in Perl arithmetic.
    fn cdf(&self, x: f64) -> f64 {
        let key = x.to_bits();
        for c in &self.cache {
            let (k, v) = c.get();
            if k == key {
                return v;
            }
        }
        let v = cdflib::pnbinom(x, self.r, self.p).unwrap_or(0.0);
        let i = self.next.get();
        self.cache[i].set((key, v));
        self.next.set((i + 1) % 4);
        v
    }

    /// Perl `dnbinom(x, r, p)`: `pnbinom(x) - pnbinom(x-1) + 1e-8`; for
    /// `0 < x < 1`, `pnbinom(x) - pnbinom(0) + 1e-8` (x is often fractional).
    fn d(&self, x: f64) -> f64 {
        if x - 1.0 < 0.0 {
            if x == 0.0 {
                self.cdf(x) + 1e-8
            } else {
                self.cdf(x) - self.cdf(0.0) + 1e-8
            }
        } else {
            self.cdf(x) - self.cdf(x - 1.0) + 1e-8
        }
    }
}

/// Port of `nb_join_pval`: walk the splits in `[kl, kr]` of total `ks`.
/// Returns `(total, observed_total)`.
#[allow(clippy::too_many_arguments)]
fn nb_join_pval(mut kl: f64, mut kr: f64, ks: f64, a: &Nb, b: &Nb, probs: f64, eps: f64) -> (f64, f64) {
    let joint = |k: f64| a.d(k) * b.d(ks - k);
    let mut lval = joint(kl);
    let mut rval = joint(kr);
    let mut prevlval = lval;
    let mut prevrval = rval;
    let mut total = lval + rval;
    let mut obstotal = 0.0;
    let mut step: f64 = 1.0;
    if lval <= probs {
        obstotal += lval;
    }
    if rval <= probs {
        obstotal += rval;
    }
    while kl < kr {
        let do_left = if (prevrval - rval).abs() / prevrval > 0.01 {
            true
        } else if (prevlval - lval).abs() / prevlval > 0.01 {
            false
        } else {
            lval > rval
        };
        if do_left {
            prevlval = lval;
            if kl + step > kr {
                step = kr - kl;
            }
            kl += step;
            lval = joint(kl);
            if step == 1.0 {
                total += lval;
            } else {
                total += lval.min(prevlval) * step;
            }
            if lval <= probs {
                if step == 1.0 {
                    obstotal += lval;
                } else if prevlval <= probs {
                    obstotal += lval.max(prevlval) * step;
                } else {
                    obstotal += lval.max(prevlval) * step * ((probs - lval) / (prevlval - lval)).abs();
                }
            }
            if (prevlval - lval).abs() / prevlval < eps {
                step = (step + 1.0).max(step * 1.5);
            }
        } else {
            prevrval = rval;
            if kr - step < kl {
                step = kr - kl;
            }
            kr -= step;
            rval = joint(kr);
            if step == 1.0 {
                total += rval;
            } else {
                total += rval.min(prevrval) * step;
            }
            if rval <= probs {
                if step == 1.0 {
                    obstotal += rval;
                } else if prevrval <= probs {
                    obstotal += rval.max(prevrval) * step;
                } else {
                    // (sic) no abs() on the denominator here, unlike the left branch
                    obstotal += rval.max(prevrval) * step * (probs - rval).abs() / (prevrval - rval);
                }
            }
            if (prevrval - rval).abs() / prevrval < eps {
                step = (step + 1.0).max(step * 1.5);
            }
        }
    }
    if obstotal > total {
        return (total, total);
    }
    (total, obstotal)
}

/// Port of `nb_pval(ka, kb, μ1, var1, μ2, var2, eps)`.
pub fn nb_pval(ka: f64, kb: f64, mu1: f64, var1: f64, mu2: f64, var2: f64, eps: f64) -> f64 {
    let ks = ka + kb;
    let (r1, p1) = nb_r_p(mu1, var1, eps);
    let (r2, p2) = nb_r_p(mu2, var2, eps);
    let a = Nb::new(r1, p1);
    let b = Nb::new(r2, p2);
    let kexp = ks * mu1 / (mu1 + mu2);
    let p = a.d(ka) * b.d(kb);
    let (ldenom, lnume) = nb_join_pval(0.0, kexp, ks, &a, &b, p, eps);
    let kl2 = if kexp + 1.0 > ks { ks } else { kexp + 1.0 };
    let (denominator, numerator) = nb_join_pval(kl2, ks, ks, &a, &b, p, eps);
    if numerator + lnume == 0.0 {
        return p / (denominator + ldenom);
    }
    (numerator + lnume) / (denominator + ldenom)
}

/// `raw_sum_mean`: expected raw group sum under H0.
fn raw_sum_mean(m1: f64, m2: f64, norm: &[f64]) -> f64 {
    norm.iter().sum::<f64>() * (m1 + m2) / 2.0
}

/// `raw_sum_var`: variance of the raw group sum (see module docs).
fn raw_sum_var(base_var: f64, m1: f64, m2: f64, norm: &[f64], raw_mean: f64, eps: f64) -> f64 {
    let base_mean = (m1 + m2) / 2.0;
    let z = norm.iter().map(|f| 1.0 / f).sum::<f64>() / norm.len() as f64;
    let s: f64 = norm.iter().map(|f| f * f).sum();
    let mut v = base_var - z * base_mean;
    if v < eps * base_mean {
        v = eps * base_mean;
    }
    v * s + raw_mean
}

/// Direction of change.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Event {
    Up,
    Down,
    /// Both groups have zero counts (cannot pass the count cutoff).
    Nons,
}

impl Event {
    pub fn as_str(self) -> &'static str {
        match self {
            Event::Up => "Up",
            Event::Down => "Down",
            Event::Nons => "Nons",
        }
    }
}

/// Result of `calc_stat` for a window or a region.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stat {
    /// Normalized means of each group.
    pub tr_mean: f64,
    pub co_mean: f64,
    /// log2(tr_mean / co_mean); +Inf if co_mean = 0, -Inf if tr_mean = 0.
    pub log2fc: f64,
    pub event: Event,
    pub pval: f64,
}

/// Port of `calc_stat` for the `nb` method. `raw_tr`/`raw_co` are raw
/// read counts per replicate; `norm_*` the factors they are divided by.
pub fn calc_stat(raw_tr: &[f64], raw_co: &[f64], norm_tr: &[f64], norm_co: &[f64]) -> Stat {
    let tr: Vec<f64> = raw_tr.iter().zip(norm_tr).map(|(c, f)| c / f).collect();
    let co: Vec<f64> = raw_co.iter().zip(norm_co).map(|(c, f)| c / f).collect();
    let tr_m = mean(&tr);
    let co_m = mean(&co);
    if tr_m == 0.0 && co_m == 0.0 {
        return Stat { tr_mean: 0.0, co_mean: 0.0, log2fc: 0.0, event: Event::Nons, pval: 1.0 };
    }
    // Perl: Inf if co_m == 0; log2(0) = -Inf if tr_m == 0.
    let log2fc = if co_m == 0.0 { f64::INFINITY } else if tr_m == 0.0 { f64::NEG_INFINITY } else { (tr_m / co_m).ln() / 2f64.ln() };
    let event = if log2fc > 0.0 { Event::Up } else { Event::Down };
    let tr_raw_m = raw_sum_mean(tr_m, co_m, norm_tr);
    let co_raw_m = raw_sum_mean(tr_m, co_m, norm_co);
    let tr_raw_v = raw_sum_var(var(&tr), tr_m, co_m, norm_tr, tr_raw_m, EPS);
    let co_raw_v = raw_sum_var(var(&co), tr_m, co_m, norm_co, co_raw_m, EPS);
    let ka: f64 = raw_tr.iter().sum();
    let kb: f64 = raw_co.iter().sum();
    let pval = nb_pval(ka, kb, tr_raw_m, tr_raw_v, co_raw_m, co_raw_v, EPS);
    Stat { tr_mean: tr_m, co_mean: co_m, log2fc, event, pval }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Against the R port of NBTest.pm (tests/data/nb_reference.R), which
    /// uses R's pbeta (TOMS 708, like Perl's Math::CDF).
    #[test]
    fn matches_r_reference() {
        let text = include_str!("../tests/data/nb_reference.tsv");
        let nums = |s: &str| s.split(';').map(|x| x.parse::<f64>().unwrap()).collect::<Vec<_>>();
        let mut worst: f64 = 0.0;
        let mut n = 0;
        for line in text.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            let (rt, rc, nt, nc) = (nums(f[0]), nums(f[1]), nums(f[2]), nums(f[3]));
            let want_fc: f64 = f[4].parse().unwrap();
            let want_p: f64 = f[5].parse().unwrap();
            let s = calc_stat(&rt, &rc, &nt, &nc);
            assert!(s.log2fc == want_fc || (s.log2fc - want_fc).abs() <= 1e-12 * want_fc.abs(), "{line}");
            let rel = (s.pval - want_p).abs() / want_p;
            worst = worst.max(rel);
            assert!(rel < 1e-9, "{line}: got {:e}", s.pval);
            n += 1;
        }
        eprintln!("{n} cases, worst relative p-value difference {worst:e}");
    }

    #[test]
    fn zero_counts() {
        let s = calc_stat(&[0.0, 0.0], &[0.0, 0.0], &[1.0, 1.0], &[1.0, 1.0]);
        assert_eq!((s.event, s.pval), (Event::Nons, 1.0));
        let s = calc_stat(&[5.0, 6.0], &[0.0, 0.0], &[1.0, 1.0], &[1.0, 1.0]);
        assert_eq!((s.event, s.log2fc), (Event::Up, f64::INFINITY));
        let s = calc_stat(&[0.0, 0.0], &[5.0, 6.0], &[1.0, 1.0], &[1.0, 1.0]);
        assert_eq!((s.event, s.log2fc), (Event::Down, f64::NEG_INFINITY));
    }
}
