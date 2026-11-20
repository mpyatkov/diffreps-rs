//! Incomplete beta function from DCDFLIB (Brown, Lovato & Russell, 1994;
//! public domain), i.e. the `bratio` algorithm of Didonato & Morris,
//! "Significant digit computation of the incomplete beta function ratios",
//! ACM TOMS 708 (1992).
//!
//! # Why a line-by-line port instead of a crate
//! diffReps computes negative binomial probabilities as differences of CDF
//! values, `pnbinom(x) - pnbinom(x-1)`, via Perl's Math::CDF, which wraps
//! this exact C code. In the test's flat tails those differences are pure
//! rounding noise around 1e-16, and the step size of the p-value walk depends
//! on them. Measured on the reference reports: differences of TOMS 708 CDFs
//! reproduce 1496/1500 diffReps p-values within 1e-12, an exact pmf only ~85%.
//! So we need *the same* rounding, not just an accurate value.
//!
//! The port keeps the C arithmetic in the same order (Rust does not
//! reassociate or fuse floating-point operations), uses `f64::exp/ln/powf`,
//! which call the platform libm like the C does, and replaces `goto`s with
//! structured control flow. C `int` conversions truncate, as `as i32` does.
//! Function names are the C names.

/// `spmpar(1)`: machine epsilon `2^(1-53)`.
const EPS_MACH: f64 = f64::EPSILON;

#[inline]
fn fifdmax1(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

#[inline]
fn fifdmin1(a: f64, b: f64) -> f64 {
    if a < b { a } else { b }
}

/// `exparg(l)`: largest (l = 0) or most negative (l != 0) `w` with `exp(w)` representable.
fn exparg(l: i32) -> f64 {
    let lnb = 0.69314718055995e0;
    if l == 0 {
        let m = 1024; // ipmpar(10): emax
        0.99999e0 * (m as f64 * lnb)
    } else {
        let m = -1021 - 1; // ipmpar(9) - 1: emin - 1
        0.99999e0 * (m as f64 * lnb)
    }
}

/// `algdiv(a, b)`: ln(Γ(b)/Γ(a+b)) for b >= 8.
fn algdiv(a: f64, b: f64) -> f64 {
    const C0: f64 = 0.833333333333333e-01;
    const C1: f64 = -0.277777777760991e-02;
    const C2: f64 = 0.793650666825390e-03;
    const C3: f64 = -0.595202931351870e-03;
    const C4: f64 = 0.837308034031215e-03;
    const C5: f64 = -0.165322962780713e-02;
    let (c, x, d);
    if a <= b {
        let h = a / b;
        c = h / (1.0e0 + h);
        x = 1.0e0 / (1.0e0 + h);
        d = b + (a - 0.5e0);
    } else {
        let h = b / a;
        c = 1.0e0 / (1.0e0 + h);
        x = h / (1.0e0 + h);
        d = a + (b - 0.5e0);
    }
    let x2 = x * x;
    let s3 = 1.0e0 + (x + x2);
    let s5 = 1.0e0 + (x + x2 * s3);
    let s7 = 1.0e0 + (x + x2 * s5);
    let s9 = 1.0e0 + (x + x2 * s7);
    let s11 = 1.0e0 + (x + x2 * s9);
    let t = (1.0e0 / b).powf(2.0);
    let mut w = ((((C5 * s11 * t + C4 * s9) * t + C3 * s7) * t + C2 * s5) * t + C1 * s3) * t + C0;
    w *= c / b;
    let t1 = a / b;
    let u = d * alnrel(t1);
    let v = a * (b.ln() - 1.0e0);
    if u <= v { w - u - v } else { w - v - u }
}

/// `alnrel(a)`: ln(1 + a).
fn alnrel(a: f64) -> f64 {
    const P1: f64 = -0.129418923021993e+01;
    const P2: f64 = 0.405303492862024e+00;
    const P3: f64 = -0.178874546012214e-01;
    const Q1: f64 = -0.162752256355323e+01;
    const Q2: f64 = 0.747811014037616e+00;
    const Q3: f64 = -0.845104217945565e-01;
    if a.abs() > 0.375e0 {
        let x = 1.0e0 + a;
        return x.ln();
    }
    let t = a / (a + 2.0e0);
    let t2 = t * t;
    let w = (((P3 * t2 + P2) * t2 + P1) * t2 + 1.0e0) / (((Q3 * t2 + Q2) * t2 + Q1) * t2 + 1.0e0);
    2.0e0 * t * w
}

/// `apser`: I_x(a, b) for a <= min(eps, eps*b), b*x <= 1, x <= 0.5.
fn apser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    const G: f64 = 0.577215664901533e0;
    let bx = b * x;
    let mut t = x - bx;
    let c = if b * eps > 2.0e-2 { bx.ln() + G + t } else { x.ln() + psi(b) + G + t };
    let tol = 5.0e0 * eps * c.abs();
    let mut j = 1.0e0;
    let mut s = 0.0e0;
    loop {
        j += 1.0e0;
        t *= x - bx / j;
        let aj = t / j;
        s += aj;
        if !(aj.abs() > tol) {
            break;
        }
    }
    -(a * (c + s))
}

/// `basym`: asymptotic expansion of I_x(a, b) for large a and b.
fn basym(a: f64, b: f64, lambda: f64, eps: f64) -> f64 {
    const E0: f64 = 1.12837916709551e0;
    const E1: f64 = 0.353553390593274e0;
    const NUM: usize = 20;
    let mut a0 = [0.0f64; 21];
    let mut b0 = [0.0f64; 21];
    let mut c = [0.0f64; 21];
    let mut d = [0.0f64; 21];
    let (h, r0, r1, w0);
    if a >= b {
        h = b / a;
        r0 = 1.0e0 / (1.0e0 + h);
        r1 = (b - a) / a;
        w0 = 1.0e0 / (b * (1.0e0 + h)).sqrt();
    } else {
        h = a / b;
        r0 = 1.0e0 / (1.0e0 + h);
        r1 = (b - a) / b;
        w0 = 1.0e0 / (a * (1.0e0 + h)).sqrt();
    }
    let t1 = -(lambda / a);
    let t2 = lambda / b;
    let f = a * rlog1(t1) + b * rlog1(t2);
    let t = (-f).exp();
    if t == 0.0e0 {
        return 0.0e0;
    }
    let z0 = f.sqrt();
    let z = 0.5e0 * (z0 / E1);
    let z2 = f + f;
    a0[0] = 2.0e0 / 3.0e0 * r1;
    c[0] = -(0.5e0 * a0[0]);
    d[0] = -c[0];
    let mut j0 = 0.5e0 / E0 * erfc1(1, z0);
    let mut j1 = E1;
    let mut sum = j0 + d[0] * w0 * j1;
    let mut s = 1.0e0;
    let h2 = h * h;
    let mut hn = 1.0e0;
    let mut w = w0;
    let mut znm1 = z;
    let mut zn = z2;
    let mut n = 2;
    while n <= NUM {
        hn = h2 * hn;
        a0[n - 1] = 2.0e0 * r0 * (1.0e0 + h * hn) / (n as f64 + 2.0e0);
        let np1 = n + 1;
        s += hn;
        a0[np1 - 1] = 2.0e0 * r1 * s / (n as f64 + 3.0e0);
        for i in n..=np1 {
            let r = -(0.5e0 * (i as f64 + 1.0e0));
            b0[0] = r * a0[0];
            for m in 2..=i {
                let mut bsum = 0.0e0;
                let mm1 = m - 1;
                for j in 1..=mm1 {
                    let mmj = m - j;
                    bsum += (j as f64 * r - mmj as f64) * a0[j - 1] * b0[mmj - 1];
                }
                b0[m - 1] = r * a0[m - 1] + bsum / m as f64;
            }
            c[i - 1] = b0[i - 1] / (i as f64 + 1.0e0);
            let mut dsum = 0.0e0;
            let im1 = i - 1;
            for j in 1..=im1 {
                let imj = i - j;
                dsum += d[imj - 1] * c[j - 1];
            }
            d[i - 1] = -(dsum + c[i - 1]);
        }
        j0 = E1 * znm1 + (n as f64 - 1.0e0) * j0;
        j1 = E1 * zn + n as f64 * j1;
        znm1 = z2 * znm1;
        zn = z2 * zn;
        w = w0 * w;
        let t0 = d[n - 1] * w * j0;
        w = w0 * w;
        let t1 = d[np1 - 1] * w * j1;
        sum += t0 + t1;
        if t0.abs() + t1.abs() <= eps * sum {
            break;
        }
        n += 2;
    }
    let u = (-bcorr(a, b)).exp();
    E0 * t * u * sum
}

/// `bcorr(a0, b0)`: del(a) + del(b) - del(a+b), where ln Γ(a) = (a-½)ln a - a + ½ln 2π + del(a).
fn bcorr(a0: f64, b0: f64) -> f64 {
    const C0: f64 = 0.833333333333333e-01;
    const C1: f64 = -0.277777777760991e-02;
    const C2: f64 = 0.793650666825390e-03;
    const C3: f64 = -0.595202931351870e-03;
    const C4: f64 = 0.837308034031215e-03;
    const C5: f64 = -0.165322962780713e-02;
    let a = fifdmin1(a0, b0);
    let b = fifdmax1(a0, b0);
    let h = a / b;
    let c = h / (1.0e0 + h);
    let x = 1.0e0 / (1.0e0 + h);
    let x2 = x * x;
    let s3 = 1.0e0 + (x + x2);
    let s5 = 1.0e0 + (x + x2 * s3);
    let s7 = 1.0e0 + (x + x2 * s5);
    let s9 = 1.0e0 + (x + x2 * s7);
    let s11 = 1.0e0 + (x + x2 * s9);
    let mut t = (1.0e0 / b).powf(2.0);
    let mut w = ((((C5 * s11 * t + C4 * s9) * t + C3 * s7) * t + C2 * s5) * t + C1 * s3) * t + C0;
    w *= c / b;
    t = (1.0e0 / a).powf(2.0);
    (((((C5 * t + C4) * t + C3) * t + C2) * t + C1) * t + C0) / a + w
}

/// `betaln(a0, b0)`: ln B(a0, b0).
fn betaln(a0: f64, b0: f64) -> f64 {
    const E: f64 = 0.918938533204673e0;
    let mut a = fifdmin1(a0, b0);
    let b = fifdmax1(a0, b0);
    // C label S60, reached with w = 0 (from S30) or w = ln(...) (from S40).
    let s60 = |a: f64, mut b: f64, w: f64| -> f64 {
        let n = (b - 1.0e0) as i32;
        let mut z = 1.0e0;
        for _ in 1..=n {
            b -= 1.0e0;
            z *= b / (a + b);
        }
        w + z.ln() + (gamln(a) + (gamln(b) - gsumln(a, b)))
    };
    if a >= 8.0e0 {
        // S100
        let w = bcorr(a, b);
        let h = a / b;
        let c = h / (1.0e0 + h);
        let u = -((a - 0.5e0) * c.ln());
        let v = b * alnrel(h);
        return if u <= v { -(0.5e0 * b.ln()) + E + w - u - v } else { -(0.5e0 * b.ln()) + E + w - v - u };
    }
    if a < 1.0e0 {
        if b >= 8.0e0 {
            return gamln(a) + algdiv(a, b);
        }
        let t1 = a + b;
        return gamln(a) + (gamln(b) - gamln(t1));
    }
    // 1 <= a < 8
    if a <= 2.0e0 {
        if b <= 2.0e0 {
            return gamln(a) + gamln(b) - gsumln(a, b);
        }
        // S30
        let w = 0.0e0;
        if b < 8.0e0 {
            return s60(a, b, w);
        }
        return gamln(a) + algdiv(a, b);
    }
    // S40: 2 < a < 8
    if b > 1000.0e0 {
        // S80
        let n = (a - 1.0e0) as i32;
        let mut w = 1.0e0;
        for _ in 1..=n {
            a -= 1.0e0;
            w *= a / (1.0e0 + a / b);
        }
        return w.ln() - n as f64 * b.ln() + (gamln(a) + algdiv(a, b));
    }
    let n = (a - 1.0e0) as i32;
    let mut w = 1.0e0;
    for _ in 1..=n {
        a -= 1.0e0;
        let h = a / b;
        w *= h / (1.0e0 + h);
    }
    w = w.ln();
    if b < 8.0e0 {
        return s60(a, b, w);
    }
    w + gamln(a) + algdiv(a, b)
}

/// `bfrac`: continued fraction expansion for I_x(a, b) when a, b > 1.
fn bfrac(a: f64, b: f64, x: f64, y: f64, lambda: f64, eps: f64) -> f64 {
    let mut bfrac = brcomp(a, b, x, y);
    if bfrac == 0.0e0 {
        return bfrac;
    }
    let c = 1.0e0 + lambda;
    let c0 = b / a;
    let c1 = 1.0e0 + 1.0e0 / a;
    let yp1 = y + 1.0e0;
    let mut n = 0.0e0;
    let mut p = 1.0e0;
    let mut s = a + 1.0e0;
    let mut an = 0.0e0;
    let mut bn = 1.0e0;
    let mut anp1 = 1.0e0;
    let mut bnp1 = c / c1;
    let mut r = c1 / c;
    loop {
        n += 1.0e0;
        let mut t = n / a;
        let w = n * (b - n) * x;
        let mut e = a / s;
        let alpha = p * (p + c0) * e * e * (w * x);
        e = (1.0e0 + t) / (c1 + t + t);
        let beta = n + w / s + e * (c + n * yp1);
        p = 1.0e0 + t;
        s += 2.0e0;
        t = alpha * an + beta * anp1;
        an = anp1;
        anp1 = t;
        t = alpha * bn + beta * bnp1;
        bn = bnp1;
        bnp1 = t;
        let r0 = r;
        r = anp1 / bnp1;
        if (r - r0).abs() <= eps * r {
            break;
        }
        an /= bnp1;
        bn /= bnp1;
        anp1 = r;
        bnp1 = 1.0e0;
    }
    bfrac *= r;
    bfrac
}

/// `bgrat`: asymptotic expansion for I_x(a, b) when a is larger than b.
/// Adds the result to `w` (as the C does through its pointer argument).
fn bgrat(a: f64, b: f64, x: f64, y: f64, w: &mut f64, eps: f64) {
    let mut c = [0.0f64; 30];
    let mut d = [0.0f64; 30];
    let bm1 = b - 0.5e0 - 0.5e0;
    let nu = a + 0.5e0 * bm1;
    let lnx = if y > 0.375e0 { x.ln() } else { alnrel(-y) };
    let z = -(nu * lnx);
    if b * z == 0.0e0 {
        return;
    }
    let mut r = b * (1.0e0 + gam1(b)) * (b * z.ln()).exp();
    r *= (a * lnx).exp() * (0.5e0 * bm1 * lnx).exp();
    let mut u = algdiv(b, a) + b * nu.ln();
    u = r * (-u).exp();
    if u == 0.0e0 {
        return;
    }
    let (_p, q) = grat1(b, z, r, eps);
    let v = 0.25e0 * (1.0e0 / nu).powf(2.0);
    let t2 = 0.25e0 * lnx * lnx;
    let l = *w / u;
    let mut j = q / r;
    let mut sum = j;
    let mut t = 1.0e0;
    let mut cn = 1.0e0;
    let mut n2 = 0.0e0;
    for n in 1..=30usize {
        let bp2n = b + n2;
        j = (bp2n * (bp2n + 1.0e0) * j + (z + bp2n + 1.0e0) * t) * v;
        n2 += 2.0e0;
        t *= t2;
        cn /= n2 * (n2 + 1.0e0);
        c[n - 1] = cn;
        let mut s = 0.0e0;
        if n != 1 {
            let nm1 = n - 1;
            let mut coef = b - n as f64;
            for i in 1..=nm1 {
                s += coef * c[i - 1] * d[n - i - 1];
                coef += b;
            }
        }
        d[n - 1] = bm1 * cn + s / n as f64;
        let dj = d[n - 1] * j;
        sum += dj;
        if sum <= 0.0e0 {
            return; // ierr = 1, w unchanged
        }
        if dj.abs() <= eps * (sum + l) {
            break;
        }
    }
    *w += u * sum;
}

/// `bpser`: power series for I_x(a, b) when b <= 1 or b*x <= 0.7.
fn bpser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    let mut bpser = 0.0e0;
    if x == 0.0e0 {
        return bpser;
    }
    let a0 = fifdmin1(a, b);
    if a0 >= 1.0e0 {
        let z = a * x.ln() - betaln(a, b);
        bpser = z.exp() / a;
    } else {
        let mut b0 = fifdmax1(a, b);
        if b0 >= 8.0e0 {
            // S90
            let u = gamln1(a0) + algdiv(a0, b0);
            let z = a * x.ln() - u;
            bpser = a0 / a * z.exp();
        } else if b0 > 1.0e0 {
            // S40
            let mut u = gamln1(a0);
            let m = (b0 - 1.0e0) as i32;
            if m >= 1 {
                let mut c = 1.0e0;
                for _ in 1..=m {
                    b0 -= 1.0e0;
                    c *= b0 / (a0 + b0);
                }
                u = c.ln() + u;
            }
            let z = a * x.ln() - u;
            b0 -= 1.0e0;
            let apb = a0 + b0;
            let t = if apb > 1.0e0 {
                let u = a0 + b0 - 1.0e0;
                (1.0e0 + gam1(u)) / apb
            } else {
                1.0e0 + gam1(apb)
            };
            bpser = z.exp() * (a0 / a) * (1.0e0 + gam1(b0)) / t;
        } else {
            bpser = x.powf(a);
            if bpser == 0.0e0 {
                return bpser;
            }
            let apb = a + b;
            let z = if apb > 1.0e0 {
                let u = a + b - 1.0e0;
                (1.0e0 + gam1(u)) / apb
            } else {
                1.0e0 + gam1(apb)
            };
            let c = (1.0e0 + gam1(a)) * (1.0e0 + gam1(b)) / z;
            bpser *= c * (b / apb);
        }
    }
    // S100
    if bpser == 0.0e0 || a <= 0.1e0 * eps {
        return bpser;
    }
    let mut sum = 0.0e0;
    let mut n = 0.0e0;
    let mut c = 1.0e0;
    let tol = eps / a;
    loop {
        n += 1.0e0;
        c *= (0.5e0 + (0.5e0 - b / n)) * x;
        let w = c / (a + n);
        sum += w;
        if !(w.abs() > tol) {
            break;
        }
    }
    bpser *= 1.0e0 + a * sum;
    bpser
}

/// `bratio(a, b, x, y)`: returns `(I_x(a,b), 1 - I_x(a,b))`; `(0, 0)` on invalid input.
pub fn bratio(a: f64, b: f64, x: f64, y: f64) -> (f64, f64) {
    enum Go {
        S90,
        S100,
        S110,
        S120,
        S130,
        S140,
        S150,
        S160,
        S200,
    }
    let mut eps = EPS_MACH;
    let mut w = 0.0e0;
    let mut w1 = 0.0e0;
    if a < 0.0e0 || b < 0.0e0 || (a == 0.0e0 && b == 0.0e0) {
        return (w, w1);
    }
    if !(0.0e0..=1.0e0).contains(&x) || !(0.0e0..=1.0e0).contains(&y) {
        return (w, w1);
    }
    let z = x + y - 0.5e0 - 0.5e0;
    if z.abs() > 3.0e0 * eps {
        return (w, w1);
    }
    if x == 0.0e0 {
        return if a == 0.0e0 { (w, w1) } else { (0.0e0, 1.0e0) };
    }
    if y == 0.0e0 {
        return if b == 0.0e0 { (w, w1) } else { (1.0e0, 0.0e0) };
    }
    if a == 0.0e0 {
        return (1.0e0, 0.0e0);
    }
    if b == 0.0e0 {
        return (0.0e0, 1.0e0);
    }
    eps = fifdmax1(eps, 1.0e-15);
    if fifdmax1(a, b) < 1.0e-3 * eps {
        return (b / (a + b), a / (a + b));
    }
    let mut ind = false;
    let mut a0 = a;
    let mut b0 = b;
    let mut x0 = x;
    let mut y0 = y;
    let mut lambda = 0.0e0;
    let mut n: i32 = 0;
    let go = if fifdmin1(a0, b0) > 1.0e0 {
        // S40
        lambda = if a > b { (a + b) * y - b } else { a - (a + b) * x };
        if lambda < 0.0e0 {
            ind = true;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
            lambda = lambda.abs();
        }
        // S70
        if b0 < 40.0e0 && b0 * x0 <= 0.7e0 {
            Go::S110
        } else if b0 < 40.0e0 {
            Go::S160
        } else if a0 > b0 {
            // S80
            if b0 <= 100.0e0 || lambda > 0.03e0 * b0 { Go::S130 } else { Go::S200 }
        } else if a0 <= 100.0e0 || lambda > 0.03e0 * a0 {
            Go::S130
        } else {
            Go::S200
        }
    } else {
        if x > 0.5e0 {
            ind = true;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
        }
        // S10
        if b0 < fifdmin1(eps, eps * a0) {
            Go::S90
        } else if a0 < fifdmin1(eps, eps * b0) && b0 * x0 <= 1.0e0 {
            Go::S100
        } else if fifdmax1(a0, b0) > 1.0e0 {
            // S20
            if b0 <= 1.0e0 {
                Go::S110
            } else if x0 >= 0.3e0 {
                Go::S120
            } else if x0 < 0.1e0 && (x0 * b0).powf(a0) <= 0.7e0 {
                Go::S110
            } else if b0 > 15.0e0 {
                // S30
                Go::S150
            } else {
                n = 20;
                Go::S140
            }
        } else if a0 >= fifdmin1(0.2e0, b0) || x0.powf(a0) <= 0.9e0 {
            Go::S110
        } else if x0 >= 0.3e0 {
            Go::S120
        } else {
            n = 20;
            Go::S140
        }
    };
    match go {
        Go::S90 => {
            w = fpser(a0, b0, x0, eps);
            w1 = 0.5e0 + (0.5e0 - w);
        }
        Go::S100 => {
            w1 = apser(a0, b0, x0, eps);
            w = 0.5e0 + (0.5e0 - w1);
        }
        Go::S110 => {
            w = bpser(a0, b0, x0, eps);
            w1 = 0.5e0 + (0.5e0 - w);
        }
        Go::S120 => {
            w1 = bpser(b0, a0, y0, eps);
            w = 0.5e0 + (0.5e0 - w1);
        }
        Go::S130 => {
            w = bfrac(a0, b0, x0, y0, lambda, 15.0e0 * eps);
            w1 = 0.5e0 + (0.5e0 - w);
        }
        Go::S140 | Go::S150 => {
            if let Go::S140 = go {
                w1 = bup(b0, a0, y0, x0, n, eps);
                b0 += n as f64;
            }
            bgrat(b0, a0, y0, x0, &mut w1, 15.0e0 * eps);
            w = 0.5e0 + (0.5e0 - w1);
        }
        Go::S160 => {
            n = b0 as i32;
            b0 -= n as f64;
            if b0 == 0.0e0 {
                n -= 1;
                b0 = 1.0e0;
            }
            w = bup(b0, a0, y0, x0, n, eps);
            if x0 <= 0.7e0 {
                w += bpser(a0, b0, x0, eps);
                w1 = 0.5e0 + (0.5e0 - w);
            } else {
                if a0 <= 15.0e0 {
                    n = 20;
                    w += bup(a0, b0, x0, y0, n, eps);
                    a0 += n as f64;
                }
                bgrat(a0, b0, x0, y0, &mut w, 15.0e0 * eps);
                w1 = 0.5e0 + (0.5e0 - w);
            }
        }
        Go::S200 => {
            w = basym(a0, b0, lambda, 100.0e0 * eps);
            w1 = 0.5e0 + (0.5e0 - w);
        }
    }
    if ind { (w1, w) } else { (w, w1) }
}

/// `brcmp1(mu, a, b, x, y)`: exp(mu) * x^a * y^b / B(a, b).
fn brcmp1(mu: i32, a: f64, b: f64, x: f64, y: f64) -> f64 {
    const CONST: f64 = 0.398942280401433e0;
    let a0 = fifdmin1(a, b);
    if a0 >= 8.0e0 {
        return brc_large(a, b, x, y, |t| esum(mu, t), CONST);
    }
    let (lnx, lny) = ln_xy(x, y);
    let mut z = a * lnx + b * lny;
    if a0 >= 1.0e0 {
        z -= betaln(a, b);
        return esum(mu, z);
    }
    let mut b0 = fifdmax1(a, b);
    if b0 >= 8.0e0 {
        let u = gamln1(a0) + algdiv(a0, b0);
        let t3 = z - u;
        return a0 * esum(mu, t3);
    }
    if b0 <= 1.0e0 {
        let mut brcmp1 = esum(mu, z);
        if brcmp1 == 0.0e0 {
            return brcmp1;
        }
        let apb = a + b;
        let z = if apb > 1.0e0 {
            let u = a + b - 1.0e0;
            (1.0e0 + gam1(u)) / apb
        } else {
            1.0e0 + gam1(apb)
        };
        let c = (1.0e0 + gam1(a)) * (1.0e0 + gam1(b)) / z;
        brcmp1 = brcmp1 * (a0 * c) / (1.0e0 + a0 / b0);
        return brcmp1;
    }
    // S70
    let mut u = gamln1(a0);
    let n = (b0 - 1.0e0) as i32;
    if n >= 1 {
        let mut c = 1.0e0;
        for _ in 1..=n {
            b0 -= 1.0e0;
            c *= b0 / (a0 + b0);
        }
        u = c.ln() + u;
    }
    z -= u;
    b0 -= 1.0e0;
    let apb = a0 + b0;
    let t = if apb > 1.0e0 {
        let u = a0 + b0 - 1.0e0;
        (1.0e0 + gam1(u)) / apb
    } else {
        1.0e0 + gam1(apb)
    };
    a0 * esum(mu, z) * (1.0e0 + gam1(b0)) / t
}

/// `brcomp(a, b, x, y)`: x^a * y^b / B(a, b).
fn brcomp(a: f64, b: f64, x: f64, y: f64) -> f64 {
    const CONST: f64 = 0.398942280401433e0;
    if x == 0.0e0 || y == 0.0e0 {
        return 0.0e0;
    }
    let a0 = fifdmin1(a, b);
    if a0 >= 8.0e0 {
        return brc_large(a, b, x, y, f64::exp, CONST);
    }
    let (lnx, lny) = ln_xy(x, y);
    let mut z = a * lnx + b * lny;
    if a0 >= 1.0e0 {
        z -= betaln(a, b);
        return z.exp();
    }
    let mut b0 = fifdmax1(a, b);
    if b0 >= 8.0e0 {
        let u = gamln1(a0) + algdiv(a0, b0);
        return a0 * (z - u).exp();
    }
    if b0 <= 1.0e0 {
        let mut brcomp = z.exp();
        if brcomp == 0.0e0 {
            return brcomp;
        }
        let apb = a + b;
        let z = if apb > 1.0e0 {
            let u = a + b - 1.0e0;
            (1.0e0 + gam1(u)) / apb
        } else {
            1.0e0 + gam1(apb)
        };
        let c = (1.0e0 + gam1(a)) * (1.0e0 + gam1(b)) / z;
        brcomp = brcomp * (a0 * c) / (1.0e0 + a0 / b0);
        return brcomp;
    }
    let mut u = gamln1(a0);
    let n = (b0 - 1.0e0) as i32;
    if n >= 1 {
        let mut c = 1.0e0;
        for _ in 1..=n {
            b0 -= 1.0e0;
            c *= b0 / (a0 + b0);
        }
        u = c.ln() + u;
    }
    z -= u;
    b0 -= 1.0e0;
    let apb = a0 + b0;
    let t = if apb > 1.0e0 {
        let u = a0 + b0 - 1.0e0;
        (1.0e0 + gam1(u)) / apb
    } else {
        1.0e0 + gam1(apb)
    };
    a0 * z.exp() * (1.0e0 + gam1(b0)) / t
}

/// Shared start of `brcomp`/`brcmp1`: (ln x, ln y) choosing `alnrel` for small x or y.
fn ln_xy(x: f64, y: f64) -> (f64, f64) {
    if x <= 0.375e0 {
        (x.ln(), alnrel(-x))
    } else if y <= 0.375e0 {
        (alnrel(-y), y.ln())
    } else {
        (x.ln(), y.ln())
    }
}

/// Shared tail (label S130) of `brcomp`/`brcmp1` for min(a, b) >= 8.
/// `e` is `exp` for brcomp and `esum(mu, .)` for brcmp1.
fn brc_large(a: f64, b: f64, x: f64, y: f64, e: impl Fn(f64) -> f64, cnst: f64) -> f64 {
    let (x0, lambda);
    if a > b {
        let h = b / a;
        x0 = 1.0e0 / (1.0e0 + h);
        lambda = (a + b) * y - b;
    } else {
        let h = a / b;
        x0 = h / (1.0e0 + h);
        lambda = a - (a + b) * x;
    }
    let y0 = if a > b {
        let h = b / a;
        h / (1.0e0 + h)
    } else {
        let h = a / b;
        1.0e0 / (1.0e0 + h)
    };
    let mut ee = -(lambda / a);
    let u = if ee.abs() > 0.6e0 { ee - (x / x0).ln() } else { rlog1(ee) };
    ee = lambda / b;
    let v = if ee.abs() > 0.6e0 { ee - (y / y0).ln() } else { rlog1(ee) };
    let z = e(-(a * u + b * v));
    cnst * (b * x0).sqrt() * z * (-bcorr(a, b)).exp()
}

/// `bup`: I_x(a, b) - I_x(a+n, b) for positive integer n.
fn bup(a: f64, b: f64, x: f64, y: f64, n: i32, eps: f64) -> f64 {
    let apb = a + b;
    let ap1 = a + 1.0e0;
    let mut mu: i32 = 0;
    let mut d = 1.0e0;
    if !(n == 1 || a < 1.0e0) && !(apb < 1.1e0 * ap1) {
        mu = exparg(1).abs() as i32;
        let k = exparg(0) as i32;
        if k < mu {
            mu = k;
        }
        let t = mu as f64;
        d = (-t).exp();
    }
    let mut bup = brcmp1(mu, a, b, x, y) / a;
    if n == 1 || bup == 0.0e0 {
        return bup;
    }
    let nm1 = n - 1;
    let mut w = d;
    let mut k: i32 = 0;
    let mut do_s30 = false;
    if b > 1.0e0 {
        if y <= 1.0e-4 {
            k = nm1;
            do_s30 = true;
        } else {
            let r = (b - 1.0e0) * x / y - a;
            if r >= 1.0e0 {
                let t = nm1 as f64;
                k = nm1;
                if r < t {
                    k = r as i32;
                }
                do_s30 = true;
            }
        }
    }
    if do_s30 {
        for i in 1..=k {
            let l = (i - 1) as f64;
            d = (apb + l) / (ap1 + l) * x * d;
            w += d;
        }
        if k == nm1 {
            bup *= w;
            return bup;
        }
    }
    // S50
    for i in (k + 1)..=nm1 {
        let l = (i - 1) as f64;
        d = (apb + l) / (ap1 + l) * x * d;
        w += d;
        if d <= eps * w {
            break;
        }
    }
    bup *= w;
    bup
}

/// `erf1(x)`: the error function.
fn erf1(x: f64) -> f64 {
    const C: f64 = 0.564189583547756e0;
    let ax = x.abs();
    if ax <= 0.5e0 {
        let t = x * x;
        let top = (((ERF_A[0] * t + ERF_A[1]) * t + ERF_A[2]) * t + ERF_A[3]) * t + ERF_A[4] + 1.0e0;
        let bot = ((ERF_B[0] * t + ERF_B[1]) * t + ERF_B[2]) * t + 1.0e0;
        return x * (top / bot);
    }
    if ax <= 4.0e0 {
        let (top, bot) = erf_pq(ax);
        let mut erf1 = 0.5e0 + (0.5e0 - (-(x * x)).exp() * top / bot);
        if x < 0.0e0 {
            erf1 = -erf1;
        }
        return erf1;
    }
    if ax >= 5.8e0 {
        return if x < 0.0e0 { -1.0e0 } else { 1.0e0 };
    }
    let x2 = x * x;
    let t = 1.0e0 / x2;
    let (top, bot) = erf_rs(t);
    let mut erf1 = (C - top / (x2 * bot)) / ax;
    erf1 = 0.5e0 + (0.5e0 - (-x2).exp() * erf1);
    if x < 0.0e0 {
        erf1 = -erf1;
    }
    erf1
}

const ERF_A: [f64; 5] = [
    0.771058495001320e-04, -0.133733772997339e-02, 0.323076579225834e-01, 0.479137145607681e-01, 0.128379167095513e+00,
];
const ERF_B: [f64; 3] = [0.301048631703895e-02, 0.538971687740286e-01, 0.375795757275549e+00];
const ERF_P: [f64; 8] = [
    -1.36864857382717e-07, 5.64195517478974e-01, 7.21175825088309e+00, 4.31622272220567e+01,
    1.52989285046940e+02, 3.39320816734344e+02, 4.51918953711873e+02, 3.00459261020162e+02,
];
const ERF_Q: [f64; 8] = [
    1.00000000000000e+00, 1.27827273196294e+01, 7.70001529352295e+01, 2.77585444743988e+02,
    6.38980264465631e+02, 9.31354094850610e+02, 7.90950925327898e+02, 3.00459260956983e+02,
];
const ERF_R: [f64; 5] = [
    2.10144126479064e+00, 2.62370141675169e+01, 2.13688200555087e+01, 4.65807828718470e+00, 2.82094791773523e-01,
];
const ERF_S: [f64; 4] = [9.41537750555460e+01, 1.87114811799590e+02, 9.90191814623914e+01, 1.80124575948747e+01];

fn erf_pq(ax: f64) -> (f64, f64) {
    let p = &ERF_P;
    let q = &ERF_Q;
    let top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax + p[6]) * ax + p[7];
    let bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax + q[6]) * ax + q[7];
    (top, bot)
}

fn erf_rs(t: f64) -> (f64, f64) {
    let r = &ERF_R;
    let s = &ERF_S;
    let top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
    let bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0e0;
    (top, bot)
}

/// `erfc1(ind, x)`: erfc(x) if ind == 0, exp(x²)·erfc(x) otherwise.
fn erfc1(ind: i32, x: f64) -> f64 {
    const C: f64 = 0.564189583547756e0;
    let ax = x.abs();
    if ax <= 0.5e0 {
        let t = x * x;
        let top = (((ERF_A[0] * t + ERF_A[1]) * t + ERF_A[2]) * t + ERF_A[3]) * t + ERF_A[4] + 1.0e0;
        let bot = ((ERF_B[0] * t + ERF_B[1]) * t + ERF_B[2]) * t + 1.0e0;
        let mut erfc1 = 0.5e0 + (0.5e0 - x * (top / bot));
        if ind != 0 {
            erfc1 = t.exp() * erfc1;
        }
        return erfc1;
    }
    let mut erfc1;
    if ax <= 4.0e0 {
        let (top, bot) = erf_pq(ax);
        erfc1 = top / bot;
    } else {
        // S20
        if x <= -5.6e0 {
            // S60
            erfc1 = 2.0e0;
            if ind != 0 {
                erfc1 = 2.0e0 * (x * x).exp();
            }
            return erfc1;
        }
        if ind == 0 && (x > 100.0e0 || x * x > -exparg(1)) {
            return 0.0e0; // S70
        }
        // S30
        let t = (1.0e0 / x).powf(2.0);
        let (top, bot) = erf_rs(t);
        erfc1 = (C - t * top / bot) / ax;
    }
    // S40
    if ind != 0 {
        if x < 0.0e0 {
            erfc1 = 2.0e0 * (x * x).exp() - erfc1;
        }
        return erfc1;
    }
    // S50
    let w = x * x;
    let t = w;
    let e = w - t;
    erfc1 = (0.5e0 + (0.5e0 - e)) * (-t).exp() * erfc1;
    if x < 0.0e0 {
        erfc1 = 2.0e0 - erfc1;
    }
    erfc1
}

/// `esum(mu, x)`: exp(mu + x).
fn esum(mu: i32, x: f64) -> f64 {
    if x > 0.0e0 {
        if mu <= 0 {
            let w = mu as f64 + x;
            if w >= 0.0e0 {
                return w.exp();
            }
        }
    } else if mu >= 0 {
        let w = mu as f64 + x;
        if w <= 0.0e0 {
            return w.exp();
        }
    }
    let w = mu as f64;
    w.exp() * x.exp()
}

/// `fpser`: I_x(a, b) for b < min(eps, eps*a) and x <= 0.5.
fn fpser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    let mut fpser = 1.0e0;
    if !(a <= 1.0e-3 * eps) {
        fpser = 0.0e0;
        let t = a * x.ln();
        if t < exparg(1) {
            return fpser;
        }
        fpser = t.exp();
    }
    fpser = b / a * fpser;
    let tol = eps / a;
    let mut an = a + 1.0e0;
    let mut t = x;
    let mut s = t / an;
    loop {
        an += 1.0e0;
        t = x * t;
        let c = t / an;
        s += c;
        if !(c.abs() > tol) {
            break;
        }
    }
    fpser *= 1.0e0 + a * s;
    fpser
}

/// `gam1(a)`: 1/Γ(a+1) - 1 for -0.5 <= a <= 1.5.
fn gam1(a: f64) -> f64 {
    const S1: f64 = 0.273076135303957e+00;
    const S2: f64 = 0.559398236957378e-01;
    const P: [f64; 7] = [
        0.577215664901533e+00, -0.409078193005776e+00, -0.230975380857675e+00, 0.597275330452234e-01,
        0.766968181649490e-02, -0.514889771323592e-02, 0.589597428611429e-03,
    ];
    const Q: [f64; 5] = [
        0.100000000000000e+01, 0.427569613095214e+00, 0.158451672430138e+00, 0.261132021441447e-01,
        0.423244297896961e-02,
    ];
    const R: [f64; 9] = [
        -0.422784335098468e+00, -0.771330383816272e+00, -0.244757765222226e+00, 0.118378989872749e+00,
        0.930357293360349e-03, -0.118290993445146e-01, 0.223047661158249e-02, 0.266505979058923e-03,
        -0.132674909766242e-03,
    ];
    let mut t = a;
    let d = a - 0.5e0;
    if d > 0.0e0 {
        t = d - 0.5e0;
    }
    if t < 0.0 {
        let top = (((((((R[8] * t + R[7]) * t + R[6]) * t + R[5]) * t + R[4]) * t + R[3]) * t + R[2]) * t + R[1]) * t
            + R[0];
        let bot = (S2 * t + S1) * t + 1.0e0;
        let w = top / bot;
        if d > 0.0e0 { t * w / a } else { a * (w + 0.5e0 + 0.5e0) }
    } else if t == 0.0 {
        0.0e0
    } else {
        let top = (((((P[6] * t + P[5]) * t + P[4]) * t + P[3]) * t + P[2]) * t + P[1]) * t + P[0];
        let bot = (((Q[4] * t + Q[3]) * t + Q[2]) * t + Q[1]) * t + 1.0e0;
        let w = top / bot;
        if d > 0.0e0 { t / a * (w - 0.5e0 - 0.5e0) } else { a * w }
    }
}

/// `gamln(a)`: ln Γ(a) for a > 0.
fn gamln(a: f64) -> f64 {
    const C0: f64 = 0.833333333333333e-01;
    const C1: f64 = -0.277777777760991e-02;
    const C2: f64 = 0.793650666825390e-03;
    const C3: f64 = -0.595202931351870e-03;
    const C4: f64 = 0.837308034031215e-03;
    const C5: f64 = -0.165322962780713e-02;
    const D: f64 = 0.418938533204673e0;
    if a <= 0.8e0 {
        return gamln1(a) - a.ln();
    }
    if a <= 2.25e0 {
        let t = a - 0.5e0 - 0.5e0;
        return gamln1(t);
    }
    if a < 10.0e0 {
        let n = (a - 1.25e0) as i32;
        let mut t = a;
        let mut w = 1.0e0;
        for _ in 1..=n {
            t -= 1.0e0;
            w = t * w;
        }
        let t1 = t - 1.0e0;
        return gamln1(t1) + w.ln();
    }
    let t = (1.0e0 / a).powf(2.0);
    let w = (((((C5 * t + C4) * t + C3) * t + C2) * t + C1) * t + C0) / a;
    D + w + (a - 0.5e0) * (a.ln() - 1.0e0)
}

/// `gamln1(a)`: ln Γ(1 + a) for -0.2 <= a <= 1.25.
fn gamln1(a: f64) -> f64 {
    const P0: f64 = 0.577215664901533e+00;
    const P1: f64 = 0.844203922187225e+00;
    const P2: f64 = -0.168860593646662e+00;
    const P3: f64 = -0.780427615533591e+00;
    const P4: f64 = -0.402055799310489e+00;
    const P5: f64 = -0.673562214325671e-01;
    const P6: f64 = -0.271935708322958e-02;
    const Q1: f64 = 0.288743195473681e+01;
    const Q2: f64 = 0.312755088914843e+01;
    const Q3: f64 = 0.156875193295039e+01;
    const Q4: f64 = 0.361951990101499e+00;
    const Q5: f64 = 0.325038868253937e-01;
    const Q6: f64 = 0.667465618796164e-03;
    const R0: f64 = 0.422784335098467e+00;
    const R1: f64 = 0.848044614534529e+00;
    const R2: f64 = 0.565221050691933e+00;
    const R3: f64 = 0.156513060486551e+00;
    const R4: f64 = 0.170502484022650e-01;
    const R5: f64 = 0.497958207639485e-03;
    const S1: f64 = 0.124313399877507e+01;
    const S2: f64 = 0.548042109832463e+00;
    const S3: f64 = 0.101552187439830e+00;
    const S4: f64 = 0.713309612391000e-02;
    const S5: f64 = 0.116165475989616e-03;
    if a < 0.6e0 {
        let w = ((((((P6 * a + P5) * a + P4) * a + P3) * a + P2) * a + P1) * a + P0)
            / ((((((Q6 * a + Q5) * a + Q4) * a + Q3) * a + Q2) * a + Q1) * a + 1.0e0);
        return -(a * w);
    }
    let x = a - 0.5e0 - 0.5e0;
    let w = (((((R5 * x + R4) * x + R3) * x + R2) * x + R1) * x + R0)
        / (((((S5 * x + S4) * x + S3) * x + S2) * x + S1) * x + 1.0e0);
    x * w
}

/// `grat1(a, x, r, eps)`: incomplete gamma ratios (P(a,x), Q(a,x)) for a <= 1;
/// `r` = exp(-x) x^a / Γ(a).
fn grat1(a: f64, x: f64, r: f64, eps: f64) -> (f64, f64) {
    if a * x == 0.0e0 {
        return if x <= a { (0.0e0, 1.0e0) } else { (1.0e0, 0.0e0) };
    }
    if a == 0.5e0 {
        if x < 0.25e0 {
            let p = erf1(x.sqrt());
            return (p, 0.5e0 + (0.5e0 - p));
        }
        let q = erfc1(0, x.sqrt());
        return (0.5e0 + (0.5e0 - q), q);
    }
    if x < 1.1e0 {
        // S10: Taylor series for P(a,x)/x^a
        let mut an = 3.0e0;
        let mut c = x;
        let mut sum = x / (a + 3.0e0);
        let tol = 0.1e0 * eps / (a + 1.0e0);
        loop {
            an += 1.0e0;
            c = -(c * (x / an));
            let t = c / (a + an);
            sum += t;
            if !(t.abs() > tol) {
                break;
            }
        }
        let j = a * x * ((sum / 6.0e0 - 0.5e0 / (a + 2.0e0)) * x + 1.0e0 / (a + 1.0e0));
        let z = a * x.ln();
        let h = gam1(a);
        let g = 1.0e0 + h;
        let s40 = if x < 0.25e0 { !(z > -0.13394e0) } else { !(a < x / 2.59e0) };
        if s40 {
            let w = z.exp();
            let p = w * g * (0.5e0 + (0.5e0 - j));
            return (p, 0.5e0 + (0.5e0 - p));
        }
        // S50
        let l = rexp(z);
        let w = 0.5e0 + (0.5e0 + l);
        let q = (w * j - l) * g - h;
        if q < 0.0e0 {
            return (1.0e0, 0.0e0);
        }
        return (0.5e0 + (0.5e0 - q), q);
    }
    // S60: continued fraction
    let mut a2nm1 = 1.0e0;
    let mut a2n = 1.0e0;
    let mut b2nm1 = x;
    let mut b2n = x + (1.0e0 - a);
    let mut c = 1.0e0;
    loop {
        a2nm1 = x * a2n + c * a2nm1;
        b2nm1 = x * b2n + c * b2nm1;
        let am0 = a2nm1 / b2nm1;
        c += 1.0e0;
        let cma = c - a;
        a2n = a2nm1 + cma * a2n;
        b2n = b2nm1 + cma * b2n;
        let an0 = a2n / b2n;
        if !((an0 - am0).abs() >= eps * an0) {
            let q = r * an0;
            return (0.5e0 + (0.5e0 - q), q);
        }
    }
}

/// `gsumln(a, b)`: ln Γ(a + b) for 1 <= a, b <= 2.
fn gsumln(a: f64, b: f64) -> f64 {
    let x = a + b - 2.0e0;
    if x <= 0.25e0 {
        let t1 = 1.0e0 + x;
        return gamln1(t1);
    }
    if x <= 1.25e0 {
        return gamln1(x) + alnrel(x);
    }
    let t2 = x - 1.0e0;
    gamln1(t2) + (x * (1.0e0 + x)).ln()
}

/// `psi(x)`: the digamma function.
fn psi(xx: f64) -> f64 {
    const DX0: f64 = 1.461632144968362341262659542325721325e0;
    const PIOV4: f64 = 0.785398163397448e0;
    const P1: [f64; 7] = [
        0.895385022981970e-02, 0.477762828042627e+01, 0.142441585084029e+03, 0.118645200713425e+04,
        0.363351846806499e+04, 0.413810161269013e+04, 0.130560269827897e+04,
    ];
    const P2: [f64; 4] = [-0.212940445131011e+01, -0.701677227766759e+01, -0.448616543918019e+01, -0.648157123766197e+00];
    const Q1: [f64; 6] = [
        0.448452573429826e+02, 0.520752771467162e+03, 0.221000799247830e+04, 0.364127349079381e+04,
        0.190831076596300e+04, 0.691091682714533e-05,
    ];
    const Q2: [f64; 4] = [0.322703493791143e+02, 0.892920700481861e+02, 0.546117738103215e+02, 0.777788548522962e+01];
    // xmax1 = min(ipmpar(3), 1/spmpar(1)) = min(2147483647, 2^52)
    let xmax1 = fifdmin1(2147483647.0, 1.0e0 / EPS_MACH);
    let xsmall = 1.0e-9;
    let mut x = xx;
    let mut aug = 0.0e0;
    if x < 0.5e0 {
        if x.abs() <= xsmall {
            if x == 0.0e0 {
                return 0.0e0;
            }
            aug = -(1.0e0 / x);
        } else {
            let mut w = -x;
            let mut sgn = PIOV4;
            if !(w > 0.0e0) {
                w = -w;
                sgn = -sgn;
            }
            if w >= xmax1 {
                return 0.0e0;
            }
            let fifidint = |v: f64| -> i64 { if v < 1.0 { 0 } else { v as i64 } };
            let mut nq = fifidint(w);
            w -= nq as f64;
            nq = fifidint(w * 4.0e0);
            w = 4.0e0 * (w - nq as f64 * 0.25e0);
            let mut n = nq / 2;
            if n + n != nq {
                w = 1.0e0 - w;
            }
            let z = PIOV4 * w;
            let mut m = n / 2;
            if m + m != n {
                sgn = -sgn;
            }
            n = (nq + 1) / 2;
            m = n / 2;
            m += m;
            if m != n {
                aug = sgn * (z.sin() / z.cos() * 4.0e0);
            } else {
                if z == 0.0e0 {
                    return 0.0e0;
                }
                aug = sgn * (z.cos() / z.sin() * 4.0e0);
            }
        }
        x = 1.0e0 - x;
    }
    if x <= 3.0e0 {
        let mut den = x;
        let mut upper = P1[0] * x;
        for i in 1..=5 {
            den = (den + Q1[i - 1]) * x;
            upper = (upper + P1[i]) * x;
        }
        den = (upper + P1[6]) / (den + Q1[5]);
        let xmx0 = x - DX0;
        return den * xmx0 + aug;
    }
    if x < xmax1 {
        let w = 1.0e0 / (x * x);
        let mut den = w;
        let mut upper = P2[0] * w;
        for i in 1..=3 {
            den = (den + Q2[i - 1]) * w;
            upper = (upper + P2[i]) * w;
        }
        aug = upper / (den + Q2[3]) - 0.5e0 / x + aug;
    }
    aug + x.ln()
}

/// `rexp(x)`: exp(x) - 1.
fn rexp(x: f64) -> f64 {
    const P1: f64 = 0.914041914819518e-09;
    const P2: f64 = 0.238082361044469e-01;
    const Q1: f64 = -0.499999999085958e+00;
    const Q2: f64 = 0.107141568980644e+00;
    const Q3: f64 = -0.119041179760821e-01;
    const Q4: f64 = 0.595130811860248e-03;
    if x.abs() <= 0.15e0 {
        return x * (((P2 * x + P1) * x + 1.0e0) / ((((Q4 * x + Q3) * x + Q2) * x + Q1) * x + 1.0e0));
    }
    let w = x.exp();
    if x > 0.0e0 { w * (0.5e0 + (0.5e0 - 1.0e0 / w)) } else { w - 0.5e0 - 0.5e0 }
}

/// `rlog1(x)`: x - ln(1 + x).
fn rlog1(x: f64) -> f64 {
    const A: f64 = 0.566749439387324e-01;
    const B: f64 = 0.456512608815524e-01;
    const P0: f64 = 0.333333333333333e+00;
    const P1: f64 = -0.224696413112536e+00;
    const P2: f64 = 0.620886815375787e-02;
    const Q1: f64 = -0.127408923933623e+01;
    const Q2: f64 = 0.354508718369557e+00;
    if x < -0.39e0 || x > 0.57e0 {
        let w = x + 0.5e0 + 0.5e0;
        return x - w.ln();
    }
    let (h, w1);
    if x < -0.18e0 {
        let h0 = x + 0.3e0;
        h = h0 / 0.7e0;
        w1 = A - h * 0.3e0;
    } else if x > 0.18e0 {
        h = 0.75e0 * x - 0.25e0;
        w1 = B + h / 3.0e0;
    } else {
        h = x;
        w1 = 0.0e0;
    }
    let r = h / (h + 2.0e0);
    let t = r * r;
    let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0e0);
    2.0e0 * t * (1.0e0 / (1.0e0 - r) - r * w) + w1
}

/// `cumbet(x, y, a, b)`: (I_x(a, b), 1 - I_x(a, b)) with y = 1 - x.
pub fn cumbet(x: f64, y: f64, a: f64, b: f64) -> (f64, f64) {
    if x <= 0.0e0 {
        return (0.0e0, 1.0e0);
    }
    if y <= 0.0e0 {
        return (1.0e0, 0.0e0);
    }
    bratio(a, b, x, y)
}

/// Perl `Math::CDF::pnbinom(s, xn, pr)`: P(S <= s) for the negative binomial
/// with `xn` successes and success probability `pr` (`s` may be fractional),
/// i.e. DCDFLIB `cdfnbn` with `which = 1`, `ompr = 1 - pr`.
/// `None` where Math::CDF returns `undef` (invalid arguments).
pub fn pnbinom(s: f64, xn: f64, pr: f64) -> Option<f64> {
    let ompr = 1.0 - pr;
    if s < 0.0 || xn < 0.0 || !(0.0..=1.0).contains(&pr) || !(0.0..=1.0).contains(&ompr) {
        return None;
    }
    let prompr = pr + ompr;
    if (prompr - 0.5e0 - 0.5e0).abs() > 3.0e0 * EPS_MACH {
        return None;
    }
    // cumnbn: cumbet(pr, ompr, xn, s + 1)
    let t1 = s + 1.0e0;
    Some(cumbet(pr, ompr, xn, t1).0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Inputs and outputs captured from the original DCDFLIB C code
    /// (gcc -O2 and -O0 give identical bits). `'b'`: cumbet(x, 1-x, a, b)
    /// with (a, b, x); `'n'`: pnbinom(s, r, p) with (s, r, p).
    /// The full comparison (230,000 cases, 0 mismatches) is described in
    /// docs/porting-notes.md.
    #[rustfmt::skip]
    const CASES: &[(char, f64, f64, f64, u64, u64)] = &[
        ('b', 0.023155711504954515, 0.16291451531358916, 5.28281099924894e-12, 0x3fdee1e4dad4c8e1, 0x3fe08f0d92959b90),
        ('b', 714809.7775381657, 2.0666466289868834, 0.9999998555463996, 0x3fefdf17a1f68fce, 0x3f70742f04b81940),
        ('b', 381.9091738526451, 67.43901106548445, 1.2656254251051295e-10, 0x0000000000000000, 0x3ff0000000000000),
        ('b', 1.7577216846839924e-05, 1.3528942841491634, 0.5392060406928092, 0x3feffff2a117e3f4, 0x3edabdd038193515),
        ('b', 0.001151389325851597, 0.7849045088370943, 0.029433143160312758, 0x3fefdae7695a51d7, 0x3f728c4b52d71483),
        ('b', 39556617.51833366, 8.384011987955448, 0.9850104010744446, 0x0000000000000000, 0x3ff0000000000000),
        ('b', 106794.21925287542, 311.5, 0.9999999999994632, 0x3ff0000000000000, 0x0000000000000000),
        ('b', 112.66169310936297, 0.0629470214853934, 1.6952269373052144e-12, 0x0000000000000000, 0x3ff0000000000000),
        ('b', 41.780228348729224, 17.0, 0.2758479752910258, 0x3d82a6ca0ba2ae63, 0x3fefffffffffb565),
        ('b', 2989.681015770249, 154.38982854310618, 0.6087329700161047, 0x0000000000000000, 0x3ff0000000000000),
        ('b', 824931493.6561302, 6844.955943918139, 0.8700979066197199, 0x0000000000000000, 0x3ff0000000000000),
        ('b', 1.398915994598982, 0.003106786762511497, 5.246593490825822e-06, 0x3dd91d0a7faa9815, 0x3feffffffff3717b),
        ('n', 0.0, 0.789268239068116, 0.02638645283952469, 0x3fad0fb6278f44c7, 0),
        ('n', 2392.3245884827015, 1410.2811754638271, 0.36701920085336853, 0x3fd442da31f92371, 0),
        ('n', 1.5, 249.35707714297763, 0.991452749757994, 0x3fe0410210c01ab2, 0),
        ('n', 15.5, 4492.252391752839, 0.9985253671243911, 0x3feff8b4e08d7a04, 0),
        ('n', 2477.5718333327823, 35194831.291730404, 0.9999281346181244, 0x3fc3a15c2868dbb7, 0),
        ('n', 208.2638113212223, 251.0864472297689, 0.5630466109141767, 0x3fe8d3057041ddec, 0),
        ('n', 732.0306639146253, 40105.97712936726, 0.9803982591260759, 0x3f7cd69cadd21c60, 0),
        ('n', 1449.15106383, 54855247.40079099, 0.9999730073380241, 0x3fcae40a7f25af80, 0),
        ('n', 9.12575757004583, 60804306.75040437, 0.9999997073309925, 0x3f933b43fe7c3d43, 0),
        ('n', 3169.803076372741, 6497.100330269635, 0.6808996921638816, 0x3fef00f24394797e, 0),
        ('n', 1042.0, 29.82739074997099, 0.033269352292443204, 0x3feb8c3efd4c053c, 0),
        ('n', 631.202577880587, 59251712.478668824, 0.9999901726024781, 0x3fef50f430b31be5, 0),
        ('n', 58.84893617, 39922.71452799382, 0.998591879458331, 0x3fe544e31e0c427a, 0),
        ('n', 5.577514568598561, 39825736.4818869, 0.9999999273623099, 0x3fee9f1881f0b15b, 0),
        ('n', 1666.3274058776476, 10623.268066528908, 0.8672465956244356, 0x3fea71b3982932fc, 0),
        ('n', 84.30048443643108, 235.81420980167266, 0.7267851225657609, 0x3fd7d8405a6f9909, 0),
        ('n', 74.5, 19.25830932382194, 0.392997393317921, 0x3fefffbcb175cf98, 0),
        ('n', 24.84893617, 12.992151802965841, 0.1665257416219696, 0x3f79f1accfd81b92, 0),
    ];

    #[test]
    fn bit_identical_to_dcdflib_c() {
        for &(kind, u, v, w, cum_bits, ccum_bits) in CASES {
            if kind == 'b' {
                let (cum, ccum) = cumbet(w, 1.0 - w, u, v);
                assert_eq!(cum.to_bits(), cum_bits, "cumbet a={u} b={v} x={w}: {cum:e}");
                assert_eq!(ccum.to_bits(), ccum_bits, "cumbet a={u} b={v} x={w}: {ccum:e}");
            } else {
                let p = pnbinom(u, v, w).unwrap();
                assert_eq!(p.to_bits(), cum_bits, "pnbinom s={u} r={v} p={w}: {p:e}");
            }
        }
    }

    #[test]
    fn invalid_arguments_are_undef() {
        assert_eq!(pnbinom(-1.0, 2.0, 0.5), None);
        assert_eq!(pnbinom(1.0, -2.0, 0.5), None);
        assert_eq!(pnbinom(1.0, 2.0, 1.5), None);
        assert!(pnbinom(0.0, 2.0, 0.5).is_some());
    }
}
