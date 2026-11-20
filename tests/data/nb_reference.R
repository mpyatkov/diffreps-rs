# Reference values for the negative binomial test.
#
# A line-by-line R port of MyBioinfo/NBTest.pm and calc_stat (ChromaModSite.pm)
# from diffReps. Perl calls Math::CDF::pnbinom(x, r, p), which is DCDFLIB's
# cdfnbn -> cumnbn -> cumbet(p, 1-p, r, x+1): the regularized incomplete beta
# I_p(r, x+1) for REAL x (no rounding). R's pbeta() uses the same TOMS 708
# algorithm, so pbeta(p, r, x + 1) reproduces Math::CDF::pnbinom.
#
# Usage: Rscript tests/data/nb_reference.R > tests/data/nb_reference.tsv

pnb <- function(x, r, p) pbeta(p, r, x + 1)

dnb <- function(x, r, p) {
  if (x - 1 < 0) {
    if (x == 0) return(pnb(x, r, p) + 1e-08)
    return(pnb(x, r, p) - pnb(0, r, p) + 1e-08)
  }
  pnb(x, r, p) - pnb(x - 1, r, p) + 1e-08
}

nb_r_p <- function(m, v, eps) {
  if (m == 0) return(c(eps, eps))
  if (v < m) v <- m + eps
  r <- m^2 / (v - m)
  c(r, r / (r + m))
}

nb_join_pval <- function(kl, kr, kS, sizeA, probA, sizeB, probB, probs, eps) {
  lval <- dnb(kl, sizeA, probA) * dnb(kS - kl, sizeB, probB)
  rval <- dnb(kr, sizeA, probA) * dnb(kS - kr, sizeB, probB)
  prevlval <- lval; prevrval <- rval
  total <- lval + rval
  obstotal <- 0; step <- 1
  if (lval <= probs) obstotal <- obstotal + lval
  if (rval <= probs) obstotal <- obstotal + rval
  while (kl < kr) {
    if (abs(prevrval - rval) / prevrval > .01) {
      do_left <- TRUE
    } else if (abs(prevlval - lval) / prevlval > .01) {
      do_left <- FALSE
    } else {
      do_left <- lval > rval
    }
    if (do_left) {
      prevlval <- lval
      if (kl + step > kr) step <- kr - kl
      kl <- kl + step
      lval <- dnb(kl, sizeA, probA) * dnb(kS - kl, sizeB, probB)
      if (step == 1) total <- total + lval else total <- total + min(lval, prevlval) * step
      if (lval <= probs) {
        if (step == 1) obstotal <- obstotal + lval
        else if (prevlval <= probs) obstotal <- obstotal + max(lval, prevlval) * step
        else obstotal <- obstotal + max(lval, prevlval) * step * abs((probs - lval) / (prevlval - lval))
      }
      if (abs(prevlval - lval) / prevlval < eps) step <- max(step + 1, step * 1.5)
    } else {
      prevrval <- rval
      if (kr - step < kl) step <- kr - kl
      kr <- kr - step
      rval <- dnb(kr, sizeA, probA) * dnb(kS - kr, sizeB, probB)
      if (step == 1) total <- total + rval else total <- total + min(rval, prevrval) * step
      if (rval <= probs) {
        if (step == 1) obstotal <- obstotal + rval
        else if (prevrval <= probs) obstotal <- obstotal + max(rval, prevrval) * step
        else obstotal <- obstotal + max(rval, prevrval) * step * abs(probs - rval) / (prevrval - rval)
      }
      if (abs(prevrval - rval) / prevrval < eps) step <- max(step + 1, step * 1.5)
    }
  }
  if (obstotal > total) return(c(total, total))
  c(total, obstotal)
}

nb_pval <- function(ka, kb, miu1, var1, miu2, var2, eps) {
  rp1 <- nb_r_p(miu1, var1, eps); rp2 <- nb_r_p(miu2, var2, eps)
  kexp <- (ka + kb) * miu1 / (miu1 + miu2)
  p <- dnb(ka, rp1[1], rp1[2]) * dnb(kb, rp2[1], rp2[2])
  l <- nb_join_pval(0, kexp, ka + kb, rp1[1], rp1[2], rp2[1], rp2[2], p, eps)
  kl2 <- kexp + 1
  if (kexp + 1 > ka + kb) kl2 <- ka + kb
  rr <- nb_join_pval(kl2, ka + kb, ka + kb, rp1[1], rp1[2], rp2[1], rp2[2], p, eps)
  if (rr[2] + l[2] == 0) return(p / (rr[1] + l[1]))
  (rr[2] + l[2]) / (rr[1] + l[1])
}

raw_sum_var <- function(base_var, m1, m2, norm, raw_mean, eps) {
  base_mean <- (m1 + m2) / 2
  z <- mean(1 / norm); s <- sum(norm^2)
  v <- base_var - z * base_mean
  if (v < eps * base_mean) v <- eps * base_mean
  v * s + raw_mean
}

# calc_stat for the 'nb' method. Returns c(log2FC, pval).
calc_stat <- function(raw_tr, raw_co, norm_tr, norm_co, eps = 1e-8) {
  tr <- raw_tr / norm_tr; co <- raw_co / norm_co
  tr_m <- mean(tr); co_m <- mean(co)
  if (tr_m == 0 && co_m == 0) return(c(0, 1))
  lfc <- if (co_m == 0) Inf else log2(tr_m / co_m)
  tr_raw_m <- sum(norm_tr) * (tr_m + co_m) / 2
  co_raw_m <- sum(norm_co) * (tr_m + co_m) / 2
  tr_raw_v <- raw_sum_var(var(tr), tr_m, co_m, norm_tr, tr_raw_m, eps)
  co_raw_v <- raw_sum_var(var(co), tr_m, co_m, norm_co, co_raw_m, eps)
  c(lfc, nb_pval(sum(raw_tr), sum(raw_co), tr_raw_m, tr_raw_v, co_raw_m, co_raw_v, eps))
}

set.seed(1)
norm_tr <- c(3.14, 0.19, 0.34, 1.64, 1.15)
norm_co <- c(0.75, 0.33, 0.45)
cat("raw_tr\traw_co\tnorm_tr\tnorm_co\tlog2fc\tpval\n")
emit <- function(rt, rc, nt, nc) {
  s <- calc_stat(rt, rc, nt, nc)
  cat(paste(rt, collapse = ";"), paste(rc, collapse = ";"),
      paste(nt, collapse = ";"), paste(nc, collapse = ";"),
      sprintf("%.17g", s[1]), sprintf("%.17g", s[2]), sep = "\t")
  cat("\n")
}
# Realistic windows: a base level, a fold change, NB noise of varying dispersion.
for (i in 1:400) {
  base <- exp(runif(1, log(2), log(300)))
  fc <- 2^runif(1, -3, 3)
  size <- exp(runif(1, log(0.5), log(1e4)))
  rt <- rnbinom(5, size = size, mu = base * fc * norm_tr)
  rc <- rnbinom(3, size = size, mu = base * norm_co)
  if (sum(rt) + sum(rc) == 0) next
  emit(rt, rc, norm_tr, norm_co)
}
# Under-dispersed replicates: variance floor, huge r (Poisson limit).
for (i in 1:60) {
  base <- exp(runif(1, log(1), log(200)))
  emit(round(base * 1.7 * norm_tr), round(base * norm_co), norm_tr, norm_co)
}
# Very small and very large totals, equal norms.
for (k in c(1, 2, 3, 5, 8)) emit(c(k, k + 1, k), c(0, 1, 0), c(1, 1, 1), c(1, 1, 1))
for (k in c(5000, 20000)) emit(c(k, k * 1.1, k * 0.9, k, k), c(k / 3, k / 2, k / 4), norm_tr, norm_co)
