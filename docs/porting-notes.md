# Porting notes: diffReps.pl → diffreps-rs

One section per milestone. Each one says what the Perl code does (with file:line
references into `../diffreps`), how the Rust version is shaped, and where the two
intentionally differ.

## M1: command line, chromosome sizes, norm file

**Perl:** `bin/diffReps.pl:36-140` parses options with `Getopt::Long`,
`choose_winstep` (`diffReps.pl:540`) fills in the window/step, and `--nsd`
`broad`/`sharp` become 2/20. `read_chrlen_ordered` (`MyShortRead.pm:234`) reads
chromosome sizes and sorts them with plain string comparison (`order_chr`,
`MyShortRead.pm:72`), so the order is `chr1, chr10, ..., chr19, chr2, ...`.
`read_norm2` (`Common.pm:495`) reads the `treatment`/`control` lines.

**Rust:** `src/cli.rs` (clap derive; `infer_long_args` reproduces Getopt::Long's
prefix matching, so `--frag` still means `--fragsize`), `src/inputs.rs`.
`Args` is what the user typed; `Config` is the validated result.

**Differences:**
- `--meth` accepts only `nb`. No annotation, hotspots or background samples.
- `--nsd` accepts any non-negative number (Perl: integers only, silently falls back to 2).
- New: `--seed`, `--sample-size` (reproducible background estimate) and
  `--bg-from-report` (developer option for comparing with a diffReps report).
- Duplicate chromosomes in `--chrlen` are an error (Perl would process them twice).

## M2: reading BED files and binning

**Perl:** `sep_chrom_bed` (`MyShortRead.pm:259`) writes one temporary file per
chromosome per sample. `bin_genome_count` (`MyShortRead.pm:142`) reduces each
read to one point, `start + frag/2` on the `+` strand and `end - frag/2`
otherwise. **A `.` strand counts as minus**, so with `--frag 0` an ATAC fragment
is counted at its end coordinate. Points go into `ceil(len/bin)` bins; points
past the last bin are ignored. Window counts are running sums of
`window/step` consecutive bins (`create_window_count`, `ChromBed.pm:123`), and
region counts sum the bins inside the region (`get_region_count`, `SRBed.pm:240`).

**Rust:** `src/bed.rs`. Each file is read once (one rayon task per file) and
each read is stored as `2*loc` in a `u32` (doubling keeps odd fragment sizes
exact). `BinCounts` builds bins of any size for one chromosome on demand and
keeps **prefix sums**, so a window or region count is one subtraction. This
replaces both the temporary files and Perl's separate window-count vectors.

**Checked:** `tests/parity_bins.rs` recomputes the raw counts of all 26,338
regions in `RIPPM_1000` and `RIPPM_200` (`round(printed * norm)`) and gets
exact agreement for every sample. Reading all eight files takes ~7.5 s and
600 MB (diffReps: 80 s for its split step).

## M3: background level and normalization factors

**Perl:** `mean_std_norm` (`diffReps.pl:594`). All of it uses bins of size
**window**, not step. Robust mode (the default) keeps each bin with probability
`100000 / total_bins` using unseeded `rand()`, then takes a right-trimmed mean
(`mean_r`, `alpha` = 0.05) and the MAD around it (×1.4826). `--std` uses the
plain mean and population SD over all bins. Normalization factors are a
median of ratios: over bins where every sample has `count >= mean + nsd*dev`
and all counts > 0, take `count / geomean(bin)`, and use the median per sample.
A factor is a **divisor** (`normalized = raw / factor`).

**Rust:** `src/background.rs`, `src/stats.rs` (`mean_r`, `mad`, `median`,
`geomean`, `var` with the same definitions). The sampling uses a seeded
SplitMix64, so runs are reproducible; `--sample-size 0` uses every bin.
`--bg-from-report` reads the mean and deviation from a diffReps report header,
which removes the randomness for parity tests. With `--norm`, the factors are
not estimated at all (diffReps estimates them and then ignores them).

**Note:** the cutoff here is `count < cutoff → skip` (so `>=` passes), while
the window scan uses `count <= cutoff → fail` (so `>` passes). Both are kept.

**Checked:** `tests/parity_background.rs`. With diffReps' header mean and
deviation, our factors print identically for `DIFFREPS_1000` (4.27 0.29 0.52
2.32 1.73 / 1.18 0.50 0.69) and `DIFFREPS_200`. Our sampled estimates are
within diffReps' own run-to-run noise.

## M4: the exact negative binomial test

**Perl:** `NBTest.pm` (the exact test of DESeq 2010, adapted) and `calc_stat`
(`ChromaModSite.pm:518`). Group sums of raw counts get NB means and variances
from `raw_sum_mean`/`raw_sum_var` (pooled normalized mean; replicate variance
minus the Poisson part, floored at `1e-8*mean`). The p-value sums
`P(a)P(kS-a)` over all splits at most as likely as the observed one, walking
in from both ends toward `kexp`. Two details make this hard to reproduce:

1. `kexp = kS*m_A/(m_A+m_B)` is **not rounded**, so `pnbinom` is called at
   fractional counts. Math::CDF passes them straight to DCDFLIB, which
   computes `I_p(r, x+1)`, a smooth function of `x`.
2. Every probability is `pnbinom(x) - pnbinom(x-1) + 1e-8`, and the walk's
   step grows when neighbouring values differ by less than a relative 1e-8.
   In flat tails that decision depends on the **rounding noise** of the
   difference of two CDF values near 1.

**What we tried first:** an exact pmf (Loader's saddle-point method, as in R's
`dbinom_raw`). It is more accurate, but only ~85% of diffReps' p-values were
reproduced, some off by 1000×. An R port using `pbeta` (TOMS 708, like DCDFLIB)
reproduced 99.7%.

**Rust:** `src/cdflib.rs` is a line-by-line port of DCDFLIB's `cumbet → bratio`
(TOMS 708) and `cdfnbn` (which=1). It is **bit-identical** to the C on
230,000 random inputs (see `tools/dcdflib-check/`). `src/nbtest.rs` ports
`nb_r_p`, `dnbinom`, `nb_join_pval`, `nb_pval` and `calc_stat` literally,
with a 4-entry cache of CDF values (the walk asks for `F(x)` and `F(x-1)` of
neighbouring points).

**Checked:** `tests/parity_nbtest.rs` recomputes the p-value, log2FC and Event
of all 26,336 regions in `RIPPM_1000` and `RIPPM_200` from their counts:
**100% match** (max relative difference 1e-14, the 15 digits diffReps prints).
`tests/data/nb_reference.R` (R port, TOMS 708 via `pbeta`) agrees within 1.4e-11.
