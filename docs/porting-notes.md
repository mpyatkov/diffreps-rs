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

## M5: sliding-window scan and region merging

**Perl:** the per-chromosome loop in `diffReps.pl` (~394-480), `isAboveCutoff`,
`isWithinGap`, `init`/`expand` (`ChromaModSite.pm`) and the `WindowCoverage`
class at the end of `diffReps.pl`. Scanning uses bins of size **step**; window
`w` covers `[w*step+1, w*step+window]`.

- A window is tested if **all** treatment replicates, or **all** control
  replicates, have a nonzero raw count `> mean + nsd*dev`.
- A significant window (`p < --pval`) extends the open region if it starts
  within `--gap` bp of the region's end **and** has the same direction as the
  region's first window; otherwise the open region is closed and a new one starts.
- A region made of one window keeps the window's statistics; a longer one
  is re-tested on its summed counts (its Event can differ from its windows').
- The best window is the one with the lowest p (the first one on ties).
- Window coverage counts the bp of passing windows once even though windows
  overlap. It gives the BH `N = floor(coverage / window)`.

**Rust:** `src/scan.rs`: one function per chromosome, run in parallel with
rayon. The open region is an `Option<Open>`; Perl's "pval undefined → needs
re-test" flag becomes `stat: Option<...>`.

**Checked:** `tests/parity_scan.rs` scans all four runs (background from the
report header; norm from `norm.txt` or recomputed). Window coverage, the number
of regions and every region's coordinates, normalized counts, Event, log2FC,
p-value, best window and window p-value match diffReps exactly:

| Run | Coverage | Regions |
|---|---|---|
| RIPPM_1000 | 222,074,700 | 6,622 |
| RIPPM_200 | 175,240,480 | 19,714 |
| DIFFREPS_1000 | 221,172,000 | 4,543 |
| DIFFREPS_200 | 175,240,480 | 13,922 |

## M6: Benjamini-Hochberg and the report table

**Perl:** `adjPval` (`ChromaModSite.pm:671`) calls `padjBH` (`Common.pm:47`)
twice, for region p-values and best-window p-values, with
`N = floor(window coverage / window)`. `print_reg` (`ChromaModSite.pm:114`)
writes the 19 columns: counts, means and fold changes with `%.2f`, p-values
with Perl's default number format (`%.15g`), `NA` for the enrichment columns
(no background samples), and `winEnd = winSta + window - 1`.

`padjBH` is the usual step-up BH except for one quirk, kept for parity: if
the smallest p-value also has the smallest `p*N/rank`, it is overwritten by
the next copy-down (standard BH would keep it). See `src/padj.rs`.

**Rust:** `src/padj.rs`, `report::fmt_g15` (Rust has no `%g`), `report::write_region`,
and the full pipeline in `src/main.rs`.

**Checked:** `tests/parity_report.rs` (formerly `parity_scan.rs`) formats
every row and compares it with diffReps column by column. All non-p-value
columns are identical in all four runs. Of the four p-value columns
(pval, padj, winP, winQ), 99.7% of fields are identical as printed and the
rest differ only in the last digits (worst relative difference 4.4e-14):

| Run | p-value fields identical | Worst rel. difference |
|---|---|---|
| RIPPM_1000 | 26,430 / 26,488 | 8.7e-15 |
| RIPPM_200 | 78,783 / 78,856 | 1.0e-14 |
| DIFFREPS_1000 | 18,113 / 18,172 | 4.4e-14 |
| DIFFREPS_200 | 55,563 / 55,688 | 1.1e-14 |

**Where the last-digit differences come from (not resolved):** not the
DCDFLIB port (bit-identical to the C, and Math::CDF's copy differs from
netlib's only by explicit `(long)` casts), not Perl's literal or `norm.txt`
parsing (checked), and the same with diffReps' own computed norm. The likely
cause is the system libm on the machine that ran diffReps: DCDFLIB calls
`exp`/`log`/`pow`, and glibc picks different implementations (e.g. FMA
variants) by version and CPU, which can round differently in rare cases. The
walk then amplifies one ulp into a few. It changes nothing in practice.

## M7: performance

Parallelism is in place since M2/M5: BED files are read in parallel, and
chromosomes are scanned in parallel (rayon, `--nproc` threads). Memory stays
low because reads are kept as one `u32` each and bins are built one
chromosome at a time (prefix sums, no window-count vectors).

Same data (5 vs 3 ATAC samples, 145 M fragments, mm10), default settings
(diffReps computes the background and norm factors itself):

| | diffReps 1.55.6 (16 cores) | diffreps-rs (4 threads, laptop) |
|---|---|---|
| window 1000: read / statistics / scan | 79 s / 192 s / 219 s | 7 s / <1 s / ~35 s |
| window 200: read / statistics / scan | 93 s / 309 s / 421 s | 6 s / ~3 s / 56 s |
| window 200: total, peak memory | ~14 min | 65 s, 1.8 GB |

The scan is dominated by the NB test: about 8.8 million windows are tested
with the 200 bp window, and the exact test walks over all splits of the total
count, evaluating two incomplete beta functions per step. Making this faster
would mean changing the algorithm (e.g. the step rule), which would change
the p-values, so it is left as is.

## M8: comparison tool, README

`tools/compare_reports.py` (Python standard library only) compares two
reports: header values (timings, paths and versions ignored), regions present
in only one report, and per-column differences for shared regions, with a
relative tolerance for the p-value columns (`--ptol`, default 1e-13).

On the full default 200 bp run (our own seeded background sample, no
`--bg-from-report`), it reports identical regions and identical columns; only
the sampled mean/deviation in the header differ in the second decimal.
