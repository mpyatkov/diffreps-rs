# diffreps-rs

A Rust port of [diffReps](https://github.com/shenlab-sinai/diffreps) 1.55.6
(Li Shen, 2013): detection of differential sites in ChIP-seq/ATAC-seq data
between two groups of samples with biological replicates, using a sliding
window and an exact negative binomial test.

Only the part of diffReps needed for the negative binomial workflow is
ported. Annotation (`region_analysis.pl`), hotspots (`findHotspots.pl`),
background samples (`--btr/--bco/--bkg`) and the T-test, G-test and
chi-square test are not.

## Download

Static Linux binaries (x86-64, no dependencies, any glibc version) are
attached to each [GitHub release](https://github.com/mpyatkov/diffreps-rs/releases):

```bash
tar xzf diffreps-vX.Y.Z-x86_64-unknown-linux-musl.tar.gz
./diffreps-vX.Y.Z-x86_64-unknown-linux-musl/diffreps --help
```

## Build

```bash
cargo build --release          # binary: target/release/diffreps
```

No runtime dependencies. Needs Rust 1.88 or newer.

### Static binary for servers with an old glibc

A normal build links against the build machine's glibc and fails on older
systems (`version GLIBC_2.29 not found`). Build a fully static binary with
musl instead (rustup ships musl; no system packages needed):

```bash
rustup target add x86_64-unknown-linux-musl
cargo build --release --target x86_64-unknown-linux-musl
# binary: target/x86_64-unknown-linux-musl/release/diffreps (runs on any x86-64 Linux)
```

musl has its own `exp`/`log`/`pow`, which can round differently from glibc in
the last bit. In our comparisons with diffReps this made no difference
(identical regions; p-values agree to the same precision as with glibc).

## Usage

Same option names as `diffReps.pl`, so existing commands work with the
script name changed:

```bash
diffreps \
  --treatment T1.bed T2.bed T3.bed T4.bed T5.bed \
  --control   C1.bed C2.bed C3.bed \
  --report    result.txt \
  --chrlen    mm10.chrom.sizes \
  --window 1000 --nsd broad --frag 0 --nproc 16 \
  --norm norm.txt                 # optional; omit to estimate factors from the data
```

| Option | Default | Meaning |
|---|---|---|
| `--treatment`, `--control` | required | BED6 files, one per replicate (≥ 2 per group) |
| `--report` | required | output file |
| `--chrlen` | required | chromosome sizes (`name length` per line); replaces `--gname` |
| `--window`, `--step`, `--mode` | peak: 1000 / 100 | window size; step = window/10 |
| `--gap` | 0 | bp allowed between significant windows of one region |
| `--nsd` | broad | count cutoff: `broad` = 2, `sharp` = 20, or a number of deviations |
| `--std`, `--alpha` | robust, 0.05 | classic mean/SD instead of trimmed mean/MAD; trim fraction |
| `--pval` | 1e-4 | window p-value cutoff |
| `--norm` | estimate | `treatment f1 f2 ...` / `control f1 f2 ...`; counts are **divided** by these |
| `--nrpass` | 1 | estimate norm factors only on bins above the `--nsd` cutoff |
| `--frag`, `--fragsize` | 100 | read shift: `start + frag/2` on `+`, `end - frag/2` otherwise |
| `--nproc` | 1 | threads |
| `--seed`, `--sample-size` | fixed, 100000 | random bin sample for the robust background; `--sample-size 0` uses all bins |
| `--bg-from-report FILE` | | take mean/deviation from a diffReps report header (for comparisons) |

The report has the same header lines and 19 columns as diffReps.

**Note on strands:** like diffReps, a read whose strand is not `+` (including
`.`) is counted at its **end** coordinate (minus `frag/2`).

## How close is it to diffReps?

Checked on four reference runs (5 vs 3 ATAC samples, mm10, windows of 1000
and 200 bp, with and without `--norm`; data not included). Using diffReps'
background values (`--bg-from-report`), every region, coordinate, count,
fold change and event is identical, and 99.7% of p-values are identical to
the 15 printed digits (the rest differ by < 5e-14 relative). Without it, the
random background sample differs from diffReps' (diffReps does not seed its
sample), which in these runs still gave identical regions.

To compare a run with a diffReps report:

```bash
python3 tools/compare_reports.py diffReps_report.txt diffreps_rs_report.txt
```

Speed on the same data: a 200 bp run takes ~1 minute on 4 threads (diffReps:
~14 minutes on 16 cores).

## Differences from diffReps

- `--chrlen` is required (no built-in genomes); duplicate chromosomes are an error.
- `--nsd` accepts fractional values; diffReps silently replaced them with 2.
- The background sample is seeded (`--seed`), so results are reproducible.
- With `--norm`, normalization factors are not estimated (diffReps estimated and then ignored them).
- Warnings about skipped BED lines are summarized per file instead of one per line.

## Tests and releases

```bash
cargo test            # unit tests, no data needed
```

GitHub Actions runs clippy and the unit tests (glibc and musl) on every push
and pull request (`.github/workflows/ci.yml`). Pushing a tag `vX.Y.Z` that
matches `version` in `Cargo.toml` builds the static musl binary and publishes
it as a release (`.github/workflows/release.yml`):

```bash
git tag v0.1.0 && git push origin v0.1.0
```

## For readers

The code was written to be read next to the Perl original:
`docs/porting-notes.md` has one section per milestone (CLI, BED and bins,
background and normalization, the NB test, the scan, BH and the report,
performance) with file:line references into diffReps and the reasons for
every deviation. The most involved part is `src/cdflib.rs`, a bit-exact port
of DCDFLIB's incomplete beta function, explained in milestone 4.

## License

diffReps is distributed under CC BY-NC-SA 4.0 (non-commercial); this port
follows its algorithms closely and should be treated as a derivative work
under the same license. DCDFLIB (`src/cdflib.rs`) is public domain.
