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
