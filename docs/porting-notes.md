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
