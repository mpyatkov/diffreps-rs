# Checking `src/cdflib.rs` against the original DCDFLIB C code

`src/cdflib.rs` is a line-by-line port of DCDFLIB's incomplete beta
(`cumbet` → `bratio`), the code inside Perl's Math::CDF. It must give the
**same bits**, because diffReps' p-values depend on the rounding of
`pnbinom(x) - pnbinom(x-1)`.

Result on 2026-11-20: 230,000 random cases, **0 mismatches** (C built with
gcc -O2 and -O0, identical). A sample of 30 cases is a unit test in
`src/cdflib.rs` (`bit_identical_to_dcdflib_c`).

## Re-running the full comparison by hand

1. Download DCDFLIB (public domain) from netlib: `random/dcdflib.c.tar.gz`,
   and unpack it somewhere outside this repository.
2. Build the C harness:
   `gcc -O2 -I<dcdflib>/src -o harness harness.c <dcdflib>/src/dcdflib.c <dcdflib>/src/ipmpar.c -lm`
3. Make inputs: `python3 gen_inputs.py > inputs.txt`
4. `./harness < inputs.txt > c.txt`
5. Run the Rust side on the same inputs: a tiny crate whose `main.rs` includes
   the module with `#[path = ".../src/cdflib.rs"] mod cdflib;`, reads the same
   lines and prints the results with C's `%a` hex float format (see `rs_main.rs`).
6. `cmp c.txt rs.txt`

Input lines: `b a b x` → `cumbet(x, 1-x, a, b)` (prints cum and ccum);
`n s r p` → `pnbinom(s, r, p)` (prints the CDF or `undef`).
