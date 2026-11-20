//! Reading BED files and counting reads in bins.
//!
//! # What diffReps does
//! `sep_chrom_bed` (MyShortRead.pm:259) splits every BED file into one
//! temporary file per chromosome. `bin_genome_count` (MyShortRead.pm:142)
//! then reduces each read to a single point and counts points per bin:
//!
//! ```text
//! loc = start + frag/2   if strand eq '+'
//! loc = end   - frag/2   otherwise            <- '.' counts as minus strand!
//! bin = floor(loc / bin_size)
//! ```
//!
//! The bin vector has `ceil(chrom_len / bin_size)` entries. Reads that fall
//! past the last bin grow Perl's array but are never read, so they are
//! effectively ignored. We drop them as well.
//!
//! # How we do it
//! Each file is read once. For every read we keep only `2 * loc` as a `u32`
//! (doubling keeps odd fragment sizes exact: `frag/2` may be `x.5`). Bins of
//! any size can then be computed on demand, per chromosome, without temp files.
//! Memory: 4 bytes per read instead of a bin vector per sample per genome.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use crate::inputs::Chrom;

/// All reads of one BED file, reduced to doubled positions, grouped by chromosome.
#[derive(Debug)]
pub struct Sample {
    pub path: PathBuf,
    /// `pos2[c]` holds `2 * loc` for each read on chromosome `c` (index into the chrom list).
    pub pos2: Vec<Vec<u32>>,
    /// Reads on known chromosomes (Perl's `$nread`).
    pub nread: u64,
    /// Lines skipped because the chromosome is not in `--chrlen`.
    pub unknown_chrom: u64,
    /// Lines skipped because they have fewer than 6 fields.
    pub malformed: u64,
    /// Reads whose shifted position is negative (cannot happen with `--frag 0`).
    pub negative: u64,
}

/// Read one BED file.
pub fn read_sample(path: &Path, chroms: &[Chrom], fragsize: u64) -> Result<Sample> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let reader = BufReader::with_capacity(1 << 20, file);
    let mut s = parse_bed(reader, chroms, fragsize)
        .with_context(|| format!("while reading {}", path.display()))?;
    s.path = path.to_path_buf();
    Ok(s)
}

/// Parse BED6 text. Separate from [`read_sample`] so tests can feed strings.
pub fn parse_bed<R: BufRead>(mut reader: R, chroms: &[Chrom], fragsize: u64) -> Result<Sample> {
    // Doubled positions must fit in u32.
    if let Some(c) = chroms.iter().find(|c| c.len >= (u32::MAX as u64) / 2 - fragsize) {
        bail!("chromosome {} is too long ({} bp) for 32-bit positions", c.name, c.len);
    }
    let index: HashMap<&[u8], usize> =
        chroms.iter().enumerate().map(|(i, c)| (c.name.as_bytes(), i)).collect();
    let mut s = Sample {
        path: PathBuf::new(),
        pos2: vec![Vec::new(); chroms.len()],
        nread: 0,
        unknown_chrom: 0,
        malformed: 0,
        negative: 0,
    };
    let frag = fragsize as i64;
    let mut line = Vec::with_capacity(256);
    let mut lineno = 0u64;
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        lineno += 1;
        let text = trim_eol(&line);
        if text.is_empty() {
            continue;
        }
        // Fields are separated by tabs or spaces (Perl's `split` on whitespace).
        let mut fields = text.split(|b| *b == b'\t' || *b == b' ').filter(|f| !f.is_empty());
        let (Some(chrom), Some(start), Some(end), Some(_name), Some(_score), Some(strand)) = (
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
        ) else {
            s.malformed += 1;
            continue;
        };
        let Some(&ci) = index.get(chrom) else {
            s.unknown_chrom += 1;
            continue;
        };
        s.nread += 1;
        let start = parse_int(start).with_context(|| format!("bad start at line {lineno}"))?;
        let end = parse_int(end).with_context(|| format!("bad end at line {lineno}"))?;
        // 2*loc: '+' -> 2*start + frag, anything else -> 2*end - frag.
        let loc2 = if strand == b"+" { 2 * start + frag } else { 2 * end - frag };
        if loc2 < 0 {
            s.negative += 1;
            continue;
        }
        s.pos2[ci].push(loc2 as u32);
    }
    Ok(s)
}

fn trim_eol(line: &[u8]) -> &[u8] {
    let mut l = line;
    while let Some((&last, rest)) = l.split_last() {
        if last == b'\n' || last == b'\r' {
            l = rest;
        } else {
            break;
        }
    }
    l
}

fn parse_int(b: &[u8]) -> Result<i64> {
    if b.is_empty() || b.len() > 18 || !b.iter().all(u8::is_ascii_digit) {
        bail!("not a non-negative integer: {:?}", String::from_utf8_lossy(b));
    }
    Ok(b.iter().fold(0i64, |acc, d| acc * 10 + (d - b'0') as i64))
}

/// Read counts of one sample on one chromosome, stored as prefix sums so that
/// any run of consecutive bins can be summed in O(1).
///
/// `prefix[i]` = number of reads in bins `0..i`; there are `n_bins + 1` entries.
#[derive(Debug, Clone)]
pub struct BinCounts {
    pub bin_size: u64,
    prefix: Vec<u32>,
}

impl BinCounts {
    /// Count reads of `pos2` in bins of `bin_size`, for a chromosome of
    /// `chrom_len` bp (`ceil(chrom_len / bin_size)` bins, like Perl).
    pub fn new(pos2: &[u32], bin_size: u64, chrom_len: u64) -> Self {
        let n_bins = chrom_len.div_ceil(bin_size) as usize;
        let mut counts = vec![0u32; n_bins];
        let width2 = 2 * bin_size;
        for &p in pos2 {
            let b = (p as u64 / width2) as usize;
            if b < n_bins {
                counts[b] += 1;
            }
        }
        let mut prefix = Vec::with_capacity(n_bins + 1);
        let mut acc = 0u32;
        prefix.push(0);
        for c in counts {
            acc += c;
            prefix.push(acc);
        }
        BinCounts { bin_size, prefix }
    }

    pub fn n_bins(&self) -> usize {
        self.prefix.len() - 1
    }

    /// Reads in bin `i`.
    #[inline]
    pub fn bin(&self, i: usize) -> u32 {
        self.prefix[i + 1] - self.prefix[i]
    }

    /// Reads in bins `first..=last`.
    #[inline]
    pub fn range(&self, first: usize, last: usize) -> u32 {
        self.prefix[last + 1] - self.prefix[first]
    }

    /// Number of sliding windows of `bins_per_window` bins
    /// (`create_window_count`, ChromBed.pm:123): `n_bins - bpw + 1`,
    /// or `None` if the chromosome is shorter than one window.
    pub fn n_windows(&self, bins_per_window: usize) -> Option<usize> {
        (self.n_bins() >= bins_per_window).then(|| self.n_bins() - bins_per_window + 1)
    }

    /// Reads in window `w` (bins `w..w+bpw`).
    #[inline]
    pub fn window(&self, w: usize, bins_per_window: usize) -> u32 {
        self.range(w, w + bins_per_window - 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chroms() -> Vec<Chrom> {
        vec![Chrom { name: "chr1".into(), len: 1000 }, Chrom { name: "chr2".into(), len: 450 }]
    }

    const BED: &str = "chr1\t10\t150\t.\t0\t.\n\
                       chr1\t0\t99\tr\t0\t+\n\
                       chr2\t300\t449\t.\t0\t-\n\
                       chrUn\t5\t10\t.\t0\t.\n\
                       chr1\t5\t10\n\
                       chr2 400 1000 x 0 .\n";

    #[test]
    fn positions_follow_strand_rules() {
        // frag = 0: '+' -> start, '.'/'-' -> end.
        let s = parse_bed(BED.as_bytes(), &chroms(), 0).unwrap();
        assert_eq!(s.pos2[0], [300, 0]);
        assert_eq!(s.pos2[1], [898, 2000]);
        assert_eq!((s.nread, s.unknown_chrom, s.malformed), (4, 1, 1));
        // frag = 101: '+' -> start + 50.5, else end - 50.5 (doubled: +-101).
        let s = parse_bed(BED.as_bytes(), &chroms(), 101).unwrap();
        assert_eq!(s.pos2[0], [199, 101]);
    }

    #[test]
    fn bins_windows_and_ranges() {
        let s = parse_bed(BED.as_bytes(), &chroms(), 0).unwrap();
        // chr1, 100 bp bins: loc 150 -> bin 1, loc 0 -> bin 0. 10 bins.
        let b = BinCounts::new(&s.pos2[0], 100, 1000);
        assert_eq!(b.n_bins(), 10);
        assert_eq!((b.bin(0), b.bin(1), b.bin(2)), (1, 1, 0));
        assert_eq!(b.range(0, 9), 2);
        assert_eq!(b.n_windows(10), Some(1));
        assert_eq!(b.n_windows(11), None);
        assert_eq!(b.window(0, 2), 2);
        assert_eq!(b.window(1, 2), 1);
        // chr2 has 450 bp -> 5 bins; loc 449 -> bin 4; loc 1000 is past the end and dropped.
        let b = BinCounts::new(&s.pos2[1], 100, 450);
        assert_eq!(b.n_bins(), 5);
        assert_eq!(b.bin(4), 1);
        assert_eq!(b.range(0, 4), 1);
    }
}
