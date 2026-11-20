//! Small text inputs: the chromosome length file and the normalization file.

use std::collections::HashSet;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result, bail};

/// One chromosome from the `--chrlen` file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chrom {
    pub name: String,
    pub len: u64,
}

/// Port of `read_chrlen_ordered` (MyShortRead.pm).
///
/// Empty lines are skipped; malformed lines are skipped with a warning, as in
/// Perl. Chromosomes are returned in diffReps' processing order, which is a
/// plain string comparison (`chr1, chr10, chr11, ..., chr19, chr2, ..., chrX`).
/// That order matters: it is the order of regions in the report, and the
/// order in which bins are visited when sampling for the robust estimate.
pub fn read_chrlen(path: &Path) -> Result<Vec<Chrom>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("cannot open chromosome length file {}", path.display()))?;
    parse_chrlen(&text, &path.display().to_string())
}

pub fn parse_chrlen(text: &str, label: &str) -> Result<Vec<Chrom>> {
    let mut chroms = Vec::new();
    let mut seen = HashSet::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split_whitespace().collect();
        if cells.len() < 2 {
            eprintln!("Insufficient chromosome information at {label}, line:{}. Skip!", i + 1);
            continue;
        }
        let (name, len) = (cells[0], cells[1]);
        // Perl accepts the length only if it matches /^[0-9]+$/.
        let digits_only = len.bytes().all(|b| b.is_ascii_digit());
        match len.parse::<u64>() {
            Ok(len) if digits_only && len > 0 => {
                if !seen.insert(name.to_string()) {
                    bail!("chromosome {name} appears twice in {label}");
                }
                chroms.push(Chrom { name: name.to_string(), len });
            }
            _ => eprintln!("Format error at {label}, line:{}. Skip!", i + 1),
        }
    }
    if chroms.is_empty() {
        bail!("no chromosomes found in {label}");
    }
    chroms.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(chroms)
}

/// Port of `read_norm2` (MyBioinfo/Common.pm).
///
/// Only the first two lines are read. Each must start with `treatment` or
/// `control` followed by one factor per replicate, in the same order as the
/// files on the command line. Both lines are required.
pub fn read_norm(path: &Path, n_tr: usize, n_co: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("error in reading {}", path.display()))?;
    parse_norm(&text, n_tr, n_co)
}

pub fn parse_norm(text: &str, n_tr: usize, n_co: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut tr: Option<Vec<f64>> = None;
    let mut co: Option<Vec<f64>> = None;
    for line in text.lines().take(2) {
        let mut cells = line.split_whitespace();
        let target = match cells.next() {
            Some("treatment") => &mut tr,
            Some("control") => &mut co,
            _ => continue,
        };
        let values = cells
            .map(|c| c.parse::<f64>().with_context(|| format!("bad normalization constant '{c}'")))
            .collect::<Result<Vec<f64>>>()?;
        *target = Some(values);
    }
    let (Some(tr), Some(co)) = (tr, co) else {
        bail!("Read normalization constants failed. Check syntax.");
    };
    if tr.len() != n_tr || co.len() != n_co {
        bail!(
            "Number of normalization constants does not equal number of samples \
             (treatment {} vs {n_tr}, control {} vs {n_co}).",
            tr.len(),
            co.len()
        );
    }
    if tr.iter().chain(&co).any(|&f| !(f > 0.0 && f.is_finite())) {
        bail!("normalization constants must be positive numbers");
    }
    Ok((tr, co))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chrlen_is_sorted_like_perl_cmp() {
        let c = parse_chrlen("chr1\t100\nchr2 200\n\nchrX\t50\nchr10\t30\nbad\n", "t").unwrap();
        let names: Vec<_> = c.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, ["chr1", "chr10", "chr2", "chrX"]);
        assert_eq!(c[2].len, 200);
    }

    #[test]
    fn chrlen_rejects_duplicates() {
        assert!(parse_chrlen("chr1 1\nchr1 2\n", "t").is_err());
    }

    #[test]
    #[allow(clippy::approx_constant)] // 3.14 is a normalization factor from the test data
    fn norm_file_like_test_data() {
        let (t, c) = parse_norm("treatment 3.14 0.19 0.34 1.64 1.15\ncontrol 0.75 0.33 0.45\n", 5, 3).unwrap();
        assert_eq!(t, [3.14, 0.19, 0.34, 1.64, 1.15]);
        assert_eq!(c, [0.75, 0.33, 0.45]);
        // Lines may come in either order.
        let (t, c) = parse_norm("control 1 2\ntreatment 3 4\n", 2, 2).unwrap();
        assert_eq!((t, c), (vec![3.0, 4.0], vec![1.0, 2.0]));
    }

    #[test]
    fn norm_file_errors() {
        assert!(parse_norm("treatment 1 2\n", 2, 2).is_err());
        assert!(parse_norm("treatment 1 2\ncontrol 1\n", 2, 2).is_err());
        assert!(parse_norm("treatment 1 0\ncontrol 1 1\n", 2, 2).is_err());
    }
}
