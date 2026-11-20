//! Helpers for the parity tests against the diffReps reference runs in `../TEST`.
//!
//! Layout expected (override the root with `DIFFREPS_TEST_DIR`):
//! ```text
//! TEST/bed/*.bed               input BED files
//! TEST/mm10.chrom.sizes
//! TEST/<RUN>/2_<RUN>_*         diffReps report
//! TEST/<RUN>/norm.txt          only for RIPPM_* runs
//! ```
//! These tests are `#[ignore]`d because they need ~4 GB of data:
//! `cargo test --release -- --ignored`.
#![allow(dead_code)]

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;

use diffreps::bed::{self, Sample};
use diffreps::inputs::{self, Chrom};
use rayon::prelude::*;

pub fn test_dir() -> PathBuf {
    std::env::var_os("DIFFREPS_TEST_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../TEST"))
}

/// A parsed diffReps report.
pub struct Report {
    /// Header `# key: values` lines, values split on tabs.
    pub header: HashMap<String, Vec<String>>,
    pub rows: Vec<Row>,
}

#[derive(Debug, Clone)]
pub struct Row {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub tr_cnt: Vec<f64>,
    pub co_cnt: Vec<f64>,
    pub event: String,
    pub log2fc: String,
    pub pval: f64,
    pub padj: f64,
    pub win_sta: u64,
    pub win_end: u64,
    pub win_fc: String,
    pub win_p: f64,
    pub win_q: f64,
    /// The raw text line.
    pub line: String,
}

impl Report {
    pub fn load(run: &str) -> Report {
        let dir = test_dir().join(run);
        let path = fs::read_dir(&dir)
            .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
            .map(|e| e.unwrap().path())
            .find(|p| p.file_name().unwrap().to_string_lossy().starts_with("2_"))
            .expect("report file 2_* not found");
        Report::parse(&fs::read_to_string(path).unwrap())
    }

    pub fn parse(text: &str) -> Report {
        let mut header = HashMap::new();
        let mut rows = Vec::new();
        for line in text.lines() {
            if let Some(h) = line.strip_prefix('#') {
                if let Some((k, v)) = h.split_once(": ") {
                    let key = k.trim_start_matches('#').trim().to_string();
                    header.insert(key, v.split('\t').map(str::to_string).collect());
                }
                continue;
            }
            if line.starts_with("Chrom\t") || line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let nums = |s: &str| s.split(';').map(|x| x.parse().unwrap()).collect();
            rows.push(Row {
                chrom: f[0].into(),
                start: f[1].parse().unwrap(),
                end: f[2].parse().unwrap(),
                tr_cnt: nums(f[4]),
                co_cnt: nums(f[5]),
                event: f[10].into(),
                log2fc: f[11].into(),
                pval: f[12].parse().unwrap(),
                padj: f[13].parse().unwrap(),
                win_sta: f[14].parse().unwrap(),
                win_end: f[15].parse().unwrap(),
                win_fc: f[16].into(),
                win_p: f[17].parse().unwrap(),
                win_q: f[18].parse().unwrap(),
                line: line.to_string(),
            });
        }
        Report { header, rows }
    }

    pub fn floats(&self, key: &str) -> Vec<f64> {
        self.header[key].iter().map(|v| v.parse().unwrap()).collect()
    }

    pub fn number(&self, key: &str) -> f64 {
        self.header[key][0].parse().unwrap()
    }

    pub fn files(&self, key: &str) -> Vec<PathBuf> {
        self.header[key].iter().map(|f| test_dir().join("bed").join(f)).collect()
    }
}

pub fn chroms() -> Vec<Chrom> {
    inputs::read_chrlen(&test_dir().join("mm10.chrom.sizes")).unwrap()
}

/// All 8 samples, read once per test binary, in the reference order
/// (treatment files then control files, as in every TEST report).
pub fn samples() -> &'static (Vec<Sample>, usize) {
    static SAMPLES: OnceLock<(Vec<Sample>, usize)> = OnceLock::new();
    SAMPLES.get_or_init(|| {
        let r = Report::load("RIPPM_1000");
        let tr = r.files("Treatment files");
        let co = r.files("Control files");
        let chroms = chroms();
        let all: Vec<PathBuf> = tr.iter().chain(&co).cloned().collect();
        let samples = all
            .par_iter()
            .map(|p| bed::read_sample(p, &chroms, 0).unwrap())
            .collect();
        (samples, tr.len())
    })
}

pub fn chrom_index(chroms: &[Chrom], name: &str) -> usize {
    chroms.iter().position(|c| c.name == name).unwrap()
}
