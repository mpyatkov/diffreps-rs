//! M3 parity: background level and normalization factors.

mod common;

use common::*;
use diffreps::background::{self, Background, Method, estimate_background, estimate_norm};
use diffreps::report::f2;

fn header_background(r: &Report) -> Background {
    Background {
        mean: [r.floats("Treatment count mean"), r.floats("Control count mean")].concat(),
        dev: [r.floats("Treatment count deviation"), r.floats("Control count deviation")].concat(),
    }
}

/// With diffReps' own mean/deviation (from the header), our median-of-ratios
/// factors must print the same as diffReps' (2 decimals).
fn norm_matches(run: &str) {
    let r = Report::load(run);
    let (samples, n_tr) = samples();
    let window = r.number("Window size") as u64;
    let bg = header_background(&r);
    let norm = estimate_norm(samples, &chroms(), window, &bg, Some(2.0)).unwrap();
    let ours = f2(&norm);
    let theirs = [&r.header["Treatment normalization constant"][..], &r.header["Control normalization constant"][..]].concat();
    eprintln!("{run}: ours {ours:?}\n{run}: diffReps {theirs:?}");
    assert_eq!(ours, theirs);
    assert_eq!(*n_tr, 5);
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn norm_factors_diffreps_1000() {
    norm_matches("DIFFREPS_1000");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn norm_factors_diffreps_200() {
    norm_matches("DIFFREPS_200");
}

/// The robust estimate is random in diffReps, so only closeness can be
/// checked. diffReps' own two runs at window 1000 differ by up to ~2%
/// (11.87 vs 11.96); allow 5% or 0.02 absolute.
fn background_close(run: &str) {
    let r = Report::load(run);
    let (samples, _) = samples();
    let window = r.number("Window size") as u64;
    let theirs = header_background(&r);
    for (label, method) in [
        ("sampled", Method::Robust { alpha: 0.05, sample_size: 100_000, seed: 20131204 }),
        ("all bins", Method::Robust { alpha: 0.05, sample_size: 0, seed: 0 }),
    ] {
        let ours = estimate_background(samples, &chroms(), window, method);
        eprintln!("{run} {label}: mean {:?} dev {:?}", f2(&ours.mean), f2(&ours.dev));
        for (a, b) in ours.mean.iter().chain(&ours.dev).zip(theirs.mean.iter().chain(&theirs.dev)) {
            assert!((a - b).abs() <= 0.05 * b.abs() + 0.02, "{run} {label}: {a} vs {b}");
        }
    }
    eprintln!("{run} diffReps: mean {:?} dev {:?}", f2(&theirs.mean), f2(&theirs.dev));
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn background_close_1000() {
    background_close("DIFFREPS_1000");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn background_close_200() {
    background_close("DIFFREPS_200");
}

#[test]
#[ignore = "needs the TEST data (~4 GB)"]
fn background_from_report_reads_header() {
    let dir = test_dir().join("DIFFREPS_1000");
    let path = std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().path())
        .find(|p| p.file_name().unwrap().to_string_lossy().starts_with("2_")).unwrap();
    let bg = background::background_from_report(&path, 5, 3, 1000).unwrap();
    assert_eq!(bg.mean[0], 11.96);
    assert_eq!(bg.dev[7], 2.13);
    assert!(background::background_from_report(&path, 5, 3, 200).is_err());
}
