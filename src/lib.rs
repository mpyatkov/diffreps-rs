//! diffreps-rs: a Rust port of diffReps 1.55.6 (Li Shen, 2013), limited to the
//! negative binomial test and without annotation/hotspots.
//!
//! See `docs/porting-notes.md` for how each module maps to the Perl code.

pub mod background;
pub mod bed;
pub mod cdflib;
pub mod cli;
pub mod inputs;
pub mod nbtest;
pub mod padj;
pub mod report;
pub mod scan;
pub mod stats;
