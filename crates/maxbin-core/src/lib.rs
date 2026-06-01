//! maxbin-core - the pure, host-agnostic core of maxbin-rs.
//!
//! This crate holds the parts of the algorithm that depend on nothing but Rust
//! and the data handed to them: tetranucleotide distance, k-mer profiling,
//! abundance/FASTA parsing, the fitted normal distributions, and the sort used
//! by the Spearman metric. It carries no external-tool calls, no threading
//! requirement, and (the gz and file-reading helpers aside) no OS dependency,
//! so it compiles to wasm and can run in a browser.
//!
//! The EM driver ([`emanager`]) lives here too: its `lgamma` is the pure-Rust
//! `libm` port and its rayon use is gated off on wasm32, so it runs serially in
//! the browser and in parallel natively. File parsing and result writing stay
//! in the maxbin-rs shell.

pub mod abundance;
pub mod distance;
pub mod emanager;
pub mod fasta;
pub mod kmer_map;
pub mod normal_distribution;
pub mod profiler;
pub mod quicksort;
