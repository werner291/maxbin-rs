//! maxbin-core - the pure, host-agnostic core of maxbin-rs.
//!
//! This crate holds the parts of the algorithm that depend on nothing but Rust
//! and the data handed to them: tetranucleotide distance, k-mer profiling,
//! abundance/FASTA parsing, the fitted normal distributions, and the sort used
//! by the Spearman metric. It carries no external-tool calls, no threading
//! requirement, and (the gz and file-reading helpers aside) no OS dependency,
//! so it compiles to wasm and can run in a browser.
//!
//! The EM driver itself (emanager) does not live here yet: it still links the
//! system `lgamma` and uses rayon. It moves in once those are made portable.

pub mod abundance;
pub mod distance;
pub mod fasta;
pub mod kmer_map;
pub mod normal_distribution;
pub mod profiler;
pub mod quicksort;
