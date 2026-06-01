//! maxbin-rs - the CLI and pipeline shell around [`maxbin_core`].
//!
//! The pure algorithm (distance, k-mer profiling, parsing, normal
//! distributions, sorting) lives in `maxbin-core`. This crate adds the parts
//! that touch the OS and external tools: the CLI, the pipeline orchestration,
//! the EM driver, the file IO, and the calls out to bowtie2, HMMER, and the
//! gene caller.

pub mod cli;
pub mod emanager;
pub mod external;
pub mod paths;
pub mod pipeline;
