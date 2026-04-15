pub mod abundance;
pub mod cli;
pub mod distance;
pub mod emanager;
pub mod external;
pub mod fasta;
pub mod kmer_map;
pub mod normal_distribution;
pub mod pipeline;
pub mod profiler;
pub mod quicksort;

// Always compiled — the C++ FFI is linked by build.rs regardless, and
// integration tests need access to these bindings.
pub mod original_ffi;
