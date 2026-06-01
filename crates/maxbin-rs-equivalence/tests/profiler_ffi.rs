//! FFI equivalence test for the tetranucleotide profiler.
//!
//! Compares `maxbin_rs::profiler::Profiler` against the original MaxBin2 C++
//! Profiler. Previously an in-source test in `profiler.rs`.
use maxbin_rs::kmer_map::KmerMap;
use maxbin_rs::profiler::Profiler;
use maxbin_rs_equivalence::original_ffi::{OriginalKmerMap, OriginalProfiler};

#[test]
fn ffi_equivalence() {
    let test_seqs = [
        "ACGTACGTACGTACGT",
        "AAAAAAAAAAAA",
        "ACGTTTTTGGGGCCCC",
        "ATCGATCGATCGATCGATCG",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
    ];

    let rust_kmap = KmerMap::new(4, true);
    let cpp_kmap = OriginalKmerMap::new(4, true);
    let entry_num = rust_kmap.get_entry_num();

    for seq in &test_seqs {
        let rust_prof = Profiler::new(4, seq.as_bytes(), &rust_kmap);
        let cpp_prof = OriginalProfiler::new(4, seq, &cpp_kmap);

        let rust_profile = rust_prof.get_profile();
        let cpp_profile = cpp_prof.get_profile(entry_num as i32);

        assert_eq!(rust_profile.len(), cpp_profile.len());
        for i in 0..rust_profile.len() {
            assert_eq!(
                rust_profile[i].to_bits(),
                cpp_profile[i].to_bits(),
                "profile not bit-identical at index {i} for seq '{seq}': rust={:e} cpp={:e}",
                rust_profile[i],
                cpp_profile[i]
            );
        }

        let rust_pn = rust_prof.get_percent_n();
        let cpp_pn = cpp_prof.get_percent_n();
        assert!(
            (rust_pn - cpp_pn).abs() < 1e-6,
            "percent_n mismatch for seq '{seq}': rust={rust_pn} cpp={cpp_pn}"
        );
    }
}
