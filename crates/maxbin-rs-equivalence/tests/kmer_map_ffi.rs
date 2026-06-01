//! FFI equivalence tests for the k-mer mapping table.
//!
//! Compares `maxbin_core::kmer_map::KmerMap` against the original MaxBin2 C++
//! kmerMap. Previously in-source tests in `kmer_map.rs`.
use maxbin_core::kmer_map::KmerMap;
use maxbin_rs_equivalence::original_ffi::OriginalKmerMap;

#[test]
fn ffi_equivalence_kmer4_symmetric() {
    let rust_km = KmerMap::new(4, true);
    let cpp_km = OriginalKmerMap::new(4, true);

    assert_eq!(rust_km.get_entry_num() as i32, cpp_km.entry_num());

    let cpp_table = cpp_km.get_full_table();
    let rust_table = rust_km.get_full_table();

    assert_eq!(rust_table.len(), cpp_table.len());
    for i in 0..rust_table.len() {
        assert_eq!(rust_table[i], cpp_table[i], "mapping mismatch at index {i}");
    }
}

#[test]
fn ffi_equivalence_kmer4_nonsymmetric() {
    let rust_km = KmerMap::new(4, false);
    let cpp_km = OriginalKmerMap::new(4, false);

    assert_eq!(rust_km.get_entry_num() as i32, cpp_km.entry_num());

    let cpp_table = cpp_km.get_full_table();
    let rust_table = rust_km.get_full_table();

    for i in 0..rust_table.len() {
        assert_eq!(rust_table[i], cpp_table[i]);
    }
}

#[test]
fn ffi_equivalence_kmer3_symmetric() {
    let rust_km = KmerMap::new(3, true);
    let cpp_km = OriginalKmerMap::new(3, true);

    assert_eq!(rust_km.get_entry_num() as i32, cpp_km.entry_num());

    let cpp_table = cpp_km.get_full_table();
    let rust_table = rust_km.get_full_table();

    for i in 0..rust_table.len() {
        assert_eq!(
            rust_table[i], cpp_table[i],
            "kmer3 symmetric mapping mismatch at index {i}"
        );
    }
}
