/// Property-based equivalence tests for kmer_map.rs.
///
/// Coverage:
/// - KmerMap::new() / get_entry_num() / get_full_table(): kmermap_mapping_equivalence (50 cases, kmerlen 1-5)
/// - get_mapping() / get_reverse_mapping_*(): kmermap_per_kmer_proptest (variable kmerlen 2-5)
/// - Invalid k-mers (non-ACGT): kmermap_invalid_kmer_equivalence
/// - Lowercase k-mers: kmermap_lowercase_kmer_equivalence

mod kmermap_proptest {
    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(50))]

        #[test]
        fn kmermap_mapping_equivalence(kmerlen in 1..6i32, symmetric in proptest::bool::ANY) {
            let rust_km = maxbin_rs::kmer_map::KmerMap::new(kmerlen as usize, symmetric);
            let cpp_km = maxbin_rs::original_ffi::OriginalKmerMap::new(kmerlen, symmetric);

            let rust_entry = rust_km.get_entry_num() as i32;
            let cpp_entry = cpp_km.entry_num();
            prop_assert_eq!(rust_entry, cpp_entry, "entry_num mismatch for kmerlen={} sym={}", kmerlen, symmetric);

            let cpp_table = cpp_km.get_full_table();
            let rust_table = rust_km.get_full_table();

            prop_assert_eq!(rust_table.len(), cpp_table.len(), "table length mismatch");

            for i in 0..rust_table.len() {
                prop_assert_eq!(rust_table[i], cpp_table[i], "mapping mismatch at index {} kmerlen={}", i, kmerlen);
            }
        }
    }
}

/// Per-kmer lookup tests with variable kmerlen (2-5), including invalid and lowercase kmers.
mod kmermap_per_kmer_proptest {
    use proptest::prelude::*;

    /// Generate a (kmerlen, symmetric, kmer_string) triple with matching lengths.
    fn arb_kmerlen_and_kmer() -> impl Strategy<Value = (usize, bool, String)> {
        (2usize..6, proptest::bool::ANY).prop_flat_map(|(kmerlen, symmetric)| {
            let kmer = prop::collection::vec(
                prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
                kmerlen,
            )
            .prop_map(|chars| chars.into_iter().collect::<String>());
            kmer.prop_map(move |k| (kmerlen, symmetric, k))
        })
    }

    /// Generate a (kmerlen, symmetric, index) triple with valid index range.
    fn arb_kmerlen_and_index() -> impl Strategy<Value = (usize, bool, u32)> {
        (2usize..6, proptest::bool::ANY).prop_flat_map(|(kmerlen, symmetric)| {
            let max_idx = 4u32.pow(kmerlen as u32);
            (0u32..max_idx).prop_map(move |idx| (kmerlen, symmetric, idx))
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(200))]

        /// get_mapping() with variable kmerlen — all results should be non-negative for valid ACGT kmers.
        #[test]
        fn kmermap_get_mapping(args in arb_kmerlen_and_kmer()) {
            let (kmerlen, symmetric, kmer) = args;
            let rust_km = maxbin_rs::kmer_map::KmerMap::new(kmerlen, symmetric);
            let cpp_km = maxbin_rs::original_ffi::OriginalKmerMap::new(kmerlen as i32, symmetric);

            let rust_idx = rust_km.get_mapping(kmer.as_bytes());
            let cpp_idx = cpp_km.get_mapping(&kmer);

            prop_assert!(rust_idx >= 0, "rust returned -1 for valid kmer '{}' kmerlen={}", kmer, kmerlen);
            prop_assert_eq!(
                rust_idx, cpp_idx,
                "get_mapping mismatch for kmer='{}' kmerlen={} symmetric={}", kmer, kmerlen, symmetric
            );
        }

        /// get_reverse_mapping_str() with variable kmerlen.
        #[test]
        fn kmermap_get_reverse_mapping_str(args in arb_kmerlen_and_kmer()) {
            let (kmerlen, symmetric, kmer) = args;
            let rust_km = maxbin_rs::kmer_map::KmerMap::new(kmerlen, symmetric);
            let cpp_km = maxbin_rs::original_ffi::OriginalKmerMap::new(kmerlen as i32, symmetric);

            let rust_idx = rust_km.get_reverse_mapping_str(kmer.as_bytes());
            let cpp_idx = cpp_km.get_reverse_mapping_str(&kmer);

            prop_assert_eq!(
                rust_idx, cpp_idx,
                "get_reverse_mapping_str mismatch for kmer='{}' kmerlen={} symmetric={}", kmer, kmerlen, symmetric
            );
        }

        /// get_reverse_mapping_idx() with variable kmerlen and valid index range.
        #[test]
        fn kmermap_get_reverse_mapping_idx(args in arb_kmerlen_and_index()) {
            let (kmerlen, symmetric, index) = args;
            let rust_km = maxbin_rs::kmer_map::KmerMap::new(kmerlen, symmetric);
            let cpp_km = maxbin_rs::original_ffi::OriginalKmerMap::new(kmerlen as i32, symmetric);

            let rust_result = rust_km.get_reverse_mapping_idx(index);
            let cpp_result = cpp_km.get_reverse_mapping_idx(index as i32);

            prop_assert_eq!(
                rust_result, cpp_result,
                "get_reverse_mapping_idx mismatch for index={} kmerlen={} symmetric={}", index, kmerlen, symmetric
            );
        }

        /// Invalid k-mers (non-ACGT characters) must return the same value from both implementations.
        #[test]
        fn kmermap_invalid_kmer_equivalence(
            symmetric in proptest::bool::ANY,
            kmer in "[A-Za-z0-9NRYWSKM]{4}"
        ) {
            let rust_km = maxbin_rs::kmer_map::KmerMap::new(4, symmetric);
            let cpp_km = maxbin_rs::original_ffi::OriginalKmerMap::new(4, symmetric);

            let rust_idx = rust_km.get_mapping(kmer.as_bytes());
            let cpp_idx = cpp_km.get_mapping(&kmer);

            prop_assert_eq!(
                rust_idx, cpp_idx,
                "get_mapping mismatch for possibly-invalid kmer='{}' symmetric={}", kmer, symmetric
            );
        }

        /// Lowercase k-mers must be handled identically by both implementations.
        #[test]
        fn kmermap_lowercase_kmer_equivalence(
            symmetric in proptest::bool::ANY,
            kmer in "[acgtACGT]{4}"
        ) {
            let rust_km = maxbin_rs::kmer_map::KmerMap::new(4, symmetric);
            let cpp_km = maxbin_rs::original_ffi::OriginalKmerMap::new(4, symmetric);

            let rust_idx = rust_km.get_mapping(kmer.as_bytes());
            let cpp_idx = cpp_km.get_mapping(&kmer);

            prop_assert_eq!(
                rust_idx, cpp_idx,
                "get_mapping mismatch for lowercase kmer='{}' symmetric={}", kmer, symmetric
            );
        }

        /// Symmetric reverse-complement invariant with variable kmerlen.
        #[test]
        fn kmermap_symmetric_revcomp_invariant(args in arb_kmerlen_and_kmer()) {
            let (kmerlen, _, kmer) = args;
            let rust_km = maxbin_rs::kmer_map::KmerMap::new(kmerlen, true);

            let fwd = rust_km.get_mapping(kmer.as_bytes());
            let rev = rust_km.get_reverse_mapping_str(kmer.as_bytes());

            prop_assert_eq!(
                fwd, rev,
                "symmetric map: get_mapping('{}')={} != get_reverse_mapping_str('{}')={} kmerlen={}",
                kmer, fwd, kmer, rev, kmerlen
            );
        }
    }
}
