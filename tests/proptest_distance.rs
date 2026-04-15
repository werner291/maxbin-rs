/// Property-based equivalence tests for distance.rs.
///
/// Coverage:
/// - euc_dist_seq() / spearman_dist_seq(): FFI equivalence (50 cases)
/// - euc_dist_profiles() / spearman_dist_profiles(): FFI equivalence (100 cases)
/// - Algebraic: symmetry, self-distance = 0
///
/// Generators exercise: short sequences (< kmerlen, producing all-zero profiles),
/// normal sequences (20-200bp), long sequences (500-2000bp), N-containing sequences.

mod distance_proptest {
    use proptest::prelude::*;

    fn arb_dna_seq() -> impl Strategy<Value = String> {
        prop_oneof![
            // Short sequences (1-4 bases, includes sub-kmerlen producing NaN/zero profiles)
            1 => prop::collection::vec(
                prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
                1..5
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // Normal range with some N characters
            6 => prop::collection::vec(
                prop_oneof![
                    9 => prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
                    1 => Just('N'),
                ],
                20..200
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // Long sequences
            2 => prop::collection::vec(
                prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
                500..2000
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // All-N (degenerate: all-zero profile, tests rank tie-breaking; must be >= kmerlen)
            1 => (4usize..30).prop_map(|len| "N".repeat(len)),
        ]
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(50))]

        #[test]
        fn euc_dist_equivalence(seq1 in arb_dna_seq(), seq2 in arb_dna_seq()) {
            let ctx = maxbin_rs::distance::DistanceContext::new(4);
            let cpp = maxbin_rs::original_ffi::OriginalEucDist::new(4);

            let rust_d = maxbin_rs::distance::euc_dist_seq(&ctx, seq1.as_bytes(), seq2.as_bytes());
            let cpp_d = cpp.get_dist_seq(&seq1, &seq2);

            let diff = (rust_d - cpp_d).abs();
            prop_assert!(
                diff < 1e-12,
                "EucDist mismatch: rust={rust_d} cpp={cpp_d} seq1_len={} seq2_len={}",
                seq1.len(), seq2.len()
            );
        }

        #[test]
        fn spearman_dist_equivalence(seq1 in arb_dna_seq(), seq2 in arb_dna_seq()) {
            let ctx = maxbin_rs::distance::DistanceContext::new(4);
            let cpp = maxbin_rs::original_ffi::OriginalSpearmanDist::new(4);

            let rust_d = maxbin_rs::distance::spearman_dist_seq(&ctx, seq1.as_bytes(), seq2.as_bytes());
            let cpp_d = cpp.get_dist_seq(&seq1, &seq2);

            let diff = (rust_d - cpp_d).abs();
            prop_assert!(
                diff < 1e-10,
                "SpearmanDist mismatch: rust={rust_d} cpp={cpp_d} seq1_len={} seq2_len={}",
                seq1.len(), seq2.len()
            );
        }
    }
}

/// Profile-based distance equivalence tests with diverse input profiles.
mod distance_profiles_proptest {
    use proptest::prelude::*;

    fn arb_dna_seq() -> impl Strategy<Value = String> {
        prop_oneof![
            // Sub-kmerlen (all-zero profile)
            1 => prop::collection::vec(
                prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
                1..4
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // Normal with N
            7 => prop::collection::vec(
                prop_oneof![
                    9 => prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
                    1 => Just('N'),
                ],
                30..200
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // Long
            1 => prop::collection::vec(
                prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
                500..2000
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // All-N
            1 => (5usize..30).prop_map(|len| "N".repeat(len)),
        ]
    }

    /// Generate a profile from a random DNA sequence.
    fn arb_profile() -> impl Strategy<Value = Vec<f64>> {
        arb_dna_seq().prop_map(|seq| {
            let kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
            let prof = maxbin_rs::profiler::Profiler::new(4, seq.as_bytes(), &kmap);
            prof.get_profile().to_vec()
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(100))]

        #[test]
        fn euc_dist_profiles_equivalence(
            seq1 in arb_dna_seq(),
            seq2 in arb_dna_seq()
        ) {
            let kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
            let prof1 = maxbin_rs::profiler::Profiler::new(4, seq1.as_bytes(), &kmap);
            let prof2 = maxbin_rs::profiler::Profiler::new(4, seq2.as_bytes(), &kmap);
            let p1 = prof1.get_profile().to_vec();
            let p2 = prof2.get_profile().to_vec();

            let rust_d = maxbin_rs::distance::euc_dist_profiles(&p1, &p2);
            let cpp = maxbin_rs::original_ffi::OriginalEucDist::new(4);
            let cpp_d = cpp.get_dist_profile(&p1, &p2);

            let diff = (rust_d - cpp_d).abs();
            prop_assert!(
                diff < 1e-12,
                "euc_dist_profiles mismatch: rust={rust_d} cpp={cpp_d} diff={diff}"
            );
        }

        #[test]
        fn spearman_dist_profiles_equivalence(
            seq1 in arb_dna_seq(),
            seq2 in arb_dna_seq()
        ) {
            let kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
            let prof1 = maxbin_rs::profiler::Profiler::new(4, seq1.as_bytes(), &kmap);
            let prof2 = maxbin_rs::profiler::Profiler::new(4, seq2.as_bytes(), &kmap);
            let p1 = prof1.get_profile().to_vec();
            let p2 = prof2.get_profile().to_vec();

            let ctx = maxbin_rs::distance::DistanceContext::new(4);
            let rust_d = maxbin_rs::distance::spearman_dist_profiles(&ctx, &p1, &p2);
            let cpp = maxbin_rs::original_ffi::OriginalSpearmanDist::new(4);
            let cpp_d = cpp.get_dist_profile(&p1, &p2);

            let diff = (rust_d - cpp_d).abs();
            prop_assert!(
                diff < 1e-10,
                "spearman_dist_profiles mismatch: rust={rust_d} cpp={cpp_d} diff={diff}"
            );
        }

        #[test]
        fn euc_dist_profiles_symmetry(p1 in arb_profile(), p2 in arb_profile()) {
            let d12 = maxbin_rs::distance::euc_dist_profiles(&p1, &p2);
            let d21 = maxbin_rs::distance::euc_dist_profiles(&p2, &p1);
            let diff = (d12 - d21).abs();
            prop_assert!(diff < 1e-14, "euc_dist not symmetric: d12={d12} d21={d21}");
        }

        #[test]
        fn spearman_dist_profiles_symmetry(p1 in arb_profile(), p2 in arb_profile()) {
            let ctx = maxbin_rs::distance::DistanceContext::new(4);
            let d12 = maxbin_rs::distance::spearman_dist_profiles(&ctx, &p1, &p2);
            let d21 = maxbin_rs::distance::spearman_dist_profiles(&ctx, &p2, &p1);
            let diff = (d12 - d21).abs();
            prop_assert!(diff < 1e-10, "spearman_dist not symmetric: d12={d12} d21={d21}");
        }

        #[test]
        fn euc_dist_profiles_self_distance(p in arb_profile()) {
            let d = maxbin_rs::distance::euc_dist_profiles(&p, &p);
            prop_assert!(d.abs() < 1e-14, "self euc_dist should be 0, got {d}");
        }

        #[test]
        fn spearman_dist_profiles_self_distance(p in arb_profile()) {
            let ctx = maxbin_rs::distance::DistanceContext::new(4);
            let d = maxbin_rs::distance::spearman_dist_profiles(&ctx, &p, &p);
            prop_assert!(d.abs() < 1e-10, "self spearman_dist should be 0, got {d}");
        }
    }
}
