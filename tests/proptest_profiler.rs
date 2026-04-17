/// Property-based equivalence tests for profiler.rs.
///
/// Coverage:
/// - Profiler::new() / get_profile() / get_percent_n(): profiler_equivalence (100 cases)
///   Generators exercise: ACGT bases, N characters (~10%), mixed sequences,
///   short sequences (1-3 bases, below kmerlen), long sequences (500-2000bp).
/// - reset() / add_profile() / calc_profile(): algebraic invariant tests
///   (C++ FFI does not expose addProfile/calcProfile directly)
mod profiler_proptest {
    use proptest::prelude::*;

    fn arb_dna_seq() -> impl Strategy<Value = String> {
        prop_oneof![
            // Normal: ACGT + ~10% N characters, length 10-200
            6 => prop::collection::vec(
                prop_oneof![
                    9 => prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
                    1 => Just('N'),
                ],
                10..200,
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // Short sequences (below kmerlen=4): both C++ and Rust produce NaN for percent_n
            1 => prop::collection::vec(
                prop_oneof![Just('A'), Just('C'), Just('G'), Just('T'), Just('N')],
                1..4,
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // Exactly kmerlen (boundary: 1 window)
            1 => prop::collection::vec(
                prop_oneof![Just('A'), Just('C'), Just('G'), Just('T'), Just('N')],
                4..=4,
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // Long sequences
            1 => prop::collection::vec(
                prop_oneof![
                    9 => prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
                    1 => Just('N'),
                ],
                500..2000,
            ).prop_map(|chars| chars.into_iter().collect::<String>()),
            // All-N sequences (>= kmerlen)
            1 => (4usize..50).prop_map(|len| "N".repeat(len)),
        ]
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(100))]

        #[test]
        fn profiler_equivalence(seq in arb_dna_seq()) {
            let rust_kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
            let cpp_kmap = maxbin_rs::original_ffi::OriginalKmerMap::new(4, true);
            let entry_num = rust_kmap.get_entry_num();

            let rust_prof = maxbin_rs::profiler::Profiler::new(4, seq.as_bytes(), &rust_kmap);
            let cpp_prof = maxbin_rs::original_ffi::OriginalProfiler::new(4, &seq, &cpp_kmap);

            let rust_profile = rust_prof.get_profile();
            let cpp_profile = cpp_prof.get_profile(entry_num as i32);

            for i in 0..rust_profile.len() {
                let diff = (rust_profile[i] - cpp_profile[i]).abs();
                prop_assert!(
                    diff < 1e-12,
                    "profile mismatch at index {i}: rust={} cpp={}",
                    rust_profile[i], cpp_profile[i]
                );
            }

            let rust_pn = rust_prof.get_percent_n();
            let cpp_pn = cpp_prof.get_percent_n();
            if rust_pn.is_nan() && cpp_pn.is_nan() {
                // Both NaN — matches (sub-kmerlen sequences produce 0/0 = NaN)
            } else {
                prop_assert!(
                    !rust_pn.is_nan() && !cpp_pn.is_nan(),
                    "percent_n NaN divergence: rust={} cpp={} seq_len={}",
                    rust_pn, cpp_pn, seq.len()
                );
                let pn_diff = (rust_pn - cpp_pn).abs();
                prop_assert!(pn_diff < 1e-6, "percent_n mismatch: rust={} cpp={}", rust_pn, cpp_pn);
            }
        }
    }
}

/// Tests for Profiler::reset(), add_profile(), and calc_profile().
///
/// These methods implement the M-step weighted average in the EM algorithm.
/// The C++ OriginalProfiler FFI does not expose addProfile/calcProfile directly,
/// so we test via algebraic invariants:
///
/// 1. reset() zeroes the profile.
/// 2. add_profile() is linear: accumulate(w1*p1 + w2*p2) == w1*p1 + w2*p2 elementwise.
/// 3. calc_profile(w) divides by w: after add_profile(p, w) then calc_profile(w),
///    the result equals p (for a single-item accumulation).
/// 4. Weighted average of two identical profiles equals that profile (any weights).
/// 5. The result of add+calc agrees with manually computing the weighted average
///    using only the already-FFI-validated get_profile() values.
mod profiler_weighted_avg_proptest {
    use proptest::prelude::*;

    fn arb_dna_seq() -> impl Strategy<Value = String> {
        prop::collection::vec(
            prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
            50..300,
        )
        .prop_map(|chars| chars.into_iter().collect::<String>())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(100))]

        /// add_profile then calc_profile with a single contributor and weight=1
        /// must reproduce the original profile exactly.
        #[test]
        fn add_calc_single_profile_identity(seq in arb_dna_seq()) {
            let kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
            let source = maxbin_rs::profiler::Profiler::new(4, seq.as_bytes(), &kmap);
            let source_profile = source.get_profile().to_vec();

            // Create accumulator: fresh Profiler from a dummy seq then reset
            let mut accum = maxbin_rs::profiler::Profiler::new(4, seq.as_bytes(), &kmap);
            accum.reset();

            // Add the source with weight=1
            accum.add_profile(&source, 1.0);
            accum.calc_profile(1.0);

            let result = accum.get_profile();
            prop_assert_eq!(result.len(), source_profile.len());
            for i in 0..result.len() {
                let diff = (result[i] - source_profile[i]).abs();
                prop_assert!(
                    diff < 1e-14,
                    "add_calc identity failed at index {i}: got {} expected {}",
                    result[i], source_profile[i]
                );
            }
        }

        /// add_profile then calc_profile with two profiles and equal weights must
        /// equal their arithmetic mean elementwise.
        #[test]
        fn add_calc_two_profiles_equal_weight(seq1 in arb_dna_seq(), seq2 in arb_dna_seq()) {
            let kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
            let prof1 = maxbin_rs::profiler::Profiler::new(4, seq1.as_bytes(), &kmap);
            let prof2 = maxbin_rs::profiler::Profiler::new(4, seq2.as_bytes(), &kmap);

            let p1 = prof1.get_profile().to_vec();
            let p2 = prof2.get_profile().to_vec();

            // Compute expected mean manually
            let expected: Vec<f64> = p1.iter().zip(p2.iter()).map(|(a, b)| (a + b) / 2.0).collect();

            // Accumulate with equal weights
            let mut accum = maxbin_rs::profiler::Profiler::new(4, seq1.as_bytes(), &kmap);
            accum.reset();
            accum.add_profile(&prof1, 1.0);
            accum.add_profile(&prof2, 1.0);
            accum.calc_profile(2.0);

            let result = accum.get_profile();
            prop_assert_eq!(result.len(), expected.len());
            for i in 0..result.len() {
                let diff = (result[i] - expected[i]).abs();
                prop_assert!(
                    diff < 1e-14,
                    "equal-weight mean failed at index {i}: got {} expected {}",
                    result[i], expected[i]
                );
            }
        }

        /// Weighted average with arbitrary positive weights.
        /// result[i] = (w1*p1[i] + w2*p2[i]) / (w1 + w2)
        #[test]
        fn add_calc_weighted_average(
            seq1 in arb_dna_seq(),
            seq2 in arb_dna_seq(),
            w1 in 0.1f64..100.0,
            w2 in 0.1f64..100.0
        ) {
            let kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
            let prof1 = maxbin_rs::profiler::Profiler::new(4, seq1.as_bytes(), &kmap);
            let prof2 = maxbin_rs::profiler::Profiler::new(4, seq2.as_bytes(), &kmap);

            let p1 = prof1.get_profile().to_vec();
            let p2 = prof2.get_profile().to_vec();

            let total = w1 + w2;
            let expected: Vec<f64> = p1.iter().zip(p2.iter())
                .map(|(a, b)| (w1 * a + w2 * b) / total)
                .collect();

            let mut accum = maxbin_rs::profiler::Profiler::new(4, seq1.as_bytes(), &kmap);
            accum.reset();
            accum.add_profile(&prof1, w1);
            accum.add_profile(&prof2, w2);
            accum.calc_profile(total);

            let result = accum.get_profile();
            prop_assert_eq!(result.len(), expected.len());
            for i in 0..result.len() {
                let diff = (result[i] - expected[i]).abs();
                prop_assert!(
                    diff < 1e-12,
                    "weighted average failed at index {i}: got {} expected {} (w1={} w2={})",
                    result[i], expected[i], w1, w2
                );
            }
        }

        /// calc_profile(0.0) must zero the profile (not produce NaN/Inf).
        #[test]
        fn calc_profile_zero_weight_produces_zeros(seq in arb_dna_seq()) {
            let kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
            let prof = maxbin_rs::profiler::Profiler::new(4, seq.as_bytes(), &kmap);
            let mut accum = maxbin_rs::profiler::Profiler::new(4, seq.as_bytes(), &kmap);
            accum.reset();
            accum.add_profile(&prof, 1.0);
            accum.calc_profile(0.0); // should zero everything

            for v in accum.get_profile() {
                prop_assert_eq!(*v, 0.0, "calc_profile(0) should zero all entries, got {}", v);
            }
        }

        /// reset() must zero all entries regardless of prior state.
        #[test]
        fn reset_zeroes_all_entries(seq in arb_dna_seq()) {
            let kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
            let source = maxbin_rs::profiler::Profiler::new(4, seq.as_bytes(), &kmap);
            let mut accum = maxbin_rs::profiler::Profiler::new(4, seq.as_bytes(), &kmap);

            // Accumulate some state
            accum.add_profile(&source, 5.0);
            // Then reset
            accum.reset();

            for v in accum.get_profile() {
                prop_assert_eq!(*v, 0.0, "reset() should zero all entries, got {}", v);
            }
        }
    }
}
