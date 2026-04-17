/// Property-based equivalence tests for normal_distribution.rs.
///
/// - NormalDistribution::new() / prob(): covered by normal_dist_equivalence (proptest, 200 cases)
///   Inline FFI test: normal_distribution.rs::tests::ffi_equivalence (hardcoded)
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    #[test]
    fn normal_dist_equivalence(
        mean in -100.0f64..100.0,
        std_dev in 0.01f64..50.0,
        input in -200.0f64..200.0
    ) {
        let rust_nd = maxbin_rs::normal_distribution::NormalDistribution::new(mean, std_dev);
        let cpp_nd = maxbin_rs::original_ffi::OriginalNormalDistribution::new(mean, std_dev);

        let rust_p = rust_nd.prob(input);
        let cpp_p = cpp_nd.prob(input);

        let diff = (rust_p - cpp_p).abs();
        // Use relative error for very small values
        let tol = if cpp_p.abs() > 1e-300 {
            (diff / cpp_p.abs()).max(diff)
        } else {
            diff
        };
        prop_assert!(
            tol < 1e-12,
            "NormalDist mismatch: mean={mean} std={std_dev} input={input} \
             rust={rust_p} cpp={cpp_p} diff={diff}"
        );
    }
}
