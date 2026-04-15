/// Property-based tests for emanager.rs probability functions.
///
/// Coverage:
/// - get_prob_dist(): validated against C++ FFI NormalDistribution on distance range 0.0..2.0.
///   Calls the real source function directly.
///
/// - get_prob_abund(): tested via algebraic properties (non-negativity, zero inputs,
///   Poisson peak). Calls the real source function directly. C++ FFI comparison is not
///   available for this function individually, but the full EM pipeline integration test
///   (tests/emanager_equivalence.rs) compares outputs byte-for-byte.
///
/// - init_em() / run_em() / classify() / write_results() / run_pipeline(): covered by
///   tests/emanager_equivalence.rs (integration tests, including synthetic full-pipeline test)

use proptest::prelude::*;

/// Compute prob_dist via C++ FFI NormalDistribution for cross-checking.
fn prob_dist_via_cpp_ffi(distance: f64) -> f64 {
    let intra = maxbin_rs::original_ffi::OriginalNormalDistribution::new(0.0, 0.01037897 / 2.0);
    let inter = maxbin_rs::original_ffi::OriginalNormalDistribution::new(0.0676654, 0.03419337);
    let d_intra = intra.prob(distance);
    let d_inter = inter.prob(distance);
    d_intra / (d_inter + d_intra)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(300))]

    /// get_prob_dist: real source function must match C++ FFI NormalDist.
    #[test]
    fn prob_dist_matches_cpp_normaldist(distance in 0.0f64..2.0) {
        let intra = maxbin_rs::emanager::make_intra_normal();
        let inter = maxbin_rs::emanager::make_inter_normal();
        let rust_val = maxbin_rs::emanager::get_prob_dist(distance, &intra, &inter);
        let cpp_val = prob_dist_via_cpp_ffi(distance);

        // At very large distances, both normal PDFs underflow to 0, producing 0/0 = NaN.
        if rust_val.is_nan() && cpp_val.is_nan() {
            return Ok(());
        }
        prop_assert!(
            !rust_val.is_nan() && !cpp_val.is_nan(),
            "NaN divergence at distance={distance}: rust={rust_val} cpp={cpp_val}"
        );

        let diff = (rust_val - cpp_val).abs();
        prop_assert!(
            diff < 1e-12,
            "get_prob_dist mismatch at distance={distance}: rust={rust_val} cpp={cpp_val}"
        );
    }

    /// get_prob_dist must be in [0, 1] for all non-negative distances.
    #[test]
    fn prob_dist_range(distance in 0.0f64..1.0) {
        let intra = maxbin_rs::emanager::make_intra_normal();
        let inter = maxbin_rs::emanager::make_inter_normal();
        let p = maxbin_rs::emanager::get_prob_dist(distance, &intra, &inter);
        prop_assert!(p >= 0.0, "prob_dist negative at distance={distance}: {p}");
        prop_assert!(p <= 1.0, "prob_dist > 1 at distance={distance}: {p}");
    }

    /// get_prob_dist is monotonically decreasing: at d=0 it's near 1,
    /// at d=0.5 it's near 0 (intra dominates at small distances).
    #[test]
    fn prob_dist_monotone_behavior(
        d_small in 0.0f64..0.02,
        d_large in 0.1f64..0.5
    ) {
        let intra = maxbin_rs::emanager::make_intra_normal();
        let inter = maxbin_rs::emanager::make_inter_normal();
        let p_small = maxbin_rs::emanager::get_prob_dist(d_small, &intra, &inter);
        let p_large = maxbin_rs::emanager::get_prob_dist(d_large, &intra, &inter);
        prop_assert!(
            p_small > p_large,
            "prob_dist should decrease with distance: p({d_small})={p_small} > p({d_large})={p_large} failed"
        );
    }

    /// get_prob_abund: zero inputs return 0.
    #[test]
    fn prob_abund_zero_inputs(val in 1.0f64..100.0) {
        prop_assert_eq!(maxbin_rs::emanager::get_prob_abund(0.0, val), 0.0);
        prop_assert_eq!(maxbin_rs::emanager::get_prob_abund(val, 0.0), 0.0);
    }

    /// get_prob_abund: positive inputs return non-negative finite values.
    #[test]
    fn prob_abund_non_negative(
        curr in 0.001f64..5000.0,
        lambda in 0.001f64..5000.0
    ) {
        let p = maxbin_rs::emanager::get_prob_abund(curr, lambda);
        prop_assert!(p >= 0.0, "prob_abund({curr}, {lambda}) = {p} < 0");
        prop_assert!(p.is_finite(), "prob_abund({curr}, {lambda}) = {p} is not finite");
    }

    /// get_prob_abund peaks near curr == lambda (Poisson mode property).
    #[test]
    fn prob_abund_peaks_near_lambda(lambda in 2.0f64..50.0) {
        let at_peak = maxbin_rs::emanager::get_prob_abund(lambda, lambda);
        let above = maxbin_rs::emanager::get_prob_abund(lambda * 2.0, lambda);
        let below = if lambda > 1.0 {
            maxbin_rs::emanager::get_prob_abund(0.1, lambda)
        } else {
            0.0
        };

        prop_assert!(
            at_peak >= above,
            "prob_abund({lambda}, {lambda})={at_peak} < prob_abund({}, {lambda})={above}",
            lambda * 2.0
        );
        if lambda > 1.0 {
            prop_assert!(
                at_peak >= below,
                "prob_abund({lambda}, {lambda})={at_peak} < prob_abund(0.1, {lambda})={below}"
            );
        }
    }
}
