//! FFI equivalence test for the normal distribution PDF.
//!
//! Compares `maxbin_rs::normal_distribution::NormalDistribution` against the
//! original MaxBin2 C++ NormalDistribution. Previously an in-source test in
//! `normal_distribution.rs`.
use maxbin_rs::normal_distribution::NormalDistribution;
use maxbin_rs_equivalence::original_ffi::OriginalNormalDistribution;

#[test]
fn ffi_equivalence() {
    let test_cases = [
        (0.0, 1.0, 0.0),
        (0.0, 1.0, 1.0),
        (0.0, 1.0, -1.0),
        (5.0, 2.0, 5.0),
        (5.0, 2.0, 3.0),
        (5.0, 2.0, 7.0),
        (10.0, 0.5, 10.0),
        (10.0, 0.5, 9.5),
        (-3.0, 1.5, -3.0),
        (-3.0, 1.5, 0.0),
        (0.0, 1.0, 3.5),
        (100.0, 10.0, 95.0),
    ];

    for &(mean, std, input) in &test_cases {
        let rust_nd = NormalDistribution::new(mean, std);
        let cpp_nd = OriginalNormalDistribution::new(mean, std);

        let rust_p = rust_nd.prob(input);
        let cpp_p = cpp_nd.prob(input);

        assert_eq!(
            rust_p.to_bits(),
            cpp_p.to_bits(),
            "NormalDist not bit-identical for mean={mean}, std={std}, input={input}: \
             rust={rust_p:e} cpp={cpp_p:e}"
        );
    }
}
