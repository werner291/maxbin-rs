//! Exact reimplementation of MaxBin2's NormalDistribution.
//!
//! Uses pi = 3.1415926 (NOT std::f64::consts::PI) to match the original.

/// The truncated pi value used in the original MaxBin2 code.
/// Matches NormalDistribution.cpp:5: `pi = 3.1415926;`
#[allow(clippy::approx_constant)]
const MAXBIN2_PI: f64 = 3.1415926;

pub struct NormalDistribution {
    mean: f64,
    std: f64,
}

impl NormalDistribution {
    /// Matches NormalDistribution.cpp:3-8 (constructor): store mean and std.
    pub fn new(mean: f64, std: f64) -> Self {
        Self { mean, std }
    }

    /// Compute the probability density function at `input`.
    /// Matches NormalDistribution.cpp:10-16 (prob()):
    /// `(1 / (std * sqrt(2 * pi))) * exp(-0.5 * pow((input - mean) / std, 2))`
    pub fn prob(&self, input: f64) -> f64 {
        let exponent = -0.5 * ((input - self.mean) / self.std).powi(2);
        (1.0 / (self.std * (2.0 * MAXBIN2_PI).sqrt())) * exponent.exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_normal_at_zero() {
        let nd = NormalDistribution::new(0.0, 1.0);
        let p = nd.prob(0.0);
        // 1 / sqrt(2 * 3.1415926) = 0.3989422816...
        let expected = 1.0 / (2.0 * MAXBIN2_PI).sqrt();
        assert!((p - expected).abs() < 1e-15);
    }

    #[test]
    fn nonzero_mean() {
        let nd = NormalDistribution::new(5.0, 2.0);
        let p = nd.prob(5.0);
        // At the mean, PDF = 1 / (std * sqrt(2*pi))
        let expected = 1.0 / (2.0 * (2.0 * MAXBIN2_PI).sqrt());
        assert!((p - expected).abs() < 1e-15);
    }

    #[test]
    fn symmetry() {
        let nd = NormalDistribution::new(0.0, 1.0);
        let p1 = nd.prob(1.0);
        let p2 = nd.prob(-1.0);
        assert!((p1 - p2).abs() < 1e-15);
    }

    #[test]
    fn uses_truncated_pi() {
        // Verify we use 3.1415926, not std::f64::consts::PI
        let nd = NormalDistribution::new(0.0, 1.0);
        let our_val = nd.prob(0.0);
        let std_pi_val = 1.0 / (2.0 * std::f64::consts::PI).sqrt();
        // They should differ slightly
        assert!(our_val != std_pi_val);
    }

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
            let cpp_nd = crate::original_ffi::OriginalNormalDistribution::new(mean, std);

            let rust_p = rust_nd.prob(input);
            let cpp_p = cpp_nd.prob(input);

            let diff = (rust_p - cpp_p).abs();
            assert!(
                diff < 1e-15,
                "NormalDist mismatch for mean={mean}, std={std}, input={input}: \
                 rust={rust_p} cpp={cpp_p} diff={diff}"
            );
        }
    }
}
