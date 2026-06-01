//! Component-level FFI equivalence tests for distance metrics.
//!
//! These compare the Rust distance functions in `maxbin_core::distance` against
//! the original MaxBin2 C++ EucDist / SpearmanDist via FFI. They previously
//! lived as in-source `#[cfg(test)]` tests inside `distance.rs`; they moved
//! here when the C++ FFI was extracted from the core crate.
use maxbin_core::distance::{
    DistanceContext, euc_dist_profiles, euc_dist_seq, spearman_dist_profiles, spearman_dist_seq,
};
use maxbin_core::kmer_map::KmerMap;
use maxbin_core::profiler::Profiler;
use maxbin_rs_equivalence::original_ffi::{OriginalEucDist, OriginalSpearmanDist};

#[test]
fn ffi_euc_dist_seq_equivalence() {
    let test_pairs = [
        ("ACGTACGTACGTACGT", "ACGTACGTACGTACGT"),
        ("AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC"),
        ("ATCGATCGATCGATCG", "TAGCTAGCTAGCTAGC"),
        ("ACGTACGTACGTACGTACGTACGT", "TTTTAAAACCCCGGGG"),
        ("AACCGGTTAACCGGTT", "GGTTAACCGGTTAACC"),
    ];

    let ctx = DistanceContext::new(4);
    let cpp = OriginalEucDist::new(4);

    for (seq1, seq2) in &test_pairs {
        let rust_d = euc_dist_seq(&ctx, seq1.as_bytes(), seq2.as_bytes());
        let cpp_d = cpp.get_dist_seq(seq1, seq2);

        assert_eq!(
            rust_d.to_bits(),
            cpp_d.to_bits(),
            "EucDist seq not bit-identical for '{seq1}' vs '{seq2}': rust={rust_d:e} cpp={cpp_d:e}"
        );
    }
}

#[test]
fn ffi_euc_dist_profile_equivalence() {
    // Build profiles from known sequences, then compare profile-based distance
    let cpp = OriginalEucDist::new(4);

    let seqs = ["ACGTACGTACGTACGT", "AAAAAAAAAAAAAAAA", "ATCGATCGATCGATCG"];

    let rust_kmap = KmerMap::new(4, true);
    let profiles: Vec<Vec<f64>> = seqs
        .iter()
        .map(|s| {
            Profiler::new(4, s.as_bytes(), &rust_kmap)
                .get_profile()
                .to_vec()
        })
        .collect();

    for i in 0..profiles.len() {
        for j in (i + 1)..profiles.len() {
            let rust_d = euc_dist_profiles(&profiles[i], &profiles[j]);
            let cpp_d = cpp.get_dist_profile(&profiles[i], &profiles[j]);

            assert_eq!(
                rust_d.to_bits(),
                cpp_d.to_bits(),
                "EucDist profile not bit-identical: rust={rust_d:e} cpp={cpp_d:e}"
            );
        }
    }
}

#[test]
fn ffi_spearman_dist_seq_equivalence() {
    let test_pairs = [
        ("ACGTACGTACGTACGT", "ACGTACGTACGTACGT"),
        ("AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC"),
        ("ATCGATCGATCGATCG", "TAGCTAGCTAGCTAGC"),
        ("ACGTACGTACGTACGTACGTACGT", "TTTTAAAACCCCGGGG"),
        ("AACCGGTTAACCGGTT", "GGTTAACCGGTTAACC"),
    ];

    let ctx = DistanceContext::new(4);
    let cpp = OriginalSpearmanDist::new(4);

    for (seq1, seq2) in &test_pairs {
        let rust_d = spearman_dist_seq(&ctx, seq1.as_bytes(), seq2.as_bytes());
        let cpp_d = cpp.get_dist_seq(seq1, seq2);

        let diff = (rust_d - cpp_d).abs();
        assert!(
            diff < 1e-10,
            "SpearmanDist seq mismatch for '{seq1}' vs '{seq2}': rust={rust_d} cpp={cpp_d}"
        );
    }
}

#[test]
fn ffi_spearman_dist_profile_equivalence() {
    let ctx = DistanceContext::new(4);
    let cpp = OriginalSpearmanDist::new(4);

    let rust_kmap = KmerMap::new(4, true);
    let seqs = ["ACGTACGTACGTACGT", "AAAAAAAAAAAAAAAA", "ATCGATCGATCGATCG"];

    let profiles: Vec<Vec<f64>> = seqs
        .iter()
        .map(|s| {
            Profiler::new(4, s.as_bytes(), &rust_kmap)
                .get_profile()
                .to_vec()
        })
        .collect();

    for i in 0..profiles.len() {
        for j in (i + 1)..profiles.len() {
            let rust_d = spearman_dist_profiles(&ctx, &profiles[i], &profiles[j]);
            let cpp_d = cpp.get_dist_profile(&profiles[i], &profiles[j]);

            let diff = (rust_d - cpp_d).abs();
            assert!(
                diff < 1e-10,
                "SpearmanDist profile mismatch: rust={rust_d} cpp={cpp_d}"
            );
        }
    }
}

#[test]
fn ffi_spearman_normalized_equivalence() {
    let mut ctx = DistanceContext::new(4);
    ctx.set_normalization(true);
    let cpp = OriginalSpearmanDist::new(4);
    cpp.set_normalization(true);

    let seq1 = "ACGTACGTACGTACGT";
    let seq2 = "AAAAAAAAAAAAAAAA";

    let rust_d = spearman_dist_seq(&ctx, seq1.as_bytes(), seq2.as_bytes());
    let cpp_d = cpp.get_dist_seq(seq1, seq2);

    let diff = (rust_d - cpp_d).abs();
    assert!(
        diff < 1e-10,
        "SpearmanDist normalized mismatch: rust={rust_d} cpp={cpp_d}"
    );
}
