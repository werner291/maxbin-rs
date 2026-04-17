/// Exact reimplementation of MaxBin2's EucDist and SpearmanDist.
///
/// Both depend on KmerMap, Profiler, and quickSort.
/// Only EucDist and SpearmanDist are used by EManager.
use crate::kmer_map::KmerMap;
use crate::profiler::Profiler;
use crate::quicksort;

/// Shared state for distance computations, mirroring AbstractDist.
/// Matches AbstractDist.cpp:3-7 / AbstractDist.cpp:26-34 (constructor / init()):
/// kmap = new kmerMap(kmerlen, true), entrynum = kmap->getEntryNum(), normalize = false.
pub struct DistanceContext {
    kmerlen: usize,
    entrynum: usize,
    kmap: KmerMap,
    normalize: bool,
}

impl DistanceContext {
    pub fn new(kmerlen: usize) -> Self {
        // Matches AbstractDist.cpp:28: kmap = new kmerMap(kmerlen, true)  (always symmetric)
        let kmap = KmerMap::new(kmerlen, true);
        let entrynum = kmap.get_entry_num();
        Self {
            kmerlen,
            entrynum,
            kmap,
            normalize: false,
        }
    }

    /// Matches AbstractDist.cpp:36-39 (setNormalization())
    pub fn set_normalization(&mut self, normalize: bool) {
        self.normalize = normalize;
    }

    pub fn entry_num(&self) -> usize {
        self.entrynum
    }

    pub fn kmap(&self) -> &KmerMap {
        &self.kmap
    }

    /// Compute rank array from a profile. Sorts profile descending, then
    /// assigns rank[tmprank[i]] = i. Note: this modifies the profile in-place
    /// (via quickSort).
    /// Matches AbstractDist.cpp:54-68 (computeRank())
    fn compute_rank(&self, profile: &mut [f64]) -> Vec<i32> {
        // Matches AbstractDist.cpp:57-61: initialize tmprank[i] = i, then sort with quickSort
        let mut tmprank: Vec<i32> = (0..self.entrynum as i32).collect();
        quicksort::sort_descending(profile, Some(&mut tmprank));
        // Matches AbstractDist.cpp:63-66: rank[tmprank[i]] = i
        let mut rank = vec![0i32; self.entrynum];
        for i in 0..self.entrynum {
            rank[tmprank[i] as usize] = i as i32;
        }
        rank
    }
}

// --- Euclidean Distance ---

/// Compute Euclidean distance between two DNA sequences.
/// Matches EucDist.cpp:4-35 (getDist(const char*, const char*)):
/// build Profiler for each sequence, then compute sum of squared differences.
pub fn euc_dist_seq(ctx: &DistanceContext, seq1: &[u8], seq2: &[u8]) -> f64 {
    let pro1 = Profiler::new(ctx.kmerlen, seq1, &ctx.kmap);
    let pro2 = Profiler::new(ctx.kmerlen, seq2, &ctx.kmap);
    euc_dist_profiles(pro1.get_profile(), pro2.get_profile())
}

/// Compute Euclidean distance between two profile arrays.
/// Matches EucDist.cpp:37-60 (getDist(double*, double*)) and EucDist.cpp:24-30:
/// sum of (p1[i] - p2[i])^2, then sqrt.
pub fn euc_dist_profiles(pro1: &[f64], pro2: &[f64]) -> f64 {
    let mut f = 0.0;
    // Matches EucDist.cpp:26-29 / EucDist.cpp:51-54: f += pow(p1[i]-p2[i], 2)
    for i in 0..pro1.len() {
        f += (pro1[i] - pro2[i]).powi(2);
    }
    // Matches EucDist.cpp:30 / EucDist.cpp:55: f = sqrt(f)
    f.sqrt()
}

// --- Spearman Footrule Distance ---

/// Compute Spearman footrule distance between two DNA sequences.
/// Matches SpearmanDist.cpp:3-51 (getDist(const char*, const char*)):
/// build profiles, compute ranks via quickSort, sum absolute rank differences.
pub fn spearman_dist_seq(ctx: &DistanceContext, seq1: &[u8], seq2: &[u8]) -> f64 {
    let pro1 = Profiler::new(ctx.kmerlen, seq1, &ctx.kmap);
    let pro2 = Profiler::new(ctx.kmerlen, seq2, &ctx.kmap);

    // computeRank modifies the profile via quickSort, so we need mutable copies
    // Matches SpearmanDist.cpp:14 / SpearmanDist.cpp:24: computeRank(profile, rank)
    let mut p1 = pro1.get_profile().to_vec();
    let mut p2 = pro2.get_profile().to_vec();

    let rank1 = ctx.compute_rank(&mut p1);
    let rank2 = ctx.compute_rank(&mut p2);

    // Matches SpearmanDist.cpp:30-41: sum |rank1[i] - rank2[i]|
    let mut j = 0.0f64;
    for i in 0..ctx.entrynum {
        j += (rank1[i] - rank2[i]).unsigned_abs() as f64;
    }

    // Matches SpearmanDist.cpp:43-46: normalize by entrynum * (entrynum + 1)
    if ctx.normalize {
        j /= (ctx.entrynum * (ctx.entrynum + 1)) as f64;
    }

    j
}

/// Compute Spearman footrule distance between two profile arrays.
/// Makes copies of the input arrays since computeRank modifies them via quickSort.
/// Matches SpearmanDist.cpp:53-93 (getDist(double*, double*))
pub fn spearman_dist_profiles(ctx: &DistanceContext, pro1: &[f64], pro2: &[f64]) -> f64 {
    // Matches SpearmanDist.cpp:58-62: copy input profiles before modifying via computeRank
    let mut p1 = pro1.to_vec();
    let mut p2 = pro2.to_vec();

    // Matches SpearmanDist.cpp:66-68: compute ranks for both profiles
    let rank1 = ctx.compute_rank(&mut p1);
    let rank2 = ctx.compute_rank(&mut p2);

    // Matches SpearmanDist.cpp:72-82: sum |rank1[i] - rank2[i]|
    let mut j = 0.0f64;
    for i in 0..ctx.entrynum {
        j += (rank1[i] - rank2[i]).unsigned_abs() as f64;
    }

    // Matches SpearmanDist.cpp:84-87: normalize by entrynum * (entrynum + 1)
    if ctx.normalize {
        j /= (ctx.entrynum * (ctx.entrynum + 1)) as f64;
    }

    j
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euc_dist_identical_sequences() {
        let ctx = DistanceContext::new(4);
        let seq = b"ACGTACGTACGTACGT";
        let d = euc_dist_seq(&ctx, seq, seq);
        assert!(
            d.abs() < 1e-12,
            "identical sequences should have distance 0, got {d}"
        );
    }

    #[test]
    fn euc_dist_different_sequences() {
        let ctx = DistanceContext::new(4);
        let d = euc_dist_seq(&ctx, b"AAAAAAAAAAAAAAAA", b"CCCCCCCCCCCCCCCC");
        assert!(d > 0.0);
    }

    #[test]
    fn spearman_dist_identical_sequences() {
        let ctx = DistanceContext::new(4);
        let seq = b"ACGTACGTACGTACGT";
        let d = spearman_dist_seq(&ctx, seq, seq);
        assert!(
            d.abs() < 1e-12,
            "identical sequences should have distance 0, got {d}"
        );
    }

    #[test]
    fn euc_dist_from_profiles() {
        let p1 = vec![0.1, 0.2, 0.3, 0.4];
        let p2 = vec![0.4, 0.3, 0.2, 0.1];
        let d = euc_dist_profiles(&p1, &p2);
        // (0.1-0.4)^2 + (0.2-0.3)^2 + (0.3-0.2)^2 + (0.4-0.1)^2 = 0.2
        let expected = 0.2f64.sqrt();
        assert!((d - expected).abs() < 1e-12, "got {d}, expected {expected}");
    }

    #[test]
    fn spearman_normalized() {
        let mut ctx = DistanceContext::new(4);
        ctx.set_normalization(true);
        let d = spearman_dist_seq(&ctx, b"AAAAAAAAAAAAAAAA", b"CCCCCCCCCCCCCCCC");
        // normalized by entrynum * (entrynum + 1)
        assert!(d >= 0.0);
        assert!(d <= 1.0, "normalized distance should be <= 1, got {d}");
    }

    // --- FFI equivalence tests ---

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
        let cpp = crate::original_ffi::OriginalEucDist::new(4);

        for (seq1, seq2) in &test_pairs {
            let rust_d = euc_dist_seq(&ctx, seq1.as_bytes(), seq2.as_bytes());
            let cpp_d = cpp.get_dist_seq(seq1, seq2);

            let diff = (rust_d - cpp_d).abs();
            assert!(
                diff < 1e-12,
                "EucDist seq mismatch for '{seq1}' vs '{seq2}': rust={rust_d} cpp={cpp_d}"
            );
        }
    }

    #[test]
    fn ffi_euc_dist_profile_equivalence() {
        // Build profiles from known sequences, then compare profile-based distance
        let cpp = crate::original_ffi::OriginalEucDist::new(4);

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

                let diff = (rust_d - cpp_d).abs();
                assert!(
                    diff < 1e-12,
                    "EucDist profile mismatch: rust={rust_d} cpp={cpp_d}"
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
        let cpp = crate::original_ffi::OriginalSpearmanDist::new(4);

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
        let cpp = crate::original_ffi::OriginalSpearmanDist::new(4);

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
        let cpp = crate::original_ffi::OriginalSpearmanDist::new(4);
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
}
