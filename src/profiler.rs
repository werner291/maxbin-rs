/// Exact reimplementation of MaxBin2's Profiler.
///
/// Computes tetranucleotide frequency profiles for DNA sequences.
/// Normalizes counts by dividing by total valid k-mers found.
use crate::kmer_map::KmerMap;

pub struct Profiler {
    profile: Vec<f64>,
    percent_n: f32,
    entry_num: usize,
}

impl Profiler {
    /// Create a new Profiler for the given sequence using the provided KmerMap.
    /// Matches Profiler.cpp:3-13 (constructor): allocate profile array and call computeProfile.
    pub fn new(kmerlen: usize, seq: &[u8], kmap: &KmerMap) -> Self {
        let entry_num = kmap.get_entry_num();
        let mut profile = vec![0.0f64; entry_num];
        let percent_n = compute_profile(kmerlen, seq, kmap, &mut profile);
        Self {
            profile,
            percent_n,
            entry_num,
        }
    }

    pub fn get_profile(&self) -> &[f64] {
        &self.profile
    }

    pub fn get_percent_n(&self) -> f32 {
        self.percent_n
    }

    /// Reset profile to zero with total_weight = 0 (for addProfile/calcProfile workflow).
    /// Matches Profiler.cpp:36-40 (reset()): zero the profile array and set total_weight = 0.
    pub fn reset(&mut self) {
        self.profile.fill(0.0);
    }

    /// Add another profile weighted by `weight`.
    /// Matches Profiler.cpp:42-50 (addProfile()): profile[i] += other->profile[i] * weight
    pub fn add_profile(&mut self, other: &Profiler, weight: f64) {
        for i in 0..self.entry_num {
            self.profile[i] += other.profile[i] * weight;
        }
    }

    /// Normalize accumulated profiles by total weight.
    /// Matches Profiler.cpp:52-71 (calcProfile()): divide by total_weight; zero out if weight==0.
    pub fn calc_profile(&mut self, total_weight: f64) {
        if total_weight == 0.0 {
            self.profile.fill(0.0);
        } else {
            for v in self.profile.iter_mut() {
                *v /= total_weight;
            }
        }
    }
}

/// Matches Profiler.cpp:73-119 (computeProfile())
fn compute_profile(kmerlen: usize, seq: &[u8], kmap: &KmerMap, profile: &mut [f64]) -> f32 {
    let len = seq.len();
    // Matches Profiler.cpp:79: loop bound is len - kmerlen + 1.
    //
    // KNOWN C++ BUG: when len < kmerlen, the C++ computes window_count =
    // (len - kmerlen + 1) as signed int. The loop doesn't execute (i < negative),
    // N stays 0, total stays 0, then percent_N = (float)0 / (float)window_count.
    // For len=3 kmerlen=4: window_count=0, so 0/0 = NaN.
    // For len=1 kmerlen=4: window_count=-2, so 0/-2 = -0.0.
    // We reproduce this bug-for-bug with the same signed integer arithmetic.
    // In practice this never fires — filter_contigs() enforces min_contig_length
    // (default 1000) before profiling.
    if len < kmerlen {
        let window_count = len as i32 - kmerlen as i32 + 1; // signed, matches C++ int
        return 0.0f32 / window_count as f32;
    }

    let mut n_count = 0;
    let mut total = 0;
    let window_count = len - kmerlen + 1;

    // Matches Profiler.cpp:79-95: slide window, count N at position i, map k-mer to index
    for i in 0..window_count {
        // Matches Profiler.cpp:81-83: count 'N' bases at the start of each window
        if seq[i] == b'N' {
            n_count += 1;
        }
        let kmer = &seq[i..i + kmerlen];
        let j = kmap.get_mapping(kmer);
        // Matches Profiler.cpp:86-94: only count valid k-mers (j != -1)
        if j != -1 {
            total += 1;
            profile[j as usize] += 1.0;
        }
    }

    // Matches Profiler.cpp:96-97: percent_N = N_count / window_count
    let percent_n = n_count as f32 / window_count as f32;
    let entry_num = kmap.get_entry_num();

    // Matches Profiler.cpp:98-118: normalize; zero out if all-N or no valid k-mers
    if percent_n == 1.0 || total == 0 {
        for i in 0..entry_num {
            profile[i] = 0.0;
        }
    } else {
        // Matches Profiler.cpp:113-117: divide each count by total valid k-mers
        for i in 0..entry_num {
            profile[i] /= total as f64;
        }
    }

    percent_n
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_profile() {
        let kmap = KmerMap::new(4, true);
        let seq = b"ACGTACGTACGT";
        let prof = Profiler::new(4, seq, &kmap);
        let profile = prof.get_profile();

        // Sum of all profile values should be ~1.0 (normalized)
        let sum: f64 = profile.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "profile sum = {sum}");
    }

    #[test]
    fn all_n_sequence() {
        let kmap = KmerMap::new(4, true);
        let seq = b"NNNNNNNNNN";
        let prof = Profiler::new(4, seq, &kmap);
        assert_eq!(prof.get_percent_n(), 1.0);
        // All zeros
        assert!(prof.get_profile().iter().all(|&v| v == 0.0));
    }

    #[test]
    fn short_sequence() {
        let kmap = KmerMap::new(4, true);
        let seq = b"ACG"; // too short for kmerlen=4
        let prof = Profiler::new(4, seq, &kmap);
        assert!(prof.get_profile().iter().all(|&v| v == 0.0));
    }

    #[test]
    fn ffi_equivalence() {
        let test_seqs = [
            "ACGTACGTACGTACGT",
            "AAAAAAAAAAAA",
            "ACGTTTTTGGGGCCCC",
            "ATCGATCGATCGATCGATCG",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        ];

        let rust_kmap = KmerMap::new(4, true);
        let cpp_kmap = crate::original_ffi::OriginalKmerMap::new(4, true);
        let entry_num = rust_kmap.get_entry_num();

        for seq in &test_seqs {
            let rust_prof = Profiler::new(4, seq.as_bytes(), &rust_kmap);
            let cpp_prof = crate::original_ffi::OriginalProfiler::new(4, seq, &cpp_kmap);

            let rust_profile = rust_prof.get_profile();
            let cpp_profile = cpp_prof.get_profile(entry_num as i32);

            assert_eq!(rust_profile.len(), cpp_profile.len());
            for i in 0..rust_profile.len() {
                let diff = (rust_profile[i] - cpp_profile[i]).abs();
                assert!(
                    diff < 1e-12,
                    "profile mismatch at index {i} for seq '{seq}': rust={} cpp={}",
                    rust_profile[i],
                    cpp_profile[i]
                );
            }

            let rust_pn = rust_prof.get_percent_n();
            let cpp_pn = cpp_prof.get_percent_n();
            assert!(
                (rust_pn - cpp_pn).abs() < 1e-6,
                "percent_n mismatch for seq '{seq}': rust={rust_pn} cpp={cpp_pn}"
            );
        }
    }
}
