/// Exact reimplementation of MaxBin2's kmerMap.
///
/// Encodes bases as A=0, T=1, C=2, G=3.
/// In symmetric mode, reverse-complement k-mers map to the same index.

/// Encode a single base character to its numeric value.
/// Returns `None` for non-ACGT characters.
/// Matches kmerMap.cpp:219-244 (calNum() switch cases)
fn base_to_num(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'T' | b't' => Some(1),
        b'C' | b'c' => Some(2),
        b'G' | b'g' => Some(3),
        _ => None,
    }
}

/// Convert a numeric index back to a k-mer string using the MaxBin2 encoding.
/// Matches kmerMap.cpp:255-285 (numToString())
fn num_to_string(mut num: u32, kmerlen: usize) -> Vec<u8> {
    let mut result = vec![b'A'; kmerlen];
    // Matches kmerMap.cpp:259-284: iterate from last position to first, extract base by mod 4
    for j in (0..kmerlen).rev() {
        let k = num % 4;
        num /= 4;
        result[j] = match k {
            0 => b'A',
            1 => b'T',
            2 => b'C',
            3 => b'G',
            _ => unreachable!(),
        };
    }
    result
}

/// Compute the reverse complement of a sequence.
/// Matches kmerMap.cpp:176-207 (revSeq())
fn rev_complement(seq: &[u8]) -> Vec<u8> {
    let len = seq.len();
    let mut result = vec![0u8; len];
    // Matches kmerMap.cpp:181-206: complement each base and store in reverse order
    for i in 0..len {
        result[len - i - 1] = match seq[i] {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other, // shouldn't happen for valid input
        };
    }
    result
}

/// Compute the numeric value of a k-mer string.
/// Returns `Some(value)` on success, `None` if any character is not ACGT.
/// This matches `calNum` which sets `err=true` on success, `err=false` on
/// encountering a non-ACGT character.
/// Matches kmerMap.cpp:209-253 (calNum())
fn cal_num(kmer: &[u8]) -> Option<u32> {
    let mut num: u64 = 0;
    // Matches kmerMap.cpp:216-250: accumulate base values, return None on unknown base
    for &b in kmer {
        num *= 4;
        match base_to_num(b) {
            Some(v) => num += v as u64,
            None => return None,
        }
    }
    Some(num as u32)
}

pub struct KmerMap {
    kmerlen: usize,
    symmetric: bool,
    kmer_map: Vec<i32>,
    entry_num: usize, // 4^kmerlen (total raw k-mer count)
}

impl KmerMap {
    /// Matches kmerMap.cpp:3-9 (constructor): compute entry_num and build mapping table.
    pub fn new(kmerlen: usize, symmetric: bool) -> Self {
        // Matches kmerMap.cpp:6: entry_num = pow(4, len_mer)
        let entry_num = 4usize.pow(kmerlen as u32);
        let kmer_map = if symmetric {
            build_symmetric_table(kmerlen, entry_num)
        } else {
            // Matches kmerMap.cpp:163-173: non-symmetric is identity mapping
            (0..entry_num as i32).collect()
        };
        Self {
            kmerlen,
            symmetric,
            kmer_map,
            entry_num,
        }
    }

    /// Number of distinct entries in the mapping (collapsed count for symmetric).
    /// Matches kmerMap.cpp:80-97 (getEntryNum())
    pub fn get_entry_num(&self) -> usize {
        if self.symmetric {
            // Matches kmerMap.cpp:83-91: even vs odd kmer length formula
            if self.kmerlen % 2 == 0 {
                (self.entry_num + 4usize.pow((self.kmerlen / 2) as u32)) / 2
            } else {
                self.entry_num / 2
            }
        } else {
            self.entry_num
        }
    }

    /// Map a k-mer string to its index. Returns -1 for invalid k-mers.
    /// Matches kmerMap.cpp:24-36 (getMapping(char*))
    pub fn get_mapping(&self, kmer: &[u8]) -> i32 {
        match cal_num(kmer) {
            // Matches kmerMap.cpp:27-31: return kmer_map[i] if err is true (valid)
            Some(i) => self.kmer_map[i as usize],
            // Matches kmerMap.cpp:33-35: return -1 if err is false (invalid base)
            None => -1,
        }
    }

    /// Map the reverse complement of a k-mer string to its index.
    /// Matches kmerMap.cpp:38-55 (getReverseMapping(char*))
    pub fn get_reverse_mapping_str(&self, kmer: &[u8]) -> i32 {
        let rev = rev_complement(kmer);
        match cal_num(&rev) {
            Some(i) => self.kmer_map[i as usize],
            None => -1,
        }
    }

    /// Map the reverse complement of a k-mer given by numeric index.
    /// Matches kmerMap.cpp:57-78 (getReverseMapping(int index))
    pub fn get_reverse_mapping_idx(&self, index: u32) -> i32 {
        let kmer = num_to_string(index, self.kmerlen);
        let rev = rev_complement(&kmer);
        match cal_num(&rev) {
            Some(i) => self.kmer_map[i as usize],
            None => -1,
        }
    }

    /// Get the full mapping table (for testing).
    pub fn get_full_table(&self) -> &[i32] {
        &self.kmer_map
    }
}

/// Build the symmetric k-mer mapping table.
/// Matches kmerMap.cpp:99-174 (buildKmerMappingTable() symmetric branch)
fn build_symmetric_table(kmerlen: usize, entry_num: usize) -> Vec<i32> {
    // Matches kmerMap.cpp:103-112: compute size of kmer_small/kmer_large arrays
    let max_pairs = if kmerlen % 2 == 0 {
        (entry_num + 4usize.pow((kmerlen / 2) as u32)) / 2
    } else {
        entry_num / 2
    };

    let mut kmer_small = vec![0i32; max_pairs];
    let mut kmer_large = vec![0i32; max_pairs];
    let mut kmer_num = 0usize;

    // Matches kmerMap.cpp:122-151: enumerate all k-mers, group canonical pairs
    for i in 0..entry_num {
        let kmer = num_to_string(i as u32, kmerlen);
        let rev = rev_complement(&kmer);
        let j = cal_num(&rev).unwrap() as usize;

        // Matches kmerMap.cpp:127-136: assign k (smaller) and p (larger) of the pair
        let (k, p) = if i <= j { (i, j) } else { (j, i) };

        // Matches kmerMap.cpp:137-144: linear scan to check if pair already recorded
        let mut found = false;
        for idx in 0..kmer_num {
            if kmer_small[idx] == k as i32 {
                found = true;
                break;
            }
        }
        // Matches kmerMap.cpp:145-150: record the new pair
        if !found {
            kmer_small[kmer_num] = k as i32;
            kmer_large[kmer_num] = p as i32;
            kmer_num += 1;
        }
    }

    // Matches kmerMap.cpp:152-159: build final kmer_map from the pair lists
    let mut kmer_map = vec![0i32; entry_num];
    for i in 0..kmer_num {
        kmer_map[kmer_small[i] as usize] = i as i32;
        kmer_map[kmer_large[i] as usize] = i as i32;
    }

    kmer_map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_symmetric_identity() {
        let km = KmerMap::new(4, false);
        assert_eq!(km.get_entry_num(), 256);
        // Identity mapping: each index maps to itself
        for i in 0..256 {
            assert_eq!(km.get_full_table()[i], i as i32);
        }
    }

    #[test]
    fn symmetric_entry_count_kmer4() {
        let km = KmerMap::new(4, true);
        // (256 + 16) / 2 = 136
        assert_eq!(km.get_entry_num(), 136);
    }

    #[test]
    fn symmetric_reverse_complement_collapse() {
        let km = KmerMap::new(4, true);
        // AAAT and its reverse complement ATTT should map to the same index
        assert_eq!(km.get_mapping(b"AAAT"), km.get_mapping(b"ATTT"));
        // ACGT and its reverse complement ACGT (palindrome) should be valid
        let idx = km.get_mapping(b"ACGT");
        assert!(idx >= 0);
    }

    #[test]
    fn invalid_kmer_returns_negative() {
        let km = KmerMap::new(4, true);
        assert_eq!(km.get_mapping(b"AAAX"), -1);
        assert_eq!(km.get_mapping(b"NACT"), -1);
    }

    #[test]
    fn base_encoding() {
        // Verify A=0, T=1, C=2, G=3
        let km = KmerMap::new(1, false);
        assert_eq!(km.get_mapping(b"A"), 0);
        assert_eq!(km.get_mapping(b"T"), 1);
        assert_eq!(km.get_mapping(b"C"), 2);
        assert_eq!(km.get_mapping(b"G"), 3);
    }

    #[test]
    fn reverse_mapping_by_index() {
        let km = KmerMap::new(4, true);
        // getReverseMapping(0) should give the mapping for the reverse complement of AAAA
        // revcomp(AAAA) = TTTT. calNum(TTTT) = 1+4+16+64 = 85. kmer_map[85].
        let rev = km.get_reverse_mapping_idx(0);
        assert_eq!(rev, km.get_mapping(b"TTTT"));
    }

    #[test]
    fn ffi_equivalence_kmer4_symmetric() {
        let rust_km = KmerMap::new(4, true);
        let cpp_km = crate::original_ffi::OriginalKmerMap::new(4, true);

        assert_eq!(rust_km.get_entry_num() as i32, cpp_km.entry_num());

        let cpp_table = cpp_km.get_full_table();
        let rust_table = rust_km.get_full_table();

        assert_eq!(rust_table.len(), cpp_table.len());
        for i in 0..rust_table.len() {
            assert_eq!(
                rust_table[i], cpp_table[i],
                "mapping mismatch at index {i}"
            );
        }
    }

    #[test]
    fn ffi_equivalence_kmer4_nonsymmetric() {
        let rust_km = KmerMap::new(4, false);
        let cpp_km = crate::original_ffi::OriginalKmerMap::new(4, false);

        assert_eq!(rust_km.get_entry_num() as i32, cpp_km.entry_num());

        let cpp_table = cpp_km.get_full_table();
        let rust_table = rust_km.get_full_table();

        for i in 0..rust_table.len() {
            assert_eq!(rust_table[i], cpp_table[i]);
        }
    }

    #[test]
    fn ffi_equivalence_kmer3_symmetric() {
        let rust_km = KmerMap::new(3, true);
        let cpp_km = crate::original_ffi::OriginalKmerMap::new(3, true);

        assert_eq!(rust_km.get_entry_num() as i32, cpp_km.entry_num());

        let cpp_table = cpp_km.get_full_table();
        let rust_table = rust_km.get_full_table();

        for i in 0..rust_table.len() {
            assert_eq!(
                rust_table[i], cpp_table[i],
                "kmer3 symmetric mapping mismatch at index {i}"
            );
        }
    }
}
