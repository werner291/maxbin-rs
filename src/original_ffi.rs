/// FFI bindings to the original MaxBin2 C++ code.
/// Used only for equivalence testing — not compiled into release builds.
use std::ffi::{CStr, CString};
use std::os::raw::c_char;

// --- fastaReader ---

unsafe extern "C" {
    fn fasta_reader_new(path: *const c_char) -> *mut std::ffi::c_void;
    fn fasta_reader_free(reader: *mut std::ffi::c_void);
    fn fasta_reader_get_num(reader: *mut std::ffi::c_void) -> u32;
    fn fasta_reader_get_header(reader: *mut std::ffi::c_void, i: u32) -> *const c_char;
    fn fasta_reader_get_seq(reader: *mut std::ffi::c_void, i: u32) -> *const c_char;
    fn fasta_reader_get_seq_len(reader: *mut std::ffi::c_void, i: u32) -> u32;
}

pub struct OriginalFastaReader {
    ptr: *mut std::ffi::c_void,
}

impl OriginalFastaReader {
    pub fn new(path: &std::path::Path) -> Self {
        let c_path = CString::new(path.to_str().unwrap()).unwrap();
        let ptr = unsafe { fasta_reader_new(c_path.as_ptr()) };
        Self { ptr }
    }

    pub fn num_records(&self) -> u32 {
        unsafe { fasta_reader_get_num(self.ptr) }
    }

    pub fn header(&self, i: u32) -> String {
        unsafe {
            let ptr = fasta_reader_get_header(self.ptr, i);
            CStr::from_ptr(ptr).to_string_lossy().into_owned()
        }
    }

    pub fn seq(&self, i: u32) -> Vec<u8> {
        unsafe {
            let ptr = fasta_reader_get_seq(self.ptr, i);
            let len = fasta_reader_get_seq_len(self.ptr, i) as usize;
            std::slice::from_raw_parts(ptr as *const u8, len).to_vec()
        }
    }

    pub fn seq_len(&self, i: u32) -> u32 {
        unsafe { fasta_reader_get_seq_len(self.ptr, i) }
    }
}

impl Drop for OriginalFastaReader {
    fn drop(&mut self) {
        unsafe { fasta_reader_free(self.ptr) }
    }
}

// --- AbundanceLoader ---

unsafe extern "C" {
    fn abundance_loader_new(path: *const c_char) -> *mut std::ffi::c_void;
    fn abundance_loader_free(loader: *mut std::ffi::c_void);
    fn abundance_loader_get_num(loader: *mut std::ffi::c_void) -> i32;
    fn abundance_loader_get_abundance_by_index(loader: *mut std::ffi::c_void, index: i32) -> f64;
    fn abundance_loader_get_abundance_by_header(
        loader: *mut std::ffi::c_void,
        header: *const c_char,
    ) -> f64;
    fn abundance_loader_is_parse_success(loader: *mut std::ffi::c_void) -> i32;
}

pub struct OriginalAbundanceLoader {
    ptr: *mut std::ffi::c_void,
}

impl OriginalAbundanceLoader {
    pub fn new(path: &std::path::Path) -> Self {
        let c_path = CString::new(path.to_str().unwrap()).unwrap();
        let ptr = unsafe { abundance_loader_new(c_path.as_ptr()) };
        Self { ptr }
    }

    pub fn num_records(&self) -> i32 {
        unsafe { abundance_loader_get_num(self.ptr) }
    }

    pub fn abundance_by_index(&self, i: i32) -> f64 {
        unsafe { abundance_loader_get_abundance_by_index(self.ptr, i) }
    }

    pub fn abundance_by_header(&self, header: &str) -> f64 {
        let c_header = CString::new(header).unwrap();
        unsafe { abundance_loader_get_abundance_by_header(self.ptr, c_header.as_ptr()) }
    }

    pub fn is_parse_success(&self) -> bool {
        unsafe { abundance_loader_is_parse_success(self.ptr) != 0 }
    }
}

impl Drop for OriginalAbundanceLoader {
    fn drop(&mut self) {
        unsafe { abundance_loader_free(self.ptr) }
    }
}

// --- quickSort ---

unsafe extern "C" {
    fn quicksort_sort_doubles(arr: *mut f64, rank: *mut i32, num: i32);
}

/// Sort doubles in descending order using the original MaxBin2 quickSort.
/// Optionally reorders a parallel rank array.
pub fn original_quicksort(arr: &mut [f64], rank: Option<&mut [i32]>) {
    let num = arr.len() as i32;
    let rank_ptr = match rank {
        Some(r) => {
            assert_eq!(r.len(), arr.len());
            r.as_mut_ptr()
        }
        None => std::ptr::null_mut(),
    };
    unsafe {
        quicksort_sort_doubles(arr.as_mut_ptr(), rank_ptr, num);
    }
}

// --- kmerMap ---

unsafe extern "C" {
    fn kmermap_new(kmerlen: i32, is_symmetric: i32) -> *mut std::ffi::c_void;
    fn kmermap_free(ptr: *mut std::ffi::c_void);
    fn kmermap_get_entry_num(ptr: *mut std::ffi::c_void) -> i32;
    fn kmermap_get_mapping(ptr: *mut std::ffi::c_void, kmer: *const c_char) -> i32;
    fn kmermap_get_reverse_mapping_str(ptr: *mut std::ffi::c_void, kmer: *const c_char) -> i32;
    fn kmermap_get_reverse_mapping_idx(ptr: *mut std::ffi::c_void, index: i32) -> i32;
    fn kmermap_get_full_table(ptr: *mut std::ffi::c_void, out: *mut i32, kmerlen: i32);
}

pub struct OriginalKmerMap {
    ptr: *mut std::ffi::c_void,
    kmerlen: i32,
}

impl OriginalKmerMap {
    pub fn new(kmerlen: i32, symmetric: bool) -> Self {
        let ptr = unsafe { kmermap_new(kmerlen, if symmetric { 1 } else { 0 }) };
        Self { ptr, kmerlen }
    }

    pub fn entry_num(&self) -> i32 {
        unsafe { kmermap_get_entry_num(self.ptr) }
    }

    pub fn get_mapping(&self, kmer: &str) -> i32 {
        let c_kmer = CString::new(kmer).unwrap();
        unsafe { kmermap_get_mapping(self.ptr, c_kmer.as_ptr()) }
    }

    pub fn get_reverse_mapping_str(&self, kmer: &str) -> i32 {
        let c_kmer = CString::new(kmer).unwrap();
        unsafe { kmermap_get_reverse_mapping_str(self.ptr, c_kmer.as_ptr()) }
    }

    pub fn get_reverse_mapping_idx(&self, index: i32) -> i32 {
        unsafe { kmermap_get_reverse_mapping_idx(self.ptr, index) }
    }

    pub fn get_full_table(&self) -> Vec<i32> {
        let total = 4i32.pow(self.kmerlen as u32);
        let mut out = vec![0i32; total as usize];
        unsafe { kmermap_get_full_table(self.ptr, out.as_mut_ptr(), self.kmerlen) };
        out
    }

    /// Return the raw pointer for use by other FFI wrappers (e.g. Profiler).
    pub fn as_ptr(&self) -> *mut std::ffi::c_void {
        self.ptr
    }
}

impl Drop for OriginalKmerMap {
    fn drop(&mut self) {
        unsafe { kmermap_free(self.ptr) }
    }
}

// --- Profiler ---

unsafe extern "C" {
    fn profiler_new(
        kmerlen: i32,
        seq: *const c_char,
        kmap_ptr: *mut std::ffi::c_void,
    ) -> *mut std::ffi::c_void;
    fn profiler_free(ptr: *mut std::ffi::c_void);
    fn profiler_get_profile(ptr: *mut std::ffi::c_void, out: *mut f64, entry_num: i32);
    fn profiler_get_percent_n(ptr: *mut std::ffi::c_void) -> f32;
}

pub struct OriginalProfiler {
    ptr: *mut std::ffi::c_void,
}

impl OriginalProfiler {
    /// Create a new Profiler. The kmerMap must outlive this Profiler.
    pub fn new(kmerlen: i32, seq: &str, kmap: &OriginalKmerMap) -> Self {
        let c_seq = CString::new(seq).unwrap();
        let ptr = unsafe { profiler_new(kmerlen, c_seq.as_ptr(), kmap.as_ptr()) };
        Self { ptr }
    }

    pub fn get_profile(&self, entry_num: i32) -> Vec<f64> {
        let mut out = vec![0.0f64; entry_num as usize];
        unsafe { profiler_get_profile(self.ptr, out.as_mut_ptr(), entry_num) };
        out
    }

    pub fn get_percent_n(&self) -> f32 {
        unsafe { profiler_get_percent_n(self.ptr) }
    }
}

impl Drop for OriginalProfiler {
    fn drop(&mut self) {
        unsafe { profiler_free(self.ptr) }
    }
}

// --- NormalDistribution ---

unsafe extern "C" {
    fn normaldist_new(mean: f64, std_dev: f64) -> *mut std::ffi::c_void;
    fn normaldist_free(ptr: *mut std::ffi::c_void);
    fn normaldist_prob(ptr: *mut std::ffi::c_void, input: f64) -> f64;
}

pub struct OriginalNormalDistribution {
    ptr: *mut std::ffi::c_void,
}

impl OriginalNormalDistribution {
    pub fn new(mean: f64, std_dev: f64) -> Self {
        let ptr = unsafe { normaldist_new(mean, std_dev) };
        Self { ptr }
    }

    pub fn prob(&self, input: f64) -> f64 {
        unsafe { normaldist_prob(self.ptr, input) }
    }
}

impl Drop for OriginalNormalDistribution {
    fn drop(&mut self) {
        unsafe { normaldist_free(self.ptr) }
    }
}

// --- EucDist ---

unsafe extern "C" {
    fn eucdist_new(kmerlen: i32) -> *mut std::ffi::c_void;
    fn eucdist_free(ptr: *mut std::ffi::c_void);
    fn eucdist_get_dist_seq(
        ptr: *mut std::ffi::c_void,
        seq1: *const c_char,
        seq2: *const c_char,
    ) -> f64;
    fn eucdist_get_dist_profile(ptr: *mut std::ffi::c_void, pro1: *mut f64, pro2: *mut f64) -> f64;
}

pub struct OriginalEucDist {
    ptr: *mut std::ffi::c_void,
}

impl OriginalEucDist {
    pub fn new(kmerlen: i32) -> Self {
        let ptr = unsafe { eucdist_new(kmerlen) };
        Self { ptr }
    }

    pub fn get_dist_seq(&self, seq1: &str, seq2: &str) -> f64 {
        let c1 = CString::new(seq1).unwrap();
        let c2 = CString::new(seq2).unwrap();
        unsafe { eucdist_get_dist_seq(self.ptr, c1.as_ptr(), c2.as_ptr()) }
    }

    pub fn get_dist_profile(&self, pro1: &[f64], pro2: &[f64]) -> f64 {
        let mut p1 = pro1.to_vec();
        let mut p2 = pro2.to_vec();
        unsafe { eucdist_get_dist_profile(self.ptr, p1.as_mut_ptr(), p2.as_mut_ptr()) }
    }
}

impl Drop for OriginalEucDist {
    fn drop(&mut self) {
        unsafe { eucdist_free(self.ptr) }
    }
}

// --- SpearmanDist ---

unsafe extern "C" {
    fn spearmandist_new(kmerlen: i32) -> *mut std::ffi::c_void;
    fn spearmandist_free(ptr: *mut std::ffi::c_void);
    fn spearmandist_get_dist_seq(
        ptr: *mut std::ffi::c_void,
        seq1: *const c_char,
        seq2: *const c_char,
    ) -> f64;
    fn spearmandist_get_dist_profile(
        ptr: *mut std::ffi::c_void,
        pro1: *mut f64,
        pro2: *mut f64,
    ) -> f64;
    fn spearmandist_set_normalization(ptr: *mut std::ffi::c_void, normalize: i32);
}

pub struct OriginalSpearmanDist {
    ptr: *mut std::ffi::c_void,
}

impl OriginalSpearmanDist {
    pub fn new(kmerlen: i32) -> Self {
        let ptr = unsafe { spearmandist_new(kmerlen) };
        Self { ptr }
    }

    pub fn get_dist_seq(&self, seq1: &str, seq2: &str) -> f64 {
        let c1 = CString::new(seq1).unwrap();
        let c2 = CString::new(seq2).unwrap();
        unsafe { spearmandist_get_dist_seq(self.ptr, c1.as_ptr(), c2.as_ptr()) }
    }

    pub fn get_dist_profile(&self, pro1: &[f64], pro2: &[f64]) -> f64 {
        let mut p1 = pro1.to_vec();
        let mut p2 = pro2.to_vec();
        unsafe { spearmandist_get_dist_profile(self.ptr, p1.as_mut_ptr(), p2.as_mut_ptr()) }
    }

    pub fn set_normalization(&self, normalize: bool) {
        unsafe { spearmandist_set_normalization(self.ptr, if normalize { 1 } else { 0 }) }
    }
}

impl Drop for OriginalSpearmanDist {
    fn drop(&mut self) {
        unsafe { spearmandist_free(self.ptr) }
    }
}

// --- EManager ---

unsafe extern "C" {
    fn emanager_new(
        fasta_path: *const c_char,
        abund_path: *const c_char,
        output_prefix: *const c_char,
    ) -> *mut std::ffi::c_void;
    fn emanager_run(ptr: *mut std::ffi::c_void, seedfile: *const c_char) -> i32;
    fn emanager_set_thread_num(ptr: *mut std::ffi::c_void, num: i32);
    fn emanager_free(ptr: *mut std::ffi::c_void);
}

pub struct OriginalEManager {
    ptr: *mut std::ffi::c_void,
}

impl OriginalEManager {
    pub fn new(
        fasta_path: &std::path::Path,
        abund_path: &std::path::Path,
        output_prefix: &str,
    ) -> Self {
        let c_fasta = CString::new(fasta_path.to_str().unwrap()).unwrap();
        let c_abund = CString::new(abund_path.to_str().unwrap()).unwrap();
        let c_output = CString::new(output_prefix).unwrap();
        let ptr = unsafe { emanager_new(c_fasta.as_ptr(), c_abund.as_ptr(), c_output.as_ptr()) };
        Self { ptr }
    }

    pub fn set_thread_num(&self, num: i32) {
        unsafe { emanager_set_thread_num(self.ptr, num) }
    }

    pub fn run(&self, seedfile: &std::path::Path) -> i32 {
        let c_seed = CString::new(seedfile.to_str().unwrap()).unwrap();
        unsafe { emanager_run(self.ptr, c_seed.as_ptr()) }
    }
}

impl Drop for OriginalEManager {
    fn drop(&mut self) {
        unsafe { emanager_free(self.ptr) }
    }
}
