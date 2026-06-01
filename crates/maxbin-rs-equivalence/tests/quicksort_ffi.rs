//! FFI equivalence test for the descending quicksort.
//!
//! Compares `maxbin_core::quicksort::sort_descending` against the original
//! MaxBin2 C++ quickSort. Previously an in-source test in `quicksort.rs`.
use maxbin_core::quicksort::sort_descending;
use maxbin_rs_equivalence::original_ffi::original_quicksort;

#[test]
fn equivalence_with_original() {
    // Test a variety of inputs against the C++ quickSort via FFI.
    let test_cases: Vec<Vec<f64>> = vec![
        vec![1.0, 5.0, 3.0, 2.0, 4.0],
        vec![5.0, 4.0, 3.0, 2.0, 1.0],
        vec![1.0, 1.0, 1.0],
        vec![3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0],
        vec![42.0],
        vec![1.0, 2.0],
        vec![2.0, 1.0],
        vec![0.001, 0.01, 0.1, 1.0, 10.0, 100.0],
        vec![-1.0, 0.0, 1.0, -0.5, 0.5],
    ];

    for (case_idx, case) in test_cases.iter().enumerate() {
        // Rust version
        let mut rust_arr = case.clone();
        let mut rust_idx: Vec<i32> = (0..case.len() as i32).collect();
        sort_descending(&mut rust_arr, Some(&mut rust_idx));

        // Original C++ version
        let mut orig_arr = case.clone();
        let mut orig_idx: Vec<i32> = (0..case.len() as i32).collect();
        original_quicksort(&mut orig_arr, Some(&mut orig_idx));

        assert_eq!(
            rust_arr, orig_arr,
            "array mismatch for case {case_idx}: {case:?}"
        );
        assert_eq!(
            rust_idx, orig_idx,
            "index mismatch for case {case_idx}: {case:?}"
        );
    }
}
