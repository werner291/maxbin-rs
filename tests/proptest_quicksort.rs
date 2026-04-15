/// Property-based equivalence tests for quicksort.rs.
///
/// - sort_descending(): covered by quicksort_equivalence (proptest, 500 cases + duplicates)
///   Inline FFI test: quicksort.rs::tests::equivalence_with_original (hardcoded)

use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(500))]

    #[test]
    fn quicksort_equivalence(input in prop::collection::vec(-1e6f64..1e6, 0..200)) {
        // Rust
        let mut rust_arr = input.clone();
        let mut rust_idx: Vec<i32> = (0..input.len() as i32).collect();
        maxbin_rs::quicksort::sort_descending(&mut rust_arr, Some(&mut rust_idx));

        // C++
        let mut orig_arr = input.clone();
        let mut orig_idx: Vec<i32> = (0..input.len() as i32).collect();
        maxbin_rs::original_ffi::original_quicksort(&mut orig_arr, Some(&mut orig_idx));

        prop_assert_eq!(&rust_arr, &orig_arr, "array mismatch");
        prop_assert_eq!(&rust_idx, &orig_idx, "index mismatch");
    }

    /// sort_descending with indices: None — the C++ also supports null rank pointer.
    #[test]
    fn quicksort_no_indices(input in prop::collection::vec(-1e6f64..1e6, 0..200)) {
        let mut rust_arr = input.clone();
        maxbin_rs::quicksort::sort_descending(&mut rust_arr, None);

        let mut orig_arr = input.clone();
        maxbin_rs::original_ffi::original_quicksort(&mut orig_arr, None);

        prop_assert_eq!(&rust_arr, &orig_arr, "array mismatch (no indices)");
    }

    #[test]
    fn quicksort_with_duplicates(
        base in prop::collection::vec(-100.0f64..100.0, 1..50),
        repeats in 1usize..5
    ) {
        // Create input with many duplicate values to stress tie-breaking
        let input: Vec<f64> = base.iter().cycle().take(base.len() * repeats).copied().collect();

        let mut rust_arr = input.clone();
        let mut rust_idx: Vec<i32> = (0..input.len() as i32).collect();
        maxbin_rs::quicksort::sort_descending(&mut rust_arr, Some(&mut rust_idx));

        let mut orig_arr = input.clone();
        let mut orig_idx: Vec<i32> = (0..input.len() as i32).collect();
        maxbin_rs::original_ffi::original_quicksort(&mut orig_arr, Some(&mut orig_idx));

        prop_assert_eq!(&rust_arr, &orig_arr, "array mismatch with duplicates");
        prop_assert_eq!(&rust_idx, &orig_idx, "index mismatch with duplicates");
    }
}
