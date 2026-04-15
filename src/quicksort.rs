/// Exact reimplementation of MaxBin2's quickSort.
/// Sorts in DESCENDING order. First-element pivot, unstable.
///
/// This must match the original's behavior exactly, including for equal
/// elements, because the Spearman distance metric depends on the sort order
/// for rank assignment.
///
/// Optionally reorders a parallel index array (used for rank tracking).
///
/// Matches quickSort.cpp:60-111 (sort_all() / sort() methods)
pub fn sort_descending(arr: &mut [f64], indices: Option<&mut [i32]>) {
    let n = arr.len();
    if n <= 1 {
        return;
    }
    // Matches quickSort.cpp:62: sort(0, num - 1)
    sort_recursive(arr, indices, 0, n as i32 - 1, n as i32);
}

// Matches quickSort.cpp:65-111 (sort(int low, int high))
fn sort_recursive(arr: &mut [f64], mut indices: Option<&mut [i32]>, low: i32, high: i32, num: i32) {
    if low < high {
        // Matches quickSort.cpp:70-72: i = low, j = high+1, pivot = arr[i]
        let mut i = low;
        let mut j = high + 1;
        let pivot = arr[low as usize];

        loop {
            // Matches quickSort.cpp:74-84: inner while loops advance i and j
            i += 1;
            while i < num && arr[i as usize] > pivot {
                i += 1;
            }
            j -= 1;
            while j > 0 && arr[j as usize] < pivot {
                j -= 1;
            }
            if i < j {
                // Matches quickSort.cpp:85-96: swap arr[i] and arr[j], plus rank[]
                arr.swap(i as usize, j as usize);
                if let Some(ref mut idx) = indices.as_deref_mut() {
                    idx.swap(i as usize, j as usize);
                }
            } else {
                break;
            }
        }

        // Matches quickSort.cpp:98-106: place pivot at arr[j], swap rank[low] and rank[j]
        arr.swap(low as usize, j as usize);
        if let Some(ref mut idx) = indices.as_deref_mut() {
            idx.swap(low as usize, j as usize);
        }

        // Matches quickSort.cpp:108-109: recurse on left and right partitions
        sort_recursive(arr, indices.as_deref_mut().map(|s| &mut *s), low, j - 1, num);
        sort_recursive(arr, indices, j + 1, high, num);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorts_descending() {
        let mut arr = vec![1.0, 5.0, 3.0, 2.0, 4.0];
        sort_descending(&mut arr, None);
        assert_eq!(arr, vec![5.0, 4.0, 3.0, 2.0, 1.0]);
    }

    #[test]
    fn preserves_parallel_indices() {
        let mut arr = vec![1.0, 5.0, 3.0, 2.0, 4.0];
        let mut idx: Vec<i32> = vec![0, 1, 2, 3, 4];
        sort_descending(&mut arr, Some(&mut idx));
        assert_eq!(arr, vec![5.0, 4.0, 3.0, 2.0, 1.0]);
        assert_eq!(idx, vec![1, 4, 2, 3, 0]);
    }

    #[test]
    fn single_element() {
        let mut arr = vec![42.0];
        sort_descending(&mut arr, None);
        assert_eq!(arr, vec![42.0]);
    }

    #[test]
    fn empty() {
        let mut arr: Vec<f64> = vec![];
        sort_descending(&mut arr, None);
        assert!(arr.is_empty());
    }

    #[test]
    fn already_sorted() {
        let mut arr = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        sort_descending(&mut arr, None);
        assert_eq!(arr, vec![5.0, 4.0, 3.0, 2.0, 1.0]);
    }

    #[test]
    fn with_duplicates() {
        let mut arr = vec![3.0, 1.0, 3.0, 2.0, 1.0];
        sort_descending(&mut arr, None);
        assert_eq!(arr, vec![3.0, 3.0, 2.0, 1.0, 1.0]);
    }

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
            crate::original_ffi::original_quicksort(&mut orig_arr, Some(&mut orig_idx));

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
}
