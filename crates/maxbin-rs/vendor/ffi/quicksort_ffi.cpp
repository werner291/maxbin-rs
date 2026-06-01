// C wrapper around MaxBin2's quickSort for FFI equivalence testing.
#include "quickSort.h"
#include <stdlib.h>
#include <string.h>

extern "C" {

// Sort an array of doubles in-place (descending, matching original behavior).
// Also reorders a parallel rank array if provided.
void quicksort_sort_doubles(double* arr, int* rank, int num) {
    quickSort qs;
    if (rank != NULL) {
        qs.input(arr, rank, num);
    } else {
        qs.input(arr, num);
    }
    qs.sort_all();
}

} // extern "C"
