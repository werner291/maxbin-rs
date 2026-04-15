// C wrapper around MaxBin2's kmerMap for FFI equivalence testing.
#include "kmerMap.h"
#include <stdlib.h>
#include <string.h>

extern "C" {

void* kmermap_new(int kmerlen, int is_symmetric) {
    return new kmerMap(kmerlen, is_symmetric != 0);
}

void kmermap_free(void* ptr) {
    delete static_cast<kmerMap*>(ptr);
}

int kmermap_get_entry_num(void* ptr) {
    return static_cast<kmerMap*>(ptr)->getEntryNum();
}

int kmermap_get_mapping(void* ptr, const char* kmer) {
    return static_cast<kmerMap*>(ptr)->getMapping(const_cast<char*>(kmer));
}

int kmermap_get_reverse_mapping_str(void* ptr, const char* kmer) {
    return static_cast<kmerMap*>(ptr)->getReverseMapping(const_cast<char*>(kmer));
}

int kmermap_get_reverse_mapping_idx(void* ptr, int index) {
    return static_cast<kmerMap*>(ptr)->getReverseMapping(index);
}

// Get the entire mapping table for a kmerMap. Caller must provide a buffer
// of size 4^kmerlen ints.
void kmermap_get_full_table(void* ptr, int* out, int kmerlen) {
    int total = 1;
    for (int i = 0; i < kmerlen; i++) total *= 4;
    for (int i = 0; i < total; i++) {
        // numToString + getMapping for each index
        // We can use getMapping with a constructed kmer string
        char buf[32];
        memset(buf, 0, sizeof(buf));
        int n = i;
        for (int j = kmerlen - 1; j >= 0; j--) {
            int k = n % 4;
            n = n / 4;
            switch (k) {
                case 0: buf[j] = 'A'; break;
                case 1: buf[j] = 'T'; break;
                case 2: buf[j] = 'C'; break;
                case 3: buf[j] = 'G'; break;
            }
        }
        out[i] = static_cast<kmerMap*>(ptr)->getMapping(buf);
    }
}

} // extern "C"
