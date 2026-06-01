// C wrapper around MaxBin2's Profiler for FFI equivalence testing.
#include "Profiler.h"
#include <stdlib.h>
#include <string.h>

extern "C" {

void* profiler_new(int kmerlen, const char* seq, void* kmap_ptr) {
    return new Profiler(kmerlen, seq, static_cast<kmerMap*>(kmap_ptr));
}

void profiler_free(void* ptr) {
    delete static_cast<Profiler*>(ptr);
}

// Copy profile data into caller-provided buffer of size entry_num doubles.
void profiler_get_profile(void* ptr, double* out, int entry_num) {
    double* profile = static_cast<Profiler*>(ptr)->getProfile();
    memcpy(out, profile, sizeof(double) * entry_num);
}

float profiler_get_percent_n(void* ptr) {
    return static_cast<Profiler*>(ptr)->getPercentN();
}

} // extern "C"
