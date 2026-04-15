// C wrapper around MaxBin2's EucDist and SpearmanDist for FFI equivalence testing.
#include "EucDist.h"
#include "SpearmanDist.h"
#include <stdlib.h>
#include <string.h>

extern "C" {

// --- EucDist ---

void* eucdist_new(int kmerlen) {
    return new EucDist(kmerlen);
}

void eucdist_free(void* ptr) {
    delete static_cast<EucDist*>(ptr);
}

double eucdist_get_dist_seq(void* ptr, const char* seq1, const char* seq2) {
    return static_cast<EucDist*>(ptr)->getDist(seq1, seq2);
}

double eucdist_get_dist_profile(void* ptr, double* pro1, double* pro2) {
    return static_cast<EucDist*>(ptr)->getDist(pro1, pro2);
}

// --- SpearmanDist ---

void* spearmandist_new(int kmerlen) {
    return new SpearmanDist(kmerlen);
}

void spearmandist_free(void* ptr) {
    delete static_cast<SpearmanDist*>(ptr);
}

double spearmandist_get_dist_seq(void* ptr, const char* seq1, const char* seq2) {
    return static_cast<SpearmanDist*>(ptr)->getDist(seq1, seq2);
}

double spearmandist_get_dist_profile(void* ptr, double* pro1, double* pro2) {
    return static_cast<SpearmanDist*>(ptr)->getDist(pro1, pro2);
}

void spearmandist_set_normalization(void* ptr, int normalize) {
    static_cast<SpearmanDist*>(ptr)->setNormalization(normalize != 0);
}

} // extern "C"
