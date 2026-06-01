// C wrapper around MaxBin2's NormalDistribution for FFI equivalence testing.
#include "NormalDistribution.h"

extern "C" {

void* normaldist_new(double mean, double std_dev) {
    return new NormalDistribution(mean, std_dev);
}

void normaldist_free(void* ptr) {
    delete static_cast<NormalDistribution*>(ptr);
}

double normaldist_prob(void* ptr, double input) {
    return static_cast<NormalDistribution*>(ptr)->prob(input);
}

} // extern "C"
