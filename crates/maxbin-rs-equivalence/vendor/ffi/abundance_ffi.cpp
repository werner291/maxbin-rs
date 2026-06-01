// C wrapper around MaxBin2's AbundanceLoader for FFI equivalence testing.
#include "AbundanceLoader.h"

extern "C" {

void* abundance_loader_new(const char* path) {
    return new AbundanceLoader(const_cast<char*>(path));
}

void abundance_loader_free(void* loader) {
    delete static_cast<AbundanceLoader*>(loader);
}

int abundance_loader_get_num(void* loader) {
    return static_cast<AbundanceLoader*>(loader)->getNum();
}

double abundance_loader_get_abundance_by_index(void* loader, int index) {
    return static_cast<AbundanceLoader*>(loader)->getAbundance(index);
}

double abundance_loader_get_abundance_by_header(void* loader, const char* header) {
    return static_cast<AbundanceLoader*>(loader)->getAbundance(const_cast<char*>(header));
}

int abundance_loader_is_parse_success(void* loader) {
    return static_cast<AbundanceLoader*>(loader)->is_parse_success() ? 1 : 0;
}

} // extern "C"
