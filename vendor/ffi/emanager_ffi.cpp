// FFI wrapper for EManager — exposes the full EM pipeline to Rust.
// For equivalence testing only.

#include "EManager.h"
#include <cstring>

extern "C" {

void* emanager_new(const char* fasta_path, const char* abund_path, const char* output_prefix) {
    // EManager takes char* (non-const) but doesn't modify them.
    char* fs = strdup(fasta_path);
    char* ab = strdup(abund_path);
    char* out = strdup(output_prefix);
    EManager* em = new EManager(fs, ab, out);
    // Note: EManager stores these pointers, so we can't free them here.
    // They'll leak, but this is only for testing.
    return (void*)em;
}

int emanager_run(void* ptr, const char* seedfile) {
    EManager* em = (EManager*)ptr;
    char* sf = strdup(seedfile);
    int ret = em->run(sf);
    free(sf);
    return ret;
}

void emanager_set_thread_num(void* ptr, int num) {
    EManager* em = (EManager*)ptr;
    em->setThreadNum(num);
}

void emanager_free(void* ptr) {
    EManager* em = (EManager*)ptr;
    delete em;
}

} // extern "C"
