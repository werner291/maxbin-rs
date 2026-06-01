// Thin C wrapper around MaxBin2's fastaReader for FFI equivalence testing.
#include "fastaReader.h"

extern "C" {

void* fasta_reader_new(const char* path) {
    return new fastaReader(const_cast<char*>(path));
}

void fasta_reader_free(void* reader) {
    delete static_cast<fastaReader*>(reader);
}

unsigned int fasta_reader_get_num(void* reader) {
    return static_cast<fastaReader*>(reader)->getNum();
}

const char* fasta_reader_get_header(void* reader, unsigned int i) {
    return static_cast<fastaReader*>(reader)->getHeaderByNum(i);
}

const char* fasta_reader_get_seq(void* reader, unsigned int i) {
    return static_cast<fastaReader*>(reader)->getSeqByNum(i);
}

unsigned int fasta_reader_get_seq_len(void* reader, unsigned int i) {
    return static_cast<fastaReader*>(reader)->getSeqLenByNum(i);
}

} // extern "C"
