# Reference Output Hashes

MaxBin2 2.2.7 run on nf-core Bacteroides fragilis test data.
Determinism verified: two identical runs produce byte-for-byte identical output.

## Invocation

```bash
run_MaxBin.pl -contig input_contigs.fa.gz -reads input_reads1.fastq.gz -thread 2 -out test1
```

Input files (from nf-core/test-datasets, modules branch):
- `test1.contigs.fa.gz` — contigs
- `test1_1.fastq.gz` — single reads file

## Output Hashes (SHA-256)

```
93fd86ff1d4d982a9f91cc7d1547346bd5eb83e1140acd5af4d2dfbc61851d2f  test1.001.fasta
53ed17404bc03ec6926598fa055aa67f28472fc8cf2b8f54bbb8f43fd4fb8a3d  test1.002.fasta
581462b4539f078a11ec79fece482efa32dbf54bb8f34539924fe3f135f8363d  test1.summary
feccc12ccf65515d156be5d1d569e6c600d79046a6ba67d9ac3713683acd544d  test1.marker
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  test1.noclass (empty)
efe680fd36316464ef803504e9bf46446adebe822b97831b330999e1fd063067  test1.tooshort
```

## Run Details

- 2 bins produced, 7 EM iterations
- Bin 1: abundance 4.09, 14.0% completeness, 911KB genome, 43.1% GC
- Bin 2: abundance 3.11, 14.0% completeness, 796KB genome, 52.1% GC
- 0 unclassified contigs, 523 lines of too-short contigs
- `prob_threshold` used: 0.50 (not 0.9 as help text claims — known bug)
- Elapsed: ~4 seconds with 2 threads

## Files NOT hashed (unstable content)

- `test1.log` — contains timestamps and absolute paths
- `test1.abund1` — derived from Bowtie2 mapping (should be stable but less interesting)
- `test1.marker_of_each_bin.tar.gz` — tarball ordering may vary
