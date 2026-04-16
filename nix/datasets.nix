# Test datasets for equivalence testing.
#
# Each entry declares URLs and SHA-256 hashes for the input files (contigs
# and sequencing reads). Nix downloads these once and caches them forever
# — re-running a test does not re-download. The hashes ensure that the
# exact same data is used on every machine.
#
# Three datasets of increasing size:
#   - B. fragilis:  nf-core/modules MaxBin2 test data (~1 MB, seconds)
#   - minigut:      nf-core/mag test profile (~10 MB, minutes)
#   - CAPES_S7:     real metagenome, 25K contigs (~2.5 GB, ~10 minutes)
#   - CAMI I High:  CAMI challenge benchmark, 39K contigs (~779 MB, ~1-3 hours)
#   - MetaHIT:      human gut, 60K contigs (~275 MB, ~3 hours confirmed)

{ fetchurl }:

let
  # --- B. fragilis (nf-core/modules test data) ---
  nfcoreBase = "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina";

  # --- minigut (nf-core/mag test profile) ---
  minigutBase = "https://raw.githubusercontent.com/nf-core/test-datasets/mag";
in {
  bfragilis = {
    contigs = fetchurl {
      url = "${nfcoreBase}/fasta/test1.contigs.fa.gz";
      hash = "sha256-Gg0r41zpq3mtceQ9rgAaVPp/fGwdcQLGSP5NvgHDzWo=";
    };
    reads1 = fetchurl {
      url = "${nfcoreBase}/fastq/test1_1.fastq.gz";
      hash = "sha256-qQ9NSKMVSIbnKlIhuN1rEVnqCQ5dviBjTpZbTugBFqs=";
    };
    reads2 = fetchurl {
      url = "${nfcoreBase}/fastq/test2_1.fastq.gz";
      hash = "sha256-M/3gs2mW3kwoogFAI6LPDJS/k29msziF17STHW2iCrE=";
    };
  };

  minigut = {
    contigs = fetchurl {
      url = "${minigutBase}/assemblies/MEGAHIT-test_minigut.contigs.fa.gz";
      hash = "sha256-RwB+UtWpnzSty4pgMnMYCUWFbSaj/KmtEt3SmasiHQY=";
    };
    reads1 = fetchurl {
      url = "${minigutBase}/test_data/test_minigut_R1.fastq.gz";
      hash = "sha256-qQ9NSKMVSIbnKlIhuN1rEVnqCQ5dviBjTpZbTugBFqs=";
    };
    reads2 = fetchurl {
      url = "${minigutBase}/test_data/test_minigut_R2.fastq.gz";
      hash = "sha256-LxSaPk04SvydNbN3glnBzSTX6i/p+AhUi+JyiyOzUGg=";
    };
  };

  capes-s7 = {
    # Real metagenome from nf-core/mag full test profile.
    # 25,244 contigs assembled by MEGAHIT, paired-end reads from ENA.
    contigs = fetchurl {
      url = "https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/mag/results-5dabb0159ac0104885e09f301db22126e8fcb394/Assembly/MEGAHIT/MEGAHIT-CAPES_S7.contigs.fa.gz";
      hash = "sha256-o7t3OksBUNrBG9HZHAJGU5NNdrH3RWzWD9KHh9IWiLM=";
    };
    reads1 = fetchurl {
      url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR320/004/ERR3201914/ERR3201914_1.fastq.gz";
      hash = "sha256-ierz51lrbqa/ZbGuh+sYGsU84BGxRsusDGYhvjWDzfs=";
    };
    reads2 = fetchurl {
      url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR320/004/ERR3201914/ERR3201914_2.fastq.gz";
      hash = "sha256-206yNZY/3rtjcYbx1vgzlTYyVvCqS5B9MAXoR95gG5Y=";
    };
  };

  # --- CAMI I High-Complexity Gold Standard Assembly ---
  # The standard metagenome binning benchmark. 39,140 contigs, 2.80 Gbp.
  # Sczyrba et al., Nature Methods 2017, doi:10.1038/nmeth.4458
  # Pre-computed abundance depth file available; no raw reads needed.
  cami-i-high = {
    contigs = fetchurl {
      url = "https://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/CAMI/CAMI_high_GoldStandardAssembly.fasta.gz";
      hash = "sha256-FMSCUBvdmKNPqU98YPClmIqPnMDVI4thhzy733zoaU4=";
    };
    depth = fetchurl {
      url = "https://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/CAMI/depth-high.txt";
      hash = "sha256-t9g39U/Qu4f8TVICE7WcU77Y07JzGP0rpa1JGlRMMj4=";
    };
  };

  # --- MetaHIT Human Gut (≥2.5 kb contigs) ---
  # 60,619 contigs from 262 human gut metagenome samples (MetaHIT project).
  # Confirmed MaxBin v1 runtime: 3h 24min on 32 threads.
  # Kang et al., PeerJ 2015, doi:10.7717/peerj.1165
  # Pre-computed MaxBin-format depth file available.
  metahit = {
    contigs = fetchurl {
      url = "https://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/MetaHIT/assembly.fa.gz";
      hash = "sha256-kGHyYT2fGHsENIHCgBk1l+4QGKpoTi/hvWoiBqXJpR4=";
    };
    depth = fetchurl {
      url = "https://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Files/MetaHIT/depth_maxbin.txt.gz";
      hash = "sha256-ohMC0ksacqLgayccbcUSxQyiU7BUf57ERzDr1iYHQKE=";
    };
  };
}
