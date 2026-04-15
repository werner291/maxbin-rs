# Pre-computed pipeline intermediates for equivalence testing.
#
# The MaxBin2 pipeline includes slow steps (Bowtie2 read alignment, HMMER
# marker gene search) that take minutes per dataset. Rather than re-running
# these every time we test, we run the original MaxBin2 pipeline once per
# dataset and cache the intermediate files as a Nix derivation. "Cached as a
# Nix derivation" means: the output is stored in the Nix store (a content-
# addressed cache on disk) and never recomputed unless the inputs change.
#
# The pipeline stage tests (tests/pipeline-stages.sh) then compare the Rust
# reimplementation against these cached outputs, stage by stage.
#
# Two variants:
#   mkPipelineIntermediates       — from raw reads (runs Bowtie2 + HMMER)
#   mkPipelineIntermediatesFromAbund — from pre-computed abundance (HMMER only)
#
# Intermediate files produced (reads variant):
#   sam      — Bowtie2 SAM alignment
#   hmmout   — HMMER marker gene hits
#   abund    — per-contig abundance table
#   tooshort — contigs below min_length (1000bp)
#   seed     — initial EM seeds (regenerated from HMM output)
#
# Intermediate files produced (abund variant — no sam):
#   hmmout   — HMMER marker gene hits
#   abund    — per-contig abundance table (MaxBin format)
#   tooshort — contigs below min_length (1000bp)
#   seed     — initial EM seeds (regenerated from HMM output)

{ runCommand, perl, gawk, maxbin2 }:

let
  mkPipelineIntermediates = { name, contigs, reads }: runCommand
    "maxbin2-intermediates-${name}" {
      nativeBuildInputs = [ maxbin2 perl ];
    } ''
      mkdir -p work
      cp ${contigs} work/contigs.fa.gz
      cp ${reads} work/reads.fq.gz
      chmod u+w work/*
      cd work
      run_MaxBin.pl -contig contigs.fa.gz -reads reads.fq.gz \
        -thread 4 -out gen -preserve_intermediate 2>&1 | tail -5

      mkdir -p $out

      # These survive -preserve_intermediate:
      cp gen.sam0 $out/sam
      cp gen.contig.tmp.hmmout $out/hmmout
      cp gen.contig.tmp.reads.abund1 $out/abund
      cp gen.tooshort $out/tooshort
      cp gen.contig.tmp.frag.faa $out/genes.faa 2>/dev/null || true

      # gen.contig.tmp and gen.seed are deleted by run_MaxBin.pl.
      # We regenerate them from the HMM output using the Perl _getmarker.pl.
      zcat contigs.fa.gz > plain_contigs.fa 2>/dev/null || cp contigs.fa.gz plain_contigs.fa

      # Filter contigs >= 1000bp (matching Perl's checkContig)
      awk 'BEGIN{seq="";name=""} /^>/{if(name!="" && length(seq)>=1000) print name"\n"seq; name=$0;seq=""} !/^>/{seq=seq$0} END{if(name!="" && length(seq)>=1000) print name"\n"seq}' plain_contigs.fa > filtered.fa

      # Regenerate seed file from HMM output using the original Perl
      MARKER_PL="${maxbin2}/libexec/maxbin2/_getmarker.pl"
      perl -e "require \"$MARKER_PL\"; my \$r = gethmmmarker(\"$out/hmmout\", \"filtered.fa\", 1000, \"$out/seed\"); if (\$r == -1) { \$r = gethmmmarker(\"$out/hmmout\", \"filtered.fa\", 1000, \"$out/seed\", 1); }" 2>/dev/null || true

      echo "Intermediates for ${name}:"
      ls -lh $out/
    '';

  # Variant for datasets that ship pre-computed abundance depth files.
  # Skips Bowtie2 entirely — passes the abundance file directly to MaxBin2
  # via -abund. The depth file format determines whether conversion is needed:
  #
  #   depthFormat = "maxbin"  — two-column TSV: contig_name<TAB>depth
  #                             Used as-is.
  #   depthFormat = "metabat" — multi-column TSV from MetaBAT's jgi_summarize_bam_depths:
  #                             contigName, contigLen, totalAvgDepth, ...
  #                             Converted to MaxBin format by extracting columns 1 and 3.
  #
  mkPipelineIntermediatesFromAbund = { name, contigs, depth, depthFormat ? "maxbin" }:
    runCommand "maxbin2-intermediates-${name}" {
      nativeBuildInputs = [ maxbin2 perl gawk ];
    } ''
      mkdir -p work
      cp ${contigs} work/contigs.fa.gz
      cp ${depth} work/depth_raw
      chmod u+w work/*
      cd work

      # Decompress depth file if gzipped (check magic bytes, not `file` command).
      if gzip -t depth_raw 2>/dev/null; then
        mv depth_raw depth_raw.gz
        gzip -d depth_raw.gz
      fi

      # Convert to MaxBin format if needed.
      ${if depthFormat == "metabat" then ''
        # MetaBAT format: contigName<TAB>contigLen<TAB>totalAvgDepth<TAB>...
        # MaxBin format:  contigName<TAB>depth
        # Skip the header line (starts with "contigName") and extract columns 1,3.
        awk -F'\t' 'NR > 1 { print $1 "\t" $3 }' depth_raw > depth.txt
      '' else ''
        cp depth_raw depth.txt
      ''}

      run_MaxBin.pl -contig contigs.fa.gz -abund depth.txt \
        -thread 4 -out gen -preserve_intermediate 2>&1 | tail -5

      mkdir -p $out

      # No SAM file when using pre-computed abundance.
      cp gen.contig.tmp.hmmout $out/hmmout
      # Use the original depth file we fed to MaxBin2 (already in MaxBin format).
      cp depth.txt $out/abund
      cp gen.tooshort $out/tooshort
      cp gen.contig.tmp.frag.faa $out/genes.faa 2>/dev/null || true

      # Regenerate seed file (same as the reads variant).
      zcat contigs.fa.gz > plain_contigs.fa 2>/dev/null || cp contigs.fa.gz plain_contigs.fa
      awk 'BEGIN{seq="";name=""} /^>/{if(name!="" && length(seq)>=1000) print name"\n"seq; name=$0;seq=""} !/^>/{seq=seq$0} END{if(name!="" && length(seq)>=1000) print name"\n"seq}' plain_contigs.fa > filtered.fa
      MARKER_PL="${maxbin2}/libexec/maxbin2/_getmarker.pl"
      perl -e "require \"$MARKER_PL\"; my \$r = gethmmmarker(\"$out/hmmout\", \"filtered.fa\", 1000, \"$out/seed\"); if (\$r == -1) { \$r = gethmmmarker(\"$out/hmmout\", \"filtered.fa\", 1000, \"$out/seed\", 1); }" 2>/dev/null || true

      echo "Intermediates for ${name}:"
      ls -lh $out/
    '';
in
  { inherit mkPipelineIntermediates mkPipelineIntermediatesFromAbund; }
