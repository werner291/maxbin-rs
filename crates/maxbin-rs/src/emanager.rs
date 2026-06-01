//! Native, file-facing shell around [`maxbin_core::emanager`].
//!
//! The pure EM algorithm (state, E-step, M-step, convergence, classification)
//! lives in `maxbin-core` so it can compile to wasm. This module adds the parts
//! that touch the filesystem: parsing the input files into the in-memory state
//! ([`init_em`]), writing the bin FASTAs and summary ([`write_results`]), and
//! the [`run_pipeline`] convenience that strings them together. The pure items
//! are re-exported so existing `emanager::*` call sites keep working.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use maxbin_core::abundance;
use maxbin_core::fasta;

pub use maxbin_core::emanager::*;

/// Initialize EM state from input files and seed contig names.
///
/// Parses the FASTA and abundance files, then hands the parsed data to
/// [`maxbin_core::emanager::build_state`] for the pure setup.
/// Matches the parsing portion of EManager.cpp:312-422 (init_EM()).
pub fn init_em(
    fasta_path: &Path,
    abundance_paths: &[&Path],
    seed_names: &[String],
    params: &EmParams,
) -> EmState {
    // Matches EManager.cpp:18-19 (constructor): parse FASTA, get seqnum
    let records = fasta::parse_file(fasta_path).expect("failed to parse FASTA");
    let seqnum = records.len();

    // Matches EManager.cpp:214-251 (addAbund()): load abundance per contig,
    // mapping each abundance file onto the contigs in FASTA order.
    let ab_num = abundance_paths.len();
    let mut seq_abundance = Vec::with_capacity(ab_num);
    for abund_path in abundance_paths {
        let abund_records =
            abundance::parse_file(abund_path).expect("failed to parse abundance file");
        let abund_map: HashMap<&str, f64> = abund_records
            .iter()
            .map(|r| (r.header.as_str(), r.abundance))
            .collect();

        let mut abund_vec = vec![0.0f64; seqnum];
        for (i, rec) in records.iter().enumerate() {
            abund_vec[i] = *abund_map
                .get(rec.header.as_str())
                .unwrap_or_else(|| panic!("abundance not found for contig: {}", rec.header));
        }
        seq_abundance.push(abund_vec);
    }

    build_state(records, seq_abundance, seed_names, params)
}

/// Write output files: .NNN.fasta bins, .noclass, .summary
/// Matches EManager.cpp:867-991 (write_result())
pub fn write_results(state: &EmState, output_prefix: &str, params: &EmParams) {
    let seed_num = state.seed_profiles.len();
    let seqnum = state.records.len();

    // Matches EManager.cpp:872-890: name bins with 4-digit zero-padding, skip empty bins
    // Matches EManager.cpp:883: sprintf(bin_name[i], "%s.%04d.fasta", outputfile, j)
    let mut bin_names: Vec<Option<String>> = vec![None; seed_num];
    let mut j = 1;
    for (bn, &count) in bin_names.iter_mut().zip(state.seed_count.iter()) {
        if count > 0 {
            *bn = Some(format!("{}.{:04}.fasta", output_prefix, j));
            j += 1;
        }
    }

    // Matches EManager.cpp:894-944 (write_result() summary section)
    let summary_path = format!("{}.summary", output_prefix);
    let mut summary = std::io::BufWriter::new(
        std::fs::File::create(&summary_path).expect("failed to create summary"),
    );

    #[allow(clippy::needless_range_loop)] // i used as value and multi-array index
    for i in 0..seed_num {
        if state.seed_count[i] > 0 {
            // Matches EManager.cpp:901-907: "Bin [name]\tabundance..."
            write!(summary, "Bin [{}]", bin_names[i].as_ref().unwrap()).unwrap();
            for k in 0..state.ab_num {
                write!(summary, "\t{:.4}", state.seed_abundance[i][k]).unwrap();
            }
            writeln!(summary).unwrap();
            for j in 0..seqnum {
                if state.seq_bin[j] == i as i32 {
                    write!(summary, "\t{}", state.records[j].header).unwrap();
                    for k in 0..state.ab_num {
                        write!(summary, "\t{:.4}", state.seq_abundance[k][j]).unwrap();
                    }
                    writeln!(summary).unwrap();
                }
            }
            writeln!(summary).unwrap();
        }
    }

    // Matches EManager.cpp:926-943: "Bins without any sequences:" section
    write!(summary, "\nBins without any sequences:\n").unwrap();
    let mut empty_j = 1;
    for i in 0..seed_num {
        if state.seed_count[i] == 0 {
            write!(summary, "{}:", empty_j).unwrap();
            for k in 0..state.ab_num {
                write!(summary, " ({:.4})", state.seed_abundance[i][k]).unwrap();
            }
            writeln!(summary).unwrap();
            empty_j += 1;
        }
    }
    drop(summary);

    // Matches EManager.cpp:947-987: open per-bin fasta files and noclass file, write sequences
    let mut bin_files: Vec<Option<std::io::BufWriter<std::fs::File>>> = (0..seed_num)
        .map(|i| {
            bin_names[i].as_ref().map(|name| {
                std::io::BufWriter::new(
                    std::fs::File::create(name).expect("failed to create bin file"),
                )
            })
        })
        .collect();

    // Matches EManager.cpp:956-957: open .noclass file
    let noclass_path = format!("{}.noclass", output_prefix);
    let mut noclass = std::io::BufWriter::new(
        std::fs::File::create(&noclass_path).expect("failed to create noclass"),
    );

    for i in 0..seqnum {
        let header_line = format!(">{}\n", state.records[i].header);
        let seq_str = std::str::from_utf8(&state.records[i].seq).unwrap();

        if state.seq_bin[i] != -1 {
            let bin_idx = state.seq_bin[i] as usize;
            let f = bin_files[bin_idx].as_mut().unwrap();
            f.write_all(header_line.as_bytes()).unwrap();
            write_fasta_seq(f, seq_str, params.fasta_line);
        } else {
            noclass.write_all(header_line.as_bytes()).unwrap();
            write_fasta_seq(&mut noclass, seq_str, params.fasta_line);
        }
    }
}

/// Write a FASTA sequence with line wrapping at `line_width` characters.
/// Matches EManager.cpp:993-1015 (write_fasta()): write FASTA_LINE chars at a time.
fn write_fasta_seq(f: &mut impl Write, seq: &str, line_width: usize) {
    let bytes = seq.as_bytes();
    let len = bytes.len();
    let mut pos = 0;
    while pos < len {
        let end = std::cmp::min(pos + line_width, len);
        f.write_all(&bytes[pos..end]).unwrap();
        f.write_all(b"\n").unwrap();
        pos = end;
    }
}

/// Full pipeline: load data, run EM, classify, write results.
pub fn run_pipeline(
    fasta_path: &Path,
    abundance_paths: &[&Path],
    seed_names: &[String],
    output_prefix: &str,
    params: &EmParams,
) {
    let mut state = init_em(fasta_path, abundance_paths, seed_names, params);
    run_em(&mut state, params);
    classify(&mut state, params);
    write_results(&state, output_prefix, params);
}
