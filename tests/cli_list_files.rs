use maxbin_rs::cli::{Command, PipelineArgs};
/// Tests for CLI options that lack equivalence test coverage:
/// -reads_list, -abund_list, -min_contig_length, -markerset.
///
/// These test CLI parsing and file expansion only — not pipeline behavior.
/// For pipeline equivalence, see tests/pipeline-stages.sh.
use std::io::Write;

fn parse(args: &[&str]) -> PipelineArgs {
    let mut full: Vec<&str> = vec!["--contig", "c.fa", "--out", "out", "--abund", "/dev/null"];
    full.extend_from_slice(args);
    let strings: Vec<String> = full.iter().map(|s| s.to_string()).collect();
    match maxbin_rs::cli::parse_from(&strings) {
        Command::Pipeline(p) => p,
        other => panic!("Expected Pipeline, got {:?}", other),
    }
}

// --- reads_list ---

#[test]
fn reads_list_expands_to_paths() {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    writeln!(f, "/data/sample1.fq.gz").unwrap();
    writeln!(f, "/data/sample2.fq.gz").unwrap();
    writeln!(f, "  /data/sample3.fq.gz  ").unwrap(); // whitespace
    writeln!(f).unwrap(); // blank line
    let cli = parse(&["--reads-list", f.path().to_str().unwrap()]);
    let files = cli.all_reads_files().unwrap();
    assert_eq!(files.len(), 3);
    assert_eq!(files[0].to_str().unwrap(), "/data/sample1.fq.gz");
    assert_eq!(files[2].to_str().unwrap(), "/data/sample3.fq.gz");
}

#[test]
fn abund_list_expands_to_paths() {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    writeln!(f, "/data/abund1.txt").unwrap();
    writeln!(f, "/data/abund2.txt").unwrap();
    let cli = parse(&["--abund-list", f.path().to_str().unwrap()]);
    let files = cli.all_abund_files().unwrap();
    // 1 from the helper's -abund /dev/null + 2 from the list file
    assert_eq!(files.len(), 3);
    assert_eq!(files[1].to_str().unwrap(), "/data/abund1.txt");
    assert_eq!(files[2].to_str().unwrap(), "/data/abund2.txt");
}

#[test]
fn reads_list_merges_with_direct_flags() {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    writeln!(f, "/data/list_sample.fq.gz").unwrap();
    let cli = parse(&[
        "--reads",
        "/data/direct.fq.gz",
        "--reads-list",
        f.path().to_str().unwrap(),
    ]);
    let files = cli.all_reads_files().unwrap();
    assert_eq!(files.len(), 2);
    assert_eq!(files[0].to_str().unwrap(), "/data/direct.fq.gz");
    assert_eq!(files[1].to_str().unwrap(), "/data/list_sample.fq.gz");
}

#[test]
fn empty_list_produces_no_extra_files() {
    let f = tempfile::NamedTempFile::new().unwrap();
    // empty file
    let cli = parse(&["--reads-list", f.path().to_str().unwrap()]);
    let files = cli.all_reads_files().unwrap();
    assert_eq!(files.len(), 0);
}

// --- min_contig_length ---

#[test]
fn min_contig_length_default() {
    let cli = parse(&[]);
    assert_eq!(cli.min_contig_length, 1000);
}

#[test]
fn min_contig_length_custom() {
    let cli = parse(&["--min-contig-length", "2500"]);
    assert_eq!(cli.min_contig_length, 2500);
}

// --- markerset ---

#[test]
fn markerset_default_uses_107() {
    let cli = parse(&[]);
    assert_eq!(cli.markerset, 107);
    assert_eq!(cli.marker_hmm_filename(), "marker.hmm");
}

#[test]
fn markerset_40_uses_bacar() {
    let cli = parse(&["--markerset", "40"]);
    assert_eq!(cli.markerset, 40);
    assert_eq!(cli.marker_hmm_filename(), "bacar_marker.hmm");
}

// --- double-dash normalization ---

#[test]
fn double_dash_reads_list() {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    writeln!(f, "/data/sample.fq.gz").unwrap();
    let cli = parse(&["--reads-list", f.path().to_str().unwrap()]);
    let files = cli.all_reads_files().unwrap();
    assert_eq!(files.len(), 1);
}
