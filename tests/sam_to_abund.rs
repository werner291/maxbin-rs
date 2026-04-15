/// Test helper: calls compute_abundance_from_sam and writes the result.
/// Usage: cargo test --test sam_to_abund -- --ignored --nocapture SAM_PATH OUT_PATH
///
/// This is invoked by tests/pipeline-stages.sh to compare our SAM→abundance
/// computation against the Perl's getsam() on the same SAM file.

#[test]
#[ignore]
fn sam_to_abund() {
    let sam = std::env::var("SAM_INPUT").expect("SAM_INPUT env var required");
    let out = std::env::var("ABUND_OUTPUT").expect("ABUND_OUTPUT env var required");

    maxbin_rs::pipeline::compute_abundance_from_sam(
        std::path::Path::new(&sam),
        std::path::Path::new(&out),
    )
    .expect("compute_abundance_from_sam failed");

    eprintln!("Wrote abundance to {out}");
}
