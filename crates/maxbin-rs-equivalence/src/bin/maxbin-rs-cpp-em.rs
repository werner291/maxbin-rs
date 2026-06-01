//! `maxbin-rs-cpp-em`: run the original MaxBin2 C++ EM via FFI.
//!
//! This binary exists only for equivalence testing and benchmarking against
//! the original implementation. The core `maxbin-rs` crate no longer carries
//! the C++ FFI, so this lives in the equivalence crate.
//!
//! The C++ EManager handles its own data loading, EM, classification, and
//! file writing internally, so the whole thing is timed as one block.
use std::path::PathBuf;

use clap::Parser;
use maxbin_rs_equivalence::original_ffi::OriginalEManager;

#[derive(Parser)]
#[command(
    name = "maxbin-rs-cpp-em",
    about = "Run the original MaxBin2 C++ EM via FFI (equivalence testing)."
)]
struct Args {
    #[arg(long)]
    contig: PathBuf,
    #[arg(long)]
    abund: PathBuf,
    #[arg(long)]
    seed: PathBuf,
    #[arg(long)]
    out: String,
    #[arg(long, default_value_t = 1)]
    thread: usize,
}

fn main() {
    let args = Args::parse();

    eprintln!("Running original C++ EM via FFI...");
    eprintln!("  contigs: {}", args.contig.display());
    eprintln!("  abund:   {}", args.abund.display());
    eprintln!("  seed:    {}", args.seed.display());
    eprintln!("  threads: {}", args.thread);

    let t = std::time::Instant::now();
    let em = OriginalEManager::new(&args.contig, &args.abund, &args.out);
    em.set_thread_num(args.thread as i32);
    let result = em.run(&args.seed);
    let secs = t.elapsed().as_secs_f64();

    eprintln!("C++ EM completed in {secs:.1}s (returned: {result})");
}
