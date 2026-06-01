fn main() {
    let cmd = maxbin_rs::cli::parse();

    let result = match cmd {
        maxbin_rs::cli::Command::Filter(ref args) => maxbin_rs::pipeline::run_filter(args),
        maxbin_rs::cli::Command::Seeds(ref args) => maxbin_rs::pipeline::run_seeds(args),
        maxbin_rs::cli::Command::Em(ref args) => maxbin_rs::pipeline::run_em(args),
        maxbin_rs::cli::Command::CppEm(ref args) => maxbin_rs::pipeline::run_cpp_em(args),
        maxbin_rs::cli::Command::SamToAbund(ref args) => {
            maxbin_rs::pipeline::run_sam_to_abund(args)
        }
        maxbin_rs::cli::Command::Pipeline(ref args) => maxbin_rs::pipeline::run_pipeline(args),
    };

    if let Err(e) = result {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
