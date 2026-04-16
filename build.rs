use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Upstream C++ source files from MaxBin2 that we need for FFI.
const UPSTREAM_CPP: &[&str] = &[
    "fastaReader.cpp",
    "fastaReader.h",
    "AbundanceLoader.cpp",
    "AbundanceLoader.h",
    "quickSort.cpp",
    "quickSort.h",
    "kmerMap.cpp",
    "kmerMap.h",
    "Profiler.cpp",
    "Profiler.h",
    "NormalDistribution.cpp",
    "NormalDistribution.h",
    "AbstractDist.cpp",
    "AbstractDist.h",
    "EucDist.cpp",
    "EucDist.h",
    "SpearmanDist.cpp",
    "SpearmanDist.h",
    "logger.cpp",
    "logger.h",
    "ThreadPool.cpp",
    "ThreadPool.h",
    "EManager.cpp",
    "EManager.h",
    "global_inc.h",
];

/// Project-specific FFI wrapper files (thin C wrappers for Rust tests).
const FFI_WRAPPERS: &[&str] = &[
    "fasta_ffi.cpp",
    "abundance_ffi.cpp",
    "quicksort_ffi.cpp",
    "kmermap_ffi.cpp",
    "profiler_ffi.cpp",
    "normaldist_ffi.cpp",
    "distance_ffi.cpp",
    "emanager_ffi.cpp",
];

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let cpp_dir = out_dir.join("maxbin2-cpp");
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    // Rebuild if patch or FFI wrappers change.
    println!("cargo:rerun-if-changed=nix/maxbin2-cpp-ffi.patch");
    println!("cargo:rerun-if-changed=vendor/ffi");
    println!("cargo:rerun-if-env-changed=MAXBIN2_SRC_TARBALL");

    // Step 1: Locate the upstream MaxBin2 source tarball.
    //
    // The nix devshell sets MAXBIN2_SRC_TARBALL to the nix store path.
    // Outside nix, the user must set it manually or place the tarball
    // at a known location.
    let tarball = env::var("MAXBIN2_SRC_TARBALL").unwrap_or_else(|_| {
        // Fallback: try to find it in the nix store by hash (the fetchurl
        // output from nix/maxbin2.nix).
        let nix_path = "/nix/store/mrff0ngxcj46v14sbv27c4rxph0miqhl-MaxBin-2.2.7.tar.gz";
        if Path::new(nix_path).exists() {
            return nix_path.to_string();
        }
        panic!(
            "Cannot find MaxBin2 source tarball.\n\
             Set MAXBIN2_SRC_TARBALL=/path/to/MaxBin-2.2.7.tar.gz\n\
             or enter the nix devshell (nix develop / direnv allow)."
        );
    });

    // Step 2: Extract only the C++ source files we need.
    if cpp_dir.exists() {
        fs::remove_dir_all(&cpp_dir).unwrap();
    }
    fs::create_dir_all(&cpp_dir).unwrap();

    let tar_output = Command::new("tar")
        .args([
            "xzf",
            &tarball,
            "-C",
            cpp_dir.to_str().unwrap(),
            "--strip-components=1",
            "MaxBin-2.2.7/src/",
        ])
        .output()
        .expect("failed to run tar");
    if !tar_output.status.success() {
        panic!(
            "tar extraction failed: {}",
            String::from_utf8_lossy(&tar_output.stderr)
        );
    }

    // The files are now at cpp_dir/src/. Move them up one level for convenience
    // (matching the flat layout build.rs expects).
    let src_subdir = cpp_dir.join("src");
    for entry in fs::read_dir(&src_subdir).unwrap() {
        let entry = entry.unwrap();
        let dest = cpp_dir.join(entry.file_name());
        fs::rename(entry.path(), dest).unwrap();
    }
    fs::remove_dir(&src_subdir).unwrap();

    // Step 3: Apply the patch.
    let patch_path = manifest_dir.join("nix/maxbin2-cpp-ffi.patch");
    // The patch has paths like a/src/EManager.cpp, so we need to adjust.
    // We'll apply with -p2 to strip "a/src/" down to the filename.
    let patch_output = Command::new("patch")
        .args(["-p2", "--no-backup-if-mismatch", "-d"])
        .arg(cpp_dir.to_str().unwrap())
        .arg("-i")
        .arg(patch_path.to_str().unwrap())
        .output()
        .expect("failed to run patch");
    if !patch_output.status.success() {
        panic!(
            "patch failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&patch_output.stdout),
            String::from_utf8_lossy(&patch_output.stderr)
        );
    }

    // Step 4: Copy FFI wrapper files alongside the upstream source.
    let ffi_dir = manifest_dir.join("vendor/ffi");
    for wrapper in FFI_WRAPPERS {
        let src = ffi_dir.join(wrapper);
        let dest = cpp_dir.join(wrapper);
        fs::copy(&src, &dest).unwrap_or_else(|e| {
            panic!("failed to copy FFI wrapper {}: {}", wrapper, e);
        });
    }

    // Step 5: Compile everything into a static library.
    // Optional: unity build — compile all C++ into one translation unit.
    // This simulates what LTO would achieve (cross-function inlining) without
    // requiring toolchain support for cross-language LTO.
    // Set MAXBIN2_CPP_LTO=1 to enable.
    println!("cargo:rerun-if-env-changed=MAXBIN2_CPP_LTO");
    let unity = env::var("MAXBIN2_CPP_LTO").as_deref() == Ok("1");

    let mut build = cc::Build::new();
    build.cpp(true).std("c++11").warnings(false).include(&cpp_dir);

    if unity {
        eprintln!("build.rs: unity build enabled (all C++ in one translation unit)");
        let unity_path = cpp_dir.join("unity.cpp");
        let mut contents = String::new();
        for file in UPSTREAM_CPP {
            if file.ends_with(".cpp") {
                contents.push_str(&format!("#include \"{}\"\n", file));
            }
        }
        for wrapper in FFI_WRAPPERS {
            contents.push_str(&format!("#include \"{}\"\n", wrapper));
        }
        fs::write(&unity_path, &contents).expect("failed to write unity.cpp");
        build.file(&unity_path);
    } else {
        for file in UPSTREAM_CPP {
            if file.ends_with(".cpp") {
                build.file(cpp_dir.join(file));
            }
        }
        for wrapper in FFI_WRAPPERS {
            build.file(cpp_dir.join(wrapper));
        }
    }

    build.compile("maxbin2_ffi");
}
