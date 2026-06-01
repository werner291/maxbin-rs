//! Equivalence-testing harness for maxbin-rs.
//!
//! This crate links the original MaxBin2 C++ code as a static library and
//! exposes Rust FFI bindings to it (`original_ffi`). It exists so the core
//! `maxbin-rs` crate can stay free of any C/C++ toolchain dependency (and
//! eventually compile to wasm), while the bit-for-bit equivalence tests
//! against the original implementation keep running here.
//!
//! The C++ sources are extracted from the upstream MaxBin2 tarball and
//! patched for FFI at build time (see `build.rs` and `vendor/`). The
//! `MAXBIN2_SRC_TARBALL` environment variable points at that tarball; the
//! Nix devshell sets it automatically.

pub mod original_ffi;
