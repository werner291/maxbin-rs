//! maxbin-rs-vis — browser teaching/visualization frontend for the maxbin-rs EM core.
//!
//! Scaffold only. This crate is the single place that will know about WASM and
//! the browser; the EM core stays plain, host-agnostic Rust and gets wired in
//! here later. For now this is an empty shell so the workspace member and the
//! browser entry point (`web/index.html`) exist and build.

/// Placeholder so the crate compiles and has something to test.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_set() {
        assert!(!super::version().is_empty());
    }
}
