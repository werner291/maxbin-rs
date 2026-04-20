//! Path management for pipeline output and intermediate files.
//!
//! Two concerns, kept separate:
//! - `OutputDir`: where final results go (bins, summary, marker, noclass, log).
//! - `WorkDir`: where generated intermediates go (filtered contigs, SAM, indices,
//!   HMMER output, seeds, recursive binning artifacts). Auto-cleaned on success
//!   unless `--keep-intermediates` is passed.

use std::path::{Path, PathBuf};

/// Final output directory — bins, summary, marker, noclass, log.
pub struct OutputDir {
    dir: PathBuf,
}

impl OutputDir {
    pub fn new(path: Option<PathBuf>) -> Self {
        Self {
            dir: path.unwrap_or_else(|| PathBuf::from("output")),
        }
    }

    /// Create the directory if it doesn't exist. Error if it exists and is
    /// non-empty, unless `force` is true.
    pub fn ensure(&self, force: bool) -> Result<(), String> {
        if self.dir.exists() {
            if !self.dir.is_dir() {
                return Err(format!(
                    "Output path {} exists but is not a directory",
                    self.dir.display()
                ));
            }
            let is_empty = self
                .dir
                .read_dir()
                .map_err(|e| format!("Can't read output dir: {e}"))?
                .next()
                .is_none();
            if !is_empty && !force {
                return Err(format!(
                    "Output directory {} is not empty. Use --force-overwrite to write anyway.",
                    self.dir.display()
                ));
            }
        } else {
            std::fs::create_dir_all(&self.dir)
                .map_err(|e| format!("Can't create output dir {}: {e}", self.dir.display()))?;
        }
        Ok(())
    }

    pub fn bin(&self, n: usize) -> PathBuf {
        self.dir.join(format!("{n:03}.fasta"))
    }

    pub fn summary(&self) -> PathBuf {
        self.dir.join("summary")
    }

    pub fn noclass(&self) -> PathBuf {
        self.dir.join("noclass")
    }

    pub fn marker(&self) -> PathBuf {
        self.dir.join("marker")
    }

    pub fn timing_json(&self) -> PathBuf {
        self.dir.join("timing.json")
    }

    pub fn log(&self) -> PathBuf {
        self.dir.join("log")
    }

    pub fn dir(&self) -> &Path {
        &self.dir
    }
}

/// Work directory for intermediate files. Owns cleanup lifecycle.
pub struct WorkDir {
    dir: PathBuf,
    keep: bool,
    /// When Some, dropping this removes the directory automatically.
    /// We call `keep()` to prevent deletion when we want to preserve.
    _tempdir: Option<tempfile::TempDir>,
}

impl WorkDir {
    /// Create a new work directory.
    ///
    /// - `keep=false`, `parent=None`: temp dir under `$TMPDIR`, auto-cleaned
    /// - `keep=true`, `parent=None`: under `./intermediates/`, preserved
    /// - `parent=Some(p)`: under `p/`, cleanup depends on `keep`
    pub fn new(keep: bool, parent: Option<&Path>) -> Result<Self, String> {
        match parent {
            Some(p) => {
                std::fs::create_dir_all(p)
                    .map_err(|e| format!("Can't create work-dir parent {}: {e}", p.display()))?;
                let td = tempfile::Builder::new()
                    .prefix("maxbin-rs-")
                    .tempdir_in(p)
                    .map_err(|e| format!("Can't create work dir in {}: {e}", p.display()))?;
                let dir = td.path().to_path_buf();
                eprintln!("Work directory: {}", dir.display());
                Ok(Self {
                    dir,
                    keep,
                    _tempdir: Some(td),
                })
            }
            None if keep => {
                let intermediates = PathBuf::from("intermediates");
                std::fs::create_dir_all(&intermediates).map_err(|e| {
                    format!(
                        "Can't create intermediates dir {}: {e}",
                        intermediates.display()
                    )
                })?;
                let td = tempfile::Builder::new()
                    .prefix("maxbin-rs-")
                    .tempdir_in(&intermediates)
                    .map_err(|e| {
                        format!("Can't create work dir in {}: {e}", intermediates.display())
                    })?;
                let dir = td.path().to_path_buf();
                eprintln!("Work directory: {}", dir.display());
                // Immediately preserve — we know we want to keep it
                let td = td.keep();
                Ok(Self {
                    dir: td,
                    keep,
                    _tempdir: None,
                })
            }
            None => {
                let td = tempfile::Builder::new()
                    .prefix("maxbin-rs-")
                    .tempdir()
                    .map_err(|e| format!("Can't create temp work dir: {e}"))?;
                let dir = td.path().to_path_buf();
                Ok(Self {
                    dir,
                    keep,
                    _tempdir: Some(td),
                })
            }
        }
    }

    /// Join a filename onto the work directory.
    pub fn path(&self, name: &str) -> PathBuf {
        self.dir.join(name)
    }

    /// Work-dir-rooted prefix string for emanager/external tools that take
    /// a string prefix and append their own suffixes.
    pub fn prefix(&self, name: &str) -> String {
        self.dir.join(name).to_string_lossy().into_owned()
    }

    // Convenience methods for well-known intermediates

    pub fn filtered_contigs(&self) -> PathBuf {
        self.path("filtered.fasta")
    }

    pub fn tooshort(&self) -> PathBuf {
        self.path("tooshort.fasta")
    }

    pub fn bowtie2_index_prefix(&self) -> String {
        self.prefix("bt2idx")
    }

    pub fn sam(&self, i: usize) -> PathBuf {
        self.path(&format!("sam{i}"))
    }

    pub fn abund_from_reads(&self, i: usize) -> PathBuf {
        self.path(&format!("reads.abund{i}"))
    }

    pub fn fgs_output_prefix(&self) -> String {
        self.prefix("genes")
    }

    pub fn fgs_faa(&self) -> PathBuf {
        // run_fraggenescan appends ".frag" to the prefix, then FGS adds ".faa"
        self.path("genes.frag.faa")
    }

    pub fn hmmout(&self) -> PathBuf {
        self.path("markers.hmmout")
    }

    /// Clean up on success: remove the work dir unless --keep-intermediates.
    pub fn cleanup(self) {
        if self.keep {
            eprintln!("Intermediates kept at: {}", self.dir.display());
            // Disengage TempDir auto-delete
            if let Some(td) = self._tempdir {
                let _ = td.keep();
            }
        }
        // If !keep and _tempdir is Some, dropping self removes it.
        // If !keep and _tempdir is None (shouldn't happen in normal flow),
        // do an explicit remove.
        else if self._tempdir.is_none() {
            let _ = std::fs::remove_dir_all(&self.dir);
        }
    }

    /// Preserve on error: keep the work dir and print its path.
    pub fn preserve(self) -> PathBuf {
        if let Some(td) = self._tempdir {
            let path = td.keep();
            eprintln!("Work dir preserved at: {}", path.display());
            path
        } else {
            eprintln!("Work dir preserved at: {}", self.dir.display());
            self.dir
        }
    }
}
