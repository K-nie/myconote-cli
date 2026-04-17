/// IQ-TREE phylogenetic tree inference wrapper
///
/// Builds maximum-likelihood phylogenetic trees from a multiple-sequence
/// alignment using IQ-TREE 2.  Supports:
///   - Automatic model selection via ModelFinder Plus (MFP)
///   - Ultrafast bootstrap (UFBoot, -B) and standard bootstrap (-b)
///   - Partitioned analysis for multi-gene/multi-locus datasets
///   - Parallel execution via OpenMP threads
///
/// IQ-TREE must be installed and available in PATH:
///   conda install -c bioconda iqtree
///
/// Output files (written to the same directory as the input alignment):
///   <prefix>.treefile   — best-fit ML tree in Newick format
///   <prefix>.iqtree     — full IQ-TREE run report
///   <prefix>.log        — console log
///   <prefix>.ckp.gz     — checkpoint (allows resuming interrupted runs)
use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct PhylogenyConfig {
    /// Substitution model, or "MFP" for automatic ModelFinder selection
    pub model: String,
    /// Number of ultrafast bootstrap replicates (0 = no bootstrap)
    pub bootstrap: usize,
    /// Partition file for multi-gene analysis (optional)
    pub partition: Option<PathBuf>,
    /// Number of parallel threads
    pub threads: usize,
    /// Output prefix (defaults to input filename stem)
    pub prefix: Option<String>,
    /// Extra IQ-TREE arguments passed verbatim
    pub extra_args: Vec<String>,
    /// Force overwrite of an existing run with the same prefix. When false
    /// and a previous `.treefile` exists for this prefix, build_tree returns
    /// early rather than overwriting — the prior version always passed
    /// `--redo` and could silently destroy a long-running result.
    pub force: bool,
}

impl Default for PhylogenyConfig {
    fn default() -> Self {
        Self {
            model: "MFP".to_string(),
            bootstrap: 1000,
            partition: None,
            threads: 4,
            prefix: None,
            extra_args: Vec::new(),
            force: false,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Availability check
// ─────────────────────────────────────────────────────────────────────────────

pub fn iqtree_available() -> bool {
    // IQ-TREE 2 ships as "iqtree2"; IQ-TREE 1 as "iqtree". Use which::which
    // for consistent PATH resolution across macOS/Linux/BSD — the prior
    // Command::new("which") path breaks on systems where /bin/which isn't
    // installed or when IPC is sandboxed.
    which::which("iqtree2").is_ok() || which::which("iqtree").is_ok()
}

/// Return the IQ-TREE binary name that is present on PATH, preferring v2.
fn iqtree_bin() -> &'static str {
    if which::which("iqtree2").is_ok() {
        "iqtree2"
    } else {
        "iqtree"
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run IQ-TREE on `alignment` with the given `config`.
///
/// Returns the path to the best-fit ML tree (`.treefile`).
pub fn build_tree<P: AsRef<Path>>(alignment: P, config: &PhylogenyConfig) -> Result<PathBuf> {
    let alignment = alignment.as_ref();

    if !alignment.exists() {
        return Err(MycoNoteError::Io(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("Alignment file not found: {}", alignment.display()),
        )));
    }

    if !iqtree_available() {
        return Err(MycoNoteError::ExternalTool(
            "IQ-TREE not found in PATH. Install with: conda install -c bioconda iqtree".to_string(),
        ));
    }

    let bin = iqtree_bin();

    // Derive output prefix from config or input filename stem
    let prefix = match &config.prefix {
        Some(p) => p.clone(),
        None => alignment
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "iqtree_out".to_string()),
    };

    // Prefix handling: if the user passes a path-like prefix ("out/run1" or
    // absolute path), treat it as a full path (IQ-TREE accepts any prefix).
    // Otherwise — bare name like "my_tree" — join with the alignment's
    // parent dir so outputs land next to the input.
    // The previous unconditional join() turned "out/run1" into
    // "<alignment_dir>/out/run1" which then failed at IQ-TREE's log-file
    // open because the nested "out/" subdir didn't exist.
    let prefix_has_path = prefix.contains(std::path::MAIN_SEPARATOR)
        || PathBuf::from(&prefix).is_absolute();
    let prefix_path = if prefix_has_path {
        PathBuf::from(&prefix)
    } else {
        alignment
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(&prefix)
    };

    // Make sure the prefix's parent directory actually exists so IQ-TREE
    // can open its log file. Otherwise a user-supplied "out/run1" with no
    // pre-existing "out/" silently fails deep inside iqtree.
    if let Some(parent) = prefix_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(MycoNoteError::Io)?;
        }
    }

    // Guard against silently overwriting a completed run. IQ-TREE can easily
    // take hours on a multi-gene partition; re-running with a stale command
    // used to wipe the .treefile without warning.
    let existing_tree = prefix_path.with_extension("treefile");
    if existing_tree.exists() && !config.force {
        println!(
            "  ✓  {} already exists — skipping IQ-TREE (pass --force to overwrite).",
            existing_tree.display()
        );
        return Ok(existing_tree);
    }

    println!("  Running IQ-TREE ({bin}) on: {}", alignment.display());
    println!(
        "  Model: {}  Bootstrap: {}  Threads: {}",
        config.model, config.bootstrap, config.threads
    );

    let mut cmd = Command::new(bin);
    cmd.arg("-s")
        .arg(alignment)
        .arg("-m")
        .arg(&config.model)
        .arg("--prefix")
        .arg(&prefix_path)
        .arg("-T")
        .arg(config.threads.to_string());
    if config.force {
        cmd.arg("--redo");
    }

    // Partition file
    if let Some(part) = &config.partition {
        cmd.arg("-p").arg(part);
    }

    // Bootstrap
    if config.bootstrap > 0 {
        // UFBoot2 (ultrafast, recommended for large datasets)
        cmd.arg("-B").arg(config.bootstrap.to_string());
        cmd.arg("--alrt").arg("1000"); // SH-aLRT branch support alongside UFBoot
    }

    // Extra user-supplied args
    for arg in &config.extra_args {
        cmd.arg(arg);
    }

    let output = cmd.output().map_err(|e| {
        MycoNoteError::Io(std::io::Error::new(
            e.kind(),
            format!("Failed to launch {bin}: {e}"),
        ))
    })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(MycoNoteError::ExternalTool(format!(
            "IQ-TREE failed (exit {})\n{}",
            output.status.code().unwrap_or(-1),
            stderr.trim()
        )));
    }

    let treefile = prefix_path.with_extension("treefile");
    if treefile.exists() {
        println!("  ✓  Tree written to: {}", treefile.display());
    } else {
        // IQ-TREE sometimes appends `.treefile` without dropping an extension
        // from the prefix — check the literal `<prefix>.treefile` form too.
        let alt = prefix_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(format!("{}.treefile", prefix));
        if alt.exists() {
            return Ok(alt);
        }
        return Err(MycoNoteError::ExternalTool(
            "IQ-TREE finished but treefile not found".to_string(),
        ));
    }

    Ok(treefile)
}

// ─────────────────────────────────────────────────────────────────────────────
// Convenience helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Read a Newick treefile and return its contents as a String.
pub fn read_treefile<P: AsRef<Path>>(path: P) -> Result<String> {
    std::fs::read_to_string(path.as_ref()).map_err(MycoNoteError::Io)
}

/// Check what IQ-TREE version is installed and return the version string.
pub fn iqtree_version() -> Option<String> {
    let bin = if Command::new("which")
        .arg("iqtree2")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
    {
        "iqtree2"
    } else {
        "iqtree"
    };

    Command::new(bin)
        .arg("--version")
        .output()
        .ok()
        .and_then(|o| {
            let out = String::from_utf8_lossy(&o.stdout).to_string()
                + &String::from_utf8_lossy(&o.stderr);
            out.lines()
                .find(|l| {
                    l.to_lowercase().contains("iq-tree") || l.to_lowercase().contains("version")
                })
                .map(|l| l.trim().to_string())
        })
}
