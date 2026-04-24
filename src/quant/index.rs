//! Decoy-aware salmon index builder with SHA256-keyed cache.
//!
//! The index is the single most expensive artifact in a quant run and
//! the one most often shared across runs on the same genome. Caching
//! it by a content-derived key avoids rebuilding when the user re-runs
//! with the same inputs, and invalidates cleanly the moment any input
//! changes.
//!
//! Cache key: `sha256(sha256(cds_fa) || sha256(genome_fa) || salmon_version || k)`.
//! The nested hashing means we don't re-hash the genome every lookup —
//! we hash it once, then combine short strings.
//!
//! Cache directory is resolved in precedence order (first hit wins):
//!
//!   1. `--index-cache <dir>` CLI flag (caller passes as `override_dir`)
//!   2. `MYCONOTE_INDEX_CACHE` environment variable
//!   3. `$XDG_CACHE_HOME/myconote/salmon_index/`
//!   4. `~/.cache/myconote/salmon_index/` (Linux XDG default when unset)
//!   5. `~/.myconote/salmon_index/` (macOS fallback)
//!
//! Silent fallback to `/tmp` on a read-only cache is deliberately NOT
//! done — a silent location shuffle breaks the reproducibility claim
//! of the bundle. We error out with a clear message naming the path
//! that failed and suggesting `--index-cache <writable-dir>` or
//! `MYCONOTE_INDEX_CACHE=<writable-dir>`.

use crate::utils::error::{MycoNoteError, Result};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use super::bundle::sha256_file;

// ─────────────────────────────────────────────────────────────────────────────
// Cache directory resolution (pure, unit-testable)
// ─────────────────────────────────────────────────────────────────────────────

/// Inputs to the cache-dir resolver. Gather them once at run start so
/// tests can pass synthetic values without poking at real env vars or
/// the real HOME directory.
#[derive(Debug, Clone, Default)]
pub struct CacheEnv {
    /// `--index-cache` CLI flag, if the user passed one.
    pub cli_override: Option<PathBuf>,
    /// Snapshot of `MYCONOTE_INDEX_CACHE` + `XDG_CACHE_HOME` + `HOME`.
    pub env: HashMap<String, String>,
}

impl CacheEnv {
    /// Build a `CacheEnv` from the real process state.
    pub fn from_process(cli_override: Option<PathBuf>) -> Self {
        let mut env = HashMap::new();
        for key in ["MYCONOTE_INDEX_CACHE", "XDG_CACHE_HOME", "HOME"] {
            if let Ok(v) = std::env::var(key) {
                env.insert(key.to_string(), v);
            }
        }
        Self { cli_override, env }
    }
}

/// Apply the precedence chain and return the cache root path. Does
/// NOT create the directory or check writability — that's the
/// caller's job (so we can unit-test resolution without filesystem
/// effects).
pub fn resolve_cache_root(env: &CacheEnv) -> Result<PathBuf> {
    // 1. CLI override wins absolutely.
    if let Some(ref p) = env.cli_override {
        return Ok(p.clone());
    }

    // 2. MYCONOTE_INDEX_CACHE env var.
    if let Some(v) = env.env.get("MYCONOTE_INDEX_CACHE") {
        if !v.is_empty() {
            return Ok(PathBuf::from(v));
        }
    }

    // 3. $XDG_CACHE_HOME/myconote/salmon_index/
    if let Some(v) = env.env.get("XDG_CACHE_HOME") {
        if !v.is_empty() {
            return Ok(PathBuf::from(v).join("myconote").join("salmon_index"));
        }
    }

    // 4 + 5. Fall back to HOME-relative paths.
    let home = env.env.get("HOME").cloned().ok_or_else(|| {
        MycoNoteError::QuantSheet(
            "cannot resolve index cache directory: no --index-cache, \
             no MYCONOTE_INDEX_CACHE, no XDG_CACHE_HOME, no HOME"
                .to_string(),
        )
    })?;

    // On platforms that follow XDG (Linux, many HPC sites), prefer
    // the ~/.cache/myconote/ path. macOS users who prefer ~/.myconote/
    // get the legacy fallback. We don't detect OS here — both paths
    // are valid and users can force either via the env var.
    Ok(PathBuf::from(&home)
        .join(".cache")
        .join("myconote")
        .join("salmon_index"))
}

/// Ensure the cache root exists and is writable, creating it if
/// missing. Fails with a precise error message that names the path
/// and suggests the override routes — so the user is never left
/// wondering *which* directory was the problem.
pub fn ensure_cache_root(root: &Path) -> Result<()> {
    if let Err(e) = fs::create_dir_all(root) {
        return Err(MycoNoteError::QuantSheet(format!(
            "index cache directory {} is not writable ({}). \
             Set MYCONOTE_INDEX_CACHE=<writable dir> or pass \
             --index-cache <writable dir>.",
            root.display(),
            e
        )));
    }
    // Probe writability via a canary file. `create_dir_all` succeeds
    // on an already-existing read-only directory on some filesystems,
    // so a stat isn't enough. Write + delete catches the real failure
    // mode (NFS mounts with no write permission).
    let canary = root.join(".myconote_write_probe");
    match File::create(&canary) {
        Ok(_) => {
            let _ = fs::remove_file(&canary);
            Ok(())
        }
        Err(e) => Err(MycoNoteError::QuantSheet(format!(
            "index cache directory {} is not writable ({}). \
             Set MYCONOTE_INDEX_CACHE=<writable dir> or pass \
             --index-cache <writable dir>.",
            root.display(),
            e
        ))),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cache key (content hash over cds + genome + salmon version + k)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the SHA256-hex cache key for a salmon index. Hashes each
/// input FASTA once (streaming), then combines the two hashes with
/// the salmon version string and the k-mer length. Two runs produce
/// the same key iff all four components are byte-identical.
pub fn compute_cache_key(
    cds_fa: &Path,
    genome_fa: &Path,
    salmon_version: &str,
    k: usize,
) -> Result<String> {
    let cds_hash = sha256_file(cds_fa)?;
    let genome_hash = sha256_file(genome_fa)?;
    let mut hasher = Sha256::new();
    hasher.update(cds_hash.as_bytes());
    hasher.update(b"|");
    hasher.update(genome_hash.as_bytes());
    hasher.update(b"|");
    hasher.update(salmon_version.as_bytes());
    hasher.update(b"|");
    hasher.update(k.to_string().as_bytes());
    Ok(format!("{:x}", hasher.finalize()))
}

// ─────────────────────────────────────────────────────────────────────────────
// Decoys.txt (genome contig names, one per line, for salmon's --decoys)
// ─────────────────────────────────────────────────────────────────────────────

/// Read a FASTA and write a newline-separated list of contig names
/// (first whitespace-delimited token of each header line) to `output`.
/// Returns the number of contigs written. This is the input salmon
/// expects on `--decoys <file>`.
///
/// The FASTA reader is intentionally minimal: we never materialize
/// sequence data in memory, we only look at header lines starting
/// with `>`. Handles gzipped FASTA? Not here — salmon itself does not
/// accept `.fa.gz` for the decoy list file anyway, and the caller
/// supplies the same genome path that would feed `salmon index`.
pub fn write_decoys_txt(genome_fa: &Path, output: &Path) -> Result<usize> {
    let f = File::open(genome_fa).map_err(|e| MycoNoteError::QuantTool {
        tool: "salmon-index".to_string(),
        message: format!("cannot open genome {}: {}", genome_fa.display(), e),
    })?;
    let reader = BufReader::new(f);

    let mut out = File::create(output)?;
    let mut n = 0usize;

    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            if name.is_empty() {
                return Err(MycoNoteError::QuantTool {
                    tool: "salmon-index".to_string(),
                    message: format!(
                        "{}: encountered a FASTA header with no name",
                        genome_fa.display()
                    ),
                });
            }
            writeln!(out, "{}", name)?;
            n += 1;
        }
    }

    if n == 0 {
        return Err(MycoNoteError::QuantTool {
            tool: "salmon-index".to_string(),
            message: format!(
                "{}: no FASTA records found; cannot use as decoy",
                genome_fa.display()
            ),
        });
    }

    Ok(n)
}

/// Count FASTA records (header lines) in a file. Used for populating
/// the bundle's `index.targets` field from a CDS FASTA.
pub fn count_fasta_records(fa: &Path) -> Result<usize> {
    let f = File::open(fa)?;
    let reader = BufReader::new(f);
    let mut n = 0usize;
    for line in reader.lines() {
        if line?.starts_with('>') {
            n += 1;
        }
    }
    Ok(n)
}

// ─────────────────────────────────────────────────────────────────────────────
// Salmon version detection + index building (subprocess boundary)
// ─────────────────────────────────────────────────────────────────────────────

/// Ask salmon for its version string. Passed into the cache key so
/// runs that straddle a salmon upgrade rebuild rather than silently
/// reusing an incompatible index.
///
/// Salmon prints to stderr in the form `salmon 1.11.4` — we take the
/// second whitespace token. If the output shape changes in a future
/// salmon release this function fails loudly rather than returning a
/// bogus version string.
pub fn detect_salmon_version(salmon_bin: &str) -> Result<String> {
    let output = duct::cmd!(salmon_bin, "--version")
        .stderr_to_stdout()
        .read()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "salmon".to_string(),
            message: format!("running '{} --version': {}", salmon_bin, e),
        })?;
    parse_salmon_version(&output).ok_or_else(|| MycoNoteError::QuantTool {
        tool: "salmon".to_string(),
        message: format!(
            "unexpected `salmon --version` output (expected 'salmon X.Y.Z'): {}",
            output.trim()
        ),
    })
}

fn parse_salmon_version(output: &str) -> Option<String> {
    for line in output.lines() {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() >= 2 && toks[0].eq_ignore_ascii_case("salmon") {
            return Some(toks[1].to_string());
        }
    }
    None
}

/// Summary of a resolved index — either freshly built or found in
/// cache. The caller records this into the reproducibility bundle.
#[derive(Debug)]
pub struct IndexResult {
    pub path: PathBuf,
    pub cache_key: String,
    pub k: usize,
    pub targets: usize,
    pub decoys: usize,
    pub salmon_version: String,
    pub cached: bool,
}

/// Inputs to `build_or_reuse_index`. Kept in its own struct so the
/// dispatcher can assemble once and re-use across tests.
#[derive(Debug, Clone)]
pub struct IndexSpec {
    pub cds_fa: PathBuf,
    pub genome_fa: PathBuf,
    pub k: usize,
    pub threads: usize,
    /// Path to the salmon binary (from `which` lookup or `--salmon`
    /// override). Defaults to `"salmon"` when None.
    pub salmon_bin: Option<String>,
}

impl IndexSpec {
    pub fn salmon(&self) -> &str {
        self.salmon_bin.as_deref().unwrap_or("salmon")
    }
}

/// Build a decoy-aware salmon index under `cache_root/<hash>/`, or
/// reuse an existing one with the same hash. Returns the result
/// struct. Requires `salmon` on PATH (or supplied via `spec.salmon_bin`).
pub fn build_or_reuse_index(spec: &IndexSpec, cache_root: &Path) -> Result<IndexResult> {
    ensure_cache_root(cache_root)?;

    let salmon_version = detect_salmon_version(spec.salmon())?;
    let key = compute_cache_key(&spec.cds_fa, &spec.genome_fa, &salmon_version, spec.k)?;
    let index_dir = cache_root.join(&key);

    // `salmon index` writes a number of files into its output dir; we
    // treat `info.json` as the sentinel — salmon produces it only on
    // successful completion, so a partial index from an interrupted
    // run won't fool the cache.
    let sentinel = index_dir.join("info.json");
    let targets = count_fasta_records(&spec.cds_fa)?;
    let decoys_count: usize;

    if sentinel.exists() {
        // Cache hit. Decoys count is best recovered from the stored
        // decoys.txt so the bundle field is accurate even on a reused
        // index.
        decoys_count = count_lines(&index_dir.join("decoys.txt")).unwrap_or(0);
        return Ok(IndexResult {
            path: index_dir,
            cache_key: key,
            k: spec.k,
            targets,
            decoys: decoys_count,
            salmon_version,
            cached: true,
        });
    }

    // Build. Staging dir avoids leaving half-written state under the
    // canonical key if salmon aborts; rename on success.
    fs::create_dir_all(&index_dir)?;
    let staging = cache_root.join(format!("{}.staging", key));
    let _ = fs::remove_dir_all(&staging);
    fs::create_dir_all(&staging)?;

    let decoys_txt = staging.join("decoys.txt");
    let n_decoys = write_decoys_txt(&spec.genome_fa, &decoys_txt)?;

    // Concatenate cds + genome into a single target file for
    // salmon — that's how the decoy-aware index works in salmon 1.x.
    let combined_fa = staging.join("targets_with_decoys.fa");
    concat_fastas(&spec.cds_fa, &spec.genome_fa, &combined_fa)?;

    let status = duct::cmd!(
        spec.salmon(),
        "index",
        "--threads",
        spec.threads.to_string(),
        "-k",
        spec.k.to_string(),
        "-t",
        combined_fa.to_string_lossy().to_string(),
        "-d",
        decoys_txt.to_string_lossy().to_string(),
        "-i",
        staging.to_string_lossy().to_string()
    )
    .stderr_to_stdout()
    .unchecked()
    .run()
    .map_err(|e| MycoNoteError::QuantTool {
        tool: "salmon-index".to_string(),
        message: format!("salmon index spawn failed: {}", e),
    })?;

    if !status.status.success() {
        return Err(MycoNoteError::QuantTool {
            tool: "salmon-index".to_string(),
            message: format!(
                "salmon index exited {}: output:\n{}",
                status.status,
                String::from_utf8_lossy(&status.stdout)
            ),
        });
    }

    // Salmon emitted its index into `staging/`; move that into place
    // under the cache key. Retain decoys.txt alongside for the
    // cache-hit path above.
    let _ = fs::remove_dir_all(&index_dir);
    fs::rename(&staging, &index_dir).map_err(|e| {
        MycoNoteError::QuantSheet(format!(
            "could not finalize index cache entry {}: {}",
            index_dir.display(),
            e
        ))
    })?;

    Ok(IndexResult {
        path: index_dir,
        cache_key: key,
        k: spec.k,
        targets,
        decoys: n_decoys,
        salmon_version,
        cached: false,
    })
}

/// Concatenate two FASTA files into one, streaming in 64 KiB chunks.
/// No attempt at header rewriting — salmon matches decoy names via
/// the decoys.txt list, not via header transformation.
fn concat_fastas(a: &Path, b: &Path, output: &Path) -> Result<()> {
    let mut out = File::create(output)?;
    for src in [a, b] {
        let mut f = File::open(src)?;
        std::io::copy(&mut f, &mut out)?;
    }
    Ok(())
}

fn count_lines(p: &Path) -> Result<usize> {
    let f = File::open(p)?;
    let reader = BufReader::new(f);
    let mut n = 0;
    for line in reader.lines() {
        if !line?.trim().is_empty() {
            n += 1;
        }
    }
    Ok(n)
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn env_with(pairs: &[(&str, &str)]) -> CacheEnv {
        let mut env = HashMap::new();
        for (k, v) in pairs {
            env.insert(k.to_string(), v.to_string());
        }
        CacheEnv {
            cli_override: None,
            env,
        }
    }

    // ── resolve_cache_root precedence ────────────────────────────────────────

    #[test]
    fn cli_override_beats_everything() {
        let e = CacheEnv {
            cli_override: Some(PathBuf::from("/opt/cache")),
            env: [
                ("MYCONOTE_INDEX_CACHE", "/env/cache"),
                ("XDG_CACHE_HOME", "/xdg"),
                ("HOME", "/home/ben"),
            ]
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        };
        assert_eq!(resolve_cache_root(&e).unwrap(), PathBuf::from("/opt/cache"));
    }

    #[test]
    fn env_var_beats_xdg_and_home() {
        let e = env_with(&[
            ("MYCONOTE_INDEX_CACHE", "/env/cache"),
            ("XDG_CACHE_HOME", "/xdg"),
            ("HOME", "/home/ben"),
        ]);
        assert_eq!(resolve_cache_root(&e).unwrap(), PathBuf::from("/env/cache"));
    }

    #[test]
    fn xdg_used_when_set() {
        let e = env_with(&[("XDG_CACHE_HOME", "/xdg"), ("HOME", "/home/ben")]);
        assert_eq!(
            resolve_cache_root(&e).unwrap(),
            PathBuf::from("/xdg/myconote/salmon_index")
        );
    }

    #[test]
    fn home_fallback_used_when_xdg_unset() {
        let e = env_with(&[("HOME", "/home/ben")]);
        assert_eq!(
            resolve_cache_root(&e).unwrap(),
            PathBuf::from("/home/ben/.cache/myconote/salmon_index")
        );
    }

    #[test]
    fn fails_when_no_env_at_all() {
        let e = env_with(&[]);
        let err = resolve_cache_root(&e).unwrap_err();
        assert!(format!("{err}").contains("no --index-cache"));
    }

    #[test]
    fn empty_env_var_is_skipped() {
        // An empty MYCONOTE_INDEX_CACHE should fall through to XDG
        // rather than resolving to the empty path.
        let e = env_with(&[
            ("MYCONOTE_INDEX_CACHE", ""),
            ("XDG_CACHE_HOME", "/xdg"),
            ("HOME", "/h"),
        ]);
        assert_eq!(
            resolve_cache_root(&e).unwrap(),
            PathBuf::from("/xdg/myconote/salmon_index")
        );
    }

    // ── ensure_cache_root writability probe ─────────────────────────────────

    #[test]
    fn ensure_cache_root_creates_new_dir() {
        let tmp = TempDir::new().unwrap();
        let target = tmp.path().join("a").join("b").join("c");
        ensure_cache_root(&target).unwrap();
        assert!(target.is_dir());
    }

    #[cfg(unix)]
    #[test]
    fn ensure_cache_root_rejects_readonly() {
        use std::os::unix::fs::PermissionsExt;
        let tmp = TempDir::new().unwrap();
        let ro = tmp.path().join("readonly");
        fs::create_dir(&ro).unwrap();
        fs::set_permissions(&ro, fs::Permissions::from_mode(0o555)).unwrap();

        // Try to use a subdir inside the read-only dir.
        let target = ro.join("sub");
        let err = ensure_cache_root(&target).unwrap_err();
        assert!(format!("{err}").contains("not writable"), "got: {err}");

        // Restore perms so tempdir cleanup can succeed.
        fs::set_permissions(&ro, fs::Permissions::from_mode(0o755)).unwrap();
    }

    // ── cache-key determinism ────────────────────────────────────────────────

    #[test]
    fn cache_key_deterministic() {
        let tmp = TempDir::new().unwrap();
        let cds = tmp.path().join("cds.fa");
        let genome = tmp.path().join("genome.fa");
        fs::write(&cds, ">t1\nACGT\n").unwrap();
        fs::write(&genome, ">chr1\nACGTACGT\n").unwrap();

        let a = compute_cache_key(&cds, &genome, "salmon 1.11.4", 31).unwrap();
        let b = compute_cache_key(&cds, &genome, "salmon 1.11.4", 31).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.len(), 64, "should be hex-encoded SHA256");
    }

    #[test]
    fn cache_key_changes_with_k() {
        let tmp = TempDir::new().unwrap();
        let cds = tmp.path().join("cds.fa");
        let genome = tmp.path().join("genome.fa");
        fs::write(&cds, ">t1\nACGT\n").unwrap();
        fs::write(&genome, ">chr1\nACGTACGT\n").unwrap();

        let k31 = compute_cache_key(&cds, &genome, "salmon 1.11.4", 31).unwrap();
        let k25 = compute_cache_key(&cds, &genome, "salmon 1.11.4", 25).unwrap();
        assert_ne!(k31, k25);
    }

    #[test]
    fn cache_key_changes_with_salmon_version() {
        let tmp = TempDir::new().unwrap();
        let cds = tmp.path().join("cds.fa");
        let genome = tmp.path().join("genome.fa");
        fs::write(&cds, ">t1\nACGT\n").unwrap();
        fs::write(&genome, ">chr1\nACGTACGT\n").unwrap();

        let old = compute_cache_key(&cds, &genome, "salmon 1.10.2", 31).unwrap();
        let new_ = compute_cache_key(&cds, &genome, "salmon 1.11.4", 31).unwrap();
        assert_ne!(old, new_);
    }

    #[test]
    fn cache_key_changes_with_cds_content() {
        let tmp = TempDir::new().unwrap();
        let cds = tmp.path().join("cds.fa");
        let genome = tmp.path().join("genome.fa");
        fs::write(&genome, ">chr1\nACGTACGT\n").unwrap();

        fs::write(&cds, ">t1\nACGT\n").unwrap();
        let a = compute_cache_key(&cds, &genome, "salmon 1.11.4", 31).unwrap();

        fs::write(&cds, ">t1\nTTTT\n").unwrap();
        let b = compute_cache_key(&cds, &genome, "salmon 1.11.4", 31).unwrap();
        assert_ne!(a, b);
    }

    // ── decoys.txt + FASTA record counting ───────────────────────────────────

    #[test]
    fn write_decoys_extracts_contig_names() {
        let tmp = TempDir::new().unwrap();
        let genome = tmp.path().join("genome.fa");
        // Mix names with + without comments after the name.
        fs::write(
            &genome,
            ">chrI dna primary\nACGT\n\
             >chrII\nACGT\n\
             >mtDNA some description\nTTTT\n",
        )
        .unwrap();

        let decoys = tmp.path().join("decoys.txt");
        let n = write_decoys_txt(&genome, &decoys).unwrap();
        assert_eq!(n, 3);

        let body = fs::read_to_string(&decoys).unwrap();
        assert_eq!(body, "chrI\nchrII\nmtDNA\n");
    }

    #[test]
    fn write_decoys_rejects_empty_fasta() {
        let tmp = TempDir::new().unwrap();
        let genome = tmp.path().join("empty.fa");
        fs::write(&genome, "").unwrap();
        let decoys = tmp.path().join("decoys.txt");
        let err = write_decoys_txt(&genome, &decoys).unwrap_err();
        assert!(format!("{err}").contains("no FASTA records"));
    }

    #[test]
    fn count_records_matches_header_count() {
        let tmp = TempDir::new().unwrap();
        let fa = tmp.path().join("x.fa");
        fs::write(&fa, ">a\nACGT\n>b\nACGT\n>c\nACGT\n").unwrap();
        assert_eq!(count_fasta_records(&fa).unwrap(), 3);
    }

    // ── salmon version string parsing ────────────────────────────────────────

    #[test]
    fn parse_salmon_version_ok() {
        assert_eq!(
            parse_salmon_version("salmon 1.11.4\n").as_deref(),
            Some("1.11.4")
        );
        // Version sometimes has trailing build info on salmon CI builds.
        assert_eq!(
            parse_salmon_version("salmon 1.10.2-devbuild\n").as_deref(),
            Some("1.10.2-devbuild")
        );
        // Case-insensitive on the "salmon" prefix.
        assert_eq!(
            parse_salmon_version("SALMON 1.0.0\n").as_deref(),
            Some("1.0.0")
        );
    }

    #[test]
    fn parse_salmon_version_none_when_prefix_missing() {
        assert_eq!(parse_salmon_version("sockeye 1.0"), None);
        assert_eq!(parse_salmon_version(""), None);
    }

    #[test]
    fn parse_salmon_version_tolerates_extra_lines() {
        // salmon sometimes emits a banner before the version line.
        let out = "Loading…\nsalmon 1.11.4\nDone.\n";
        assert_eq!(parse_salmon_version(out).as_deref(), Some("1.11.4"));
    }

    // ── live salmon smoke test (ignored unless salmon is on PATH) ──────────
    //
    // Run with `cargo test --release -- --ignored salmon_index` after
    // `conda activate fast_myco_env` (or any env where salmon resolves
    // via `which`). Exercises the real subprocess boundary so we catch
    // version skew the moment salmon changes its flag surface.

    #[test]
    #[ignore]
    fn salmon_index_live_build_and_cache_hit() {
        if which::which("salmon").is_err() {
            eprintln!("salmon not on PATH — skipping live index test");
            return;
        }

        let tmp = TempDir::new().unwrap();
        let cds = tmp.path().join("cds.fa");
        let genome = tmp.path().join("genome.fa");
        // Two short "transcripts" + two "chromosomes". Sequences are
        // long enough for salmon k=21; we won't run quant against this
        // index here, only verify the build.
        fs::write(
            &cds,
            ">t1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
             >t2\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA\n",
        )
        .unwrap();
        fs::write(
            &genome,
            ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
             >chr2\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA\n",
        )
        .unwrap();

        let cache_root = tmp.path().join("cache");
        let spec = IndexSpec {
            cds_fa: cds,
            genome_fa: genome,
            k: 21,
            threads: 2,
            salmon_bin: None,
        };

        let r1 = build_or_reuse_index(&spec, &cache_root).unwrap();
        assert!(!r1.cached, "first build should be fresh");
        assert!(
            r1.path.join("info.json").is_file(),
            "salmon sentinel missing"
        );
        assert_eq!(r1.targets, 2);
        assert_eq!(r1.decoys, 2);

        // Second call with same inputs → cache hit, no rebuild.
        let r2 = build_or_reuse_index(&spec, &cache_root).unwrap();
        assert!(r2.cached, "second call should reuse cache");
        assert_eq!(r1.cache_key, r2.cache_key);
        assert_eq!(r1.path, r2.path);
    }
}
