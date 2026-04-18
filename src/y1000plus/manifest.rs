//! Persistent manifest for the Y1000+ cache.
//!
//! Lives at `<cache_root>/MANIFEST.toml`. Tracks which subsets are installed,
//! their byte sizes on disk, SHA-256 checksums of the source archives (for
//! re-verification), and the bundle version so upgrades can invalidate old
//! caches cleanly.

use crate::utils::error::{MycoNoteError, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

pub const MANIFEST_FILENAME: &str = "MANIFEST.toml";
/// Bump this when the on-disk layout or extraction logic changes
/// incompatibly; subcommands refuse to use caches built by older versions.
pub const BUNDLE_SCHEMA_VERSION: u32 = 1;

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Manifest {
    /// Schema version used to build the cache.
    #[serde(default)]
    pub schema_version: u32,
    /// Paper citation to keep with the cache for reproducibility.
    #[serde(default = "default_citation")]
    pub citation: String,
    /// Subsets currently installed: key -> installed info.
    #[serde(default)]
    pub installed: BTreeMap<String, InstalledSubset>,
}

fn default_citation() -> String {
    "Opulente DA et al. (2024). Genomic factors shape carbon and nitrogen metabolic niche \
     breadth across Saccharomycotina yeasts. Science 384(6694): eadj4503."
        .to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InstalledSubset {
    /// Bytes on disk after extraction (best-effort — may drift if user edits cache).
    pub on_disk_bytes: u64,
    /// ISO-8601 UTC timestamp of install.
    pub installed_at: String,
    /// For resumable/verifiable re-downloads. Empty if unverified.
    #[serde(default)]
    pub source_checksum_sha256: String,
}

/// Resolve the cache root. Priority: `MYCONOTE_Y1000_CACHE` env var, then
/// `$HOME/.myconote/y1000plus`.
pub fn cache_root() -> PathBuf {
    if let Ok(p) = std::env::var("MYCONOTE_Y1000_CACHE") {
        return PathBuf::from(p);
    }
    let home = std::env::var("HOME").unwrap_or_else(|_| ".".to_string());
    PathBuf::from(home).join(".myconote").join("y1000plus")
}

pub fn manifest_path(root: &Path) -> PathBuf {
    root.join(MANIFEST_FILENAME)
}

/// Load the manifest from disk. Returns an empty (default) manifest if none
/// exists — a fresh install is just an empty manifest that grows over time.
pub fn load(root: &Path) -> Result<Manifest> {
    let path = manifest_path(root);
    if !path.exists() {
        return Ok(Manifest {
            schema_version: BUNDLE_SCHEMA_VERSION,
            citation: default_citation(),
            installed: BTreeMap::new(),
        });
    }
    let raw = std::fs::read_to_string(&path).map_err(MycoNoteError::Io)?;
    let m: Manifest = toml::from_str(&raw)
        .map_err(|e| MycoNoteError::InvalidFormat(format!("Corrupt MANIFEST.toml: {e}")))?;
    if m.schema_version != BUNDLE_SCHEMA_VERSION {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Y1000+ cache at {} was built with schema v{} but myconote expects v{}. \
             Either `setup --y1000plus --upgrade` or delete the cache and re-install.",
            root.display(),
            m.schema_version,
            BUNDLE_SCHEMA_VERSION
        )));
    }
    Ok(m)
}

pub fn save(root: &Path, manifest: &Manifest) -> Result<()> {
    std::fs::create_dir_all(root).map_err(MycoNoteError::Io)?;
    let toml_str = toml::to_string_pretty(manifest)
        .map_err(|e| MycoNoteError::InvalidFormat(format!("Serialising MANIFEST.toml: {e}")))?;
    std::fs::write(manifest_path(root), toml_str).map_err(MycoNoteError::Io)?;
    Ok(())
}

/// Best-effort directory size in bytes. Used to track on_disk_bytes after
/// extraction; silently ignores unreadable entries.
pub fn dir_size_bytes(path: &Path) -> u64 {
    fn walk(p: &Path, acc: &mut u64) {
        if let Ok(entries) = std::fs::read_dir(p) {
            for e in entries.flatten() {
                if let Ok(meta) = e.metadata() {
                    if meta.is_file() {
                        *acc += meta.len();
                    } else if meta.is_dir() {
                        walk(&e.path(), acc);
                    }
                }
            }
        }
    }
    let mut acc = 0u64;
    walk(path, &mut acc);
    acc
}

/// Utc now as "2026-04-17T21:30:12Z" without pulling chrono in.
pub fn utc_now_iso() -> String {
    // We don't need perfect precision for a human-readable timestamp; use the
    // raw seconds-since-epoch and let it stand in if we ever need sub-second.
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Naive calendar math to avoid the chrono dependency. Good enough for
    // human-readable install timestamps; not leap-second aware.
    let days_since_epoch = (secs / 86_400) as i64;
    let (y, m, d) = civil_from_days(days_since_epoch);
    let time_of_day = secs % 86_400;
    let hh = time_of_day / 3600;
    let mm = (time_of_day % 3600) / 60;
    let ss = time_of_day % 60;
    format!("{y:04}-{m:02}-{d:02}T{hh:02}:{mm:02}:{ss:02}Z")
}

/// Howard Hinnant's civil-from-days algorithm.
fn civil_from_days(z: i64) -> (i32, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = y + if m <= 2 { 1 } else { 0 };
    (y as i32, m as u32, d as u32)
}
