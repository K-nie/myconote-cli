//! HTTP downloader for Y1000+ subsets.
//!
//! Uses the blocking `reqwest` client already vendored by myconote. For each
//! file spec we:
//!   * check whether the archive is already on disk and has the expected size
//!     (treat it as cached if so — resumable)
//!   * otherwise stream the body into `<cache_root>/_archives/<filename>` with
//!     an indicatif progress bar
//!   * compute SHA-256 on the fly for the manifest
//!
//! The download step never unpacks the archive — extraction lives in
//! `extract.rs` and runs after all files for a subset are in place. Keeping
//! the two phases separate means a mid-extraction crash doesn't cost the
//! user the download again.

use crate::utils::error::{MycoNoteError, Result};
use crate::y1000plus::subsets::{format_bytes, FileSpec};
use indicatif::{ProgressBar, ProgressStyle};
use reqwest::blocking::Client;
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

/// Subdirectory under the cache root where raw archives accumulate before
/// extraction. Kept around so re-runs can resume without a fresh download.
const ARCHIVE_SUBDIR: &str = "_archives";

pub fn archive_path(cache_root: &Path, spec: &FileSpec) -> PathBuf {
    cache_root.join(ARCHIVE_SUBDIR).join(spec.name)
}

/// Download one file spec to `<cache_root>/_archives/<spec.name>` if not
/// already present at the expected size. Returns the path + the SHA-256 hex
/// digest of the bytes on disk.
pub fn download_file(cache_root: &Path, spec: &FileSpec) -> Result<(PathBuf, String)> {
    let dest = archive_path(cache_root, spec);
    std::fs::create_dir_all(dest.parent().unwrap()).map_err(MycoNoteError::Io)?;

    // Fast-path: archive already present and size matches → hash it + return.
    if dest.exists() {
        let actual = std::fs::metadata(&dest).map_err(MycoNoteError::Io)?.len();
        if actual == spec.size_bytes {
            println!(
                "   · {} already cached ({}) — skipping download",
                spec.name,
                format_bytes(actual)
            );
            let digest = sha256_of_file(&dest)?;
            return Ok((dest, digest));
        } else {
            // Stale partial — drop it and re-fetch.
            let _ = std::fs::remove_file(&dest);
        }
    }

    let client = Client::builder()
        .timeout(Duration::from_secs(60 * 60))  // figshare can be slow; 1 h cap
        .connect_timeout(Duration::from_secs(30))
        .build()
        .map_err(|e| MycoNoteError::ExternalTool(format!("reqwest build: {e}")))?;

    println!("   ↓ downloading {} ({})", spec.name, format_bytes(spec.size_bytes));
    let mut resp = client
        .get(spec.url)
        .send()
        .map_err(|e| MycoNoteError::ExternalTool(format!("GET {}: {e}", spec.url)))?;
    if !resp.status().is_success() {
        return Err(MycoNoteError::ExternalTool(format!(
            "figshare returned HTTP {} for {}",
            resp.status(),
            spec.url
        )));
    }

    let bar = ProgressBar::new(spec.size_bytes);
    bar.set_style(
        ProgressStyle::with_template(
            "     {bar:40.cyan/blue} {bytes}/{total_bytes} ({bytes_per_sec}, eta {eta})",
        )
        .unwrap()
        .progress_chars("=>-"),
    );

    let mut file = BufWriter::new(File::create(&dest).map_err(MycoNoteError::Io)?);
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 64 * 1024];
    let mut total = 0u64;
    loop {
        let n = resp.read(&mut buf).map_err(MycoNoteError::Io)?;
        if n == 0 {
            break;
        }
        file.write_all(&buf[..n]).map_err(MycoNoteError::Io)?;
        hasher.update(&buf[..n]);
        total += n as u64;
        bar.set_position(total);
    }
    file.flush().map_err(MycoNoteError::Io)?;
    bar.finish_and_clear();

    // Size sanity. Our registry stores sizes rounded to the kB/MB reported
    // by the figshare UI, so we only alarm on big divergences (>5%); a small
    // delta just reflects rounding rather than a truncated download.
    let on_disk = std::fs::metadata(&dest).map_err(MycoNoteError::Io)?.len();
    if spec.size_bytes > 0 {
        let delta = on_disk.abs_diff(spec.size_bytes);
        let tolerance = (spec.size_bytes / 20).max(64 * 1024);
        if delta > tolerance {
            eprintln!(
                "     ⚠  expected ~{} for {}, got {} (delta {}) — kept for diagnosis",
                format_bytes(spec.size_bytes),
                spec.name,
                format_bytes(on_disk),
                format_bytes(delta),
            );
        }
    }

    let digest = hex_digest(hasher.finalize());
    Ok((dest, digest))
}

fn sha256_of_file(path: &Path) -> Result<String> {
    let mut f = File::open(path).map_err(MycoNoteError::Io)?;
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = f.read(&mut buf).map_err(MycoNoteError::Io)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hex_digest(hasher.finalize()))
}

fn hex_digest(bytes: impl AsRef<[u8]>) -> String {
    let mut out = String::with_capacity(bytes.as_ref().len() * 2);
    for b in bytes.as_ref() {
        out.push_str(&format!("{:02x}", b));
    }
    out
}
