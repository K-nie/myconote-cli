//! Extraction: take a downloaded archive and unpack/normalise it into the
//! Y1000+ cache layout. Delegates by `ExtractKind`.
//!
//! After extraction we deliberately do **not** remove the source archive —
//! keeping it around lets a subsequent `--upgrade` or corrupted-cache run
//! re-extract without re-downloading. Users can reclaim the space via
//! `setup --y1000plus --uninstall <subset>`, which will wipe both the
//! extracted tree and the archive for that subset.

use crate::utils::error::{MycoNoteError, Result};
use crate::y1000plus::subsets::{ExtractKind, FileSpec};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Unpack a single file spec into `cache_root/<spec.dest_subdir>/`.
pub fn extract_one(cache_root: &Path, archive_path: &Path, spec: &FileSpec) -> Result<()> {
    let dest_dir = cache_root.join(spec.dest_subdir);
    std::fs::create_dir_all(&dest_dir).map_err(MycoNoteError::Io)?;

    match spec.extract {
        ExtractKind::Passthrough => passthrough(archive_path, &dest_dir, spec),
        ExtractKind::TarGz | ExtractKind::TarGzThenDiamondIndex => {
            extract_tar_gz(archive_path, &dest_dir)
        }
        ExtractKind::Zip => extract_zip(archive_path, &dest_dir),
        ExtractKind::Xlsx => {
            // For xlsx we just copy the file into place — conversion to TSV
            // happens lazily the first time a subcommand reads from the
            // metabolism/phenotypes subset, so a calamine parse failure
            // doesn't brick the install step.
            passthrough(archive_path, &dest_dir, spec)
        }
    }
}

fn passthrough(src: &Path, dest_dir: &Path, spec: &FileSpec) -> Result<()> {
    let dest = dest_dir.join(spec.name);
    std::fs::copy(src, &dest).map_err(MycoNoteError::Io)?;
    println!("   · extracted {}", dest.display());
    Ok(())
}

fn extract_tar_gz(archive: &Path, dest_dir: &Path) -> Result<()> {
    let f = File::open(archive).map_err(MycoNoteError::Io)?;
    let gz = GzDecoder::new(BufReader::new(f));
    let mut tar = tar::Archive::new(gz);
    // `unpack` handles directory creation + symlinks + permissions.
    tar.unpack(dest_dir).map_err(|e| {
        MycoNoteError::ExternalTool(format!(
            "tar.gz extract failed for {}: {e}",
            archive.display()
        ))
    })?;
    println!("   · extracted → {}/", dest_dir.display());
    Ok(())
}

fn extract_zip(archive: &Path, dest_dir: &Path) -> Result<()> {
    let f = File::open(archive).map_err(MycoNoteError::Io)?;
    let mut zip = zip::ZipArchive::new(BufReader::new(f))
        .map_err(|e| MycoNoteError::ExternalTool(format!("zip open: {e}")))?;
    for i in 0..zip.len() {
        let mut entry = zip
            .by_index(i)
            .map_err(|e| MycoNoteError::ExternalTool(format!("zip entry {i}: {e}")))?;
        let Some(rel) = entry.enclosed_name() else {
            // Skip entries with absolute or traversal-laden paths — defensive
            // even though figshare archives are trusted.
            continue;
        };
        let out_path = dest_dir.join(rel);
        if entry.is_dir() {
            std::fs::create_dir_all(&out_path).map_err(MycoNoteError::Io)?;
            continue;
        }
        if let Some(parent) = out_path.parent() {
            std::fs::create_dir_all(parent).map_err(MycoNoteError::Io)?;
        }
        let mut out = File::create(&out_path).map_err(MycoNoteError::Io)?;
        std::io::copy(&mut entry, &mut out).map_err(MycoNoteError::Io)?;
    }
    println!("   · extracted → {}/", dest_dir.display());
    Ok(())
}
