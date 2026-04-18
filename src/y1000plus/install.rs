//! Install/uninstall orchestration for Y1000+ subsets.
//!
//! For each requested subset we:
//!   1. Check the manifest — if the subset is already installed, skip unless
//!      `force` is set.
//!   2. Download every file spec the subset owns (resumable via
//!      [`download::download_file`]).
//!   3. Extract per `ExtractKind` into `<cache_root>/<dest_subdir>/`.
//!   4. Measure the on-disk size and stamp an entry in MANIFEST.toml.
//!
//! `uninstall` is the reverse: for each subset we wipe every unique
//! `dest_subdir` it owns plus every archive under `_archives/`, then
//! drop the manifest entry.

use crate::utils::error::{MycoNoteError, Result};
use crate::y1000plus::download::{archive_path, download_file};
use crate::y1000plus::extract::extract_one;
use crate::y1000plus::manifest::{
    cache_root, dir_size_bytes, load as load_manifest, save as save_manifest, utc_now_iso,
    InstalledSubset, Manifest, BUNDLE_SCHEMA_VERSION,
};
use crate::y1000plus::subsets::{format_bytes, Subset};
use std::collections::BTreeSet;
use std::path::Path;

pub struct InstallOptions {
    pub force: bool,
}

impl Default for InstallOptions {
    fn default() -> Self {
        Self { force: false }
    }
}

/// Install every subset in `targets`, updating the manifest as we go so a
/// partial failure still leaves a usable cache for the completed subsets.
pub fn install_many(targets: &[Subset], opts: &InstallOptions) -> Result<()> {
    let root = cache_root();
    std::fs::create_dir_all(&root).map_err(MycoNoteError::Io)?;
    let mut manifest = load_manifest(&root).unwrap_or_default();
    manifest.schema_version = BUNDLE_SCHEMA_VERSION;

    println!(
        "Y1000+ install: {} subset(s) → {}",
        targets.len(),
        root.display()
    );
    println!();

    let mut installed_bytes = 0u64;
    for subset in targets {
        if manifest.installed.contains_key(subset.key()) && !opts.force {
            println!(
                "● {} already installed — skipping (pass --force to re-install)",
                subset.key()
            );
            continue;
        }
        println!(
            "▶ {} ({}) — {}",
            subset.key(),
            format_bytes(subset.total_bytes()),
            subset.summary()
        );
        install_one(&root, *subset, &mut manifest)?;
        save_manifest(&root, &manifest)?;
        if let Some(info) = manifest.installed.get(subset.key()) {
            installed_bytes += info.on_disk_bytes;
        }
        println!();
    }

    println!(
        "✓ done — {} of on-disk reference data available at {}",
        format_bytes(installed_bytes),
        root.display()
    );
    Ok(())
}

fn install_one(root: &Path, subset: Subset, manifest: &mut Manifest) -> Result<()> {
    let mut last_digest = String::new();
    for spec in subset.files() {
        let (archive, digest) = download_file(root, spec)?;
        extract_one(root, &archive, spec)?;
        last_digest = digest;
    }
    let dest_root = root.join(subset.files()[0].dest_subdir);
    let size = dir_size_bytes(&dest_root);
    manifest.installed.insert(
        subset.key().to_string(),
        InstalledSubset {
            on_disk_bytes: size,
            installed_at: utc_now_iso(),
            source_checksum_sha256: last_digest,
        },
    );
    Ok(())
}

pub fn uninstall_many(targets: &[Subset]) -> Result<()> {
    let root = cache_root();
    let mut manifest = load_manifest(&root).unwrap_or_default();
    if !root.exists() {
        println!(
            "No Y1000+ cache at {} — nothing to uninstall.",
            root.display()
        );
        return Ok(());
    }

    // Collect unique dest_subdirs and archive filenames so repeated subsets
    // (e.g. PhylogenyPlace bundles multiple files into one subdir) get cleaned
    // exactly once.
    let mut dirs_to_remove: BTreeSet<String> = BTreeSet::new();
    let mut archives_to_remove: BTreeSet<String> = BTreeSet::new();
    for subset in targets {
        for spec in subset.files() {
            dirs_to_remove.insert(spec.dest_subdir.to_string());
            archives_to_remove.insert(spec.name.to_string());
        }
    }

    for dir in &dirs_to_remove {
        let p = root.join(dir);
        if p.exists() {
            std::fs::remove_dir_all(&p).map_err(MycoNoteError::Io)?;
            println!("   · removed {}/", p.display());
        }
    }
    let archives_root = root.join("_archives");
    for name in &archives_to_remove {
        let p = archives_root.join(name);
        if p.exists() {
            std::fs::remove_file(&p).map_err(MycoNoteError::Io)?;
            println!("   · removed {}", p.display());
        }
    }
    for subset in targets {
        manifest.installed.remove(subset.key());
    }
    save_manifest(&root, &manifest)?;
    println!(
        "✓ uninstalled {} subset(s) from {}",
        targets.len(),
        root.display()
    );
    Ok(())
}

/// Used as a last-step safety check before starting a potentially-huge
/// download. Returns the sum of bytes for subsets that are *not* yet
/// installed (per the manifest).
pub fn pending_download_bytes(targets: &[Subset]) -> u64 {
    let root = cache_root();
    let manifest = load_manifest(&root).unwrap_or_default();
    targets
        .iter()
        .filter(|s| {
            let spec_exists = s.files().iter().all(|f| {
                let p = archive_path(&root, f);
                p.exists()
                    && std::fs::metadata(&p)
                        .map(|m| m.len() == f.size_bytes)
                        .unwrap_or(false)
            });
            !manifest.installed.contains_key(s.key()) && !spec_exists
        })
        .map(|s| s.total_bytes())
        .sum()
}
