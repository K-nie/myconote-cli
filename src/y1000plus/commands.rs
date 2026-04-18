//! User-facing commands for the Y1000+ bundle: `--list`, `--dry-run`,
//! `--include <csv>`, `--preset <name>`. Actual downloading lives in
//! a sibling `download.rs` module (arriving in the next checkpoint).
//!
//! Citation stamped on every `--list`: Opulente et al. (2024), Science.

use crate::utils::error::{MycoNoteError, Result};
use crate::y1000plus::manifest::{cache_root, dir_size_bytes, load as load_manifest};
use crate::y1000plus::presets::Preset;
use crate::y1000plus::subsets::{format_bytes, Subset};
use crate::y1000plus::CITATION;
use std::collections::BTreeSet;

/// Parsed form of the flags that target the Y1000+ bundle.
#[derive(Debug, Default, Clone)]
pub struct Y1000Args {
    pub list: bool,
    pub dry_run: bool,
    pub yes: bool,
    pub include: Vec<Subset>,
    pub preset: Option<Preset>,
    pub uninstall: Vec<Subset>,
}

/// Extract the Y1000+ flags from a setup-command args slice. Unknown tokens
/// are left for the caller's existing DB-setup flow to handle, so this is a
/// pure *additive* layer on top of `setup`.
///
/// Returns `(y_args, remaining)` where `remaining` are the args that didn't
/// match any Y1000+ flag.
pub fn parse_args(args: &[String]) -> Result<(Y1000Args, Vec<String>)> {
    let mut y = Y1000Args::default();
    let mut rest = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--list" => {
                y.list = true;
                i += 1;
            }
            "--dry-run" => {
                y.dry_run = true;
                i += 1;
            }
            "--yes" | "-y" => {
                y.yes = true;
                i += 1;
            }
            "--include" if i + 1 < args.len() => {
                for tok in args[i + 1].split(',') {
                    let tok = tok.trim();
                    if tok.is_empty() {
                        continue;
                    }
                    let s = Subset::from_str(tok).ok_or_else(|| {
                        MycoNoteError::InvalidFormat(format!(
                            "Unknown Y1000+ subset '{tok}'. Run `setup --y1000plus --list` \
                             to see available subsets.",
                        ))
                    })?;
                    y.include.push(s);
                }
                i += 2;
            }
            "--preset" if i + 1 < args.len() => {
                let name = &args[i + 1];
                y.preset = Some(Preset::from_str(name).ok_or_else(|| {
                    MycoNoteError::InvalidFormat(format!(
                        "Unknown Y1000+ preset '{name}'. Choose one of: \
                         starter | phylogeny | compare | reference | full",
                    ))
                })?);
                i += 2;
            }
            "--uninstall" if i + 1 < args.len() => {
                for tok in args[i + 1].split(',') {
                    let tok = tok.trim();
                    if tok.is_empty() {
                        continue;
                    }
                    let s = Subset::from_str(tok).ok_or_else(|| {
                        MycoNoteError::InvalidFormat(format!("Unknown subset '{tok}'"))
                    })?;
                    y.uninstall.push(s);
                }
                i += 2;
            }
            _ => {
                rest.push(args[i].clone());
                i += 1;
            }
        }
    }
    Ok((y, rest))
}

/// `setup --y1000plus --list`
///
/// Prints a table of every subset: whether it's installed locally, its size
/// on disk + source size, and what it enables. Totals at the bottom.
pub fn run_list() -> Result<()> {
    let root = cache_root();
    let manifest = load_manifest(&root).unwrap_or_default();

    println!("Y1000+ bundle  ({})", root.display());
    println!();
    println!("{}", CITATION);
    println!();

    println!(
        "{:<18} {:>10} {:>12} {:<50}",
        "Subset", "Source", "Installed", "Enables"
    );
    println!("{}", "─".repeat(92));

    let mut total_source: u64 = 0;
    let mut total_installed: u64 = 0;
    for s in Subset::ALL {
        let source_size = s.total_bytes();
        total_source += source_size;

        let installed_info = manifest.installed.get(s.key());
        let installed_str = match installed_info {
            Some(info) => {
                total_installed += info.on_disk_bytes;
                format_bytes(info.on_disk_bytes)
            }
            None => "—".to_string(),
        };

        println!(
            "{:<18} {:>10} {:>12} {:<50}",
            s.key(),
            format_bytes(source_size),
            installed_str,
            s.summary(),
        );
    }
    println!("{}", "─".repeat(92));
    println!(
        "{:<18} {:>10} {:>12}",
        "Total",
        format_bytes(total_source),
        format_bytes(total_installed),
    );
    println!();
    println!("Presets (use `--preset <name>`):");
    for p in Preset::ALL {
        println!("  {:<10}  {}", p.key(), p.description());
    }
    println!();
    println!("Examples:");
    println!("  myconote-cli setup --y1000plus --dry-run --preset starter");
    println!("  myconote-cli setup --y1000plus --include kegg,busco");
    println!("  myconote-cli setup --y1000plus --uninstall domains");
    println!();
    println!("Cache location is overridable via the MYCONOTE_Y1000_CACHE env var.");
    Ok(())
}

/// Resolve the user's request (preset + include) into a de-duplicated subset
/// list in stable order.
pub fn resolve_target_subsets(args: &Y1000Args) -> Vec<Subset> {
    let mut set: BTreeSet<&'static str> = BTreeSet::new();
    let mut ordered: Vec<Subset> = Vec::new();
    let mut push = |set: &mut BTreeSet<&'static str>, ordered: &mut Vec<Subset>, s: Subset| {
        if set.insert(s.key()) {
            ordered.push(s);
        }
    };
    if let Some(p) = args.preset {
        for s in p.subsets() {
            push(&mut set, &mut ordered, s);
        }
    }
    for &s in &args.include {
        push(&mut set, &mut ordered, s);
    }
    ordered
}

/// `setup --y1000plus --dry-run [--include … | --preset …]`
///
/// Prints exactly what *would* be downloaded and how much space it would
/// take, without touching the network. Also surfaces skips for already-
/// installed subsets so repeated runs are idempotent.
pub fn run_dry_run(args: &Y1000Args) -> Result<()> {
    let targets = resolve_target_subsets(args);
    if targets.is_empty() {
        println!(
            "No subsets selected. Pass `--preset <name>` or `--include a,b,c`.\n\
             Run `setup --y1000plus --list` to see all options."
        );
        return Ok(());
    }

    let root = cache_root();
    let manifest = load_manifest(&root).unwrap_or_default();

    println!("Y1000+ dry-run  (cache: {})", root.display());
    println!();
    println!(
        "{:<18} {:>10} {:<18} {:<40}",
        "Subset", "Download", "Status", "Enables"
    );
    println!("{}", "─".repeat(88));

    let mut to_download: u64 = 0;
    let mut already: u64 = 0;
    for s in &targets {
        let bytes = s.total_bytes();
        let (status, counts_toward_download) = if manifest.installed.contains_key(s.key()) {
            ("already installed", false)
        } else {
            ("will download", true)
        };
        if counts_toward_download {
            to_download += bytes;
        } else {
            already += bytes;
        }
        println!(
            "{:<18} {:>10} {:<18} {:<40}",
            s.key(),
            format_bytes(bytes),
            status,
            s.summary(),
        );
    }
    println!("{}", "─".repeat(88));
    println!("  to download : {}", format_bytes(to_download));
    println!("  already here: {}", format_bytes(already));
    println!();
    println!("Run without --dry-run to proceed. Add --yes to skip the confirmation prompt.");
    Ok(())
}

/// Route the parsed Y1000+ args to the right handler. Returns `Ok(true)` if
/// the command was fully handled here (so the caller should skip its
/// existing DB-setup path); `Ok(false)` if the Y1000+ layer didn't apply
/// and the caller should continue with regular database setup.
pub fn dispatch(args: &Y1000Args) -> Result<bool> {
    use crate::y1000plus::install::{
        install_many, pending_download_bytes, uninstall_many, InstallOptions,
    };

    if args.list {
        run_list()?;
        return Ok(true);
    }
    if args.dry_run {
        run_dry_run(args)?;
        return Ok(true);
    }
    if !args.uninstall.is_empty() {
        uninstall_many(&args.uninstall)?;
        return Ok(true);
    }

    let targets = resolve_target_subsets(args);
    if targets.is_empty() {
        return Ok(false); // Nothing Y1000+-ish asked for; let caller fall through.
    }

    // Confirmation prompt unless --yes was set. Tells the user the real byte
    // cost (excluding subsets already on disk) before any network hit.
    if !args.yes {
        let bytes = pending_download_bytes(&targets);
        println!(
            "About to download {} for {} subset(s). Continue? [y/N] ",
            format_bytes(bytes),
            targets.len()
        );
        let mut line = String::new();
        std::io::stdin()
            .read_line(&mut line)
            .map_err(MycoNoteError::Io)?;
        let trimmed = line.trim().to_lowercase();
        if trimmed != "y" && trimmed != "yes" {
            println!("Aborted. Re-run with `--yes` to skip this prompt.");
            return Ok(true);
        }
    }

    install_many(&targets, &InstallOptions::default())?;
    Ok(true)
}
