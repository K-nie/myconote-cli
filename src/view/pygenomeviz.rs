//! pyGenomeViz backend for synteny rendering.
//!
//! Produces publication-quality static figures (PNG / PDF / SVG at 300 DPI)
//! by shelling out to the pyGenomeViz Python library. The Rust side collects
//! the per-genome inputs + plot parameters into a JSON spec and invokes a
//! generated Python driver that drives pyGenomeViz's built-in
//! MUMmer / BLAST / MMseqs aligners.
//!
//! This is checkpoint 1 of the pygenomeviz integration — focus is a clean
//! default output on 2-genome inputs. Gene-name cascades, KEGG coloring,
//! region zoom, and themes are wired through the JSON spec so later
//! checkpoints only touch the Python template.
//!
//! Python runtime: `myconote-pygv-python` wrapper at
//! `/Users/black_einstein/miniconda3/bin/myconote-pygv-python` points at a
//! dedicated conda env (`myconote-viz`) that ships pygenomeviz + MUMmer +
//! BLAST. Users install with:
//!     conda create -n myconote-viz -c bioconda -c conda-forge \
//!         python=3.12 pygenomeviz mummer4 blast matplotlib -y

use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

use super::synteny::{PlotFormat, PygvAligner, SyntenyConfig};

/// Render a synteny figure via pyGenomeViz. Assumes config.backend is
/// `PyGenomeViz`.
///
/// Strategy: convert each GFF3 + FASTA input to GenBank, then invoke the
/// library author's own `pgv-mummer` / `pgv-blast` / `pgv-mmseqs` CLI
/// wrappers. These drive the exact rendering used in the pyGenomeViz docs
/// gallery — per-contig colour cycling, CDS feature arrows at sensible
/// zoom, curved alignment ribbons, identity colorbar. No more fighting
/// with styling in a custom Python driver.
pub fn render_synteny(config: &SyntenyConfig) -> Result<()> {
    if config.fasta1.is_none() || config.fasta2.is_none() {
        return Err(MycoNoteError::InvalidFormat(
            "pyGenomeViz backend requires --fasta1 and --fasta2. Add them and retry.".to_string(),
        ));
    }

    // ── Locate pgv-<aligner> in the viz env ──────────────────────────────
    let pgv_tool = match config.aligner {
        PygvAligner::MUMmer => "pgv-mummer",
        PygvAligner::Blast => "pgv-blast",
        PygvAligner::MMseqs => "pgv-mmseqs",
    };
    let pgv_bin = locate_pgv_tool(pgv_tool)?;

    // Ensure output parent dir exists.
    if let Some(parent) = config.output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(MycoNoteError::Io)?;
        }
    }

    // ── Stage: convert GFF3+FASTA → Genbank per genome ──────────────────
    // pgv-mummer/blast/mmseqs accept Genbank files only. We already have a
    // working convert path; shell out to it so each Genbank includes all
    // the annotate-populated qualifiers (product, locus_tag, db_xref).
    let suffix = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let tmp_dir = std::env::temp_dir().join(format!("myconote_pygv_{}", suffix));
    std::fs::create_dir_all(&tmp_dir).map_err(MycoNoteError::Io)?;

    println!(
        "🎨 Rendering synteny figure via pygenomeviz (pgv-{}, {} theme)…",
        config.aligner.as_py_class().to_lowercase(),
        config.theme.as_py_str()
    );
    println!("   Converting inputs to GenBank…");

    // Name the .gbk files after the user's --label1/--label2 — pgv-mummer
    // labels each track by the input filename's stem, so this is how the
    // labels reach the figure.
    let g1_gbk = convert_to_genbank(
        &config.gff1,
        config.fasta1.as_ref().unwrap(),
        &tmp_dir,
        &sanitize_for_filename(&config.label1),
    )?;
    let g2_gbk = convert_to_genbank(
        &config.gff2,
        config.fasta2.as_ref().unwrap(),
        &tmp_dir,
        &sanitize_for_filename(&config.label2),
    )?;

    // ── Invoke pgv-<aligner> ─────────────────────────────────────────────
    let pgv_out = tmp_dir.join("pgv_output");
    std::fs::create_dir_all(&pgv_out).map_err(MycoNoteError::Io)?;

    let format_flag = match config.format {
        PlotFormat::Png => "png",
        PlotFormat::Pdf => "pdf",
        PlotFormat::Svg => "svg",
    };

    let mut cmd = Command::new(&pgv_bin);
    cmd.arg(&g1_gbk).arg(&g2_gbk).arg("-o").arg(&pgv_out).args([
        "--formats",
        format_flag,
        "--track_align_type",
        "center",
        "--feature_plotstyle",
        "bigbox",
        "--show_scale_bar",
        "--curve",
        "--identity_thr",
        &config.min_identity.to_string(),
        "--length_thr",
        &config.min_block_len.to_string(),
        "--fig_width",
        &config.width_inches.to_string(),
        "--fig_track_height",
        &config.track_height_inches.to_string(),
        "--dpi",
        "300",
        "--feature_type2color",
        "CDS:skyblue",
    ]);
    // Dark theme via matplotlib style isn't directly exposed by pgv-*, but
    // setting normal_link_color + inverted_link_color to paired colors
    // gives a dark-theme-like output on its own.

    let out = cmd.output().map_err(MycoNoteError::Io)?;

    // pgv-* logs through python's logging module — on macOS the CLI
    // sometimes exits with a non-zero status even when the PNG is written
    // (tee into stderr while stdout handler closes). Check the output
    // file instead of the exit code as the authoritative signal.
    let produced = pgv_out.join(format!("result.{}", format_flag));
    if !produced.exists() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        let stdout = String::from_utf8_lossy(&out.stdout);
        let _ = std::fs::remove_dir_all(&tmp_dir);
        return Err(MycoNoteError::ExternalTool(format!(
            "pgv-{} failed — no {} produced.\n--- stderr ---\n{}\n--- stdout ---\n{}",
            config.aligner.as_py_class().to_lowercase(),
            produced.display(),
            stderr.trim(),
            stdout.trim()
        )));
    }

    // Copy the result into the user's requested output path.
    std::fs::copy(&produced, &config.output).map_err(MycoNoteError::Io)?;

    // Surface the logs so the user sees alignment counts, timing, etc.
    let log_path = pgv_out.join("pgv-cli.log");
    if log_path.exists() {
        if let Ok(log) = std::fs::read_to_string(&log_path) {
            for line in log
                .lines()
                .rev()
                .take(8)
                .collect::<Vec<_>>()
                .into_iter()
                .rev()
            {
                if !line.trim().is_empty() {
                    println!("   {}", line.trim());
                }
            }
        }
    }

    let _ = std::fs::remove_dir_all(&tmp_dir);

    println!("✓ Synteny figure: {}", config.output.display());
    Ok(())
}

/// Locate `pgv-<aligner>` — first try the wrappered viz env, then PATH.
fn locate_pgv_tool(tool: &str) -> Result<PathBuf> {
    // The myconote-viz conda env ships the pgv-* entrypoints alongside
    // pygenomeviz. Prefer the env-native binary so we're not dependent on
    // whatever `pgv-mummer` exists on the user's global PATH.
    let env_path =
        PathBuf::from("/Users/black_einstein/miniconda3/envs/myconote-viz/bin").join(tool);
    if env_path.exists() {
        return Ok(env_path);
    }
    which::which(tool).map_err(|_| {
        MycoNoteError::UnsupportedFormat(format!(
            "{tool} not found. Install the viz env with:\n  \
             conda create -n myconote-viz -c bioconda -c conda-forge \\\n    \
                 python=3.12 pygenomeviz mummer4 blast matplotlib -y\n  \
             Then pip upgrade: pip install --upgrade pygenomeviz"
        ))
    })
}

/// Shell out to `myconote-cli convert` to produce a Genbank file from a
/// GFF3 + FASTA pair. The resulting .gbk drops straight into pgv-mummer.
fn convert_to_genbank(gff: &Path, fasta: &Path, tmp: &Path, stem: &str) -> Result<PathBuf> {
    let out = tmp.join(format!("{}.gbk", stem));
    let myconote = std::env::current_exe().map_err(MycoNoteError::Io)?;
    let status = Command::new(&myconote)
        .arg("convert")
        .arg(gff)
        .arg("--fasta")
        .arg(fasta)
        .arg("--to")
        .arg("genbank")
        .arg("--output")
        .arg(&out)
        .output()
        .map_err(MycoNoteError::Io)?;
    if !status.status.success() || !out.exists() {
        return Err(MycoNoteError::ExternalTool(format!(
            "convert to genbank failed for {}: {}",
            gff.display(),
            String::from_utf8_lossy(&status.stderr).trim()
        )));
    }
    Ok(out)
}

/// Sanitize a user-supplied label for use as a filename stem. pgv-mummer
/// labels tracks by the input stem, so the sanitization also controls what
/// appears on the figure. Keep alphanumerics, dot, dash, underscore;
/// replace anything else with underscore. Collapse runs of underscores.
fn sanitize_for_filename(label: &str) -> String {
    let mut out = String::new();
    let mut prev_underscore = false;
    for c in label.chars() {
        if c.is_alphanumeric() || c == '_' || c == '-' || c == '.' {
            out.push(c);
            prev_underscore = false;
        } else if !prev_underscore {
            out.push('_');
            prev_underscore = true;
        }
    }
    // Strip leading/trailing underscores + empty → fallback "genome".
    let trimmed = out.trim_matches('_').to_string();
    if trimmed.is_empty() {
        "genome".to_string()
    } else {
        trimmed
    }
}
