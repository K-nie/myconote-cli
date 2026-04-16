use crate::parser::region::RegionSelector;
use crate::stats::GenomeStatistics;
use crate::utils::error::{MycoNoteError, Result};
use serde::Serialize;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Stage enum
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage {
    Sort,
    Mask,
    Train,
    Predict,
    Update,
    Annotate,
    Submit,
}

impl Stage {
    pub fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "sort" => Ok(Stage::Sort),
            "mask" => Ok(Stage::Mask),
            "train" => Ok(Stage::Train),
            "predict" => Ok(Stage::Predict),
            "update" => Ok(Stage::Update),
            "annotate" => Ok(Stage::Annotate),
            "submit" => Ok(Stage::Submit),
            _ => Err(MycoNoteError::ChatContext(format!(
                "unknown stage '{}'. Expected one of: sort, mask, train, predict, update, annotate, submit",
                s
            ))),
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            Stage::Sort => "sort",
            Stage::Mask => "mask",
            Stage::Train => "train",
            Stage::Predict => "predict",
            Stage::Update => "update",
            Stage::Annotate => "annotate",
            Stage::Submit => "submit",
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// StageContext — the structured payload fed to the rule engine and prompt
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Serialize, Debug, Clone)]
pub struct ArtifactSummary {
    pub name: String,
    pub size_bytes: u64,
    pub line_count: Option<usize>,
    /// First N bytes of the file (for small text files only).
    pub preview: Option<String>,
}

#[derive(Serialize, Debug, Clone)]
pub struct StageContext {
    pub stage: String,
    pub dir: String,
    pub artifacts: Vec<ArtifactSummary>,
    pub stats: Option<serde_json::Value>,
    pub notes: Vec<String>,
}

/// Maximum size of a text artifact to inline in the context.
const MAX_INLINE_BYTES: u64 = 50_000;

/// Maximum number of preview lines for small text files.
const MAX_PREVIEW_LINES: usize = 200;

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

pub fn build_context(stage: Stage, dir: &Path) -> Result<StageContext> {
    match stage {
        Stage::Sort => build_sort_context(dir),
        Stage::Mask => build_mask_context(dir),
        Stage::Train => build_train_context(dir),
        Stage::Predict => build_predict_context(dir),
        Stage::Update => build_update_context(dir),
        Stage::Annotate => build_annotate_context(dir),
        Stage::Submit => build_submit_context(dir),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-stage builders
// ─────────────────────────────────────────────────────────────────────────────

fn build_sort_context(dir: &Path) -> Result<StageContext> {
    let mut artifacts = Vec::new();
    let mut notes = Vec::new();

    // Look for name_map.tsv and genome FASTA
    collect_matching_files(
        dir,
        &["name_map.tsv", "rename_table.tsv"],
        &mut artifacts,
        &mut notes,
    );
    collect_matching_files_by_ext(
        dir,
        &["fa", "fas", "fasta", "fna"],
        &mut artifacts,
        &mut notes,
    );

    if artifacts.is_empty() {
        notes.push("No sort output files found in directory.".to_string());
    }

    Ok(StageContext {
        stage: "sort".to_string(),
        dir: dir.display().to_string(),
        artifacts,
        stats: None,
        notes,
    })
}

fn build_mask_context(dir: &Path) -> Result<StageContext> {
    let mut artifacts = Vec::new();
    let mut notes = Vec::new();

    collect_matching_files_by_ext(
        dir,
        &["fa", "fas", "fasta", "fna"],
        &mut artifacts,
        &mut notes,
    );
    collect_matching_files(
        dir,
        &["families.fa", "families.fasta"],
        &mut artifacts,
        &mut notes,
    );

    // Calculate masking percentage from FASTA if available
    let masked_fasta = find_file_containing(dir, "masked");
    if let Some(path) = masked_fasta {
        match compute_mask_percentage(&path) {
            Ok(pct) => notes.push(format!("Soft-masked content: {:.1}%", pct)),
            Err(e) => notes.push(format!("Could not compute mask %: {}", e)),
        }
    }

    Ok(StageContext {
        stage: "mask".to_string(),
        dir: dir.display().to_string(),
        artifacts,
        stats: None,
        notes,
    })
}

fn build_train_context(dir: &Path) -> Result<StageContext> {
    let mut artifacts = Vec::new();
    let mut notes = Vec::new();

    // Training produces small report files
    collect_matching_files_by_ext(dir, &["tsv", "txt", "log"], &mut artifacts, &mut notes);

    Ok(StageContext {
        stage: "train".to_string(),
        dir: dir.display().to_string(),
        artifacts,
        stats: None,
        notes,
    })
}

fn build_predict_context(dir: &Path) -> Result<StageContext> {
    let mut artifacts = Vec::new();
    let mut notes = Vec::new();

    // Look for predict_summary.txt
    collect_matching_files(
        dir,
        &["predict_summary.txt", "summary.txt"],
        &mut artifacts,
        &mut notes,
    );

    // Collect GFF3 files
    collect_matching_files_by_ext(dir, &["gff3", "gff"], &mut artifacts, &mut notes);

    // Try to compute GenomeStatistics from the consensus GFF3
    let stats = find_gff3_and_compute_stats(dir, &mut notes);

    Ok(StageContext {
        stage: "predict".to_string(),
        dir: dir.display().to_string(),
        artifacts,
        stats,
        notes,
    })
}

fn build_update_context(dir: &Path) -> Result<StageContext> {
    let mut artifacts = Vec::new();
    let mut notes = Vec::new();

    collect_matching_files_by_ext(dir, &["gff3", "gff"], &mut artifacts, &mut notes);

    let stats = find_gff3_and_compute_stats(dir, &mut notes);

    Ok(StageContext {
        stage: "update".to_string(),
        dir: dir.display().to_string(),
        artifacts,
        stats,
        notes,
    })
}

fn build_annotate_context(dir: &Path) -> Result<StageContext> {
    let mut artifacts = Vec::new();
    let mut notes = Vec::new();

    collect_matching_files_by_ext(dir, &["gff3", "gff"], &mut artifacts, &mut notes);
    collect_matching_files_by_ext(dir, &["tsv", "txt"], &mut artifacts, &mut notes);

    let stats = find_gff3_and_compute_stats(dir, &mut notes);

    // Count functional annotation sources that have non-empty results
    let annotation_sources = [
        "pfam",
        "eggnog",
        "cazyme",
        "merops",
        "interproscan",
        "busco",
        "mmseqs",
    ];
    let mut found_sources = Vec::new();
    if let Ok(entries) = std::fs::read_dir(dir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_lowercase();
            for src in &annotation_sources {
                if name.contains(src) {
                    let meta = entry.metadata().ok();
                    let size = meta.map(|m| m.len()).unwrap_or(0);
                    if size > 0 {
                        found_sources.push(format!("{} ({} bytes)", src, size));
                    }
                }
            }
        }
    }
    if !found_sources.is_empty() {
        notes.push(format!(
            "Annotation sources found: {}",
            found_sources.join(", ")
        ));
    }

    Ok(StageContext {
        stage: "annotate".to_string(),
        dir: dir.display().to_string(),
        artifacts,
        stats,
        notes,
    })
}

fn build_submit_context(dir: &Path) -> Result<StageContext> {
    let mut artifacts = Vec::new();
    let mut notes = Vec::new();

    // errorsummary.val is the key file for submission validation
    collect_matching_files(dir, &["errorsummary.val"], &mut artifacts, &mut notes);
    collect_matching_files_by_ext(dir, &["sqn", "tbl", "fsa"], &mut artifacts, &mut notes);

    Ok(StageContext {
        stage: "submit".to_string(),
        dir: dir.display().to_string(),
        artifacts,
        stats: None,
        notes,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Collect files matching exact names.
fn collect_matching_files(
    dir: &Path,
    names: &[&str],
    artifacts: &mut Vec<ArtifactSummary>,
    notes: &mut Vec<String>,
) {
    if let Ok(entries) = std::fs::read_dir(dir) {
        for entry in entries.flatten() {
            let fname = entry.file_name().to_string_lossy().to_string();
            if names.iter().any(|n| fname.contains(n)) {
                if let Ok(art) = summarize_file(&entry.path()) {
                    artifacts.push(art);
                }
            }
        }
    }
    let _ = notes; // suppress unused warning; notes used by callers for context
}

/// Collect files matching extensions.
fn collect_matching_files_by_ext(
    dir: &Path,
    exts: &[&str],
    artifacts: &mut Vec<ArtifactSummary>,
    _notes: &mut Vec<String>,
) {
    if let Ok(entries) = std::fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
                if exts.iter().any(|e| ext.eq_ignore_ascii_case(e)) {
                    if let Ok(art) = summarize_file(&path) {
                        artifacts.push(art);
                    }
                }
            }
        }
    }
}

/// Build an ArtifactSummary for a file. Inline small text files.
fn summarize_file(path: &Path) -> Result<ArtifactSummary> {
    let meta = std::fs::metadata(path).map_err(|e| {
        MycoNoteError::ChatContext(format!("cannot stat {}: {}", path.display(), e))
    })?;

    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();

    let size_bytes = meta.len();

    // Only inline text files under the size cap
    let (line_count, preview) = if size_bytes <= MAX_INLINE_BYTES && is_likely_text(path) {
        match std::fs::read_to_string(path) {
            Ok(text) => {
                let lines: Vec<&str> = text.lines().collect();
                let count = lines.len();
                let preview_text: String = lines
                    .iter()
                    .take(MAX_PREVIEW_LINES)
                    .cloned()
                    .collect::<Vec<&str>>()
                    .join("\n");
                (Some(count), Some(preview_text))
            }
            Err(_) => (None, None),
        }
    } else if size_bytes > MAX_INLINE_BYTES && is_likely_text(path) {
        // Count lines but don't inline
        let count = std::fs::read_to_string(path)
            .ok()
            .map(|t| t.lines().count());
        (count, None)
    } else {
        (None, None)
    };

    Ok(ArtifactSummary {
        name,
        size_bytes,
        line_count,
        preview,
    })
}

/// Heuristic: is the file likely a text file?
fn is_likely_text(path: &Path) -> bool {
    let text_exts = [
        "txt", "tsv", "csv", "log", "gff3", "gff", "gtf", "bed", "val", "tbl", "fsa", "fa", "fas",
        "fasta", "fna",
    ];
    match path.extension().and_then(|e| e.to_str()) {
        Some(ext) => text_exts.iter().any(|e| ext.eq_ignore_ascii_case(e)),
        None => false,
    }
}

/// Find the first GFF3 file in a directory and compute GenomeStatistics.
fn find_gff3_and_compute_stats(dir: &Path, notes: &mut Vec<String>) -> Option<serde_json::Value> {
    // Prefer "consensus" or "final" GFF3
    let preferred = [
        "consensus.gff3",
        "final.gff3",
        "annotated.gff3",
        "updated.gff3",
    ];
    let mut gff_path: Option<PathBuf> = None;

    if let Ok(entries) = std::fs::read_dir(dir) {
        let mut all_gff = Vec::new();
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_lowercase();
            if name.ends_with(".gff3") || name.ends_with(".gff") {
                if preferred.iter().any(|p| name.contains(p)) {
                    gff_path = Some(entry.path());
                    break;
                }
                all_gff.push(entry.path());
            }
        }
        if gff_path.is_none() {
            gff_path = all_gff.into_iter().next();
        }
    }

    let path = gff_path?;
    notes.push(format!(
        "GFF3 used for stats: {}",
        path.file_name().unwrap_or_default().to_string_lossy()
    ));

    let selector = RegionSelector::new();
    match GenomeStatistics::from_gff_with_selector(&path, &selector, false) {
        Ok(stats) => serde_json::to_value(&stats).ok(),
        Err(e) => {
            notes.push(format!("GenomeStatistics error: {}", e));
            None
        }
    }
}

/// Find a file whose name contains the given substring.
fn find_file_containing(dir: &Path, substr: &str) -> Option<PathBuf> {
    let entries = std::fs::read_dir(dir).ok()?;
    for entry in entries.flatten() {
        let name = entry.file_name().to_string_lossy().to_lowercase();
        let ext = entry
            .path()
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_lowercase();
        if name.contains(substr) && ["fa", "fas", "fasta", "fna"].contains(&ext.as_str()) {
            return Some(entry.path());
        }
    }
    None
}

/// Compute the percentage of lowercase (soft-masked) bases in a FASTA.
fn compute_mask_percentage(path: &Path) -> Result<f64> {
    let text = std::fs::read_to_string(path).map_err(|e| {
        MycoNoteError::ChatContext(format!("cannot read {}: {}", path.display(), e))
    })?;

    let mut total: u64 = 0;
    let mut masked: u64 = 0;

    for line in text.lines() {
        if line.starts_with('>') {
            continue;
        }
        for c in line.chars() {
            if c.is_alphabetic() {
                total += 1;
                if c.is_lowercase() {
                    masked += 1;
                }
            }
        }
    }

    if total == 0 {
        return Ok(0.0);
    }

    Ok((masked as f64 / total as f64) * 100.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn stage_from_str_valid() {
        assert_eq!(Stage::from_str("predict").unwrap(), Stage::Predict);
        assert_eq!(Stage::from_str("SORT").unwrap(), Stage::Sort);
        assert_eq!(Stage::from_str("Annotate").unwrap(), Stage::Annotate);
    }

    #[test]
    fn stage_from_str_invalid() {
        let err = Stage::from_str("bogus").unwrap_err();
        assert!(err.to_string().contains("unknown stage 'bogus'"));
    }

    #[test]
    fn stage_name_roundtrip() {
        for s in [
            "sort", "mask", "train", "predict", "update", "annotate", "submit",
        ] {
            let stage = Stage::from_str(s).unwrap();
            assert_eq!(stage.name(), s);
        }
    }

    #[test]
    fn summarize_small_text_file() {
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("test.txt");
        fs::write(&file, "line1\nline2\nline3\n").unwrap();

        let art = summarize_file(&file).unwrap();
        assert_eq!(art.name, "test.txt");
        assert_eq!(art.line_count, Some(3));
        assert!(art.preview.is_some());
    }

    #[test]
    fn build_context_empty_dir() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = build_context(Stage::Predict, dir.path()).unwrap();
        assert_eq!(ctx.stage, "predict");
        assert!(ctx.artifacts.is_empty());
        assert!(ctx.stats.is_none());
    }

    #[test]
    fn build_context_with_gff3() {
        let dir = tempfile::tempdir().unwrap();
        // Create a minimal GFF3
        let gff = dir.path().join("consensus.gff3");
        let content = "##gff-version 3\nscaffold_1\tmyconote\tgene\t1\t1000\t.\t+\t.\tID=gene1\n\
                        scaffold_1\tmyconote\tmRNA\t1\t1000\t.\t+\t.\tID=mrna1;Parent=gene1\n\
                        scaffold_1\tmyconote\texon\t1\t500\t.\t+\t.\tID=exon1;Parent=mrna1\n\
                        scaffold_1\tmyconote\tCDS\t1\t500\t.\t+\t0\tID=cds1;Parent=mrna1\n";
        fs::write(&gff, content).unwrap();

        let ctx = build_context(Stage::Predict, dir.path()).unwrap();
        assert_eq!(ctx.stage, "predict");
        assert!(!ctx.artifacts.is_empty());
        assert!(ctx.stats.is_some());
    }

    #[test]
    fn mask_percentage_computation() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("genome_masked.fa");
        // 5 uppercase + 5 lowercase = 50% masked
        fs::write(&fa, ">scaffold_1\nACGTNacgtn\n").unwrap();
        let pct = compute_mask_percentage(&fa).unwrap();
        // N and n are not alphabetic in Rust so they don't count
        assert!(pct > 40.0 && pct < 60.0);
    }
}
