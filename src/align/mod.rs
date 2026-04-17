//! Multiple-sequence alignment helpers used by the compare → phylogeny
//! bridge. Wraps MAFFT (one alignment per single-copy orthogroup) and
//! builds the concatenated supermatrix + partition file that IQ-TREE
//! consumes via `phylogeny`.
//!
//! Scope is intentionally narrow: this module does not replace dedicated
//! alignment pipelines like MACSE or Muscle-Hybrid. It's the minimum
//! needed to turn `compare`'s single-copy ortholog FASTA directory into a
//! species-tree-ready alignment.

use crate::utils::error::{MycoNoteError, Result};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Run MAFFT on a single orthogroup FASTA. Returns a map of
/// taxon_name → aligned-sequence-string. `taxon_names` is the canonical
/// ordering we'll enforce at concat time so empty gaps fill correctly
/// when an orthogroup happens to lack a taxon.
pub fn align_orthogroup(
    input_fa: &Path,
    taxon_names: &[String],
) -> Result<BTreeMap<String, String>> {
    let mafft = which::which("mafft").map_err(|_| {
        MycoNoteError::ExternalTool(
            "mafft not found in PATH. Install with: conda install -c bioconda mafft".to_string(),
        )
    })?;

    let output = Command::new(&mafft)
        .args(["--auto", "--quiet"])
        .arg(input_fa)
        .output()
        .map_err(MycoNoteError::Io)?;

    if !output.status.success() {
        return Err(MycoNoteError::ExternalTool(format!(
            "mafft failed on {}: {}",
            input_fa.display(),
            String::from_utf8_lossy(&output.stderr)
        )));
    }

    parse_aligned_fasta(
        &String::from_utf8_lossy(&output.stdout),
        taxon_names,
        input_fa,
    )
}

/// Parse MAFFT output into a taxon → aligned-sequence map. We match each
/// header against the known taxon names (OrthoFinder prefixes sequences
/// with the input-file stem, which is also our taxon name). Fallback: use
/// the first token of the header.
fn parse_aligned_fasta(
    aligned: &str,
    taxon_names: &[String],
    source: &Path,
) -> Result<BTreeMap<String, String>> {
    let mut by_taxon: BTreeMap<String, String> = BTreeMap::new();
    let mut current_taxon: Option<String> = None;
    let mut current_seq = String::new();

    for line in aligned.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(t) = current_taxon.take() {
                by_taxon.insert(t, std::mem::take(&mut current_seq));
            } else {
                current_seq.clear();
            }
            let matched = taxon_names
                .iter()
                .find(|t| header.starts_with(t.as_str()))
                .cloned()
                .unwrap_or_else(|| {
                    header
                        .split(|c: char| c.is_whitespace() || c == '|' || c == '_')
                        .next()
                        .unwrap_or("unknown")
                        .to_string()
                });
            current_taxon = Some(matched);
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if let Some(t) = current_taxon.take() {
        by_taxon.insert(t, current_seq);
    }

    if by_taxon.is_empty() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "MAFFT produced no aligned sequences for {}",
            source.display()
        )));
    }
    Ok(by_taxon)
}

#[derive(Debug, Clone)]
pub struct Partition {
    pub name: String,
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone)]
pub struct SpeciesTreeAlignment {
    pub supermatrix: PathBuf,
    pub partition: PathBuf,
    pub n_orthogroups: usize,
    pub n_taxa: usize,
    pub total_columns: usize,
}

/// Align every orthogroup FASTA under `sco_dir` with MAFFT, then
/// concatenate them into a single supermatrix. Writes `supermatrix.fa`
/// and a NEXUS-format `partitions.nex`. Parallelism: one MAFFT process
/// per orthogroup via rayon.
pub fn build_species_tree_alignment(
    sco_dir: &Path,
    taxon_names: &[String],
    out_dir: &Path,
    threads: usize,
) -> Result<SpeciesTreeAlignment> {
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let entries: Vec<PathBuf> = std::fs::read_dir(sco_dir)
        .map_err(MycoNoteError::Io)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| matches!(s, "fa" | "fasta" | "faa"))
                .unwrap_or(false)
        })
        .collect();

    if entries.is_empty() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "No orthogroup FASTAs found in {}",
            sco_dir.display()
        )));
    }

    println!(
        "  Aligning {} single-copy orthogroups with MAFFT ({} threads)…",
        entries.len(),
        threads
    );

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|e| MycoNoteError::ExternalTool(format!("rayon pool: {}", e)))?;

    let aligned: Vec<(String, BTreeMap<String, String>)> = pool.install(|| {
        entries
            .par_iter()
            .filter_map(|p| {
                let og = p
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("OG")
                    .to_string();
                match align_orthogroup(p, taxon_names) {
                    Ok(map) => Some((og, map)),
                    Err(e) => {
                        eprintln!("  ⚠  skipping {} — {}", p.display(), e);
                        None
                    }
                }
            })
            .collect()
    });

    if aligned.is_empty() {
        return Err(MycoNoteError::ExternalTool(
            "No orthogroups aligned successfully. Check that mafft runs.".to_string(),
        ));
    }

    let mut all_taxa: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for (_, m) in &aligned {
        for t in m.keys() {
            all_taxa.insert(t.clone());
        }
    }
    let taxa: Vec<String> = all_taxa.into_iter().collect();

    let mut concat: BTreeMap<String, String> = BTreeMap::new();
    for t in &taxa {
        concat.insert(t.clone(), String::new());
    }
    let mut partitions: Vec<Partition> = Vec::with_capacity(aligned.len());
    let mut cursor = 1usize;

    let mut sorted = aligned;
    sorted.sort_by(|a, b| a.0.cmp(&b.0));

    for (og_id, by_taxon) in &sorted {
        let width = by_taxon.values().map(|s| s.len()).next().unwrap_or(0);
        if width == 0 {
            continue;
        }
        for t in &taxa {
            let slot = concat.get_mut(t).unwrap();
            match by_taxon.get(t) {
                Some(seq) if seq.len() == width => slot.push_str(seq),
                _ => slot.push_str(&"-".repeat(width)),
            }
        }
        partitions.push(Partition {
            name: og_id.clone(),
            start: cursor,
            end: cursor + width - 1,
        });
        cursor += width;
    }

    let super_path = out_dir.join("supermatrix.fa");
    let mut f = std::fs::File::create(&super_path).map_err(MycoNoteError::Io)?;
    for t in &taxa {
        writeln!(f, ">{}", t).map_err(MycoNoteError::Io)?;
        let seq = concat.get(t).map(|s| s.as_str()).unwrap_or("");
        for chunk in seq.as_bytes().chunks(60) {
            f.write_all(chunk).map_err(MycoNoteError::Io)?;
            writeln!(f).map_err(MycoNoteError::Io)?;
        }
    }

    let part_path = out_dir.join("partitions.nex");
    let mut p = std::fs::File::create(&part_path).map_err(MycoNoteError::Io)?;
    writeln!(p, "#nexus").map_err(MycoNoteError::Io)?;
    writeln!(p, "begin sets;").map_err(MycoNoteError::Io)?;
    for part in &partitions {
        writeln!(
            p,
            "  charset {} = {}-{};",
            part.name, part.start, part.end
        )
        .map_err(MycoNoteError::Io)?;
    }
    writeln!(p, "end;").map_err(MycoNoteError::Io)?;

    let total_columns = cursor - 1;
    println!(
        "  ✓  Supermatrix: {} taxa × {} columns, {} partitions",
        taxa.len(),
        total_columns,
        partitions.len()
    );

    Ok(SpeciesTreeAlignment {
        supermatrix: super_path,
        partition: part_path,
        n_orthogroups: partitions.len(),
        n_taxa: taxa.len(),
        total_columns,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_aligned_fasta_matches_known_taxa() {
        let body = ">ct_a|gene1\nMKKLTLL--V\n>ct_b|gene2\nMKKLTLL--V\n";
        let map = parse_aligned_fasta(
            body,
            &["ct_a".into(), "ct_b".into()],
            Path::new("/tmp/og.fa"),
        )
        .unwrap();
        assert_eq!(map.len(), 2);
        assert_eq!(map.get("ct_a").unwrap(), "MKKLTLL--V");
        assert_eq!(map.get("ct_b").unwrap(), "MKKLTLL--V");
    }

    #[test]
    fn parse_empty_alignment_errors() {
        assert!(parse_aligned_fasta("", &[], Path::new("/tmp/og.fa")).is_err());
    }

    #[test]
    fn parse_handles_multiline_sequences() {
        let body = ">ct_a|gene1\nMKKL\nTLL--V\n>ct_b|gene2\nMKKL\nTLL--V\n";
        let map = parse_aligned_fasta(
            body,
            &["ct_a".into(), "ct_b".into()],
            Path::new("/tmp/og.fa"),
        )
        .unwrap();
        assert_eq!(map.get("ct_a").unwrap().len(), 10);
    }
}
