//! Multiple-sequence alignment helpers used by the compare → phylogeny
//! bridge. Wraps MAFFT (one alignment per single-copy orthogroup) and
//! builds the concatenated supermatrix + partition file that IQ-TREE
//! consumes via `phylogeny`.
//!
//! Scope is intentionally narrow: this module does not replace dedicated
//! alignment pipelines like MACSE or Muscle-Hybrid. It's the minimum
//! needed to turn `compare`'s single-copy ortholog FASTA directory into a
//! species-tree-ready alignment.

use crate::compare::orthofinder::Orthogroup;
use crate::utils::error::{MycoNoteError, Result};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Run MAFFT on a single orthogroup FASTA. Returns a map of
/// taxon_name → aligned-sequence-string.
///
/// `taxon_assignments` supplies the authoritative per-sequence taxon from
/// the ortholog table — we take each aligned record in FASTA order and
/// zip with `taxon_assignments`. This handles the pathological case where
/// two genomes share identical gene IDs (self-compare edge case), which
/// header-parsing heuristics can't distinguish.
pub fn align_orthogroup(
    input_fa: &Path,
    taxon_assignments: &[String],
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
        taxon_assignments,
        input_fa,
    )
}

/// Parse MAFFT output into a taxon → aligned-sequence map. The caller
/// supplies `taxon_assignments` — one taxon name per sequence in the INPUT
/// FASTA, in the same order OrthoFinder wrote them. Because MAFFT (and the
/// `--auto` pipeline) doesn't reorder records, we can zip positionally.
///
/// Why not match by header prefix? OrthoFinder's `Single_Copy_
/// Orthologue_Sequences/*.fa` files strip the species prefix and write the
/// bare gene ID, so header-matching against taxon names fails in the
/// general case — and fails pathologically (false-positive matches) when
/// two species have overlapping gene IDs (self-compare or shared-locus-
/// tag clades). Positional mapping is robust to both.
fn parse_aligned_fasta(
    aligned: &str,
    taxon_assignments: &[String],
    source: &Path,
) -> Result<BTreeMap<String, String>> {
    let mut seqs: Vec<String> = Vec::new();
    let mut current = String::new();
    let mut in_record = false;

    for line in aligned.lines() {
        if line.starts_with('>') {
            if in_record {
                seqs.push(std::mem::take(&mut current));
            }
            in_record = true;
        } else if in_record {
            current.push_str(line.trim());
        }
    }
    if in_record {
        seqs.push(current);
    }

    if seqs.is_empty() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "MAFFT produced no aligned sequences for {}",
            source.display()
        )));
    }

    if seqs.len() != taxon_assignments.len() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "MAFFT produced {} sequences but taxon table expects {} for {}",
            seqs.len(),
            taxon_assignments.len(),
            source.display()
        )));
    }

    let mut by_taxon: BTreeMap<String, String> = BTreeMap::new();
    for (taxon, seq) in taxon_assignments.iter().zip(seqs.into_iter()) {
        // When a single-copy orthogroup happens to collapse two identical
        // gene IDs into one (artificial self-compare), we still want both
        // entries preserved — the later overwrites, accepting the second
        // sequence as the canonical one. This matches OrthoFinder's own
        // tie-breaking.
        by_taxon.insert(taxon.clone(), seq);
    }
    Ok(by_taxon)
}

/// Build the per-orthogroup taxon ordering. For each orthogroup, we emit a
/// list the same length as the orthogroup's FASTA — OrthoFinder writes
/// single-copy FASTAs in taxon order, so we expand every genome's
/// gene-list entry into its taxon name.
fn taxon_assignments_for_orthogroup(og: &Orthogroup, taxon_names: &[String]) -> Vec<String> {
    let mut out = Vec::new();
    for taxon in taxon_names {
        if let Some(genes) = og.members.get(taxon) {
            for _ in genes {
                out.push(taxon.clone());
            }
        }
    }
    out
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
///
/// `orthogroups` supplies the authoritative per-orthogroup gene-to-taxon
/// mapping from the parsed OrthoFinder output — we use this (not header
/// parsing) to assign each aligned sequence to its taxon, which is the
/// only robust way to handle clades where two species have overlapping
/// gene IDs.
pub fn build_species_tree_alignment(
    sco_dir: &Path,
    taxon_names: &[String],
    orthogroups: &[Orthogroup],
    out_dir: &Path,
    threads: usize,
) -> Result<SpeciesTreeAlignment> {
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    // Build a map OG_id → taxon_assignments. Only single-copy orthogroups
    // make it into `sco_dir`, so we filter to those to avoid alignment on
    // paralog-containing clusters.
    let og_lookup: HashMap<String, Vec<String>> = orthogroups
        .iter()
        .filter(|og| og.is_single_copy(taxon_names.len()))
        .map(|og| {
            (
                og.id.clone(),
                taxon_assignments_for_orthogroup(og, taxon_names),
            )
        })
        .collect();

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
                let Some(assignments) = og_lookup.get(&og) else {
                    eprintln!(
                        "  ⚠  skipping {} — no single-copy taxon map (orthogroup not in table?)",
                        p.display()
                    );
                    return None;
                };
                match align_orthogroup(p, assignments) {
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
        writeln!(p, "  charset {} = {}-{};", part.name, part.start, part.end)
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
    use crate::compare::orthofinder::Orthogroup;

    #[test]
    fn positional_taxon_mapping_works_for_identical_gene_ids() {
        // Pathological self-compare case: same gene ID appears in two
        // different taxa. Header-parsing heuristics can't distinguish them;
        // positional zip with taxon_assignments does.
        let body = ">g001031\nMKKLTLL--V\n>g001031\nMKKLTII--V\n";
        let map = parse_aligned_fasta(
            body,
            &["ct_a".into(), "ct_b".into()],
            Path::new("/tmp/og.fa"),
        )
        .unwrap();
        assert_eq!(map.len(), 2);
        assert_eq!(map.get("ct_a").unwrap(), "MKKLTLL--V");
        assert_eq!(map.get("ct_b").unwrap(), "MKKLTII--V");
    }

    #[test]
    fn parse_empty_alignment_errors() {
        assert!(parse_aligned_fasta("", &[], Path::new("/tmp/og.fa")).is_err());
    }

    #[test]
    fn parse_handles_multiline_sequences() {
        let body = ">ct_a_g1\nMKKL\nTLL--V\n>ct_b_g2\nMKKL\nTLL--V\n";
        let map = parse_aligned_fasta(
            body,
            &["ct_a".into(), "ct_b".into()],
            Path::new("/tmp/og.fa"),
        )
        .unwrap();
        assert_eq!(map.get("ct_a").unwrap().len(), 10);
        assert_eq!(map.get("ct_b").unwrap().len(), 10);
    }

    #[test]
    fn mismatched_seq_count_errors() {
        // 2 sequences but we claim 3 taxa → structural mismatch, must
        // error loudly rather than silently truncate.
        let body = ">a\nMKKL\n>b\nMKKI\n";
        let err = parse_aligned_fasta(
            body,
            &["a".into(), "b".into(), "c".into()],
            Path::new("/tmp/og.fa"),
        );
        assert!(err.is_err());
    }

    #[test]
    fn taxon_assignments_handles_missing_taxon() {
        // Orthogroup missing one of the 3 expected taxa — assignments
        // reflect only present taxa so MAFFT output will be 2 sequences.
        let mut og = Orthogroup {
            id: "OG0001".into(),
            members: std::collections::HashMap::new(),
        };
        og.members.insert("ct_a".into(), vec!["g1".into()]);
        og.members.insert("ct_c".into(), vec!["g5".into()]);
        // ct_b absent
        let assignments =
            taxon_assignments_for_orthogroup(&og, &["ct_a".into(), "ct_b".into(), "ct_c".into()]);
        assert_eq!(assignments, vec!["ct_a".to_string(), "ct_c".to_string()]);
    }
}
