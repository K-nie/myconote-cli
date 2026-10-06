/// SNAP HMM self-training wrapper.
///
/// SNAP has no `-train` mode — the earlier version of this file invoked
/// `snap -train …`, which this `snap` build rejects with
/// `zoeParseOptions: unknown option (-train)`, so SNAP silently contributed
/// nothing. Real SNAP self-training is the fathom → forge → hmm-assembler.pl
/// chain that MAKER and funannotate drive, from a training GFF + genome:
///
///   1. convert the first-pass GFF3 to SNAP's ZFF annotation format (.ann),
///      with the genome copied alongside as the .dna file
///   2. `fathom training.ann training.dna -categorize 1000`   → uni.* etc.
///   3. `fathom uni.ann uni.dna -export 1000 -plus`           → export.*
///   4. `forge export.ann export.dna`                         → parameter files
///   5. `hmm-assembler.pl <name> . > snap_trained.hmm`
///
/// Those helper tools ship with the bioconda `snap` package. When they are
/// absent we skip cleanly and non-fatally so the pipeline falls back to
/// Augustus + GeneMark instead of emitting a confusing ZOE error.
use crate::parser::gff::{GFFReader, GFFRecord};
use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

pub fn snap_available() -> bool {
    tool_on_path("snap")
}

fn tool_on_path(bin: &str) -> bool {
    Command::new("which")
        .arg(bin)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Returns the list of missing companion tools, or `None` when all present.
fn missing_training_tools() -> Option<Vec<String>> {
    let missing: Vec<String> = ["fathom", "forge", "hmm-assembler.pl"]
        .iter()
        .filter(|t| !tool_on_path(t))
        .map(|t| t.to_string())
        .collect();
    if missing.is_empty() {
        None
    } else {
        Some(missing)
    }
}

pub fn train_snap(training_gff3: &Path, genome_fasta: &Path, out_dir: &Path) -> Result<PathBuf> {
    if !snap_available() {
        return Err(MycoNoteError::ExternalTool(
            "SNAP not found. Install: conda install -c bioconda snap".to_string(),
        ));
    }

    // Self-training needs SNAP's companion tools. They normally ship with the
    // bioconda `snap` package, but a stripped install can be missing them —
    // detect that up front and bail cleanly rather than run a broken chain.
    if let Some(missing) = missing_training_tools() {
        return Err(MycoNoteError::ExternalTool(format!(
            "SNAP self-training unavailable (needs {}); skipping SNAP — \
             pipeline continues on Augustus/GeneMark",
            missing.join(", ")
        )));
    }

    let snap_dir = out_dir.join("snap_training");
    std::fs::create_dir_all(&snap_dir).map_err(MycoNoteError::Io)?;

    // 1. GFF3 → ZFF (.ann), genome copied beside it as .dna. fathom resolves
    //    the .dna sequences by the names in the .ann header lines, so the two
    //    must share seqids — copying the same genome guarantees that.
    let ann_path = snap_dir.join("training.ann");
    let dna_path = snap_dir.join("training.dna");

    let records: Vec<GFFRecord> = GFFReader::from_path(training_gff3)?
        .filter_map(|r| r.ok())
        .collect();
    let zff = records_to_zff(&records);
    // Count exon rows (anything that is not a `>seq` header). No CDS means no
    // training signal, so skip cleanly instead of handing fathom an empty set.
    let exon_rows = zff.lines().filter(|l| !l.starts_with('>')).count();
    if exon_rows == 0 {
        return Err(MycoNoteError::ExternalTool(
            "SNAP self-training skipped: no CDS features in the training GFF".to_string(),
        ));
    }
    std::fs::write(&ann_path, zff).map_err(MycoNoteError::Io)?;
    std::fs::copy(genome_fasta, &dna_path).map_err(MycoNoteError::Io)?;

    // 2. Categorise genes (valid / errored / overlapping). The 1000 is the
    //    flank fathom keeps around each gene.
    run_tool(
        "fathom",
        &["training.ann", "training.dna", "-categorize", "1000"],
        &snap_dir,
    )?;

    // 3. Export the clean "unique" genes, converting everything to the plus
    //    strand (standard for SNAP training, matching MAKER's recipe).
    run_tool(
        "fathom",
        &["uni.ann", "uni.dna", "-export", "1000", "-plus"],
        &snap_dir,
    )?;

    // 4. Estimate the HMM parameters from the exported set.
    run_tool("forge", &["export.ann", "export.dna"], &snap_dir)?;

    // 5. Assemble the parameter files into a single HMM on stdout.
    let assembled = run_tool("hmm-assembler.pl", &["myconote_snap", "."], &snap_dir)?;

    let final_hmm = snap_dir.join("snap_trained.hmm");
    std::fs::write(&final_hmm, &assembled.stdout).map_err(MycoNoteError::Io)?;

    // hmm-assembler succeeds with empty output if forge produced nothing
    // usable; treat that as a clean failure rather than a zero-byte HMM.
    if std::fs::metadata(&final_hmm).map(|m| m.len()).unwrap_or(0) == 0 {
        return Err(MycoNoteError::ExternalTool(
            "SNAP self-training produced an empty HMM (no usable training genes); skipping SNAP"
                .to_string(),
        ));
    }

    Ok(final_hmm)
}

/// Run a companion tool in `cwd`, capturing its output. Returns the captured
/// `Output` on success (stdout is used by `hmm-assembler.pl`), or an error
/// carrying the tool's stderr so the caller can log a real diagnostic instead
/// of a bare "failed".
fn run_tool(bin: &str, args: &[&str], cwd: &Path) -> Result<std::process::Output> {
    let out = Command::new(bin)
        .args(args)
        .current_dir(cwd)
        .output()
        .map_err(|e| MycoNoteError::ExternalTool(format!("{} failed to launch: {}", bin, e)))?;
    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        return Err(MycoNoteError::ExternalTool(format!(
            "{} {} failed: {}",
            bin,
            args.join(" "),
            stderr.trim()
        )));
    }
    Ok(out)
}

/// Convert parsed GFF3 records to SNAP's ZFF annotation format.
///
/// ZFF lists, under each `>seqid` header, one line per coding exon:
///   `<label>\t<start>\t<end>\t<group>`
/// where `group` ties a gene's exons together, `label` is `Esngl` for a
/// single-exon gene or `Einit` / `Exon` / `Eterm` for the first / internal /
/// last exon in transcription order, and the coordinate pair runs 5'→3' — so
/// on the minus strand `start > end`. We train on CDS segments (SNAP models
/// coding structure), grouped by their transcript Parent.
fn records_to_zff(records: &[GFFRecord]) -> String {
    use std::collections::BTreeMap;

    // transcript id → (seqid, strand, CDS segments)
    let mut by_tx: BTreeMap<String, (String, char, Vec<(u64, u64)>)> = BTreeMap::new();
    for r in records {
        if r.feature_type != "CDS" {
            continue;
        }
        let parent = match r.parent() {
            Some(p) => p.clone(),
            None => continue,
        };
        let entry = by_tx
            .entry(parent)
            .or_insert_with(|| (r.seqid.clone(), r.strand, Vec::new()));
        entry.2.push((r.start, r.end));
    }

    // Bucket transcripts under their sequence, ordered by genomic start, so the
    // emitted ZFF is deterministic regardless of input row order.
    let mut by_seq: BTreeMap<String, Vec<(u64, String, char, Vec<(u64, u64)>)>> = BTreeMap::new();
    for (name, (seqid, strand, exons)) in by_tx {
        let min_start = exons.iter().map(|e| e.0).min().unwrap_or(0);
        by_seq
            .entry(seqid)
            .or_default()
            .push((min_start, name, strand, exons));
    }

    let mut out = String::new();
    for (seqid, mut txs) in by_seq {
        txs.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        out.push_str(&format!(">{}\n", seqid));
        for (_, name, strand, mut exons) in txs {
            exons.sort_by_key(|e| e.0); // ascending genomic coordinate
            if strand == '-' {
                // Transcription runs high→low on the minus strand.
                exons.reverse();
            }
            let n = exons.len();
            for (i, (lo, hi)) in exons.iter().enumerate() {
                let label = if n == 1 {
                    "Esngl"
                } else if i == 0 {
                    "Einit"
                } else if i == n - 1 {
                    "Eterm"
                } else {
                    "Exon"
                };
                // 5'→3' coordinate order: plus = (lo, hi); minus = (hi, lo).
                let (c1, c2) = if strand == '-' {
                    (*hi, *lo)
                } else {
                    (*lo, *hi)
                };
                out.push_str(&format!("{}\t{}\t{}\t{}\n", label, c1, c2, name));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn cds(seqid: &str, start: u64, end: u64, strand: char, parent: &str) -> GFFRecord {
        let mut attributes = HashMap::new();
        attributes.insert("Parent".to_string(), parent.to_string());
        GFFRecord {
            seqid: seqid.into(),
            source: "test".into(),
            feature_type: "CDS".into(),
            start,
            end,
            score: None,
            strand,
            phase: Some(0),
            attributes,
        }
    }

    #[test]
    fn single_exon_plus_gene_is_esngl() {
        let records = vec![cds("chr1", 101, 250, '+', "g1.mRNA")];
        let zff = records_to_zff(&records);
        assert_eq!(zff, ">chr1\nEsngl\t101\t250\tg1.mRNA\n");
    }

    #[test]
    fn multi_exon_plus_gene_labels_init_internal_term() {
        let records = vec![
            cds("chr1", 100, 200, '+', "g1.mRNA"),
            cds("chr1", 300, 350, '+', "g1.mRNA"),
            cds("chr1", 500, 600, '+', "g1.mRNA"),
        ];
        let zff = records_to_zff(&records);
        assert_eq!(
            zff,
            ">chr1\n\
             Einit\t100\t200\tg1.mRNA\n\
             Exon\t300\t350\tg1.mRNA\n\
             Eterm\t500\t600\tg1.mRNA\n"
        );
    }

    #[test]
    fn minus_gene_reverses_order_and_swaps_coordinates() {
        // Minus strand: transcription starts at the high coordinate, so Einit
        // is the 500-600 exon written as 600 500.
        let records = vec![
            cds("chr1", 100, 200, '-', "g2.mRNA"),
            cds("chr1", 500, 600, '-', "g2.mRNA"),
        ];
        let zff = records_to_zff(&records);
        assert_eq!(
            zff,
            ">chr1\n\
             Einit\t600\t500\tg2.mRNA\n\
             Eterm\t200\t100\tg2.mRNA\n"
        );
    }

    #[test]
    fn empty_when_no_cds() {
        let mut r = cds("chr1", 100, 200, '+', "g1.mRNA");
        r.feature_type = "exon".into();
        let zff = records_to_zff(&[r]);
        assert!(zff.is_empty());
    }
}
