//! Extract one primary-transcript protein per gene for ortholog inference.
//!
//! Alternative-splice isoforms cluster into the same orthogroup anyway, so
//! feeding them all to OrthoFinder inflates the all-vs-all similarity search
//! without adding information. The convention across comparative-genomics
//! tools is to keep one representative transcript per gene; we pick the
//! longest coding sequence (the standard "longest CDS" heuristic).

use crate::annotate::genetic_code::GeneticCode;
use crate::parser::fasta::{read_fasta_index, reverse_complement};
use crate::parser::gff::{GFFReader, GFFRecord};
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Extract one protein sequence per gene (longest-CDS isoform) from a GFF3 +
/// FASTA pair and write them to `out_fa`. Returns the number of proteins
/// written.
///
/// Headers are of the form `>{gene_id}` so OrthoFinder's per-species gene
/// identifiers round-trip to the same IDs used elsewhere in the MycoNote
/// output tree.
pub fn extract_primary_proteins(
    gff: &Path,
    fasta: &Path,
    out_fa: &Path,
    genetic_code_table: u8,
) -> Result<usize> {
    let fasta_index = read_fasta_index(fasta)?;
    let records: Vec<GFFRecord> = GFFReader::from_path(gff)?
        .filter_map(|r| r.ok())
        .collect();

    // 1. Build mRNA_id → gene_id (via Parent on mRNA row).
    let mut mrna_to_gene: HashMap<String, String> = HashMap::new();
    for rec in &records {
        if rec.feature_type != "mRNA" && rec.feature_type != "transcript" {
            continue;
        }
        if let (Some(id), Some(parent)) = (rec.id(), rec.parent()) {
            mrna_to_gene.insert(id.clone(), parent.clone());
        }
    }

    // 2. Group CDS rows by mRNA parent, summing coding length.
    let mut mrna_cds: HashMap<String, Vec<GFFRecord>> = HashMap::new();
    for rec in &records {
        if rec.feature_type != "CDS" {
            continue;
        }
        let Some(parent) = rec.parent() else { continue };
        // Parent can be a comma-list for isoform-sharing CDS; use the first.
        let key = parent.split(',').next().unwrap_or(parent).to_string();
        mrna_cds.entry(key).or_default().push(rec.clone());
    }

    // 3. Per gene, pick the mRNA whose CDSes sum to the longest total.
    let mut gene_to_best_mrna: HashMap<String, (String, u64)> = HashMap::new();
    for (mrna_id, cdss) in &mrna_cds {
        let Some(gene_id) = mrna_to_gene.get(mrna_id) else {
            // mRNA row missing — fall back to the CDS's direct parent and
            // treat it as its own gene (conservative: keeps the data).
            let total_len: u64 = cdss.iter().map(|c| c.end.saturating_sub(c.start) + 1).sum();
            let entry = gene_to_best_mrna
                .entry(mrna_id.clone())
                .or_insert((mrna_id.clone(), 0));
            if total_len > entry.1 {
                *entry = (mrna_id.clone(), total_len);
            }
            continue;
        };
        let total_len: u64 = cdss.iter().map(|c| c.end.saturating_sub(c.start) + 1).sum();
        let entry = gene_to_best_mrna
            .entry(gene_id.clone())
            .or_insert((mrna_id.clone(), 0));
        if total_len > entry.1 {
            *entry = (mrna_id.clone(), total_len);
        }
    }

    // 4. Translate each winning mRNA and write to the output FASTA.
    let gc = GeneticCode::from_table_number(genetic_code_table).unwrap_or(GeneticCode::Standard);
    let mut out = std::fs::File::create(out_fa).map_err(MycoNoteError::Io)?;
    let mut written = 0usize;

    for (gene_id, (mrna_id, _len)) in &gene_to_best_mrna {
        let Some(cdss) = mrna_cds.get(mrna_id) else { continue };
        if cdss.is_empty() { continue; }

        let mut sorted = cdss.clone();
        sorted.sort_by_key(|r| r.start);
        let strand = sorted[0].strand;

        let Some(seq_rec) = fasta_index.get(&sorted[0].seqid) else { continue };

        let mut cds_seq = String::new();
        for cds in &sorted {
            cds_seq.push_str(seq_rec.subsequence(cds.start, cds.end));
        }
        if strand == '-' {
            cds_seq = reverse_complement(&cds_seq);
        }

        let protein = gc.translate(&cds_seq);
        // Strip a trailing stop codon (`*`) but reject runts and internal-stop
        // pseudo-genes — they're not useful training material for OrthoFinder
        // and inflate the all-vs-all search.
        let trimmed = protein.trim_end_matches('*');
        if trimmed.len() < 10 {
            continue;
        }
        if trimmed.contains('*') {
            continue;
        }

        writeln!(out, ">{}", gene_id).map_err(MycoNoteError::Io)?;
        for chunk in trimmed.as_bytes().chunks(60) {
            writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                .map_err(MycoNoteError::Io)?;
        }
        written += 1;
    }

    Ok(written)
}

/// Count the number of unique gene IDs in a GFF3 — used for pre-flight tier
/// estimation before running the (potentially expensive) protein extraction.
pub fn count_genes(gff: &Path) -> Result<usize> {
    let mut n = 0usize;
    for rec in GFFReader::from_path(gff)?.filter_map(|r| r.ok()) {
        if rec.feature_type == "gene" {
            n += 1;
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write_tmp(label: &str, ext: &str, body: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!(
            "myconote_compare_primary_{}_{}.{}",
            std::process::id(),
            label,
            ext
        ));
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        p
    }

    #[test]
    fn longest_isoform_wins_when_two_mrnas_share_a_gene() {
        // Gene g1 has two mRNAs: short (9 bp CDS) and long (21 bp CDS).
        // Only the long one should survive primary-transcript selection.
        // Note: the 9-bp variant's translation would be 3 aa (<10) and thus
        // filtered anyway, so this test validates the selection logic rather
        // than the length filter.
        let gff = "\
##gff-version 3\n\
chr1\ttest\tgene\t1\t60\t.\t+\t.\tID=g1\n\
chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=g1.t1;Parent=g1\n\
chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g1.t1\n\
chr1\ttest\tmRNA\t1\t60\t.\t+\t.\tID=g1.t2;Parent=g1\n\
chr1\ttest\tCDS\t1\t60\t.\t+\t0\tParent=g1.t2\n\
";
        // 60 nt = 20 codons → 20 aa after translation (actually 19 + stop or
        // 20 sense depending on sequence). Use M followed by glycines.
        let fa = ">chr1\nATGGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGT\n";
        let gff_p = write_tmp("longest_two", "gff3", gff);
        let fa_p = write_tmp("longest_two", "fa", fa);
        let out = write_tmp("longest_two", "faa", "");

        let n = extract_primary_proteins(&gff_p, &fa_p, &out, 1).unwrap();
        assert_eq!(n, 1, "expected exactly one primary-transcript protein");

        let body = std::fs::read_to_string(&out).unwrap();
        assert!(body.starts_with(">g1\n"), "header should be gene_id:\n{}", body);
        // 60 nt / 3 = 20 aa. The actual protein depends on reading frame
        // but must be longer than the short (30-nt) variant's 10 aa.
        let seq: String = body.lines().skip(1).collect();
        assert!(seq.len() >= 15, "should pick the 60-nt isoform:\n{}", seq);

        for p in [&gff_p, &fa_p, &out] {
            let _ = std::fs::remove_file(p);
        }
    }

    #[test]
    fn internal_stop_genes_are_filtered_out() {
        // A CDS that translates to a protein with an internal stop codon is
        // a pseudogene or misannotation — should be dropped.
        let gff = "\
##gff-version 3\n\
chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g_pseudo\n\
chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=g_pseudo.t1;Parent=g_pseudo\n\
chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g_pseudo.t1\n\
";
        // ATG GGT TAA GGT GGT ... → M G * G G ...  (TAA = stop codon 3rd)
        let fa = ">chr1\nATGGGTTAAGGTGGTGGTGGTGGTGGTGGT\n";
        let gff_p = write_tmp("internal_stop", "gff3", gff);
        let fa_p = write_tmp("internal_stop", "fa", fa);
        let out = write_tmp("internal_stop", "faa", "");

        let n = extract_primary_proteins(&gff_p, &fa_p, &out, 1).unwrap();
        assert_eq!(n, 0, "pseudogene with internal stop must be filtered");

        for p in [&gff_p, &fa_p, &out] {
            let _ = std::fs::remove_file(p);
        }
    }

    #[test]
    fn count_genes_counts_only_gene_features() {
        let gff = "\
##gff-version 3\n\
chr1\ttest\tgene\t1\t60\t.\t+\t.\tID=g1\n\
chr1\ttest\tmRNA\t1\t60\t.\t+\t.\tID=g1.t1;Parent=g1\n\
chr1\ttest\tCDS\t1\t60\t.\t+\t0\tParent=g1.t1\n\
chr1\ttest\tgene\t100\t200\t.\t-\t.\tID=g2\n\
chr1\ttest\tmRNA\t100\t200\t.\t-\t.\tID=g2.t1;Parent=g2\n\
chr1\ttest\tgene\t300\t400\t.\t+\t.\tID=g3\n\
";
        let p = write_tmp("count", "gff3", gff);
        assert_eq!(count_genes(&p).unwrap(), 3);
        let _ = std::fs::remove_file(&p);
    }
}
