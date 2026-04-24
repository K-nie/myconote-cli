/// GFF3 → spliced CDS FASTA
///
/// Extracts per-mRNA spliced coding nucleotide sequences:
///   1. Collect CDS children of each `mRNA` / `transcript` parent.
///   2. Sort CDS segments by genomic start (ascending).
///   3. Concatenate segments in genomic order, then reverse-complement
///      the whole concatenation if the parent is on the `-` strand.
///      (Equivalently: for `-` strand, segments are emitted in
///      transcript order, which is genomic-descending.)
///   4. Apply the phase offset from the *first transcript-order CDS*
///      per the GFF3 spec so that the emitted CDS starts in-frame.
///   5. Emit one FASTA record per mRNA.
///
/// Fungal-specific note: partial CDS at contig ends are kept and
/// flagged with `partial=true` in the header so downstream tools
/// (and users) can distinguish real annotation artifacts from
/// genuine low-expression transcripts. Phase is only applied when
/// non-zero (zero phase is the common case and a no-op).
///
/// Author: Benjamin Narh-Madey
use crate::parser::fasta::{read_fasta_index, reverse_complement, FastaRecord};
use crate::parser::gff::{GFFReader, GFFRecord};
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Extract spliced CDS sequences from a GFF3 + FASTA pair, writing a
/// multi-FASTA to `output`.
///
/// Returns the number of transcripts written. Errors on I/O failures;
/// records with missing seqids or empty CDS sets are logged to stderr
/// and skipped, not fatal.
pub fn extract_spliced_cds(gff_path: &Path, fasta_path: &Path, output: &Path) -> Result<usize> {
    let fasta_index = read_fasta_index(fasta_path)?;

    // Index records in one pass. GFF3 does not guarantee CDS come after their
    // mRNA parent, so we buffer the whole file before emitting.
    let records: Vec<GFFRecord> = GFFReader::from_path(gff_path)?
        .filter_map(|r| r.ok())
        .collect();

    let mut mrnas: HashMap<String, GFFRecord> = HashMap::new();
    let mut mrna_to_gene: HashMap<String, String> = HashMap::new();
    let mut cds_by_parent: HashMap<String, Vec<GFFRecord>> = HashMap::new();

    for rec in &records {
        match rec.feature_type.as_str() {
            "mRNA" | "transcript" => {
                if let Some(id) = rec.id() {
                    mrnas.insert(id.clone(), rec.clone());
                    if let Some(p) = rec.parent() {
                        mrna_to_gene.insert(id.clone(), p.clone());
                    }
                }
            }
            "CDS" => {
                if let Some(p) = rec.parent() {
                    cds_by_parent
                        .entry(p.clone())
                        .or_default()
                        .push(rec.clone());
                }
            }
            _ => {}
        }
    }

    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut count = 0usize;

    // Iterate mRNAs in sorted order so output is deterministic across runs
    // regardless of HashMap iteration order.
    let mut mrna_ids: Vec<&String> = mrnas.keys().collect();
    mrna_ids.sort();

    for mrna_id in mrna_ids {
        let mrna = &mrnas[mrna_id];
        let transcript_id = mrna
            .attributes
            .get("transcript_id")
            .cloned()
            .unwrap_or_else(|| mrna_id.clone());

        let cds_list = match cds_by_parent.get(mrna_id) {
            Some(list) if !list.is_empty() => list,
            _ => {
                // Empty or missing CDS children: common for ncRNA-like mRNAs
                // or broken annotations. Log and skip rather than abort.
                eprintln!("warn: mRNA {} has no CDS children; skipping", mrna_id);
                continue;
            }
        };

        // Verify all CDS children agree with the mRNA's contig and strand.
        // Heterogeneous strand among a single mRNA's CDSs is a malformed GFF3;
        // skip loudly rather than silently concatenate wrong.
        let seqid = &mrna.seqid;
        let strand = mrna.strand;
        if cds_list
            .iter()
            .any(|c| &c.seqid != seqid || c.strand != strand)
        {
            eprintln!(
                "warn: mRNA {} has CDS children with mismatched seqid/strand; skipping",
                mrna_id
            );
            continue;
        }

        let seq_rec: &FastaRecord = match fasta_index.get(seqid) {
            Some(s) => s,
            None => {
                eprintln!(
                    "warn: mRNA {} references seqid '{}' not found in FASTA; skipping",
                    mrna_id, seqid
                );
                continue;
            }
        };

        // Sort CDS segments by genomic start ascending (canonical GFF3 order).
        let mut ordered = cds_list.clone();
        ordered.sort_by_key(|c| c.start);

        // Partial-at-contig-ends check. We treat a CDS whose start == 1
        // or whose end == contig length as a candidate partial. This is
        // a heuristic — a CDS flush against the contig edge might be
        // complete biologically — but it's the signal the spec asked for.
        let contig_len = seq_rec.len() as u64;
        let partial = ordered.iter().any(|c| c.start == 1 || c.end == contig_len);

        // Concatenate segments in genomic order into a single buffer.
        let mut cds_seq = String::new();
        for cds in &ordered {
            cds_seq.push_str(seq_rec.subsequence(cds.start, cds.end));
        }

        // Reverse-complement once for minus-strand mRNAs. After this step,
        // cds_seq is in transcript orientation (5' → 3'), and the first
        // nucleotide corresponds to the *transcript-first* CDS segment.
        if strand == '-' {
            cds_seq = reverse_complement(&cds_seq);
        }

        // Apply phase offset from the transcript-first CDS per GFF3 spec.
        // For `+` strand that is the genomic-first segment (ordered[0]);
        // for `-` strand that is the genomic-last segment (ordered.last()).
        let first_transcript_cds = if strand == '-' {
            ordered.last()
        } else {
            ordered.first()
        };
        let phase = first_transcript_cds.and_then(|c| c.phase).unwrap_or(0) as usize;
        if phase > 0 && phase < cds_seq.len() {
            cds_seq = cds_seq[phase..].to_string();
        }

        if cds_seq.is_empty() {
            eprintln!(
                "warn: mRNA {} produced an empty CDS after phase trim; skipping",
                mrna_id
            );
            continue;
        }

        let gene_id = mrna_to_gene
            .get(mrna_id)
            .map(|s| s.as_str())
            .unwrap_or(mrna_id.as_str());

        // Header format:  >transcript_id gene=<gene_id> strand=<+/-> segments=<n> length=<bp>[ partial=true]
        let partial_tag = if partial { " partial=true" } else { "" };
        writeln!(
            out,
            ">{} gene={} strand={} segments={} length={}{}",
            transcript_id,
            gene_id,
            strand,
            ordered.len(),
            cds_seq.len(),
            partial_tag,
        )
        .map_err(MycoNoteError::Io)?;
        for chunk in cds_seq.as_bytes().chunks(60) {
            writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                .map_err(MycoNoteError::Io)?;
        }

        count += 1;
    }

    if count == 0 {
        return Err(MycoNoteError::InvalidFormat(
            "No spliced CDS written: check that GFF3 seqids match FASTA headers \
             (seqid must equal the first whitespace-delimited token of the '>' line) \
             and that mRNA features have CDS children"
                .to_string(),
        ));
    }

    Ok(count)
}
