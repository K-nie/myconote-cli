/// GenBank flat-file writer
///
/// Converts GFF3 + FASTA into NCBI GenBank format (.gbk / .gb).
///
/// # Features
/// - Correct LOCUS / DEFINITION / ACCESSION / FEATURES / ORIGIN sections
/// - Multi-exon CDS written as `join(...)` locations
/// - Reverse-strand features wrapped in `complement(...)`
/// - Full gene → mRNA → CDS hierarchy with propagated `locus_tag`
/// - SeqID matching: strips FASTA header to first whitespace token
///   (fixes the notorious mismatch problem from progressiveMauve workflows)
use crate::parser::fasta::FastaRecord;
use crate::parser::gff::GFFRecord;
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;

// ─────────────────────────────────────────────────────────────────────────────
// Line-width constants (NCBI spec)
// ─────────────────────────────────────────────────────────────────────────────

const LOCUS_WIDTH: usize = 80;
const FEATURE_INDENT: &str = "     "; // 5 spaces
const QUALIFIER_INDENT: &str = "                     "; // 21 spaces
const SEQUENCE_WIDTH: usize = 60;
const SEQUENCE_CHUNK: usize = 10;

// ─────────────────────────────────────────────────────────────────────────────
// Internal data model
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct Interval {
    start: u64,
    end: u64,
    strand: char,
}

#[derive(Debug)]
struct FeatureBlock {
    kind: String,
    intervals: Vec<Interval>,
    qualifiers: Vec<(String, String)>,
}

impl FeatureBlock {
    fn new(kind: &str) -> Self {
        FeatureBlock {
            kind: kind.to_string(),
            intervals: Vec::new(),
            qualifiers: Vec::new(),
        }
    }

    fn add_qualifier(&mut self, key: &str, value: &str) {
        self.qualifiers.push((key.to_string(), value.to_string()));
    }

    fn strand(&self) -> char {
        self.intervals.first().map(|i| i.strand).unwrap_or('+')
    }

    fn location_string(&self) -> String {
        if self.intervals.is_empty() {
            return String::new();
        }

        let mut sorted = self.intervals.clone();
        sorted.sort_by_key(|i| i.start);

        let parts: Vec<String> = sorted
            .iter()
            .map(|i| format!("{}..{}", i.start, i.end))
            .collect();

        let loc = if parts.len() == 1 {
            parts[0].clone()
        } else {
            format!("join({})", parts.join(","))
        };

        if self.strand() == '-' {
            format!("complement({})", loc)
        } else {
            loc
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GFF3 → FeatureBlock conversion
// ─────────────────────────────────────────────────────────────────────────────

fn strand_char(s: char) -> char {
    if s == '-' {
        '-'
    } else {
        '+'
    }
}

fn build_feature_blocks(records: &[GFFRecord], chr: &str) -> Vec<FeatureBlock> {
    let mut id_to_gene: HashMap<String, &GFFRecord> = HashMap::new();
    let mut id_to_mrna: HashMap<String, &GFFRecord> = HashMap::new();
    let mut mrna_to_cds: HashMap<String, Vec<&GFFRecord>> = HashMap::new();
    let mut gene_to_mrnas: HashMap<String, Vec<String>> = HashMap::new();

    for rec in records.iter().filter(|r| r.seqid == chr) {
        match rec.feature_type.as_str() {
            "gene" => {
                if let Some(id) = rec.id() {
                    id_to_gene.insert(id.to_string(), rec);
                }
            }
            "mRNA" | "transcript" => {
                if let Some(id) = rec.id() {
                    id_to_mrna.insert(id.to_string(), rec);
                    if let Some(parent) = rec.parent() {
                        gene_to_mrnas
                            .entry(parent.to_string())
                            .or_default()
                            .push(id.to_string());
                    }
                }
            }
            "CDS" => {
                if let Some(parent) = rec.parent() {
                    mrna_to_cds.entry(parent.to_string()).or_default().push(rec);
                }
            }
            _ => {}
        }
    }

    let mut blocks: Vec<FeatureBlock> = Vec::new();

    let mut gene_ids: Vec<&String> = id_to_gene.keys().collect();
    gene_ids.sort_by_key(|id| id_to_gene[*id].start);

    for gene_id in gene_ids {
        let gene_rec = id_to_gene[gene_id];
        let locus_tag = gene_rec
            .id()
            .map(|s| s.clone())
            .unwrap_or_else(|| format!("gene_{}", gene_rec.start));

        // gene block
        let mut gene_block = FeatureBlock::new("gene");
        gene_block.intervals.push(Interval {
            start: gene_rec.start,
            end: gene_rec.end,
            strand: strand_char(gene_rec.strand),
        });
        gene_block.add_qualifier("locus_tag", &locus_tag);
        if let Some(name) = gene_rec.attributes.get("Name") {
            gene_block.add_qualifier("gene", name);
        }
        blocks.push(gene_block);

        let empty: Vec<String> = Vec::new();
        let mrnas = gene_to_mrnas.get(gene_id).unwrap_or(&empty);

        for mrna_id in mrnas {
            let Some(mrna_rec) = id_to_mrna.get(mrna_id.as_str()) else {
                continue;
            };

            // mRNA block
            let mut mrna_block = FeatureBlock::new("mRNA");
            mrna_block.intervals.push(Interval {
                start: mrna_rec.start,
                end: mrna_rec.end,
                strand: strand_char(mrna_rec.strand),
            });
            mrna_block.add_qualifier("locus_tag", &locus_tag);
            if let Some(name) = mrna_rec.attributes.get("Name") {
                mrna_block.add_qualifier("product", name);
            }
            mrna_block.add_qualifier("transcript_id", mrna_id);
            blocks.push(mrna_block);

            // CDS block (aggregates all CDS children)
            if let Some(cds_recs) = mrna_to_cds.get(mrna_id.as_str()) {
                let mut cds_block = FeatureBlock::new("CDS");
                for cds in cds_recs.iter() {
                    cds_block.intervals.push(Interval {
                        start: cds.start,
                        end: cds.end,
                        strand: strand_char(cds.strand),
                    });
                }
                cds_block.add_qualifier("locus_tag", &locus_tag);
                if let Some(name) = mrna_rec.attributes.get("Name") {
                    cds_block.add_qualifier("product", name);
                }
                cds_block.add_qualifier("protein_id", &format!("gnl|local|{}", mrna_id));
                blocks.push(cds_block);
            }
        }
    }

    blocks
}

// ─────────────────────────────────────────────────────────────────────────────
// Formatting helpers
// ─────────────────────────────────────────────────────────────────────────────

fn write_feature<W: Write>(w: &mut W, block: &FeatureBlock) -> std::io::Result<()> {
    let loc = block.location_string();
    let key_field = format!("{:<16}", block.kind);
    writeln!(w, "{}{}{}", FEATURE_INDENT, key_field, loc)?;

    for (key, value) in &block.qualifiers {
        let qualifier = if value.parse::<i64>().is_ok() {
            format!("/{}={}", key, value)
        } else {
            format!("/{}=\"{}\"", key, value)
        };
        write_wrapped(w, QUALIFIER_INDENT, &qualifier)?;
    }
    Ok(())
}

fn write_wrapped<W: Write>(w: &mut W, indent: &str, text: &str) -> std::io::Result<()> {
    const MAX: usize = 79;
    let avail = MAX.saturating_sub(indent.len());

    if text.len() + indent.len() <= MAX {
        return writeln!(w, "{}{}", indent, text);
    }

    let mut remaining = text;
    while !remaining.is_empty() {
        let take = remaining.len().min(avail);
        let split = if take >= remaining.len() {
            take
        } else {
            remaining[..take]
                .rfind(|c: char| c == ' ' || c == ',')
                .map(|i| i + 1)
                .unwrap_or(take)
        };
        writeln!(w, "{}{}", indent, &remaining[..split])?;
        remaining = &remaining[split..];
    }
    Ok(())
}

fn write_origin<W: Write>(w: &mut W, sequence: &str) -> std::io::Result<()> {
    writeln!(w, "ORIGIN")?;
    let seq = sequence.to_lowercase();
    let bytes = seq.as_bytes();
    let total = bytes.len();
    let mut pos = 0;
    while pos < total {
        write!(w, "{:>9}", pos + 1)?;
        let line_end = (pos + SEQUENCE_WIDTH).min(total);
        let mut chunk = pos;
        while chunk < line_end {
            let end = (chunk + SEQUENCE_CHUNK).min(line_end);
            write!(
                w,
                " {}",
                std::str::from_utf8(&bytes[chunk..end]).unwrap_or("")
            )?;
            chunk = end;
        }
        writeln!(w)?;
        pos += SEQUENCE_WIDTH;
    }
    writeln!(w, "//")?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Write one GenBank record for a single chromosome / contig.
pub fn write_genbank_record<W: Write>(
    w: &mut W,
    seqid: &str,
    records: &[GFFRecord],
    fasta: &FastaRecord,
    organism: Option<&str>,
) -> Result<()> {
    let seq_len = fasta.len();
    let organism_str = organism.unwrap_or("Unknown fungal organism");

    // LOCUS
    let locus_line = format!(
        "LOCUS       {:<16} {:>11} bp    DNA     linear   FUN",
        seqid, seq_len
    );
    writeln!(w, "{}", &locus_line[..locus_line.len().min(LOCUS_WIDTH)])
        .map_err(MycoNoteError::Io)?;

    // DEFINITION
    writeln!(w, "DEFINITION  {} chromosome {}.", organism_str, seqid).map_err(MycoNoteError::Io)?;

    // ACCESSION / VERSION
    writeln!(w, "ACCESSION   {}", seqid).map_err(MycoNoteError::Io)?;
    writeln!(w, "VERSION     {}", seqid).map_err(MycoNoteError::Io)?;

    // SOURCE / ORGANISM
    writeln!(w, "SOURCE      {}", organism_str).map_err(MycoNoteError::Io)?;
    writeln!(w, "  ORGANISM  {}", organism_str).map_err(MycoNoteError::Io)?;
    writeln!(w, "            Eukaryota; Fungi.").map_err(MycoNoteError::Io)?;

    // FEATURES
    writeln!(w, "FEATURES             Location/Qualifiers").map_err(MycoNoteError::Io)?;

    // source spanning the whole record
    let mut source = FeatureBlock::new("source");
    source.intervals.push(Interval {
        start: 1,
        end: seq_len as u64,
        strand: '+',
    });
    source.add_qualifier("organism", organism_str);
    source.add_qualifier("mol_type", "genomic DNA");
    write_feature(w, &source).map_err(MycoNoteError::Io)?;

    // annotation features
    let blocks = build_feature_blocks(records, seqid);
    for block in &blocks {
        write_feature(w, block).map_err(MycoNoteError::Io)?;
    }

    // ORIGIN
    write_origin(w, &fasta.sequence).map_err(MycoNoteError::Io)?;

    Ok(())
}

/// Write a multi-record GenBank file from all GFF records + a FASTA index.
///
/// The FASTA index must be keyed by the **bare sequence ID** (first whitespace-
/// delimited token of the `>` header line).  Chromosomes without a matching
/// FASTA entry are skipped with a warning.
pub fn write_genbank<W: Write>(
    w: &mut W,
    all_records: &[GFFRecord],
    fasta_index: &HashMap<String, FastaRecord>,
    organism: Option<&str>,
) -> Result<()> {
    use std::collections::BTreeSet;

    let seqids: BTreeSet<&str> = all_records.iter().map(|r| r.seqid.as_str()).collect();

    let mut written = 0usize;
    for seqid in seqids {
        if let Some(fasta_rec) = fasta_index.get(seqid) {
            let owned: Vec<GFFRecord> = all_records
                .iter()
                .filter(|r| r.seqid == seqid)
                .cloned()
                .collect();
            write_genbank_record(w, seqid, &owned, fasta_rec, organism)?;
            written += 1;
        } else {
            eprintln!("  ⚠  No FASTA sequence for seqid '{}' — skipping", seqid);
        }
    }

    if written == 0 {
        return Err(MycoNoteError::InvalidFormat(
            "No records written: check that GFF3 seqids match FASTA headers \
             (seqid must equal the first whitespace-delimited token of the '>' line)"
                .to_string(),
        ));
    }

    Ok(())
}
