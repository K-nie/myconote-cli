/// File format auto-detection
///
/// Determines the file format from:
/// 1. File extension (fast path)
/// 2. Content sniffing (first non-empty line) when the extension is ambiguous

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Recognised genomic file formats.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FileFormat {
    Fasta,
    Gff3,
    Gff2,
    Gtf,
    GenBank,
    Unknown,
}

impl std::fmt::Display for FileFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FileFormat::Fasta   => write!(f, "FASTA"),
            FileFormat::Gff3    => write!(f, "GFF3"),
            FileFormat::Gff2    => write!(f, "GFF2"),
            FileFormat::Gtf     => write!(f, "GTF"),
            FileFormat::GenBank => write!(f, "GenBank"),
            FileFormat::Unknown => write!(f, "Unknown"),
        }
    }
}

/// Detect the format of a file at `path`.
///
/// Returns `FileFormat::Unknown` rather than an error when detection fails,
/// so callers can decide how to handle unsupported files.
pub fn detect_format<P: AsRef<Path>>(path: P) -> FileFormat {
    let path = path.as_ref();

    // ── 1. Extension heuristic ────────────────────────────────────────────
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext.to_lowercase().as_str() {
            "fa" | "fna" | "faa" | "fasta" | "fas" => return FileFormat::Fasta,
            "gbk" | "gb" | "genbank" => return FileFormat::GenBank,
            "gff"  => {} // need content sniff – could be GFF2 or GFF3
            "gff3" => return FileFormat::Gff3,
            "gtf"  => return FileFormat::Gtf,
            _ => {}
        }
    }

    // ── 2. Content sniffing ───────────────────────────────────────────────
    sniff_content(path)
}

fn sniff_content(path: &Path) -> FileFormat {
    let Ok(file) = File::open(path) else { return FileFormat::Unknown };
    let reader = BufReader::new(file);

    for line in reader.lines().flatten() {
        let trimmed = line.trim();
        if trimmed.is_empty() { continue; }

        // FASTA
        if trimmed.starts_with('>') {
            return FileFormat::Fasta;
        }
        // GenBank
        if trimmed.starts_with("LOCUS") {
            return FileFormat::GenBank;
        }
        // GFF pragma
        if trimmed.starts_with("##gff-version") {
            let ver = trimmed.trim_start_matches("##gff-version").trim();
            return if ver.starts_with('3') { FileFormat::Gff3 } else { FileFormat::Gff2 };
        }
        // GTF has `gene_id` / `transcript_id` attributes in column 9
        if trimmed.starts_with('#') { continue; } // skip other pragma lines
        // Tab-separated data line – inspect column 9
        let cols: Vec<&str> = trimmed.splitn(9, '\t').collect();
        if cols.len() == 9 {
            let attrs = cols[8];
            if attrs.contains("gene_id \"") || attrs.contains("transcript_id \"") {
                return FileFormat::Gtf;
            }
            return FileFormat::Gff3; // default for tab-data without GTF markers
        }

        break; // first meaningful non-comment line exhausted
    }

    FileFormat::Unknown
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_extension() {
        assert_eq!(detect_format("genome.fa"),    FileFormat::Fasta);
        assert_eq!(detect_format("prot.faa"),     FileFormat::Fasta);
        assert_eq!(detect_format("seq.fasta"),    FileFormat::Fasta);
    }

    #[test]
    fn test_gff3_extension() {
        assert_eq!(detect_format("annot.gff3"),   FileFormat::Gff3);
    }

    #[test]
    fn test_genbank_extension() {
        assert_eq!(detect_format("genome.gbk"),   FileFormat::GenBank);
        assert_eq!(detect_format("genome.gb"),    FileFormat::GenBank);
    }

    #[test]
    fn test_gtf_extension() {
        assert_eq!(detect_format("genes.gtf"),    FileFormat::Gtf);
    }
}
