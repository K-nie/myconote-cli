use serde::Serialize;

// ─────────────────────────────────────────────────────────────────────────────
// Format detection
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Serialize)]
pub enum PasteFormat {
    Gff3,
    FastaHeader,
    ValidationError,
    LogOutput,
    Tsv,
    Unknown,
}

impl std::fmt::Display for PasteFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PasteFormat::Gff3 => write!(f, "GFF3"),
            PasteFormat::FastaHeader => write!(f, "FASTA header"),
            PasteFormat::ValidationError => write!(f, "validation error"),
            PasteFormat::LogOutput => write!(f, "log output"),
            PasteFormat::Tsv => write!(f, "TSV/CSV"),
            PasteFormat::Unknown => write!(f, "unknown format"),
        }
    }
}

/// Read all of stdin and detect the format.
pub fn read_stdin() -> (String, PasteFormat) {
    let mut input = String::new();
    let _ = std::io::Read::read_to_string(&mut std::io::stdin(), &mut input);
    let format = detect_format(&input);
    (input, format)
}

/// Detect the format of pasted text.
pub fn detect_format(text: &str) -> PasteFormat {
    let lines: Vec<&str> = text.lines().take(20).collect();
    if lines.is_empty() {
        return PasteFormat::Unknown;
    }

    // Check for GFF3: tab-separated, 9 fields, recognized feature types
    let gff_features = ["gene", "mrna", "transcript", "cds", "exon", "three_prime_utr", "five_prime_utr"];
    let gff_count = lines.iter().filter(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        fields.len() == 9 && gff_features.iter().any(|f| fields[2].to_lowercase() == *f)
    }).count();
    if gff_count >= 1 {
        return PasteFormat::Gff3;
    }

    // Check for FASTA header
    if lines.iter().any(|l| l.starts_with('>')) {
        return PasteFormat::FastaHeader;
    }

    // Check for NCBI validation errors
    let error_patterns = ["ERROR:", "WARNING:", "SEQ_FEAT.", "SEQ_DESCR.", "SEQ_INST."];
    let error_count = lines.iter().filter(|l| {
        error_patterns.iter().any(|p| l.contains(p))
    }).count();
    if error_count >= 1 {
        return PasteFormat::ValidationError;
    }

    // Check for log output (timestamps, tool names)
    let log_patterns = ["augustus", "snap", "glimmerhmm", "genemark", "pasa", "evm",
                        "repeatmasker", "repeatmodeler", "trinity", "busco",
                        "[INFO]", "[ERROR]", "[WARN]", "ERROR:", "WARNING:"];
    let log_count = lines.iter().filter(|l| {
        let lower = l.to_lowercase();
        log_patterns.iter().any(|p| lower.contains(&p.to_lowercase()))
    }).count();
    if log_count >= 2 {
        return PasteFormat::LogOutput;
    }

    // Check for TSV/CSV
    let tab_count = lines.iter().filter(|l| l.split('\t').count() >= 3).count();
    if tab_count >= 2 {
        return PasteFormat::Tsv;
    }
    let comma_count = lines.iter().filter(|l| l.split(',').count() >= 3).count();
    if comma_count >= 2 {
        return PasteFormat::Tsv;
    }

    PasteFormat::Unknown
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_gff3_line() {
        let text = "scaffold_1\tmyconote\tgene\t1\t1000\t.\t+\t.\tID=gene1";
        assert_eq!(detect_format(text), PasteFormat::Gff3);
    }

    #[test]
    fn detect_gff3_multiple_lines() {
        let text = "scaffold_1\tmyconote\tgene\t1\t1000\t.\t+\t.\tID=gene1\n\
                    scaffold_1\tmyconote\tmRNA\t1\t1000\t.\t+\t.\tID=mrna1;Parent=gene1";
        assert_eq!(detect_format(text), PasteFormat::Gff3);
    }

    #[test]
    fn detect_fasta_header() {
        let text = ">scaffold_1 some description\nACGTACGT";
        assert_eq!(detect_format(text), PasteFormat::FastaHeader);
    }

    #[test]
    fn detect_ncbi_error() {
        let text = "ERROR: valid [SEQ_FEAT.NoStop] No stop codon found\n\
                    WARNING: valid [SEQ_FEAT.NotSpliceConsensus] Splice site not consensus";
        assert_eq!(detect_format(text), PasteFormat::ValidationError);
    }

    #[test]
    fn detect_log_output() {
        let text = "[INFO] Running Augustus with species model saccharomyces\n\
                    [INFO] Running SNAP with hmm saccharomyces.hmm\n\
                    [INFO] EVM consensus complete";
        assert_eq!(detect_format(text), PasteFormat::LogOutput);
    }

    #[test]
    fn detect_tsv() {
        let text = "gene1\tPF00001\t1.2e-10\n\
                    gene2\tPF00002\t3.4e-20\n\
                    gene3\tPF00003\t5.6e-15";
        assert_eq!(detect_format(text), PasteFormat::Tsv);
    }

    #[test]
    fn detect_unknown() {
        let text = "This is just some random text about biology.";
        assert_eq!(detect_format(text), PasteFormat::Unknown);
    }

    #[test]
    fn detect_empty() {
        assert_eq!(detect_format(""), PasteFormat::Unknown);
    }
}
