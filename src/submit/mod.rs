/// NCBI submission preparation module
///
/// Prepares genome annotations for NCBI GenBank/RefSeq submission by:
///   1. Validating GFF3 + FASTA against NCBI requirements
///   2. Converting to NCBI feature table format (.tbl)
///   3. Running table2asn (successor to tbl2asn) for .sqn generation
///   4. Validating output with NCBI's asnval
///
/// This is myconote-cli's own NCBI submission pipeline — completely
/// independent and not derived from any other annotation tool.
///
/// Reference: https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
use crate::parser::gff::{GFFReader, GFFRecord};
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct SubmitConfig {
    /// Input GFF3 (annotated)
    pub gff: PathBuf,
    /// Genome FASTA
    pub fasta: PathBuf,
    /// Output directory
    pub out_dir: PathBuf,
    /// Organism name (required)
    pub organism: String,
    /// Strain name
    pub strain: Option<String>,
    /// NCBI BioProject accession
    pub bioproject: Option<String>,
    /// NCBI BioSample accession
    pub biosample: Option<String>,
    /// Locus tag prefix (registered with NCBI)
    pub locus_tag_prefix: String,
    /// Molecule type: "genomic DNA" (default)
    pub mol_type: String,
    /// Topology: "linear" (default)
    pub topology: String,
    /// Genetic code table number
    pub genetic_code: u8,
    /// Annotation pipeline name for source qualifier
    pub pipeline: String,
    /// Contact email
    pub email: String,
    /// Run table2asn validation
    pub validate: bool,
}

impl Default for SubmitConfig {
    fn default() -> Self {
        Self {
            gff: PathBuf::new(),
            fasta: PathBuf::new(),
            out_dir: PathBuf::from("submit_out"),
            organism: String::new(),
            strain: None,
            bioproject: None,
            biosample: None,
            locus_tag_prefix: "MYCO".to_string(),
            mol_type: "genomic DNA".to_string(),
            topology: "linear".to_string(),
            genetic_code: 1,
            pipeline: format!("myconote-cli v{}", env!("CARGO_PKG_VERSION")),
            email: String::new(),
            validate: true,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Validation results
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
pub struct ValidationResult {
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
    pub gene_count: usize,
    pub trna_count: usize,
    pub rrna_count: usize,
    pub cds_count: usize,
    pub has_locus_tags: bool,
    pub has_products: bool,
    pub duplicate_ids: Vec<String>,
    pub orphan_features: Vec<String>,
    pub coordinate_issues: Vec<String>,
}

impl ValidationResult {
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }

    pub fn print_summary(&self) {
        if self.is_valid() {
            println!("  Validation: PASSED");
        } else {
            println!("  Validation: FAILED ({} errors)", self.errors.len());
            for e in &self.errors {
                println!("    ERROR: {}", e);
            }
        }
        if !self.warnings.is_empty() {
            println!("  Warnings: {}", self.warnings.len());
            for w in self.warnings.iter().take(10) {
                println!("    WARN: {}", w);
            }
            if self.warnings.len() > 10 {
                println!("    ... and {} more", self.warnings.len() - 10);
            }
        }
        println!(
            "  Features: {} genes, {} CDS, {} tRNA, {} rRNA",
            self.gene_count, self.cds_count, self.trna_count, self.rrna_count
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry points
// ─────────────────────────────────────────────────────────────────────────────

/// Validate GFF3 + FASTA for NCBI compliance.
pub fn validate_for_ncbi(gff: &Path, fasta: &Path) -> Result<ValidationResult> {
    let mut result = ValidationResult::default();

    // Load all GFF records
    let records: Vec<GFFRecord> = GFFReader::from_path(gff)?.filter_map(|r| r.ok()).collect();

    // Check for duplicate IDs
    let mut seen_ids: HashSet<String> = HashSet::new();
    for rec in &records {
        if let Some(id) = rec.id() {
            if !seen_ids.insert(id.clone()) {
                result.duplicate_ids.push(id.clone());
            }
        }
    }
    if !result.duplicate_ids.is_empty() {
        result.errors.push(format!(
            "{} duplicate feature IDs found",
            result.duplicate_ids.len()
        ));
    }

    // Check Parent references
    let all_ids: HashSet<String> = records.iter().filter_map(|r| r.id().cloned()).collect();

    for rec in &records {
        if let Some(parent) = rec.parent() {
            if !all_ids.contains(parent.as_str()) {
                result.orphan_features.push(format!(
                    "{}:{}-{} references missing parent '{}'",
                    rec.seqid, rec.start, rec.end, parent
                ));
            }
        }
    }
    if !result.orphan_features.is_empty() {
        result.errors.push(format!(
            "{} features reference missing parents",
            result.orphan_features.len()
        ));
    }

    // Count feature types
    for rec in &records {
        match rec.feature_type.as_str() {
            "gene" => result.gene_count += 1,
            "CDS" => result.cds_count += 1,
            "tRNA" => result.trna_count += 1,
            "rRNA" => result.rrna_count += 1,
            _ => {}
        }
    }

    // Check locus_tags and products
    let genes_with_locus: usize = records
        .iter()
        .filter(|r| r.feature_type == "gene")
        .filter(|r| r.attributes.contains_key("locus_tag"))
        .count();
    result.has_locus_tags = genes_with_locus == result.gene_count;

    if !result.has_locus_tags && result.gene_count > 0 {
        result.warnings.push(format!(
            "Only {}/{} genes have locus_tag attributes",
            genes_with_locus, result.gene_count
        ));
    }

    let genes_with_product: usize = records
        .iter()
        .filter(|r| r.feature_type == "gene")
        .filter(|r| r.attributes.contains_key("product"))
        .count();
    result.has_products = genes_with_product == result.gene_count;

    if !result.has_products && result.gene_count > 0 {
        result.warnings.push(format!(
            "Only {}/{} genes have product descriptions",
            genes_with_product, result.gene_count
        ));
    }

    // Check coordinate ordering
    for rec in &records {
        if rec.start > rec.end {
            result.coordinate_issues.push(format!(
                "{}:{} start ({}) > end ({})",
                rec.seqid, rec.feature_type, rec.start, rec.end
            ));
        }
    }
    if !result.coordinate_issues.is_empty() {
        result.errors.push(format!(
            "{} features have invalid coordinates",
            result.coordinate_issues.len()
        ));
    }

    // Validate FASTA exists and has sequences matching GFF seqids
    let fasta_seqids = read_fasta_seqids(fasta)?;
    let gff_seqids: HashSet<String> = records.iter().map(|r| r.seqid.clone()).collect();

    let missing: Vec<_> = gff_seqids.difference(&fasta_seqids).collect();
    if !missing.is_empty() {
        result.errors.push(format!(
            "{} GFF seqids not found in FASTA: {}",
            missing.len(),
            missing
                .iter()
                .take(5)
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }

    Ok(result)
}

/// Convert GFF3 + FASTA to NCBI feature table (.tbl) format.
pub fn write_feature_table(gff: &Path, output_tbl: &Path, config: &SubmitConfig) -> Result<usize> {
    let records: Vec<GFFRecord> = GFFReader::from_path(gff)?.filter_map(|r| r.ok()).collect();

    let mut f = std::fs::File::create(output_tbl).map_err(MycoNoteError::Io)?;
    let mut current_seqid = String::new();
    let mut feature_count = 0usize;

    for rec in &records {
        // Write sequence header when seqid changes
        if rec.seqid != current_seqid {
            if !current_seqid.is_empty() {
                writeln!(f).map_err(MycoNoteError::Io)?;
            }
            writeln!(f, ">Feature {}", rec.seqid).map_err(MycoNoteError::Io)?;
            current_seqid = rec.seqid.clone();
        }

        match rec.feature_type.as_str() {
            "gene" => {
                let (s, e) = if rec.strand == '-' {
                    (rec.end, rec.start)
                } else {
                    (rec.start, rec.end)
                };
                writeln!(f, "{}\t{}\tgene", s, e).map_err(MycoNoteError::Io)?;

                if let Some(locus_tag) = rec.attributes.get("locus_tag") {
                    writeln!(f, "\t\t\tlocus_tag\t{}", locus_tag).map_err(MycoNoteError::Io)?;
                }
                if let Some(gene_name) = rec.attributes.get("Name") {
                    writeln!(f, "\t\t\tgene\t{}", gene_name).map_err(MycoNoteError::Io)?;
                }
                feature_count += 1;
            }
            "CDS" => {
                let (s, e) = if rec.strand == '-' {
                    (rec.end, rec.start)
                } else {
                    (rec.start, rec.end)
                };
                writeln!(f, "{}\t{}\tCDS", s, e).map_err(MycoNoteError::Io)?;

                let product = rec
                    .attributes
                    .get("product")
                    .cloned()
                    .unwrap_or_else(|| "hypothetical protein".to_string());
                writeln!(f, "\t\t\tproduct\t{}", product).map_err(MycoNoteError::Io)?;

                if let Some(ref codon) = rec.phase {
                    writeln!(f, "\t\t\tcodon_start\t{}", codon + 1).map_err(MycoNoteError::Io)?;
                }
                writeln!(f, "\t\t\ttransl_table\t{}", config.genetic_code)
                    .map_err(MycoNoteError::Io)?;
                feature_count += 1;
            }
            "tRNA" => {
                let (s, e) = if rec.strand == '-' {
                    (rec.end, rec.start)
                } else {
                    (rec.start, rec.end)
                };
                writeln!(f, "{}\t{}\ttRNA", s, e).map_err(MycoNoteError::Io)?;

                let product = rec
                    .attributes
                    .get("product")
                    .cloned()
                    .unwrap_or_else(|| "tRNA-Xxx".to_string());
                writeln!(f, "\t\t\tproduct\t{}", product).map_err(MycoNoteError::Io)?;
                feature_count += 1;
            }
            "rRNA" => {
                let (s, e) = if rec.strand == '-' {
                    (rec.end, rec.start)
                } else {
                    (rec.start, rec.end)
                };
                writeln!(f, "{}\t{}\trRNA", s, e).map_err(MycoNoteError::Io)?;

                let product = rec.attributes.get("product").cloned().unwrap_or_default();
                writeln!(f, "\t\t\tproduct\t{}", product).map_err(MycoNoteError::Io)?;
                feature_count += 1;
            }
            _ => {}
        }
    }

    Ok(feature_count)
}

/// Write NCBI submission template (.sbt) for table2asn.
pub fn write_submission_template(output: &Path, config: &SubmitConfig) -> Result<()> {
    let mut f = std::fs::File::create(output).map_err(MycoNoteError::Io)?;

    writeln!(f, "Submit-block ::= {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "  contact {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "    contact {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "      name name {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "        last \"\",").map_err(MycoNoteError::Io)?;
    writeln!(f, "        first \"\"").map_err(MycoNoteError::Io)?;
    writeln!(f, "      }},").map_err(MycoNoteError::Io)?;
    writeln!(f, "      affil std {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "        affil \"\",").map_err(MycoNoteError::Io)?;
    writeln!(f, "        email \"{}\"", config.email).map_err(MycoNoteError::Io)?;
    writeln!(f, "      }}").map_err(MycoNoteError::Io)?;
    writeln!(f, "    }}").map_err(MycoNoteError::Io)?;
    writeln!(f, "  }}").map_err(MycoNoteError::Io)?;
    writeln!(f, "}}").map_err(MycoNoteError::Io)?;

    Ok(())
}

/// Run table2asn to generate .sqn file for NCBI submission.
pub fn run_table2asn(config: &SubmitConfig) -> Result<PathBuf> {
    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;

    // First validate
    if config.validate {
        println!("  Validating GFF3 for NCBI compliance...");
        let validation = validate_for_ncbi(&config.gff, &config.fasta)?;
        validation.print_summary();
        if !validation.is_valid() {
            return Err(MycoNoteError::InvalidFormat(
                "GFF3 failed NCBI validation. Fix errors above before submitting.".to_string(),
            ));
        }
    }

    // Write feature table
    let tbl_path = config.out_dir.join("annotation.tbl");
    println!("  Writing NCBI feature table...");
    let n_features = write_feature_table(&config.gff, &tbl_path, config)?;
    println!(
        "  {} features written to {}",
        n_features,
        tbl_path.display()
    );

    // Write submission template
    let sbt_path = config.out_dir.join("template.sbt");
    write_submission_template(&sbt_path, config)?;

    // Copy FASTA to output dir (table2asn expects it alongside .tbl)
    let fsa_path = config.out_dir.join("annotation.fsa");
    std::fs::copy(&config.fasta, &fsa_path).map_err(MycoNoteError::Io)?;

    // Build source modifiers file
    let src_path = config.out_dir.join("source_modifiers.src");
    write_source_modifiers(&src_path, config)?;

    // Check for table2asn
    let table2asn_ok = Command::new("table2asn")
        .arg("--help")
        .output()
        .map(|o| o.status.success() || !o.stderr.is_empty())
        .unwrap_or(false);

    if !table2asn_ok {
        println!("  table2asn not found - .tbl and .fsa files are ready for manual submission.");
        println!("  Download table2asn: https://ftp.ncbi.nlm.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn/");
        return Ok(tbl_path);
    }

    // Run table2asn
    println!("  Running table2asn...");
    let sqn_path = config.out_dir.join("annotation.sqn");

    let mut cmd = Command::new("table2asn");
    cmd.args([
        "-indir",
        config.out_dir.to_str().unwrap_or(""),
        "-t",
        sbt_path.to_str().unwrap_or(""),
        "-outdir",
        config.out_dir.to_str().unwrap_or(""),
        "-V",
        "vb",
        "-j",
        &format!("[organism={}]", config.organism),
    ]);

    if let Some(ref strain) = config.strain {
        cmd.arg("-j").arg(format!("[strain={}]", strain));
    }

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("table2asn: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "table2asn failed. Check validation errors in the output directory.".to_string(),
        ));
    }

    println!("  Submission file: {}", sqn_path.display());
    Ok(sqn_path)
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn read_fasta_seqids(fasta: &Path) -> Result<HashSet<String>> {
    let file = std::fs::File::open(fasta).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut ids = HashSet::new();

    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();
        if let Some(header) = trimmed.strip_prefix('>') {
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            if !id.is_empty() {
                ids.insert(id);
            }
        }
    }

    Ok(ids)
}

fn write_source_modifiers(path: &Path, config: &SubmitConfig) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "organism\t{}", config.organism).map_err(MycoNoteError::Io)?;
    writeln!(f, "mol_type\t{}", config.mol_type).map_err(MycoNoteError::Io)?;
    writeln!(f, "topology\t{}", config.topology).map_err(MycoNoteError::Io)?;

    if let Some(ref strain) = config.strain {
        writeln!(f, "strain\t{}", strain).map_err(MycoNoteError::Io)?;
    }
    if let Some(ref bp) = config.bioproject {
        writeln!(f, "BioProject\t{}", bp).map_err(MycoNoteError::Io)?;
    }
    if let Some(ref bs) = config.biosample {
        writeln!(f, "BioSample\t{}", bs).map_err(MycoNoteError::Io)?;
    }

    Ok(())
}
