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
    /// Submitter first name (goes into template.sbt name.first)
    pub contact_first: String,
    /// Submitter last name (goes into template.sbt name.last)
    pub contact_last: String,
    /// Submitter institution (goes into template.sbt affil.affil)
    pub institution: String,
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
            contact_first: String::new(),
            contact_last: String::new(),
            institution: String::new(),
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
///
/// Handles NCBI's multi-segment feature syntax correctly: a CDS (or tRNA)
/// that spans multiple exons is emitted as a single feature block with one
/// interval line per segment, not as N separate features. See
/// https://www.ncbi.nlm.nih.gov/genbank/feature_table/ — "To indicate that a
/// feature is on the complementary strand, the location should be indicated
/// by putting the larger coordinate in column 1, and the smaller one in
/// column 2". For multi-segment features on the minus strand, segments are
/// also emitted in descending order (5'→3' reading order).
pub fn write_feature_table(gff: &Path, output_tbl: &Path, config: &SubmitConfig) -> Result<usize> {
    use std::collections::HashMap;

    let records: Vec<GFFRecord> = GFFReader::from_path(gff)?.filter_map(|r| r.ok()).collect();

    // Group CDS / exon-like segments by Parent (= mRNA ID). This is how a
    // multi-exon eukaryotic gene's coding segments get joined into one
    // feature in the .tbl output.
    let mut cds_by_parent: HashMap<String, Vec<&GFFRecord>> = HashMap::new();
    for rec in &records {
        if rec.feature_type == "CDS" {
            if let Some(parent) = rec.attributes.get("Parent") {
                // A GFF3 Parent can be a comma-separated list when a CDS is
                // shared between isoforms; use the first parent so each
                // joined feature block still corresponds to a single mRNA.
                let first_parent = parent.split(',').next().unwrap_or(parent).to_string();
                cds_by_parent.entry(first_parent).or_default().push(rec);
            }
        }
    }

    // Build ID → record index so we can walk the Parent chain
    // (CDS → mRNA → gene) when looking up `product`. Annotate writes the
    // product on the gene row only; without this walk every CDS in a
    // annotated.gff3 would fall through to "hypothetical protein" in the
    // emitted .tbl — exactly the bug caught by the E2E pipeline run.
    let mut by_id: HashMap<String, &GFFRecord> = HashMap::new();
    for rec in &records {
        if let Some(id) = rec.attributes.get("ID") {
            by_id.insert(id.clone(), rec);
        }
    }

    let mut f = std::fs::File::create(output_tbl).map_err(MycoNoteError::Io)?;
    let mut current_seqid = String::new();
    let mut feature_count = 0usize;
    let mut emitted_cds_parents: std::collections::HashSet<String> =
        std::collections::HashSet::new();

    for rec in &records {
        if rec.seqid != current_seqid {
            if !current_seqid.is_empty() {
                writeln!(f).map_err(MycoNoteError::Io)?;
            }
            writeln!(f, ">Feature {}", rec.seqid).map_err(MycoNoteError::Io)?;
            current_seqid = rec.seqid.clone();
        }

        match rec.feature_type.as_str() {
            "gene" => {
                let (s, e) = coords_for_strand(rec);
                writeln!(f, "{}\t{}\tgene", s, e).map_err(MycoNoteError::Io)?;

                if let Some(locus_tag) = rec.attributes.get("locus_tag") {
                    writeln!(f, "\t\t\tlocus_tag\t{}", locus_tag).map_err(MycoNoteError::Io)?;
                }
                if let Some(gene_name) = rec.attributes.get("Name") {
                    writeln!(f, "\t\t\tgene\t{}", gene_name).map_err(MycoNoteError::Io)?;
                }
                // Surface cross-references written by annotate (Dbxref=UniProtKB:...).
                // NCBI's .tbl syntax uses `db_xref` at the gene level so downstream
                // asnval / tbl2asn accept the reference.
                if let Some(dbxref) = rec.attributes.get("Dbxref") {
                    for xref in dbxref.split(',') {
                        let xref = xref.trim();
                        if !xref.is_empty() {
                            writeln!(f, "\t\t\tdb_xref\t{}", xref).map_err(MycoNoteError::Io)?;
                        }
                    }
                }
                feature_count += 1;
            }
            "mRNA" => {
                // Emit a CDS block once per mRNA using all its CDS segments.
                // We key on the mRNA's ID (which is what the child CDS rows
                // list as their Parent). Segments are ordered 5'→3'.
                let mrna_id = match rec.attributes.get("ID") {
                    Some(id) => id.clone(),
                    None => continue,
                };
                if !emitted_cds_parents.insert(mrna_id.clone()) {
                    continue;
                }
                let Some(cdss) = cds_by_parent.get(&mrna_id) else {
                    continue;
                };
                if cdss.is_empty() {
                    continue;
                }
                emit_joined_feature(&mut f, cdss, "CDS", rec, &by_id, config)?;
                feature_count += 1;
            }
            "tRNA" => {
                let (s, e) = coords_for_strand(rec);
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
                let (s, e) = coords_for_strand(rec);
                writeln!(f, "{}\t{}\trRNA", s, e).map_err(MycoNoteError::Io)?;

                let product = rec.attributes.get("product").cloned().unwrap_or_default();
                writeln!(f, "\t\t\tproduct\t{}", product).map_err(MycoNoteError::Io)?;
                feature_count += 1;
            }
            _ => {
                // Fall-through: GFF3s produced by some pipelines omit the
                // mRNA feature row, so we'd never hit the `"mRNA"` arm above
                // and CDS segments would be silently skipped. If we see a
                // bare CDS whose Parent we haven't emitted yet, emit the
                // joined block now using this record as the metadata source.
                if rec.feature_type == "CDS" {
                    let parent = match rec.attributes.get("Parent") {
                        Some(p) => p.split(',').next().unwrap_or(p).to_string(),
                        None => continue,
                    };
                    if !emitted_cds_parents.insert(parent.clone()) {
                        continue;
                    }
                    if let Some(cdss) = cds_by_parent.get(&parent) {
                        emit_joined_feature(&mut f, cdss, "CDS", rec, &by_id, config)?;
                        feature_count += 1;
                    }
                }
            }
        }
    }

    Ok(feature_count)
}

/// Return (first, second) column for an NCBI .tbl interval line, flipping
/// the order on the minus strand per NCBI spec.
fn coords_for_strand(rec: &GFFRecord) -> (u64, u64) {
    if rec.strand == '-' {
        (rec.end, rec.start)
    } else {
        (rec.start, rec.end)
    }
}

/// Walk up the GFF3 Parent chain from `start` looking for the first record
/// whose attributes contain `key`. Annotate writes functional qualifiers
/// (product, Name, locus_tag, Ontology_term, Dbxref) on the `gene` row only,
/// so a CDS → mRNA → gene walk is needed to surface them in the emitted
/// feature table.
fn attr_from_parent_chain<'a>(
    by_id: &'a std::collections::HashMap<String, &'a GFFRecord>,
    start: &'a GFFRecord,
    key: &str,
) -> Option<String> {
    if let Some(v) = start.attributes.get(key) {
        return Some(v.clone());
    }
    let mut current = start;
    for _ in 0..8 {
        let parent_attr = match current.attributes.get("Parent") {
            Some(p) => p,
            None => return None,
        };
        let parent_id = parent_attr.split(',').next().unwrap_or(parent_attr);
        let parent_rec = match by_id.get(parent_id) {
            Some(r) => *r,
            None => return None,
        };
        if let Some(v) = parent_rec.attributes.get(key) {
            return Some(v.clone());
        }
        current = parent_rec;
    }
    None
}

/// Emit a multi-segment feature block (CDS or exon-set) using NCBI's
/// implicit `join()`: header line has `start\tend\tFEATURE`, subsequent
/// interval lines are just `start\tend`, and qualifiers come after the last
/// interval. `meta_source` supplies the `product`, `codon_start`, etc. —
/// typically the mRNA row or (fallback) the first CDS row. `by_id` lets the
/// function walk up the Parent chain to find product/Name on the gene row
/// when the mRNA/CDS don't carry them.
fn emit_joined_feature<W: Write>(
    f: &mut W,
    segments: &[&GFFRecord],
    feature_key: &str,
    meta_source: &GFFRecord,
    by_id: &std::collections::HashMap<String, &GFFRecord>,
    config: &SubmitConfig,
) -> Result<()> {
    if segments.is_empty() {
        return Ok(());
    }
    let strand = segments[0].strand;

    // Order segments in 5'→3' reading direction: ascending on '+', descending
    // on '-'. Coordinates within each interval are also flipped on '-'.
    let mut sorted: Vec<&GFFRecord> = segments.to_vec();
    if strand == '-' {
        sorted.sort_by(|a, b| b.start.cmp(&a.start));
    } else {
        sorted.sort_by_key(|r| r.start);
    }

    for (i, seg) in sorted.iter().enumerate() {
        let (s, e) = coords_for_strand(seg);
        if i == 0 {
            writeln!(f, "{}\t{}\t{}", s, e, feature_key).map_err(MycoNoteError::Io)?;
        } else {
            writeln!(f, "{}\t{}", s, e).map_err(MycoNoteError::Io)?;
        }
    }

    // Qualifiers — emit once per feature, after all interval lines.
    // Look on the provided meta_source (mRNA row), then the first CDS
    // segment, then walk up the Parent chain to the gene row (where
    // annotate actually writes `product=` in annotated.gff3).
    let product = meta_source
        .attributes
        .get("product")
        .cloned()
        .or_else(|| {
            sorted
                .first()
                .and_then(|r| r.attributes.get("product").cloned())
        })
        .or_else(|| attr_from_parent_chain(by_id, meta_source, "product"))
        .or_else(|| {
            sorted
                .first()
                .and_then(|r| attr_from_parent_chain(by_id, r, "product"))
        })
        .unwrap_or_else(|| "hypothetical protein".to_string());
    writeln!(f, "\t\t\tproduct\t{}", product).map_err(MycoNoteError::Io)?;

    // codon_start defaults to 1 but can be overridden by the first segment's
    // GFF3 phase (a phase of 0 → codon_start=1, 1 → 2, 2 → 3).
    if let Some(first) = sorted.first() {
        if let Some(phase) = first.phase {
            writeln!(f, "\t\t\tcodon_start\t{}", phase + 1).map_err(MycoNoteError::Io)?;
        }
    }

    writeln!(f, "\t\t\ttransl_table\t{}", config.genetic_code).map_err(MycoNoteError::Io)?;
    Ok(())
}

/// Write NCBI submission template (.sbt) for table2asn.
///
/// The template needs a valid contact name, institution, and email or
/// table2asn rejects the `.sbt` with `Error loading template file`
/// (F7). We now wire the CLI values (`--contact-first`, `--contact-last`,
/// `--institution`, `--email`) into the ASN.1 fields and emit a warning
/// when any of them is missing so users know table2asn will fail on the
/// stub before they discover it through a cryptic error.
pub fn write_submission_template(output: &Path, config: &SubmitConfig) -> Result<()> {
    let missing: Vec<&str> = [
        ("--contact-first", &config.contact_first),
        ("--contact-last", &config.contact_last),
        ("--institution", &config.institution),
        ("--email", &config.email),
    ]
    .iter()
    .filter_map(|(flag, val)| if val.is_empty() { Some(*flag) } else { None })
    .collect();
    if !missing.is_empty() {
        eprintln!(
            "  ⚠ template.sbt: missing fields {} — table2asn will reject the stub. \
             Provide via CLI flags to make submit work end-to-end.",
            missing.join(", ")
        );
    }

    let mut f = std::fs::File::create(output).map_err(MycoNoteError::Io)?;

    writeln!(f, "Submit-block ::= {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "  contact {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "    contact {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "      name name {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "        last \"{}\",", escape_asn_string(&config.contact_last))
        .map_err(MycoNoteError::Io)?;
    writeln!(f, "        first \"{}\"", escape_asn_string(&config.contact_first))
        .map_err(MycoNoteError::Io)?;
    writeln!(f, "      }},").map_err(MycoNoteError::Io)?;
    writeln!(f, "      affil std {{").map_err(MycoNoteError::Io)?;
    writeln!(f, "        affil \"{}\",", escape_asn_string(&config.institution))
        .map_err(MycoNoteError::Io)?;
    writeln!(f, "        email \"{}\"", escape_asn_string(&config.email))
        .map_err(MycoNoteError::Io)?;
    writeln!(f, "      }}").map_err(MycoNoteError::Io)?;
    writeln!(f, "    }}").map_err(MycoNoteError::Io)?;
    writeln!(f, "  }}").map_err(MycoNoteError::Io)?;
    writeln!(f, "}}").map_err(MycoNoteError::Io)?;

    Ok(())
}

/// Escape a string for embedding inside an ASN.1 double-quoted value.
/// Doubles internal quotes, drops embedded newlines.
fn escape_asn_string(s: &str) -> String {
    s.replace('"', "\"\"").replace(['\n', '\r'], " ")
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write_gff(label: &str, body: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!(
            "myconote_submit_test_{}_{}.gff3",
            std::process::id(),
            label
        ));
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, "##gff-version 3").unwrap();
        f.write_all(body.as_bytes()).unwrap();
        p
    }

    fn tbl_path(label: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!(
            "myconote_submit_test_{}_{}.tbl",
            std::process::id(),
            label
        ));
        p
    }

    fn default_cfg() -> SubmitConfig {
        SubmitConfig {
            genetic_code: 1,
            ..SubmitConfig::default()
        }
    }

    #[test]
    fn multi_exon_plus_strand_is_joined_into_one_cds_block() {
        // One mRNA with three CDS exons. The old bug emitted three separate
        // CDS blocks each with their own qualifiers; correct output has one
        // header line + two bare interval lines + one product/codon block.
        let gff = "\
NW_1\tmaker\tgene\t100\t900\t.\t+\t.\tID=g1;Name=G1;locus_tag=MYCO_0001\n\
NW_1\tmaker\tmRNA\t100\t900\t.\t+\t.\tID=g1.mRNA;Parent=g1;product=ABC transporter\n\
NW_1\tmaker\tCDS\t100\t200\t.\t+\t0\tID=cds1a;Parent=g1.mRNA\n\
NW_1\tmaker\tCDS\t300\t500\t.\t+\t0\tID=cds1b;Parent=g1.mRNA\n\
NW_1\tmaker\tCDS\t700\t900\t.\t+\t0\tID=cds1c;Parent=g1.mRNA\n";
        let gff_path = write_gff("plus", gff);
        let tbl = tbl_path("plus");
        let n = write_feature_table(&gff_path, &tbl, &default_cfg()).unwrap();
        let out = std::fs::read_to_string(&tbl).unwrap();

        // Exactly one gene block + one CDS block.
        assert_eq!(
            n, 2,
            "expected 2 features (gene + joined CDS), got {}: {}",
            n, out
        );

        // Header line is the first CDS interval.
        assert!(
            out.contains("100\t200\tCDS"),
            "missing CDS header line:\n{}",
            out
        );
        // Subsequent intervals appear as bare start\tend lines — no CDS keyword.
        assert!(
            out.contains("\n300\t500\n"),
            "second segment not emitted as bare interval:\n{}",
            out
        );
        assert!(
            out.contains("\n700\t900\n"),
            "third segment not emitted as bare interval:\n{}",
            out
        );
        // Product qualifier appears exactly once (after last interval).
        assert_eq!(
            out.matches("\t\t\tproduct\tABC transporter").count(),
            1,
            "product qualifier should appear exactly once:\n{}",
            out
        );
        // transl_table also appears exactly once.
        assert_eq!(
            out.matches("\t\t\ttransl_table\t1").count(),
            1,
            "transl_table must be emitted once:\n{}",
            out
        );

        let _ = std::fs::remove_file(&gff_path);
        let _ = std::fs::remove_file(&tbl);
    }

    #[test]
    fn multi_exon_minus_strand_joins_in_reverse_reading_order() {
        // On the minus strand, segments are emitted 5'→3' (descending by
        // coordinate) and each interval has its larger coord in column 1.
        let gff = "\
NW_1\tmaker\tgene\t100\t900\t.\t-\t.\tID=g2;locus_tag=MYCO_0002\n\
NW_1\tmaker\tmRNA\t100\t900\t.\t-\t.\tID=g2.mRNA;Parent=g2;product=reverse gene\n\
NW_1\tmaker\tCDS\t100\t200\t.\t-\t0\tID=cds2a;Parent=g2.mRNA\n\
NW_1\tmaker\tCDS\t700\t900\t.\t-\t0\tID=cds2b;Parent=g2.mRNA\n";
        let gff_path = write_gff("minus", gff);
        let tbl = tbl_path("minus");
        write_feature_table(&gff_path, &tbl, &default_cfg()).unwrap();
        let out = std::fs::read_to_string(&tbl).unwrap();

        // Gene line uses flipped coords on '-'.
        assert!(
            out.contains("900\t100\tgene"),
            "gene line should be 'end start' on - strand:\n{}",
            out
        );
        // CDS header should be the downstream (in genomic terms, higher-coord) segment first.
        assert!(
            out.contains("900\t700\tCDS"),
            "first (5'-most on - strand) CDS segment should lead:\n{}",
            out
        );
        assert!(
            out.contains("\n200\t100\n"),
            "second CDS segment (lower genomic coords) should follow as flipped interval:\n{}",
            out
        );

        let _ = std::fs::remove_file(&gff_path);
        let _ = std::fs::remove_file(&tbl);
    }

    #[test]
    fn single_exon_cds_still_emits_one_block() {
        let gff = "\
NW_1\tmaker\tgene\t10\t100\t.\t+\t.\tID=g3;locus_tag=MYCO_0003\n\
NW_1\tmaker\tmRNA\t10\t100\t.\t+\t.\tID=g3.mRNA;Parent=g3;product=single-exon thing\n\
NW_1\tmaker\tCDS\t10\t100\t.\t+\t0\tID=cds3;Parent=g3.mRNA\n";
        let gff_path = write_gff("single", gff);
        let tbl = tbl_path("single");
        let n = write_feature_table(&gff_path, &tbl, &default_cfg()).unwrap();
        let out = std::fs::read_to_string(&tbl).unwrap();

        assert_eq!(n, 2);
        assert!(out.contains("10\t100\tCDS"));
        assert_eq!(
            out.matches("\tCDS").count(),
            1,
            "single-exon CDS should emit exactly one header:\n{}",
            out
        );

        let _ = std::fs::remove_file(&gff_path);
        let _ = std::fs::remove_file(&tbl);
    }

    #[test]
    fn cds_without_mrna_row_is_still_joined_via_parent_fallback() {
        // Some tools write a gene + CDSes but no mRNA row. The fallback path
        // in the `_ =>` arm should still group CDSes by Parent.
        let gff = "\
NW_1\tmaker\tgene\t100\t500\t.\t+\t.\tID=g4;locus_tag=MYCO_0004\n\
NW_1\tmaker\tCDS\t100\t200\t.\t+\t0\tID=cds4a;Parent=g4;product=no-mRNA gene\n\
NW_1\tmaker\tCDS\t300\t500\t.\t+\t0\tID=cds4b;Parent=g4\n";
        let gff_path = write_gff("nomrna", gff);
        let tbl = tbl_path("nomrna");
        write_feature_table(&gff_path, &tbl, &default_cfg()).unwrap();
        let out = std::fs::read_to_string(&tbl).unwrap();

        assert!(out.contains("100\t200\tCDS"));
        assert!(out.contains("\n300\t500\n"));
        // Exactly one CDS header.
        assert_eq!(
            out.matches("\tCDS").count(),
            1,
            "Parent-fallback must still emit one joined CDS:\n{}",
            out
        );

        let _ = std::fs::remove_file(&gff_path);
        let _ = std::fs::remove_file(&tbl);
    }
}
