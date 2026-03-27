/// Functional annotation pipeline
///
/// Assigns biological function to predicted genes by:
///   1. MMseqs2 homology search against UniProt/Swiss-Prot
///   2. hmmscan domain search against Pfam-A
///   3. GO term assignment from UniProt hits
///   4. BUSCO completeness assessment
///
/// Each step is optional — the pipeline degrades gracefully if a tool
/// is absent or a database has not been downloaded.
///
/// Download databases with: `myconote annotate --download-dbs`

pub mod mmseqs;
pub mod pfam;
pub mod go;
pub mod busco;
pub mod db;
pub mod interproscan;
pub mod eggnog;
pub mod cazyme;
pub mod secretome;
pub mod antismash;
pub mod merops;

use crate::utils::error::{MycoNoteError, Result};
use crate::predict::kingdom::Kingdom;
use crate::progress;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct AnnotateConfig {
    /// Predicted gene GFF3 (output of `myconote predict`)
    pub gff:          PathBuf,
    /// Genome FASTA (masked or unmasked)
    pub fasta:        PathBuf,
    /// Output directory
    pub out_dir:      PathBuf,
    /// Kingdom (affects BUSCO lineage selection)
    pub kingdom:      Kingdom,
    /// Locus tag prefix (must match the one used in `predict`)
    pub locus_prefix: String,
    /// Organism name for report header
    pub organism:     Option<String>,
    /// Parallel threads
    pub threads:      usize,
    /// Run MMseqs2 homology search
    pub run_mmseqs:   bool,
    /// Run hmmscan Pfam domain search
    pub run_pfam:     bool,
    /// Run BUSCO completeness check
    pub run_busco:    bool,
    /// Path to Swiss-Prot MMseqs2 database (auto-detected from db_dir)
    pub swissprot_db: Option<PathBuf>,
    /// Path to Pfam-A HMM database (auto-detected from db_dir)
    pub pfam_db:      Option<PathBuf>,
    /// Directory where databases are stored (default: ~/.myconote/dbs)
    pub db_dir:       PathBuf,
    /// E-value cutoff for MMseqs2 and hmmscan
    pub evalue:       f64,
    /// Minimum sequence identity for MMseqs2 hits (0.0–1.0)
    pub min_identity: f64,
    /// Run InterProScan via EBI REST API (requires internet; very thorough)
    pub run_interproscan: bool,
    /// Email address for EBI InterProScan submissions (required by EBI)
    pub interproscan_email: String,

    // ── New annotation modules ────────────────────────────────────────────
    /// Run EggNog-mapper (COG/NOG functional categories)
    pub run_eggnog:       bool,
    /// Path to EggNog-mapper database dir (default: auto-detect)
    pub eggnog_db:        Option<PathBuf>,
    /// Pre-computed emapper.annotations file (skip running emapper)
    pub eggnog_results:   Option<PathBuf>,

    /// Run CAZyme annotation (dbCAN / DIAMOND vs dbCAN database)
    pub run_cazyme:       bool,
    /// Path to dbCAN diamond database (.dmnd) for fallback DIAMOND search
    pub cazyme_db:        Option<PathBuf>,

    /// Run secretome prediction (SignalP + TMHMM)
    pub run_secretome:    bool,
    /// SignalP organism type: "euk" (default), "gram+", "gram-"
    pub signalp_organism: String,

    /// Run antiSMASH BGC cluster prediction
    pub run_antismash:    bool,
    /// Pre-computed antiSMASH output directory (skip running antiSMASH)
    pub antismash_dir:    Option<PathBuf>,
    /// antiSMASH taxon: "fungi" | "bacteria" | "plants"
    pub antismash_taxon:  String,

    /// Run MEROPS protease annotation (DIAMOND vs merops.dmnd)
    pub run_merops:       bool,
    /// Path to MEROPS DIAMOND database (auto-detected from db_dir)
    pub merops_db:        Option<PathBuf>,
}

impl Default for AnnotateConfig {
    fn default() -> Self {
        let db_dir = dirs_home().join(".myconote").join("dbs");
        Self {
            gff:          PathBuf::new(),
            fasta:        PathBuf::new(),
            out_dir:      PathBuf::from("annotate_out"),
            kingdom:      Kingdom::Fungi,
            locus_prefix: "GENE".to_string(),
            organism:     None,
            threads:      4,
            run_mmseqs:   true,
            run_pfam:     true,
            run_busco:    true,
            swissprot_db: None,
            pfam_db:      None,
            db_dir,
            evalue:              1e-5,
            run_interproscan:    false,
            interproscan_email:  String::new(),
            min_identity:        0.3,

            run_eggnog:       false,
            eggnog_db:        None,
            eggnog_results:   None,

            run_cazyme:       false,
            cazyme_db:        None,

            run_secretome:    false,
            signalp_organism: "euk".to_string(),

            run_antismash:    false,
            antismash_dir:    None,
            antismash_taxon:  "fungi".to_string(),

            run_merops:       false,
            merops_db:        None,
        }
    }
}

fn dirs_home() -> PathBuf {
    std::env::var("HOME")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/tmp"))
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-gene annotation record
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct GeneAnnotation {
    pub locus_tag:    String,
    /// Best Swiss-Prot hit description
    pub product:      Option<String>,
    /// Matched UniProt accession
    pub uniprot_acc:  Option<String>,
    /// Swiss-Prot hit identity (0–100)
    pub identity:     Option<f64>,
    /// e-value of best hit
    pub evalue:       Option<f64>,
    /// GO terms assigned (from UniProt mapping)
    pub go_terms:     Vec<String>,
    /// Pfam domain IDs found
    pub pfam_domains: Vec<String>,
    /// Whether BUSCO marked this as complete (if applicable)
    pub busco_status: Option<String>,
    /// InterPro accessions from InterProScan (IPR...)
    pub ipr_accessions: Vec<String>,
    /// Additional databases from InterProScan (TIGRFAM, Gene3D, etc.)
    pub ipr_databases:  Vec<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Annotation results container
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
pub struct AnnotationResults {
    pub genes:        HashMap<String, GeneAnnotation>,
    pub busco_summary: Option<busco::BuscoSummary>,
}

impl AnnotationResults {
    /// Annotated fraction (genes with a product description)
    pub fn annotated_fraction(&self) -> f64 {
        if self.genes.is_empty() { return 0.0; }
        let ann = self.genes.values().filter(|g| g.product.is_some()).count();
        ann as f64 / self.genes.len() as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Pipeline entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run the full functional annotation pipeline.
pub fn run_annotation(config: &AnnotateConfig) -> Result<AnnotationResults> {
    // ── Validate inputs ───────────────────────────────────────────────────────
    if !config.gff.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "GFF3 not found: {}", config.gff.display()
        )));
    }
    if !config.fasta.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "FASTA not found: {}", config.fasta.display()
        )));
    }

    std::fs::create_dir_all(&config.out_dir).map_err(MycoNoteError::Io)?;

    println!("── Functional annotation ────────────────────────────────────");
    println!("  GFF3    : {}", config.gff.display());
    println!("  FASTA   : {}", config.fasta.display());
    println!("  Output  : {}", config.out_dir.display());

    // ── Count total steps ─────────────────────────────────────────────────────
    let total_steps = 1  // protein extraction (always)
        + if config.run_mmseqs       { 1 } else { 0 }
        + if config.run_interproscan { 1 } else { 0 }
        + if config.run_pfam         { 1 } else { 0 }
        + if config.run_busco        { 1 } else { 0 };
    let mut step = 0usize;

    // ── Extract CDS protein sequences ─────────────────────────────────────────
    step += 1;
    progress::step(step, total_steps, "Extracting protein sequences…");
    let pb = progress::spinner("Translating CDS features…");
    let proteins_fa = config.out_dir.join("proteins.fa");
    let gene_ids = extract_proteins(&config.gff, &config.fasta, &proteins_fa)?;
    progress::finish_spinner(&pb, format!("{} proteins extracted", gene_ids.len()));

    let mut results = AnnotationResults::default();
    for id in &gene_ids {
        results.genes.insert(id.clone(), GeneAnnotation {
            locus_tag: id.clone(),
            ..Default::default()
        });
    }

    // ── MMseqs2 homology ──────────────────────────────────────────────────────
    if config.run_mmseqs {
        step += 1;
        progress::step(step, total_steps, "MMseqs2 homology (Swiss-Prot)…");
        let swissprot = resolve_swissprot(config);
        match swissprot {
            Some(db) => {
                let pb2 = progress::spinner("Running mmseqs easy-search…");
                let hits_tsv = config.out_dir.join("mmseqs_hits.tsv");
                match mmseqs::run(&proteins_fa, &db, &hits_tsv, config) {
                    Ok(hits) => {
                        progress::finish_spinner(&pb2, format!("{} hits found", hits.len()));
                        merge_mmseqs_hits(&mut results, hits);
                    }
                    Err(e) => progress::warn_spinner(&pb2, format!("MMseqs2 failed: {}", e)),
                }
            }
            None => {
                eprintln!("  ⚠  Swiss-Prot database not found.");
                eprintln!("     Run: myconote annotate --download-dbs");
            }
        }
    }

    // ── InterProScan (EBI REST API) ───────────────────────────────────────────
    if config.run_interproscan {
        step += 1;
        progress::step(step, total_steps, "InterProScan (EBI REST API)…");
        if config.interproscan_email.is_empty() {
            eprintln!("  ⚠  InterProScan requires --email <address> (EBI policy)");
        } else {
            let ipr_tsv = config.out_dir.join("interproscan_hits.tsv");
            match interproscan::run(&proteins_fa, &ipr_tsv, &config.interproscan_email) {
                Ok(ipr_map) => {
                    let n_with_hits = ipr_map.values().filter(|v| !v.is_empty()).count();
                    println!("  InterProScan: {}/{} proteins have hits", n_with_hits, gene_ids.len());
                    merge_interproscan_hits(&mut results, ipr_map);
                }
                Err(e) => eprintln!("  ⚠  InterProScan failed: {}", e),
            }
        }
    }

    // ── GO term assignment (from MMseqs2 UniProt hits) ────────────────────────
    {
        let go_tsv = config.out_dir.join("go_terms.tsv");
        let uniprot_accs: Vec<String> = results.genes.values()
            .filter_map(|g| g.uniprot_acc.clone())
            .collect();
        if !uniprot_accs.is_empty() {
            let pb2 = progress::spinner("Fetching GO terms from UniProt…");
            match go::assign_go_terms(&uniprot_accs, &go_tsv) {
                Ok(go_map) => {
                    progress::finish_spinner(&pb2, format!("GO terms for {} genes", go_map.len()));
                    merge_go_terms(&mut results, go_map);
                }
                Err(e) => progress::warn_spinner(&pb2, format!("GO assignment failed: {}", e)),
            }
        }
    }

    // ── Pfam domain search ────────────────────────────────────────────────────
    if config.run_pfam {
        step += 1;
        progress::step(step, total_steps, "Pfam domain search (hmmscan)…");
        let pfam_db = resolve_pfam(config);
        match pfam_db {
            Some(db) => {
                let pb2 = progress::spinner("Running hmmscan…");
                let pfam_out = config.out_dir.join("pfam_hits.tsv");
                match pfam::run(&proteins_fa, &db, &pfam_out, config) {
                    Ok(domains) => {
                        progress::finish_spinner(&pb2, format!("{} domain hits", domains.len()));
                        merge_pfam_domains(&mut results, domains);
                    }
                    Err(e) => progress::warn_spinner(&pb2, format!("Pfam search failed: {}", e)),
                }
            }
            None => {
                eprintln!("  ⚠  Pfam-A database not found.");
                eprintln!("     Run: myconote annotate --download-dbs");
            }
        }
    }

    // ── BUSCO completeness ────────────────────────────────────────────────────
    if config.run_busco {
        step += 1;
        let lineage = config.kingdom.busco_lineage();
        progress::step(step, total_steps, &format!("BUSCO completeness ({})…", lineage));
        let pb2 = progress::spinner(format!("Running BUSCO (lineage: {})…", lineage));
        let busco_dir = config.out_dir.join("busco");
        match busco::run(&proteins_fa, lineage, &busco_dir, config.threads) {
            Ok(summary) => {
                progress::finish_spinner(&pb2, format!(
                    "BUSCO: {:.1}% complete ({} single, {} dup, {} missing)",
                    summary.percent_complete(), summary.single_copy,
                    summary.duplicated, summary.missing,
                ));
                results.busco_summary = Some(summary);
            }
            Err(e) => progress::warn_spinner(&pb2, format!("BUSCO failed (non-fatal): {}", e)),
        }
    }

    // ── MEROPS protease annotation ────────────────────────────────────────────
    if config.run_merops {
        let merops_dir = config.out_dir.join("merops");
        let db_dir = config.merops_db.as_ref()
            .and_then(|p| p.parent().map(|d| d.to_path_buf()))
            .unwrap_or_else(|| config.db_dir.clone());
        let pb2 = progress::spinner("Running MEROPS protease annotation…");
        match merops::run_merops(&proteins_fa, &db_dir, &merops_dir, config.threads, config.evalue) {
            Ok(hits_map) if !hits_map.is_empty() => {
                let merops_tsv = config.out_dir.join("merops_hits.tsv");
                let n = merops::write_merops_table(&hits_map, &merops_tsv)
                    .unwrap_or(0);
                progress::finish_spinner(&pb2, format!("{} protease genes annotated", n));
                merops::print_merops_summary(&hits_map);
            }
            Ok(_) => progress::warn_spinner(&pb2, "MEROPS: no database found or no hits"),
            Err(e) => progress::warn_spinner(&pb2, format!("MEROPS failed (non-fatal): {}", e)),
        }
    }

    // ── Write annotated GFF3 ──────────────────────────────────────────────────
    let annotated_gff = config.out_dir.join("annotated.gff3");
    write_annotated_gff(&config.gff, &annotated_gff, &results)?;
    println!("  ✓  Annotated GFF3 → {}", annotated_gff.display());

    // ── Write functional annotation TSV ──────────────────────────────────────
    let annot_tsv = config.out_dir.join("annotations.tsv");
    write_annotation_tsv(&annot_tsv, &results)?;
    println!("  ✓  Annotation table → {}", annot_tsv.display());

    // ── Write summary report ──────────────────────────────────────────────────
    let report_path = config.out_dir.join("annotation_report.txt");
    write_report(&report_path, config, &results)?;
    println!("  ✓  Report → {}", report_path.display());
    println!("  Annotated: {:.1}% of genes have a functional description",
        results.annotated_fraction() * 100.0);

    Ok(results)
}

// ─────────────────────────────────────────────────────────────────────────────
// Protein extraction from GFF3 + FASTA
// ─────────────────────────────────────────────────────────────────────────────

/// Extract translated CDS sequences for all genes → proteins.fa
/// Returns list of locus_tag / gene IDs extracted.
fn extract_proteins(
    gff_path:  &Path,
    fasta_path: &Path,
    out_fa:    &Path,
) -> Result<Vec<String>> {
    use crate::parser::gff::GFFReader;
    use crate::parser::fasta::read_fasta_index;

    let fasta_index = read_fasta_index(fasta_path)?;
    let records: Vec<_> = GFFReader::from_path(gff_path)?
        .filter_map(|r| r.ok())
        .collect();

    // Collect CDS records grouped by gene (via Parent chain)
    let mut gene_cds: HashMap<String, Vec<_>> = HashMap::new();
    for rec in &records {
        if rec.feature_type != "CDS" { continue; }
        // Walk up Parent chain to find gene locus_tag
        let parent = match rec.parent() {
            Some(p) => p.clone(),
            None    => continue,
        };
        gene_cds.entry(parent).or_default().push(rec.clone());
    }

    let mut out = std::fs::File::create(out_fa).map_err(MycoNoteError::Io)?;
    let mut extracted: Vec<String> = Vec::new();

    for (mrna_id, mut cds_list) in gene_cds {
        // Find the gene parent of this mRNA
        let gene_id = records.iter()
            .find(|r| r.id().map(|id| id == &mrna_id).unwrap_or(false))
            .and_then(|r| r.parent())
            .cloned()
            .unwrap_or_else(|| mrna_id.clone());

        let locus_tag = records.iter()
            .find(|r| r.id().map(|id| id == &gene_id).unwrap_or(false))
            .and_then(|r| r.attributes.get("locus_tag"))
            .cloned()
            .unwrap_or_else(|| gene_id.clone());

        // Sort CDS by start position
        cds_list.sort_by_key(|r| r.start);

        // Get genome sequence for this seqid
        let seq_rec = match fasta_index.get(&cds_list[0].seqid) {
            Some(s) => s,
            None    => continue,
        };

        // Concatenate CDS bases (1-based inclusive)
        let mut cds_seq = String::new();
        let strand = cds_list[0].strand;
        for cds in &cds_list {
            let subseq = seq_rec.subsequence(cds.start, cds.end);
            cds_seq.push_str(subseq);
        }

        // Reverse complement if on minus strand
        if strand == '-' {
            cds_seq = crate::parser::fasta::reverse_complement(&cds_seq);
        }

        // Translate to protein
        let protein = translate_dna(&cds_seq);
        if protein.len() < 10 { continue; }  // skip very short ORFs

        writeln!(out, ">{}", locus_tag).map_err(MycoNoteError::Io)?;
        for chunk in protein.as_bytes().chunks(60) {
            writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                .map_err(MycoNoteError::Io)?;
        }

        extracted.push(locus_tag);
    }

    Ok(extracted)
}

/// Standard genetic code translation (stop = *)
fn translate_dna(dna: &str) -> String {
    let bytes = dna.as_bytes();
    let mut prot = String::with_capacity(bytes.len() / 3);
    let mut i = 0;
    while i + 2 < bytes.len() {
        let codon = [
            bytes[i].to_ascii_uppercase(),
            bytes[i+1].to_ascii_uppercase(),
            bytes[i+2].to_ascii_uppercase(),
        ];
        prot.push(codon_to_aa(&codon));
        i += 3;
    }
    // Remove trailing stop codon if present
    if prot.ends_with('*') { prot.pop(); }
    prot
}

fn codon_to_aa(c: &[u8; 3]) -> char {
    match c {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"AGT" | b"AGC" => 'S',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => 'X',
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Database path resolution
// ─────────────────────────────────────────────────────────────────────────────

fn resolve_swissprot(config: &AnnotateConfig) -> Option<PathBuf> {
    if let Some(ref p) = config.swissprot_db {
        if p.exists() { return Some(p.clone()); }
    }
    // Look in db_dir
    let candidates = [
        config.db_dir.join("swissprot").join("swissprot"),
        config.db_dir.join("uniprot_sprot.mmseqs"),
    ];
    candidates.into_iter().find(|p| p.exists())
}

fn resolve_pfam(config: &AnnotateConfig) -> Option<PathBuf> {
    if let Some(ref p) = config.pfam_db {
        if p.exists() { return Some(p.clone()); }
    }
    let candidates = [
        config.db_dir.join("pfam").join("Pfam-A.hmm"),
        config.db_dir.join("Pfam-A.hmm"),
    ];
    candidates.into_iter().find(|p| p.exists())
}

// ─────────────────────────────────────────────────────────────────────────────
// Result merging helpers
// ─────────────────────────────────────────────────────────────────────────────

fn merge_mmseqs_hits(results: &mut AnnotationResults, hits: Vec<mmseqs::MmseqsHit>) {
    for hit in hits {
        if let Some(gene) = results.genes.get_mut(&hit.query_id) {
            if gene.product.is_none() {
                gene.product      = Some(hit.description.clone());
                gene.uniprot_acc  = Some(hit.target_id.clone());
                gene.identity     = Some(hit.identity);
                gene.evalue       = Some(hit.evalue);
            }
        }
    }
}

fn merge_go_terms(results: &mut AnnotationResults, go_map: HashMap<String, Vec<String>>) {
    for gene in results.genes.values_mut() {
        if let Some(ref acc) = gene.uniprot_acc.clone() {
            if let Some(terms) = go_map.get(acc) {
                gene.go_terms = terms.clone();
            }
        }
    }
}

fn merge_pfam_domains(results: &mut AnnotationResults, domains: Vec<pfam::PfamHit>) {
    for hit in domains {
        if let Some(gene) = results.genes.get_mut(&hit.query_id) {
            if !gene.pfam_domains.contains(&hit.domain_id) {
                gene.pfam_domains.push(hit.domain_id.clone());
            }
        }
    }
}

fn merge_interproscan_hits(
    results: &mut AnnotationResults,
    ipr_map: HashMap<String, Vec<interproscan::IprScanResult>>,
) {
    for (protein_id, hits) in ipr_map {
        if let Some(gene) = results.genes.get_mut(&protein_id) {
            for hit in &hits {
                // Add InterPro accession if not already present
                if let Some(ref ipr_acc) = hit.ipr_accession {
                    if !gene.ipr_accessions.contains(ipr_acc) {
                        gene.ipr_accessions.push(ipr_acc.clone());
                    }
                }
                // Add database accession to ipr_databases
                let db_entry = format!("{}:{}", hit.database, hit.db_accession);
                if !gene.ipr_databases.contains(&db_entry) {
                    gene.ipr_databases.push(db_entry);
                }
                // Supplement GO terms from InterProScan (these are often more complete)
                for go in &hit.go_terms {
                    if !gene.go_terms.contains(go) {
                        gene.go_terms.push(go.clone());
                    }
                }
                // Add product description if not yet set
                if gene.product.is_none() && !hit.db_desc.is_empty()
                    && hit.db_desc != "Uncharacterised protein" {
                    gene.product = Some(hit.db_desc.clone());
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Output writers
// ─────────────────────────────────────────────────────────────────────────────

/// Copy input GFF3, adding product= and Dbxref= attributes to gene features.
fn write_annotated_gff(
    input_gff:  &Path,
    output_gff: &Path,
    results:    &AnnotationResults,
) -> Result<()> {
    use crate::parser::gff::GFFReader;

    let mut out = std::fs::File::create(output_gff).map_err(MycoNoteError::Io)?;
    writeln!(out, "##gff-version 3").map_err(MycoNoteError::Io)?;

    for rec_res in GFFReader::from_path(input_gff)? {
        let mut rec = match rec_res { Ok(r) => r, Err(_) => continue };

        if rec.feature_type == "gene" {
            let locus_tag = rec.attributes.get("locus_tag")
                .or_else(|| rec.attributes.get("ID"))
                .cloned()
                .unwrap_or_default();

            if let Some(ann) = results.genes.get(&locus_tag) {
                if let Some(ref product) = ann.product {
                    rec.attributes.insert("product".into(), product.clone());
                }
                // Build Dbxref from UniProt + InterPro accessions
                let mut dbxrefs: Vec<String> = Vec::new();
                if let Some(ref acc) = ann.uniprot_acc {
                    dbxrefs.push(format!("UniProtKB:{}", acc));
                }
                for ipr in &ann.ipr_accessions {
                    dbxrefs.push(format!("InterPro:{}", ipr));
                }
                if !dbxrefs.is_empty() {
                    rec.attributes.insert("Dbxref".into(), dbxrefs.join(","));
                }
                if !ann.go_terms.is_empty() {
                    rec.attributes.insert("Ontology_term".into(),
                        ann.go_terms.join(","));
                }
                // Note: Pfam domains + InterPro database entries
                let mut notes: Vec<String> = Vec::new();
                if !ann.pfam_domains.is_empty() {
                    notes.push(format!("Pfam:{}", ann.pfam_domains.join(",")));
                }
                if !ann.ipr_databases.is_empty() {
                    notes.push(format!("IPR_db:{}", ann.ipr_databases.join(",")));
                }
                if !notes.is_empty() {
                    rec.attributes.insert("Note".into(), notes.join("|"));
                }
            }
        }

        writeln!(out, "{}", rec.to_gff3_line()).map_err(MycoNoteError::Io)?;
    }

    Ok(())
}

/// Write tab-separated annotation table.
fn write_annotation_tsv(path: &Path, results: &AnnotationResults) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f,
        "locus_tag\tproduct\tuniprot_acc\tidentity\tevalue\tgo_terms\tpfam_domains\tipr_accessions\tipr_databases"
    ).map_err(MycoNoteError::Io)?;

    let mut genes: Vec<_> = results.genes.values().collect();
    genes.sort_by(|a, b| a.locus_tag.cmp(&b.locus_tag));

    for g in genes {
        writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g.locus_tag,
            g.product.as_deref().unwrap_or("hypothetical protein"),
            g.uniprot_acc.as_deref().unwrap_or(""),
            g.identity.map(|v| format!("{:.1}", v)).unwrap_or_default(),
            g.evalue.map(|v| format!("{:.2e}", v)).unwrap_or_default(),
            g.go_terms.join("|"),
            g.pfam_domains.join("|"),
            g.ipr_accessions.join("|"),
            g.ipr_databases.join("|"),
        ).map_err(MycoNoteError::Io)?;
    }

    Ok(())
}

/// Human-readable annotation report.
fn write_report(path: &Path, config: &AnnotateConfig, results: &AnnotationResults) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;

    let org = config.organism.as_deref().unwrap_or("Unknown organism");

    writeln!(f, "myconote annotate — Functional Annotation Report").map_err(MycoNoteError::Io)?;
    writeln!(f, "=================================================").map_err(MycoNoteError::Io)?;
    writeln!(f, "Organism  : {}", org).map_err(MycoNoteError::Io)?;
    writeln!(f, "Kingdom   : {}", config.kingdom.display_name()).map_err(MycoNoteError::Io)?;
    writeln!(f, "Input GFF : {}", config.gff.display()).map_err(MycoNoteError::Io)?;
    writeln!(f, "").map_err(MycoNoteError::Io)?;

    let total = results.genes.len();
    let annotated = results.genes.values().filter(|g| g.product.is_some()).count();
    let with_go   = results.genes.values().filter(|g| !g.go_terms.is_empty()).count();
    let with_pfam = results.genes.values().filter(|g| !g.pfam_domains.is_empty()).count();

    writeln!(f, "Gene totals").map_err(MycoNoteError::Io)?;
    writeln!(f, "  Total genes             : {}", total).map_err(MycoNoteError::Io)?;
    writeln!(f, "  With product description: {} ({:.1}%)",
        annotated, annotated as f64 / total.max(1) as f64 * 100.0)
        .map_err(MycoNoteError::Io)?;
    writeln!(f, "  With GO terms           : {} ({:.1}%)",
        with_go, with_go as f64 / total.max(1) as f64 * 100.0)
        .map_err(MycoNoteError::Io)?;
    writeln!(f, "  With Pfam domains       : {} ({:.1}%)",
        with_pfam, with_pfam as f64 / total.max(1) as f64 * 100.0)
        .map_err(MycoNoteError::Io)?;

    if let Some(ref b) = results.busco_summary {
        writeln!(f, "").map_err(MycoNoteError::Io)?;
        writeln!(f, "BUSCO ({} lineage)", b.lineage).map_err(MycoNoteError::Io)?;
        writeln!(f, "  Complete    : {} ({:.1}%)", b.complete(), b.percent_complete())
            .map_err(MycoNoteError::Io)?;
        writeln!(f, "    Single    : {}", b.single_copy).map_err(MycoNoteError::Io)?;
        writeln!(f, "    Duplicated: {}", b.duplicated).map_err(MycoNoteError::Io)?;
        writeln!(f, "  Fragmented  : {}", b.fragmented).map_err(MycoNoteError::Io)?;
        writeln!(f, "  Missing     : {}", b.missing).map_err(MycoNoteError::Io)?;
        writeln!(f, "  Total BUSCOs: {}", b.total).map_err(MycoNoteError::Io)?;
    }

    Ok(())
}
