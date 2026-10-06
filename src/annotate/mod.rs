pub mod antismash;
pub mod busco;
pub mod cazyme;
pub mod db;
pub mod eggnog;
pub mod genetic_code;
pub mod go;
pub mod interproscan;
pub mod merops;
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
pub mod secretome;
pub mod trnascan;

use crate::predict::kingdom::Kingdom;
use crate::progress;
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct AnnotateConfig {
    /// Predicted gene GFF3 (output of `myconote predict`)
    pub gff: PathBuf,
    /// Genome FASTA (masked or unmasked)
    pub fasta: PathBuf,
    /// Output directory
    pub out_dir: PathBuf,
    /// Kingdom (affects BUSCO lineage selection)
    pub kingdom: Kingdom,
    /// Locus tag prefix (must match the one used in `predict`)
    pub locus_prefix: String,
    /// Organism name for report header
    pub organism: Option<String>,
    /// Parallel threads
    pub threads: usize,
    /// Run MMseqs2 homology search
    pub run_mmseqs: bool,
    /// Run hmmscan Pfam domain search
    pub run_pfam: bool,
    /// Run BUSCO completeness check
    pub run_busco: bool,
    /// Path to Swiss-Prot MMseqs2 database (auto-detected from db_dir)
    pub swissprot_db: Option<PathBuf>,
    /// Path to Pfam-A HMM database (auto-detected from db_dir)
    pub pfam_db: Option<PathBuf>,
    /// Directory where databases are stored (default: ~/.myconote/dbs)
    pub db_dir: PathBuf,
    /// E-value cutoff for MMseqs2 and hmmscan
    pub evalue: f64,
    /// Minimum sequence identity for MMseqs2 hits (0.0–1.0)
    pub min_identity: f64,
    /// Run InterProScan via EBI REST API (requires internet; very thorough)
    pub run_interproscan: bool,
    /// Email address for EBI InterProScan submissions (required by EBI)
    pub interproscan_email: String,

    // ── New annotation modules ────────────────────────────────────────────
    /// Run EggNog-mapper (COG/NOG functional categories)
    pub run_eggnog: bool,
    /// Path to EggNog-mapper database dir (default: auto-detect)
    pub eggnog_db: Option<PathBuf>,
    /// Pre-computed emapper.annotations file (skip running emapper)
    pub eggnog_results: Option<PathBuf>,

    /// Run CAZyme annotation (dbCAN / DIAMOND vs dbCAN database)
    pub run_cazyme: bool,
    /// Path to dbCAN diamond database (.dmnd) for fallback DIAMOND search
    pub cazyme_db: Option<PathBuf>,

    /// Run secretome prediction (SignalP + TMHMM)
    pub run_secretome: bool,
    /// SignalP organism type: "euk" (default), "gram+", "gram-"
    pub signalp_organism: String,

    /// Run antiSMASH BGC cluster prediction
    pub run_antismash: bool,
    /// Pre-computed antiSMASH output directory (skip running antiSMASH)
    pub antismash_dir: Option<PathBuf>,
    /// antiSMASH taxon: "fungi" | "bacteria" | "plants"
    pub antismash_taxon: String,

    /// Run MEROPS protease annotation (DIAMOND vs merops.dmnd)
    pub run_merops: bool,
    /// Path to MEROPS DIAMOND database (auto-detected from db_dir)
    pub merops_db: Option<PathBuf>,

    // ── tRNA + genetic code ──────────────────────────────────────────────
    /// Run tRNAscan-SE for tRNA gene prediction
    pub run_trnascan: bool,
    /// tRNAscan mode: "eukaryotic" | "mitochondrial" | "general"
    pub trnascan_mode: String,
    /// Genetic code table (1=standard, 12=Candida CTG, etc.)
    pub genetic_code: u8,
}

impl Default for AnnotateConfig {
    fn default() -> Self {
        let db_dir = dirs_home().join(".myconote").join("dbs");
        Self {
            gff: PathBuf::new(),
            fasta: PathBuf::new(),
            out_dir: PathBuf::from("annotate_out"),
            kingdom: Kingdom::Fungi,
            locus_prefix: "GENE".to_string(),
            organism: None,
            threads: 4,
            run_mmseqs: true,
            run_pfam: true,
            run_busco: true,
            swissprot_db: None,
            pfam_db: None,
            db_dir,
            evalue: 1e-5,
            run_interproscan: false,
            interproscan_email: String::new(),
            min_identity: 0.3,

            run_eggnog: false,
            eggnog_db: None,
            eggnog_results: None,

            run_cazyme: false,
            cazyme_db: None,

            run_secretome: false,
            signalp_organism: "euk".to_string(),

            run_antismash: false,
            antismash_dir: None,
            antismash_taxon: "fungi".to_string(),

            run_merops: false,
            merops_db: None,

            run_trnascan: false,
            trnascan_mode: "eukaryotic".to_string(),
            genetic_code: 1,
        }
    }
}

fn dirs_home() -> PathBuf {
    std::env::var("HOME")
        .or_else(|_| std::env::var("USERPROFILE"))
        .map(PathBuf::from)
        .unwrap_or_else(|_| std::env::temp_dir())
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-gene annotation record
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct GeneAnnotation {
    pub locus_tag: String,
    /// Best Swiss-Prot hit description
    pub product: Option<String>,
    /// Matched UniProt accession
    pub uniprot_acc: Option<String>,
    /// Swiss-Prot hit identity (0–100)
    pub identity: Option<f64>,
    /// e-value of best hit
    pub evalue: Option<f64>,
    /// GO terms assigned (from UniProt mapping)
    pub go_terms: Vec<String>,
    /// Pfam domain IDs found
    pub pfam_domains: Vec<String>,
    /// Whether BUSCO marked this as complete (if applicable)
    pub busco_status: Option<String>,
    /// InterPro accessions from InterProScan (IPR...)
    pub ipr_accessions: Vec<String>,
    /// Additional databases from InterProScan (TIGRFAM, Gene3D, etc.)
    pub ipr_databases: Vec<String>,
    /// EggNog COG single-letter category (e.g. "J" for translation)
    pub cog_category: Option<String>,
    /// EggNog orthologous group (e.g. "COG0012@1|root")
    pub eggnog_og: Option<String>,
    /// CAZyme family assignments (e.g. ["GH5", "CBM1"]); multiple per gene possible
    pub cazyme_families: Vec<String>,
    /// Signal-peptide prediction: "SP(Sec/SPI)", "NO_SP", etc. (SignalP label)
    pub signal_peptide: Option<String>,
    /// Count of predicted transmembrane helices (from TMHMM / DeepTMHMM)
    pub tm_helices: Option<u32>,
    /// antiSMASH biosynthetic gene cluster ID containing this gene, if any
    pub bgc_cluster: Option<String>,
    /// Type of the BGC cluster (e.g. "T1PKS", "NRPS", "terpene")
    pub bgc_type: Option<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Annotation results container
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
pub struct AnnotationResults {
    pub genes: HashMap<String, GeneAnnotation>,
    pub busco_summary: Option<busco::BuscoSummary>,
}

impl AnnotationResults {
    /// Annotated fraction (genes with a product description)
    pub fn annotated_fraction(&self) -> f64 {
        if self.genes.is_empty() {
            return 0.0;
        }
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
            "GFF3 not found: {}",
            config.gff.display()
        )));
    }
    if !config.fasta.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "FASTA not found: {}",
            config.fasta.display()
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
        + if config.run_busco        { 1 } else { 0 }
        + if config.run_eggnog       { 1 } else { 0 }
        + if config.run_cazyme       { 1 } else { 0 }
        + if config.run_secretome    { 1 } else { 0 }
        + if config.run_antismash    { 1 } else { 0 };
    let mut step = 0usize;

    // ── Extract CDS protein sequences ─────────────────────────────────────────
    step += 1;
    progress::step(step, total_steps, "Extracting protein sequences…");
    let pb = progress::spinner("Translating CDS features…");
    let proteins_fa = config.out_dir.join("proteins.fa");
    let gene_ids = extract_proteins(
        &config.gff,
        &config.fasta,
        &proteins_fa,
        config.genetic_code,
    )?;
    progress::finish_spinner(&pb, format!("{} proteins extracted", gene_ids.len()));

    let mut results = AnnotationResults::default();
    for id in &gene_ids {
        results.genes.insert(
            id.clone(),
            GeneAnnotation {
                locus_tag: id.clone(),
                ..Default::default()
            },
        );
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
                    println!(
                        "  InterProScan: {}/{} proteins have hits",
                        n_with_hits,
                        gene_ids.len()
                    );
                    merge_interproscan_hits(&mut results, ipr_map);
                }
                Err(e) => eprintln!("  ⚠  InterProScan failed: {}", e),
            }
        }
    }

    // ── GO term assignment (from MMseqs2 UniProt hits) ────────────────────────
    {
        let go_tsv = config.out_dir.join("go_terms.tsv");
        let uniprot_accs: Vec<String> = results
            .genes
            .values()
            .filter_map(|g| g.uniprot_acc.clone())
            .collect();
        if !uniprot_accs.is_empty() {
            let pb2 = progress::spinner("Fetching GO terms from UniProt…");
            match go::assign_go_terms(&uniprot_accs, &go_tsv) {
                Ok(go_map) => {
                    let with_go = go_map.values().filter(|v| !v.is_empty()).count();
                    progress::finish_spinner(
                        &pb2,
                        format!(
                            "GO terms for {} / {} genes (via UniProt)",
                            with_go,
                            go_map.len()
                        ),
                    );
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
        progress::step(
            step,
            total_steps,
            &format!("BUSCO completeness ({})…", lineage),
        );
        let pb2 = progress::spinner(format!("Running BUSCO (lineage: {})…", lineage));
        let busco_dir = config.out_dir.join("busco");
        match busco::run(&proteins_fa, lineage, &busco_dir, config.threads) {
            Ok(summary) => {
                progress::finish_spinner(
                    &pb2,
                    format!(
                        "BUSCO: {:.1}% complete ({} single, {} dup, {} missing)",
                        summary.percent_complete(),
                        summary.single_copy,
                        summary.duplicated,
                        summary.missing,
                    ),
                );
                results.busco_summary = Some(summary);
            }
            Err(e) => progress::warn_spinner(&pb2, format!("BUSCO failed (non-fatal): {}", e)),
        }
    }

    // ── tRNAscan-SE tRNA prediction ──────────────────────────────────────────
    if config.run_trnascan {
        let pb_trna = progress::spinner("Running tRNAscan-SE…");
        let trna_dir = config.out_dir.join("trnascan");
        let trna_cfg = trnascan::TrnaScanConfig {
            mode: config.trnascan_mode.clone(),
            threads: config.threads,
            ..trnascan::TrnaScanConfig::default()
        };
        match trnascan::run_trnascan(&config.fasta, &trna_dir, &trna_cfg) {
            Ok(trna_result) => {
                progress::finish_spinner(
                    &pb_trna,
                    format!("tRNAscan-SE: {} tRNA genes found", trna_result.total),
                );
                trna_result.summarize();
                // Write tRNA GFF3
                let trna_gff = config.out_dir.join("trnascan.gff3");
                if let Err(e) =
                    trnascan::write_trna_gff3(&trna_result.trnas, &trna_gff, &config.locus_prefix)
                {
                    eprintln!("  Warning: failed to write tRNA GFF3: {}", e);
                }
            }
            Err(e) => {
                progress::warn_spinner(&pb_trna, format!("tRNAscan-SE failed (non-fatal): {}", e))
            }
        }
    }

    // ── MEROPS protease annotation ────────────────────────────────────────────
    if config.run_merops {
        let merops_dir = config.out_dir.join("merops");
        let db_dir = config
            .merops_db
            .as_ref()
            .and_then(|p| p.parent().map(|d| d.to_path_buf()))
            .unwrap_or_else(|| config.db_dir.clone());
        let pb2 = progress::spinner("Running MEROPS protease annotation…");
        match merops::run_merops(
            &proteins_fa,
            &db_dir,
            &merops_dir,
            config.threads,
            config.evalue,
        ) {
            Ok(hits_map) if !hits_map.is_empty() => {
                let merops_tsv = config.out_dir.join("merops_hits.tsv");
                let n = merops::write_merops_table(&hits_map, &merops_tsv).unwrap_or(0);
                progress::finish_spinner(&pb2, format!("{} protease genes annotated", n));
                merops::print_merops_summary(&hits_map);
            }
            Ok(_) => progress::warn_spinner(&pb2, "MEROPS: no database found or no hits"),
            Err(e) => progress::warn_spinner(&pb2, format!("MEROPS failed (non-fatal): {}", e)),
        }
    }

    // ── EggNog-mapper (COG/NOG orthology + KEGG + GO) ─────────────────────────
    if config.run_eggnog {
        step += 1;
        progress::step(step, total_steps, "EggNog-mapper (COG/NOG + KEGG)…");
        // Users can either point at an existing emapper.annotations file or
        // let us run emapper.py end-to-end. The pre-computed path is useful
        // on HPC where emapper takes hours and is typically run once per
        // genome via a separate submit script.
        let eggnog_dir = config.out_dir.join("eggnog");
        let pb2 = progress::spinner("Reading / running EggNog-mapper…");
        let annotations_path: Option<PathBuf> = if let Some(pre) = &config.eggnog_results {
            if pre.exists() {
                Some(pre.clone())
            } else {
                progress::warn_spinner(
                    &pb2,
                    format!("Pre-computed emapper file not found: {}", pre.display()),
                );
                None
            }
        } else if eggnog::emapper_available() {
            match eggnog::run_emapper(
                &proteins_fa,
                &eggnog_dir,
                config.eggnog_db.as_deref(),
                config.threads,
            ) {
                Ok(p) => Some(p),
                Err(e) => {
                    progress::warn_spinner(&pb2, format!("emapper failed: {}", e));
                    None
                }
            }
        } else {
            progress::warn_spinner(
                &pb2,
                "emapper.py not in PATH (install: conda install -c bioconda eggnog-mapper)",
            );
            None
        };

        if let Some(path) = annotations_path {
            match eggnog::parse_emapper_results(&path) {
                Ok(hits) => {
                    let summary_tsv = config.out_dir.join("eggnog_hits.tsv");
                    let _ = eggnog::write_eggnog_table(&hits, &summary_tsv);
                    for (gene_id, hit) in &hits {
                        if let Some(g) = results.genes.get_mut(gene_id) {
                            if !hit.cog_cat.is_empty() && hit.cog_cat != "-" {
                                g.cog_category = Some(hit.cog_cat.clone());
                            }
                            if !hit.best_og.is_empty() && hit.best_og != "-" {
                                g.eggnog_og = Some(hit.best_og.clone());
                            }
                            // Prefer an explicit EggNog product description when
                            // the gene has no Swiss-Prot hit (keeps well-curated
                            // mmseqs labels untouched).
                            if g.product.is_none()
                                && !hit.description.is_empty()
                                && hit.description != "-"
                            {
                                g.product = Some(hit.description.clone());
                            }
                            for go in &hit.go_terms {
                                if !go.is_empty() && !g.go_terms.contains(go) {
                                    g.go_terms.push(go.clone());
                                }
                            }
                        }
                    }
                    progress::finish_spinner(
                        &pb2,
                        format!("EggNog: {} proteins with COG assignments", hits.len()),
                    );
                }
                Err(e) => progress::warn_spinner(&pb2, format!("parse failed: {}", e)),
            }
        }
    }

    // ── CAZyme annotation (dbCAN preferred, DIAMOND fallback) ─────────────────
    if config.run_cazyme {
        step += 1;
        progress::step(step, total_steps, "CAZyme families (dbCAN / DIAMOND)…");
        let cazyme_dir = config.out_dir.join("cazyme");
        let pb2 = progress::spinner("Running CAZyme annotation…");
        let overview: Option<PathBuf> = if cazyme::dbcan_available() {
            match cazyme::run_dbcan(
                &proteins_fa,
                &cazyme_dir,
                Some(&config.db_dir),
                config.threads,
            ) {
                Ok(p) => Some(p),
                Err(e) => {
                    progress::warn_spinner(&pb2, format!("run_dbcan failed: {}", e));
                    None
                }
            }
        } else if let Some(db) = &config.cazyme_db {
            if db.exists() {
                match cazyme::run_diamond_cazyme(
                    &proteins_fa,
                    db,
                    &cazyme_dir,
                    config.threads,
                    config.evalue,
                ) {
                    Ok(p) => Some(p),
                    Err(e) => {
                        progress::warn_spinner(&pb2, format!("diamond CAZyme failed: {}", e));
                        None
                    }
                }
            } else {
                progress::warn_spinner(&pb2, format!("dbCAN db not found: {}", db.display()));
                None
            }
        } else {
            progress::warn_spinner(
                &pb2,
                "No dbCAN backend: install run_dbcan.py or pass --cazyme-db <dbCAN.dmnd>",
            );
            None
        };

        if let Some(path) = overview {
            // parse_dbcan_overview works for dbCAN overview.txt; for the DIAMOND
            // fallback we parse_diamond_cazyme into the same structure.
            let hits_result = if path.file_name().and_then(|s| s.to_str()) == Some("overview.txt") {
                cazyme::parse_dbcan_overview(&path)
            } else {
                cazyme::parse_diamond_cazyme(&path, 0.3)
            };
            match hits_result {
                Ok(hits) => {
                    let out_tsv = config.out_dir.join("cazyme_hits.tsv");
                    let _ = cazyme::write_cazyme_table(&hits, &out_tsv);
                    for (gene_id, families) in &hits {
                        if let Some(g) = results.genes.get_mut(gene_id) {
                            for fam in families {
                                let name = &fam.family;
                                if !g.cazyme_families.contains(name) {
                                    g.cazyme_families.push(name.clone());
                                }
                            }
                        }
                    }
                    let total: usize = hits.values().map(|v| v.len()).sum();
                    progress::finish_spinner(
                        &pb2,
                        format!("{} CAZyme annotations across {} genes", total, hits.len()),
                    );
                }
                Err(e) => progress::warn_spinner(&pb2, format!("CAZyme parse failed: {}", e)),
            }
        }
    }

    // ── Secretome (signal peptide + transmembrane) ────────────────────────────
    if config.run_secretome {
        step += 1;
        progress::step(step, total_steps, "Secretome (SignalP / DeepSig + TMHMM)…");
        let sec_dir = config.out_dir.join("secretome");
        let pb2 = progress::spinner("Running signal-peptide prediction…");
        match secretome::run_signalp(&proteins_fa, &sec_dir, &config.signalp_organism) {
            Ok(sp_path) => {
                match secretome::parse_signalp_output(&sp_path) {
                    Ok(sp_hits) => {
                        for (gene_id, sp) in &sp_hits {
                            if let Some(g) = results.genes.get_mut(gene_id) {
                                g.signal_peptide = Some(sp.prediction.clone());
                            }
                        }
                        // Fire TMHMM next; failures here are non-fatal so the
                        // signal-peptide data we just merged stays intact.
                        let tmhmm_res = secretome::run_tmhmm(&proteins_fa, &sec_dir)
                            .and_then(|tm_path| secretome::parse_tmhmm_output(&tm_path));
                        let tm_map = tmhmm_res.ok();
                        if let Some(tm_hits) = &tm_map {
                            for (gene_id, tm) in tm_hits {
                                if let Some(g) = results.genes.get_mut(gene_id) {
                                    g.tm_helices = Some(tm.tm_count as u32);
                                }
                            }
                        }
                        // A gene is "secreted" = has a signal peptide AND no
                        // downstream TM helix (so it actually exits the cell).
                        let secretome_set: std::collections::HashSet<String> = sp_hits
                            .iter()
                            .filter(|(id, sp)| {
                                if !sp.has_signal {
                                    return false;
                                }
                                let tm_count = tm_map
                                    .as_ref()
                                    .and_then(|m| m.get(*id))
                                    .map(|t| t.tm_count)
                                    .unwrap_or(0);
                                tm_count == 0
                            })
                            .map(|(id, _)| id.clone())
                            .collect();
                        let sec_tsv = config.out_dir.join("secretome_hits.tsv");
                        let _ = secretome::write_secretome_table(
                            &secretome_set,
                            &sp_hits,
                            tm_map.as_ref(),
                            &sec_tsv,
                        );
                        let secreted = sp_hits
                            .values()
                            .filter(|s| s.prediction.starts_with("SP"))
                            .count();
                        progress::finish_spinner(
                            &pb2,
                            format!(
                                "Signal peptides on {} / {} proteins",
                                secreted,
                                sp_hits.len()
                            ),
                        );
                    }
                    Err(e) => progress::warn_spinner(&pb2, format!("SignalP parse failed: {}", e)),
                }
            }
            Err(e) => {
                progress::warn_spinner(&pb2, format!("signal-peptide tool unavailable: {}", e))
            }
        }
    }

    // ── antiSMASH biosynthetic gene clusters ──────────────────────────────────
    if config.run_antismash {
        step += 1;
        progress::step(step, total_steps, "antiSMASH (biosynthetic gene clusters)…");
        let pb2 = progress::spinner("Loading / running antiSMASH…");
        // Users typically run antiSMASH separately (it's slow and stateful);
        // --antismash-dir lets them point at a pre-computed output directory.
        let as_dir: Option<PathBuf> = config.antismash_dir.clone();

        let as_dir = match as_dir {
            Some(d) if d.exists() => Some(d),
            Some(d) => {
                progress::warn_spinner(&pb2, format!("antiSMASH dir not found: {}", d.display()));
                None
            }
            None if antismash::antismash_available() => {
                // Auto-run: needs a GenBank input, which we don't build here
                // (submit writes .gbk). For now, emit a helpful message —
                // integrating the end-to-end auto-run requires feeding
                // submit's .gbk output back to annotate, which is a pipeline
                // re-order rather than a single-block fix.
                progress::warn_spinner(
                    &pb2,
                    "antiSMASH auto-run requires a .gbk input — run `myconote submit` first \
                     then re-run with --antismash-dir pointing at the antiSMASH output",
                );
                None
            }
            None => {
                progress::warn_spinner(
                    &pb2,
                    "antiSMASH not installed. Install: conda install -c bioconda antismash",
                );
                None
            }
        };

        if let Some(dir) = as_dir {
            match antismash::parse_antismash_gff(&dir) {
                Ok(clusters) => {
                    let out_tsv = config.out_dir.join("antismash_clusters.tsv");
                    let _ = antismash::write_bgc_table(&clusters, &out_tsv);
                    for cluster in &clusters {
                        for gene_id in &cluster.gene_ids {
                            if let Some(g) = results.genes.get_mut(gene_id) {
                                g.bgc_cluster = Some(cluster.cluster_id.clone());
                                g.bgc_type = Some(cluster.bgc_type.clone());
                            }
                        }
                    }
                    progress::finish_spinner(
                        &pb2,
                        format!("{} biosynthetic gene cluster(s) parsed", clusters.len()),
                    );
                }
                Err(e) => progress::warn_spinner(&pb2, format!("antiSMASH parse failed: {}", e)),
            }
        }
    }

    // ── Write annotated GFF3 ──────────────────────────────────────────────────
    let annotated_gff = config.out_dir.join("annotated.gff3");
    write_annotated_gff(&config.gff, &annotated_gff, &results)?;
    println!("  ✓  Annotated GFF3 → {}", annotated_gff.display());

    // ── Validate output GFF3 ────────────────────────────────────────────────
    match crate::utils::validation::validate_gff3(&annotated_gff) {
        Ok(validation) => {
            if !validation.is_valid() {
                println!(
                    "  ⚠  Output GFF3 has {} validation issue(s) — see annotation_report.txt",
                    validation.errors.len()
                );
            }
        }
        Err(_) => {} // validation failure is non-fatal
    }

    // ── Validate output protein FASTA ────────────────────────────────────────
    //
    // When any protein carries an internal stop codon, emit both a printed
    // summary and a per-gene TSV so users can act on the finding without
    // re-running validation. Fixes F5: the previous message only reported
    // the count, leaving users blind to *which* proteins were affected.
    match crate::utils::validation::validate_protein_fasta(&proteins_fa) {
        Ok(pval) => {
            if !pval.internal_stops.is_empty() {
                let list_path = config.out_dir.join("internal_stop_codon_genes.tsv");
                if let Ok(mut f) = std::fs::File::create(&list_path) {
                    use std::io::Write;
                    let _ = writeln!(f, "gene_id\tinternal_stop_positions");
                    for (gid, positions) in &pval.internal_stops {
                        let pos_str: Vec<String> =
                            positions.iter().map(|p| p.to_string()).collect();
                        let _ = writeln!(f, "{}\t{}", gid, pos_str.join(","));
                    }
                    println!(
                        "  ⚠  {} protein(s) have internal stop codons → {}",
                        pval.internal_stops.len(),
                        list_path.display()
                    );
                } else {
                    println!(
                        "  ⚠  {} protein(s) have internal stop codons \
                         (could not write list to {})",
                        pval.internal_stops.len(),
                        list_path.display()
                    );
                }
            }
        }
        Err(_) => {}
    }

    // ── Write functional annotation TSV ──────────────────────────────────────
    let annot_tsv = config.out_dir.join("annotations.tsv");
    write_annotation_tsv(&annot_tsv, &results)?;
    println!("  ✓  Annotation table → {}", annot_tsv.display());

    // ── Write summary report ──────────────────────────────────────────────────
    let report_path = config.out_dir.join("annotation_report.txt");
    write_report(&report_path, config, &results)?;
    println!("  ✓  Report → {}", report_path.display());
    println!(
        "  Annotated: {:.1}% of genes have a functional description",
        results.annotated_fraction() * 100.0
    );

    Ok(results)
}

// ─────────────────────────────────────────────────────────────────────────────
// Protein extraction from GFF3 + FASTA
// ─────────────────────────────────────────────────────────────────────────────

/// Extract translated CDS sequences for all genes → proteins.fa
/// Returns list of locus_tag / gene IDs extracted.
fn extract_proteins(
    gff_path: &Path,
    fasta_path: &Path,
    out_fa: &Path,
    genetic_code_table: u8,
) -> Result<Vec<String>> {
    use crate::parser::fasta::read_fasta_index;
    use crate::parser::gff::GFFReader;

    let fasta_index = read_fasta_index(fasta_path)?;
    let records: Vec<_> = GFFReader::from_path(gff_path)?
        .filter_map(|r| r.ok())
        .collect();

    // Collect CDS records grouped by gene (via Parent chain)
    let mut gene_cds: HashMap<String, Vec<_>> = HashMap::new();
    for rec in &records {
        if rec.feature_type != "CDS" {
            continue;
        }
        // Walk up Parent chain to find gene locus_tag
        let parent = match rec.parent() {
            Some(p) => p.clone(),
            None => continue,
        };
        gene_cds.entry(parent).or_default().push(rec.clone());
    }

    let mut out = std::fs::File::create(out_fa).map_err(MycoNoteError::Io)?;
    let mut extracted: Vec<String> = Vec::new();

    for (mrna_id, mut cds_list) in gene_cds {
        // Find the gene parent of this mRNA
        let gene_id = records
            .iter()
            .find(|r| r.id().map(|id| id == &mrna_id).unwrap_or(false))
            .and_then(|r| r.parent())
            .cloned()
            .unwrap_or_else(|| mrna_id.clone());

        let locus_tag = records
            .iter()
            .find(|r| r.id().map(|id| id == &gene_id).unwrap_or(false))
            .and_then(|r| r.attributes.get("locus_tag"))
            .cloned()
            .unwrap_or_else(|| gene_id.clone());

        // Sort CDS by start position
        cds_list.sort_by_key(|r| r.start);

        // Get genome sequence for this seqid
        let seq_rec = match fasta_index.get(&cds_list[0].seqid) {
            Some(s) => s,
            None => continue,
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

        // Translate to protein (genetic code-aware)
        let gc = genetic_code::GeneticCode::from_table_number(genetic_code_table)
            .unwrap_or(genetic_code::GeneticCode::STANDARD);
        let protein = gc.translate(&cds_seq);
        if protein.len() < 10 {
            continue;
        } // skip very short ORFs

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
/// Retained for backwards compatibility; prefer genetic_code::GeneticCode::translate()
#[allow(dead_code)]
fn translate_dna(dna: &str) -> String {
    let bytes = dna.as_bytes();
    let mut prot = String::with_capacity(bytes.len() / 3);
    let mut i = 0;
    while i + 2 < bytes.len() {
        let codon = [
            bytes[i].to_ascii_uppercase(),
            bytes[i + 1].to_ascii_uppercase(),
            bytes[i + 2].to_ascii_uppercase(),
        ];
        prot.push(codon_to_aa(&codon));
        i += 3;
    }
    // Remove trailing stop codon if present
    if prot.ends_with('*') {
        prot.pop();
    }
    prot
}

#[allow(dead_code)]
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
        if p.exists() {
            return Some(p.clone());
        }
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
        if p.exists() {
            return Some(p.clone());
        }
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
                gene.product = Some(hit.description.clone());
                gene.uniprot_acc = Some(hit.target_id.clone());
                gene.identity = Some(hit.identity);
                gene.evalue = Some(hit.evalue);
            }
        }
    }
}

fn merge_go_terms(results: &mut AnnotationResults, go_map: HashMap<String, Vec<String>>) {
    for gene in results.genes.values_mut() {
        if let Some(ref acc) = gene.uniprot_acc.clone() {
            if let Some(terms) = go_map.get(acc) {
                // Extend rather than replace — preserve GO terms already added
                // by InterProScan (merge_interproscan_hits runs first).
                for t in terms {
                    if !gene.go_terms.contains(t) {
                        gene.go_terms.push(t.clone());
                    }
                }
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
                if gene.product.is_none()
                    && !hit.db_desc.is_empty()
                    && hit.db_desc != "Uncharacterised protein"
                {
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
    input_gff: &Path,
    output_gff: &Path,
    results: &AnnotationResults,
) -> Result<()> {
    use crate::parser::gff::{include_stop_codon_in_cds, GFFReader};
    use std::collections::BTreeMap;

    // Buffer the whole GFF so we can fold stop codons into their terminal CDS
    // per transcript before writing. `annotate` is a file-in / file-out stage,
    // not a streaming pipe, so holding the records is fine. The fold is
    // idempotent: when the input already carries stop-inclusive CDS (the usual
    // case, since `predict` now emits it that way) nothing changes here.
    let mut records: Vec<crate::parser::gff::GFFRecord> = Vec::new();
    for rec_res in GFFReader::from_path(input_gff)? {
        if let Ok(r) = rec_res {
            records.push(r);
        }
    }

    // Group the coding/exon/stop rows by transcript (their Parent) and run the
    // fold on each group, then write the corrected coordinates back in place.
    let mut by_parent: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, r) in records.iter().enumerate() {
        if matches!(r.feature_type.as_str(), "CDS" | "exon" | "stop_codon") {
            if let Some(parent) = r.parent() {
                by_parent.entry(parent.clone()).or_default().push(i);
            }
        }
    }
    for idxs in by_parent.values() {
        let mut group: Vec<crate::parser::gff::GFFRecord> =
            idxs.iter().map(|&i| records[i].clone()).collect();
        include_stop_codon_in_cds(&mut group);
        for (k, &i) in idxs.iter().enumerate() {
            records[i] = group[k].clone();
        }
    }

    let mut out = std::fs::File::create(output_gff).map_err(MycoNoteError::Io)?;
    writeln!(out, "##gff-version 3").map_err(MycoNoteError::Io)?;

    for mut rec in records {
        if rec.feature_type == "gene" {
            let locus_tag = rec
                .attributes
                .get("locus_tag")
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
                    rec.attributes
                        .insert("Ontology_term".into(), ann.go_terms.join(","));
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
        "locus_tag\tproduct\tuniprot_acc\tidentity\tevalue\tgo_terms\tpfam_domains\tipr_accessions\tipr_databases\tcog_category\teggnog_og\tcazyme_families\tsignal_peptide\ttm_helices\tbgc_cluster\tbgc_type"
    ).map_err(MycoNoteError::Io)?;

    let mut genes: Vec<_> = results.genes.values().collect();
    genes.sort_by(|a, b| a.locus_tag.cmp(&b.locus_tag));

    for g in genes {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g.locus_tag,
            g.product.as_deref().unwrap_or("hypothetical protein"),
            g.uniprot_acc.as_deref().unwrap_or(""),
            g.identity.map(|v| format!("{:.1}", v)).unwrap_or_default(),
            g.evalue.map(|v| format!("{:.2e}", v)).unwrap_or_default(),
            g.go_terms.join("|"),
            g.pfam_domains.join("|"),
            g.ipr_accessions.join("|"),
            g.ipr_databases.join("|"),
            g.cog_category.as_deref().unwrap_or(""),
            g.eggnog_og.as_deref().unwrap_or(""),
            g.cazyme_families.join("|"),
            g.signal_peptide.as_deref().unwrap_or(""),
            g.tm_helices.map(|n| n.to_string()).unwrap_or_default(),
            g.bgc_cluster.as_deref().unwrap_or(""),
            g.bgc_type.as_deref().unwrap_or(""),
        )
        .map_err(MycoNoteError::Io)?;
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
    let annotated = results
        .genes
        .values()
        .filter(|g| g.product.is_some())
        .count();
    let with_go = results
        .genes
        .values()
        .filter(|g| !g.go_terms.is_empty())
        .count();
    let with_pfam = results
        .genes
        .values()
        .filter(|g| !g.pfam_domains.is_empty())
        .count();

    writeln!(f, "Gene totals").map_err(MycoNoteError::Io)?;
    writeln!(f, "  Total genes             : {}", total).map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "  With product description: {} ({:.1}%)",
        annotated,
        annotated as f64 / total.max(1) as f64 * 100.0
    )
    .map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "  With GO terms           : {} ({:.1}%)",
        with_go,
        with_go as f64 / total.max(1) as f64 * 100.0
    )
    .map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "  With Pfam domains       : {} ({:.1}%)",
        with_pfam,
        with_pfam as f64 / total.max(1) as f64 * 100.0
    )
    .map_err(MycoNoteError::Io)?;

    if let Some(ref b) = results.busco_summary {
        writeln!(f, "").map_err(MycoNoteError::Io)?;
        writeln!(f, "BUSCO ({} lineage)", b.lineage).map_err(MycoNoteError::Io)?;
        writeln!(
            f,
            "  Complete    : {} ({:.1}%)",
            b.complete(),
            b.percent_complete()
        )
        .map_err(MycoNoteError::Io)?;
        writeln!(f, "    Single    : {}", b.single_copy).map_err(MycoNoteError::Io)?;
        writeln!(f, "    Duplicated: {}", b.duplicated).map_err(MycoNoteError::Io)?;
        writeln!(f, "  Fragmented  : {}", b.fragmented).map_err(MycoNoteError::Io)?;
        writeln!(f, "  Missing     : {}", b.missing).map_err(MycoNoteError::Io)?;
        writeln!(f, "  Total BUSCOs: {}", b.total).map_err(MycoNoteError::Io)?;
    }

    Ok(())
}
