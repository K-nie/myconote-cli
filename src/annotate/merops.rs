/// MEROPS protease annotation
///
/// MEROPS (https://www.ebi.ac.uk/merops/) is the reference database for
/// proteolytic enzymes (proteases, peptidases) and their inhibitors.
///
/// Each protease family has a code: A (aspartic), C (cysteine), G (glutamic),
/// M (metallo), N (asparagine lyase), S (serine), T (threonine), U (unknown),
/// P (mixed catalytic type).  Clan codes group evolutionarily related families.
///
/// Annotation strategy (in priority order):
///   1. DIAMOND BLASTp vs merops_scan.lib (pre-formatted FASTA from MEROPS)
///   2. hmmscan vs MEROPS HMM profiles (merops.hmm)
///   3. Keyword rescue from product descriptions already in GeneAnnotation
///
/// The database file is downloaded by `myconote setup --download-dbs merops`.
/// Expected path: ~/.myconote/dbs/merops/merops_scan.lib
///                ~/.myconote/dbs/merops/merops.dmnd   (pre-built DIAMOND db)
///                ~/.myconote/dbs/merops/merops.hmm    (optional)

use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Data structures
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct MeropsHit {
    /// Gene / protein query ID
    pub query_id:     String,
    /// MEROPS identifier, e.g. "S08.001" (family S08, peptidase 001)
    pub merops_id:    String,
    /// Protease family code (A, C, G, M, N, S, T, U, P, I for inhibitors)
    pub family:       String,
    /// Clan code, e.g. "PA", "CA", "MA" (empty if not available)
    pub clan:         String,
    /// Substrate specificity / name, e.g. "subtilisin A"
    pub name:         String,
    /// Organism of the matched MEROPS sequence
    pub organism:     String,
    /// DIAMOND / BLAST e-value
    pub evalue:       f64,
    /// Percent identity (0–100)
    pub identity:     f64,
    /// Whether this is a protease inhibitor (true) rather than a peptidase (false)
    pub is_inhibitor: bool,
    /// Method: "diamond" or "hmmscan"
    pub method:       String,
}

// ─────────────────────────────────────────────────────────────────────────────
// Availability checks
// ─────────────────────────────────────────────────────────────────────────────

pub fn diamond_available() -> bool {
    Command::new("diamond").arg("version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

pub fn hmmscan_available() -> bool {
    Command::new("hmmscan").arg("-h")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

// ─────────────────────────────────────────────────────────────────────────────
// Find MEROPS database
// ─────────────────────────────────────────────────────────────────────────────

pub fn find_merops_db(db_dir: &Path) -> Option<PathBuf> {
    let candidates = [
        db_dir.join("merops").join("merops.dmnd"),
        db_dir.join("merops").join("merops_scan.dmnd"),
        db_dir.join("merops.dmnd"),
    ];
    for c in &candidates {
        if c.exists() { return Some(c.clone()); }
    }
    None
}

pub fn find_merops_hmm(db_dir: &Path) -> Option<PathBuf> {
    let candidates = [
        db_dir.join("merops").join("merops.hmm"),
        db_dir.join("merops.hmm"),
    ];
    for c in &candidates {
        if c.exists() { return Some(c.clone()); }
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Run MEROPS annotation
// ─────────────────────────────────────────────────────────────────────────────

/// Main entry point. Runs DIAMOND (preferred) or hmmscan, returns per-gene hits.
pub fn run_merops(
    proteins_fa: &Path,
    db_dir:      &Path,
    out_dir:     &Path,
    threads:     usize,
    evalue:      f64,
) -> Result<HashMap<String, Vec<MeropsHit>>> {
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    if let Some(dmnd) = find_merops_db(db_dir) {
        println!("  MEROPS: using DIAMOND database {}", dmnd.display());
        let tsv = out_dir.join("merops_diamond.tsv");
        run_diamond_merops(proteins_fa, &dmnd, &tsv, threads, evalue)?;
        return parse_diamond_merops(&tsv);
    }

    if let Some(hmm) = find_merops_hmm(db_dir) {
        println!("  MEROPS: using hmmscan profiles {}", hmm.display());
        let tbl = out_dir.join("merops_hmm.tbl");
        run_hmmscan_merops(proteins_fa, &hmm, &tbl, threads, evalue)?;
        return parse_hmmscan_merops(&tbl);
    }

    eprintln!(
        "  ⚠  MEROPS database not found in {}.\n  \
         Run: myconote setup --download-dbs merops", db_dir.display()
    );
    Ok(HashMap::new())
}

// ─────────────────────────────────────────────────────────────────────────────
// DIAMOND search
// ─────────────────────────────────────────────────────────────────────────────

fn run_diamond_merops(
    query:   &Path,
    db:      &Path,
    out_tsv: &Path,
    threads: usize,
    evalue:  f64,
) -> Result<()> {
    if !diamond_available() {
        return Err(MycoNoteError::ExternalTool(
            "diamond not found. Install: conda install -c bioconda diamond".to_string()
        ));
    }

    // Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
    let status = Command::new("diamond")
        .args(["blastp",
               "--query", query.to_str().unwrap_or(""),
               "--db",    db.to_str().unwrap_or(""),
               "--out",   out_tsv.to_str().unwrap_or(""),
               "--outfmt", "6",
               "qseqid", "sseqid", "pident", "evalue", "stitle",
               "--evalue", &format!("{}", evalue),
               "--max-target-seqs", "1",
               "--threads", &threads.to_string(),
               "--sensitive",
        ])
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("diamond blastp: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool("DIAMOND MEROPS search failed".to_string()));
    }

    Ok(())
}

fn parse_diamond_merops(tsv: &Path) -> Result<HashMap<String, Vec<MeropsHit>>> {
    let file   = std::fs::File::open(tsv).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, Vec<MeropsHit>> = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') { continue; }

        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 5 { continue; }

        let query_id  = cols[0].to_string();
        let target_id = cols[1].to_string();
        let identity: f64 = cols[2].parse().unwrap_or(0.0);
        let evalue:   f64 = cols[3].parse().unwrap_or(1.0);
        // stitle: e.g. "MER000001 - Subtilisin A (Bacillus amyloliquefaciens) - S08.001"
        let stitle = cols[4..].join("\t");

        let (family, clan, name, organism, merops_id, is_inhibitor) =
            parse_merops_title(&target_id, &stitle);

        let hit = MeropsHit {
            query_id:    query_id.clone(),
            merops_id,
            family,
            clan,
            name,
            organism,
            evalue,
            identity,
            is_inhibitor,
            method: "diamond".to_string(),
        };

        map.entry(query_id).or_default().push(hit);
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// hmmscan search (fallback)
// ─────────────────────────────────────────────────────────────────────────────

fn run_hmmscan_merops(
    query:   &Path,
    hmm_db:  &Path,
    out_tbl: &Path,
    threads: usize,
    evalue:  f64,
) -> Result<()> {
    if !hmmscan_available() {
        return Err(MycoNoteError::ExternalTool(
            "hmmscan not found. Install: conda install -c bioconda hmmer".to_string()
        ));
    }

    let status = Command::new("hmmscan")
        .args([
            "--tblout", out_tbl.to_str().unwrap_or(""),
            "-E",       &format!("{}", evalue),
            "--cpu",    &threads.to_string(),
            "--noali",
            hmm_db.to_str().unwrap_or(""),
            query.to_str().unwrap_or(""),
        ])
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("hmmscan: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool("hmmscan MEROPS search failed".to_string()));
    }

    Ok(())
}

fn parse_hmmscan_merops(tbl: &Path) -> Result<HashMap<String, Vec<MeropsHit>>> {
    let file   = std::fs::File::open(tbl).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, Vec<MeropsHit>> = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') { continue; }

        // hmmscan --tblout columns (space-separated, max 18 cols + description):
        // target_name acc query_name acc sequence_evalue sequence_score sequence_bias
        // best_domain_evalue best_domain_score best_domain_bias exp reg clu ov env dom rep inc description
        let cols: Vec<&str> = t.splitn(19, char::is_whitespace)
            .filter(|s| !s.is_empty())
            .collect();
        if cols.len() < 5 { continue; }

        let target_name = cols[0];  // e.g. MER000001 or family code
        let query_id    = cols[2].to_string();
        let evalue:  f64 = cols[4].parse().unwrap_or(1.0);
        let desc = if cols.len() > 18 { cols[18..].join(" ") } else { String::new() };

        let (family, clan, name, organism, merops_id, is_inhibitor) =
            parse_merops_title(target_name, &desc);

        let hit = MeropsHit {
            query_id:    query_id.clone(),
            merops_id,
            family,
            clan,
            name,
            organism,
            evalue,
            identity:    0.0,   // hmmscan doesn't report identity
            is_inhibitor,
            method: "hmmscan".to_string(),
        };

        map.entry(query_id).or_default().push(hit);
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// Title parsing
// ─────────────────────────────────────────────────────────────────────────────

/// Parse MEROPS family/clan/name/organism from a MEROPS sequence identifier
/// and title string.
///
/// MEROPS FASTA headers are structured as:
///   >MER000001 - Subtilisin A (Bacillus amyloliquefaciens) - S08.001
///   or just the MEROPS ID like "S08.001" in the HMM name field.
fn parse_merops_title(
    id_str: &str,
    title:  &str,
) -> (String, String, String, String, String, bool) {
    // Try to extract MEROPS ID (e.g. S08.001) from the id or title
    let merops_id = extract_merops_id(id_str)
        .or_else(|| extract_merops_id(title))
        .unwrap_or_else(|| id_str.to_string());

    // Family = first letter of MEROPS ID
    let family = merops_id.chars().next()
        .map(|c| c.to_string())
        .unwrap_or_default();

    // Is this an inhibitor? (I family in MEROPS)
    let is_inhibitor = family == "I";

    // Clan from MEROPS ID prefix (e.g. "S08" → clan usually in a lookup)
    // Simple heuristic: S → serine, A → aspartic, C → cysteine, M → metallo
    let clan = family_to_clan(&family);

    // Name: try to parse from title "- Subtilisin A (...) - S08.001"
    let name = extract_merops_name(title)
        .unwrap_or_else(|| title.split(" - ").next().unwrap_or("").trim().to_string());

    // Organism: in parentheses in MEROPS titles
    let organism = extract_organism(title).unwrap_or_default();

    (family, clan, name, organism, merops_id, is_inhibitor)
}

fn extract_merops_id(s: &str) -> Option<String> {
    // Match pattern like S08.001 or A01.005
    for word in s.split_whitespace() {
        let word = word.trim_matches(|c: char| !c.is_alphanumeric() && c != '.');
        if word.len() >= 3 {
            let chars: Vec<char> = word.chars().collect();
            if chars[0].is_ascii_uppercase()
                && chars[1..].iter().take(2).all(|c| c.is_ascii_digit())
                && word.contains('.')
            {
                return Some(word.to_string());
            }
        }
    }
    None
}

fn extract_merops_name(title: &str) -> Option<String> {
    // Title form: "MER000001 - Subtilisin A (Bacillus...) - S08.001"
    let parts: Vec<&str> = title.split(" - ").collect();
    if parts.len() >= 2 {
        // Second part is typically "Name (Organism)"
        let part = parts[1].trim();
        // Strip organism in parentheses
        let name = if let Some(paren) = part.find('(') {
            part[..paren].trim().to_string()
        } else {
            part.to_string()
        };
        if !name.is_empty() { return Some(name); }
    }
    None
}

fn extract_organism(title: &str) -> Option<String> {
    if let Some(start) = title.find('(') {
        if let Some(end) = title.find(')') {
            if end > start {
                return Some(title[start+1..end].to_string());
            }
        }
    }
    None
}

fn family_to_clan(family: &str) -> String {
    match family {
        "S" => "PA".to_string(),   // serine peptidases
        "C" => "CA".to_string(),   // cysteine peptidases
        "A" => "AA".to_string(),   // aspartic peptidases
        "M" => "MA".to_string(),   // metallo peptidases
        "T" => "PB".to_string(),   // threonine peptidases
        "G" => "GC".to_string(),   // glutamic peptidases
        "N" => "PD".to_string(),   // asparagine lyases
        "I" => "IA".to_string(),   // inhibitors
        _   => String::new(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Output writers
// ─────────────────────────────────────────────────────────────────────────────

/// Write per-gene MEROPS hits to TSV.
pub fn write_merops_table(
    hits_map: &HashMap<String, Vec<MeropsHit>>,
    out_path: &Path,
) -> Result<usize> {
    let mut f = std::fs::File::create(out_path).map_err(MycoNoteError::Io)?;
    writeln!(f,
        "gene_id\tmerops_id\tfamily\tclan\tname\torganism\tevalue\tidentity\tis_inhibitor\tmethod"
    ).map_err(MycoNoteError::Io)?;

    let mut count = 0usize;
    let mut genes: Vec<&String> = hits_map.keys().collect();
    genes.sort();

    for gene_id in genes {
        let hits = &hits_map[gene_id];
        // Write only the best hit (lowest evalue)
        if let Some(best) = hits.iter().min_by(|a, b|
            a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal)
        ) {
            writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}\t{}\t{}",
                best.query_id,
                best.merops_id,
                best.family,
                best.clan,
                best.name,
                best.organism,
                best.evalue,
                best.identity,
                if best.is_inhibitor { "yes" } else { "no" },
                best.method,
            ).map_err(MycoNoteError::Io)?;
            count += 1;
        }
    }

    Ok(count)
}

/// Print a summary table of MEROPS families found.
pub fn print_merops_summary(hits_map: &HashMap<String, Vec<MeropsHit>>) {
    if hits_map.is_empty() {
        println!("  No MEROPS proteases found.");
        return;
    }

    // Tally by family
    let mut family_counts: HashMap<String, (u32, u32)> = HashMap::new(); // (peptidases, inhibitors)
    for hits in hits_map.values() {
        if let Some(best) = hits.iter().min_by(|a, b|
            a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal)
        ) {
            let entry = family_counts.entry(best.family.clone()).or_insert((0, 0));
            if best.is_inhibitor { entry.1 += 1; } else { entry.0 += 1; }
        }
    }

    let mut families: Vec<(&String, &(u32, u32))> = family_counts.iter().collect();
    families.sort_by(|a, b| {
        let total_b = b.1.0 + b.1.1;
        let total_a = a.1.0 + a.1.1;
        total_b.cmp(&total_a)
    });

    let total_peptidases: u32 = families.iter().map(|(_, c)| c.0).sum();
    let total_inhibitors: u32 = families.iter().map(|(_, c)| c.1).sum();

    println!("  MEROPS summary: {} protease-related genes ({} peptidases, {} inhibitors)",
        hits_map.len(), total_peptidases, total_inhibitors);
    println!("  {:>8}  {:>12}  {:>12}  {:>20}",
        "Family", "Peptidases", "Inhibitors", "Catalytic type");

    for (fam, (pep, inh)) in &families {
        let cat = family_catalytic_type(fam);
        println!("  {:>8}  {:>12}  {:>12}  {:>20}", fam, pep, inh, cat);
    }
}

fn family_catalytic_type(family: &str) -> &'static str {
    match family {
        "S" => "Serine",
        "C" => "Cysteine",
        "A" => "Aspartic",
        "M" => "Metallo",
        "T" => "Threonine",
        "G" => "Glutamic",
        "N" => "Asparagine lyase",
        "P" => "Mixed",
        "U" => "Unknown",
        "I" => "Inhibitor",
        _   => "Unknown",
    }
}
