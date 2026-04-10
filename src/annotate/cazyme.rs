/// CAZyme annotation via dbCAN
///
/// Carbohydrate-Active enZYmes (CAZymes) are critical for fungi —
/// they degrade plant cell walls, synthesise fungal cell wall components,
/// and are key virulence determinants in plant pathogens.
///
/// This module wraps the dbCAN2 pipeline (run_dbcan.py) which uses
/// three complementary methods and takes a majority-vote:
///   1. DIAMOND blast against CAZy database (dbCAN-sub)
///   2. HMMER search against dbCAN HMM profiles
///   3. Hotpep peptide pattern search
///
/// If run_dbcan.py is not available, falls back to standalone DIAMOND
/// search against the included dbCAN.dmnd database bundled with myconote.
///
/// CAZyme families are reported as:
///   GH   — glycoside hydrolases
///   GT   — glycosyltransferases
///   PL   — polysaccharide lyases
///   CE   — carbohydrate esterases
///   AA   — auxiliary activities
///   CBM  — carbohydrate-binding modules
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Data types
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct CazymeHit {
    /// Gene / protein ID
    pub gene_id: String,
    /// CAZyme family (e.g. "GH18", "GT2")
    pub family: String,
    /// E-value of best hit
    pub evalue: f64,
    /// Coverage of the HMM profile (0.0–1.0) — only for HMMER hits
    pub coverage: f64,
    /// Detection method used
    pub method: CazymeMethod,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CazymeMethod {
    Diamond,
    HMMER,
    Hotpep,
    DbCan2, // majority vote from run_dbcan.py
}

impl std::fmt::Display for CazymeMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CazymeMethod::Diamond => write!(f, "DIAMOND"),
            CazymeMethod::HMMER => write!(f, "HMMER"),
            CazymeMethod::Hotpep => write!(f, "Hotpep"),
            CazymeMethod::DbCan2 => write!(f, "dbCAN2"),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Run dbCAN2 (run_dbcan.py) — preferred method
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_dbcan(
    protein_fasta: &Path,
    out_dir: &Path,
    db_dir: Option<&Path>,
    threads: usize,
) -> Result<PathBuf> {
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let mut cmd = Command::new("run_dbcan.py");
    cmd.arg(protein_fasta)
        .arg("protein")
        .arg("--out_dir")
        .arg(out_dir)
        .arg("--cpu")
        .arg(threads.to_string());

    if let Some(db) = db_dir {
        cmd.arg("--db_dir").arg(db);
    }

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("run_dbcan.py: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "run_dbcan.py exited with non-zero status".to_string(),
        ));
    }

    Ok(out_dir.join("overview.txt"))
}

pub fn dbcan_available() -> bool {
    Command::new("which")
        .arg("run_dbcan.py")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

// ─────────────────────────────────────────────────────────────────────────────
// Fallback: DIAMOND search against dbCAN.dmnd
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_diamond_cazyme(
    protein_fasta: &Path,
    diamond_db: &Path,
    out_dir: &Path,
    threads: usize,
    evalue_cutoff: f64,
) -> Result<PathBuf> {
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let out_tsv = out_dir.join("cazyme_diamond.tsv");

    let status = Command::new("diamond")
        .arg("blastp")
        .arg("--query")
        .arg(protein_fasta)
        .arg("--db")
        .arg(diamond_db)
        .arg("--out")
        .arg(&out_tsv)
        .arg("--outfmt")
        .arg("6")
        .arg("qseqid")
        .arg("sseqid")
        .arg("pident")
        .arg("length")
        .arg("evalue")
        .arg("bitscore")
        .arg("qcovhsp")
        .arg("--evalue")
        .arg(evalue_cutoff.to_string())
        .arg("--max-target-seqs")
        .arg("1")
        .arg("--threads")
        .arg(threads.to_string())
        .arg("--quiet")
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("diamond: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "diamond blastp (CAZyme) exited with non-zero status".to_string(),
        ));
    }

    Ok(out_tsv)
}

// ─────────────────────────────────────────────────────────────────────────────
// Parse dbCAN2 overview.txt
// ─────────────────────────────────────────────────────────────────────────────

/// Parse the `overview.txt` produced by run_dbcan.py.
/// Returns map of gene_id → Vec<CazymeHit> (one entry per family).
pub fn parse_dbcan_overview(path: &Path) -> Result<HashMap<String, Vec<CazymeHit>>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, Vec<CazymeHit>> = HashMap::new();

    // overview.txt columns: Gene ID | HMMER | Hotpep | DIAMOND | #ofTools
    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.starts_with("Gene ID") || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }

        let gene_id = fields[0].to_string();
        // Collect unique families from all three columns
        let mut families: Vec<(String, CazymeMethod)> = vec![];
        let method_fields = [
            (fields[1], CazymeMethod::HMMER),
            (fields[2], CazymeMethod::Hotpep),
            (fields[3], CazymeMethod::Diamond),
        ];
        for (field, method) in method_fields {
            if field != "-" {
                for fam in field.split('+') {
                    let fam = fam.split('(').next().unwrap_or(fam).trim().to_string();
                    if !fam.is_empty() {
                        families.push((fam, method.clone()));
                    }
                }
            }
        }

        for (family, method) in families {
            map.entry(gene_id.clone()).or_default().push(CazymeHit {
                gene_id: gene_id.clone(),
                family,
                evalue: 0.0,
                coverage: 0.0,
                method,
            });
        }
    }

    Ok(map)
}

/// Parse a DIAMOND tabular output from the fallback CAZyme search.
pub fn parse_diamond_cazyme(
    path: &Path,
    evalue_cutoff: f64,
) -> Result<HashMap<String, Vec<CazymeHit>>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, Vec<CazymeHit>> = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let gene_id = fields[0].to_string();
        let subject = fields[1]; // e.g. "GH18_Chitinase|CAZy:GH18"
        let evalue: f64 = fields[4].parse().unwrap_or(1.0);
        if evalue > evalue_cutoff {
            continue;
        }

        // Extract family name — take the first "|"-delimited token, strip after "_"
        let family = subject
            .split('|')
            .next()
            .unwrap_or(subject)
            .split('_')
            .next()
            .unwrap_or(subject)
            .to_string();

        map.entry(gene_id.clone()).or_default().push(CazymeHit {
            gene_id,
            family,
            evalue,
            coverage: 0.0,
            method: CazymeMethod::Diamond,
        });
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// Write CAZyme results
// ─────────────────────────────────────────────────────────────────────────────

pub fn write_cazyme_table(hits: &HashMap<String, Vec<CazymeHit>>, output: &Path) -> Result<usize> {
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(out, "gene_id\tcazyme_family\tevalue\tmethod").map_err(MycoNoteError::Io)?;

    let mut total = 0usize;
    let mut genes: Vec<&String> = hits.keys().collect();
    genes.sort();

    for gene_id in genes {
        for hit in &hits[gene_id] {
            writeln!(
                out,
                "{}\t{}\t{:.2e}\t{}",
                hit.gene_id, hit.family, hit.evalue, hit.method
            )
            .map_err(MycoNoteError::Io)?;
            total += 1;
        }
    }

    Ok(total)
}

/// Print a family-class breakdown (GH/GT/PL/CE/AA/CBM) to stdout.
pub fn print_cazyme_summary(hits: &HashMap<String, Vec<CazymeHit>>) {
    let mut class_counts: HashMap<&str, usize> = HashMap::new();
    for hit_list in hits.values() {
        for hit in hit_list {
            let cls = cazyme_class(&hit.family);
            *class_counts.entry(cls).or_insert(0) += 1;
        }
    }

    println!("  CAZyme family breakdown ({} genes):", hits.len());
    let mut sorted: Vec<(&&str, &usize)> = class_counts.iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(a.1));
    for (cls, count) in sorted {
        println!("    {:4}  {}", count, cazyme_class_name(cls));
    }
}

fn cazyme_class(family: &str) -> &'static str {
    if family.starts_with("GH") {
        "GH"
    } else if family.starts_with("GT") {
        "GT"
    } else if family.starts_with("PL") {
        "PL"
    } else if family.starts_with("CE") {
        "CE"
    } else if family.starts_with("AA") {
        "AA"
    } else if family.starts_with("CBM") {
        "CBM"
    } else {
        "Other"
    }
}

fn cazyme_class_name(cls: &str) -> &'static str {
    match cls {
        "GH" => "Glycoside Hydrolases (GH) — cell wall degradation",
        "GT" => "Glycosyltransferases (GT) — cell wall biosynthesis",
        "PL" => "Polysaccharide Lyases (PL) — pectin degradation",
        "CE" => "Carbohydrate Esterases (CE) — de-acetylation",
        "AA" => "Auxiliary Activities (AA) — redox enzymes (LPMO, laccase)",
        "CBM" => "Carbohydrate-Binding Modules (CBM) — substrate binding",
        _ => "Other CAZymes",
    }
}
