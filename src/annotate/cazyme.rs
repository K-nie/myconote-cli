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
// Native HMMER vs dbCAN HMM profiles (preferred — funannotate-comparable)
// ─────────────────────────────────────────────────────────────────────────────
//
// This is the same method funannotate uses: hmmsearch the predicted proteins
// against the dbCAN HMM database (dbCAN.hmm) and keep hits passing the canonical
// dbCAN cutoffs. It calls hmmsearch directly (reusing the same engine as the
// Pfam step) rather than the external run_dbcan pipeline, so it is self-contained
// and insensitive to run_dbcan's cross-version CLI / overview.txt changes.
//
// DIAMOND-blastp against the full CAZyDB sequence set (run_diamond_cazyme above)
// measures similarity to full-length multidomain proteins and massively
// over-calls CAZymes; the HMM-profile search tests family-domain membership and
// is the correct comparator for the benchmark.

/// Run hmmsearch of `protein_fasta` against the dbCAN HMM database and return
/// the path to the `--domtblout` table. A loose E-value pre-filter is applied
/// here; the strict dbCAN cutoffs are enforced in `parse_dbcan_domtblout`.
pub fn run_hmmer_cazyme(
    protein_fasta: &Path,
    dbcan_hmm: &Path,
    out_dir: &Path,
    threads: usize,
) -> Result<PathBuf> {
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    // Press the HMM database if the .h3i index is missing (dbCAN.hmm is usually
    // shipped pre-pressed by funannotate setup; this is a safety net).
    let h3i = dbcan_hmm.with_extension("hmm.h3i");
    if !h3i.exists() {
        let hmmpress = which::which("hmmpress").map_err(|_| {
            MycoNoteError::UnsupportedFormat("hmmpress not found. Install hmmer.".to_string())
        })?;
        let press_ok = Command::new(&hmmpress)
            .args(["-f", &dbcan_hmm.to_string_lossy()])
            .status()
            .map_err(MycoNoteError::Io)?;
        if !press_ok.success() {
            return Err(MycoNoteError::InvalidFormat(
                "hmmpress (dbCAN) failed.".to_string(),
            ));
        }
    }

    let hmmsearch = which::which("hmmsearch").map_err(|_| {
        MycoNoteError::UnsupportedFormat("hmmsearch not found. Install hmmer.".to_string())
    })?;

    let domtblout = out_dir.join("dbcan_hmmer.domtblout");
    let hmm_s = dbcan_hmm.to_string_lossy();
    let query_s = protein_fasta.to_string_lossy();
    let dom_s = domtblout.to_string_lossy();

    let status = Command::new(&hmmsearch)
        .args([
            "--domtblout",
            dom_s.as_ref(),
            "--cpu",
            &threads.to_string(),
            // Loose pre-filter; dbCAN E/coverage cutoffs applied in the parser.
            "-E",
            "1e-10",
            "--domE",
            "1e-10",
            "--noali",
            hmm_s.as_ref(),
            query_s.as_ref(),
        ])
        .stdout(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "hmmsearch (dbCAN) failed. Check that dbCAN.hmm is pressed (hmmpress).".to_string(),
        ));
    }

    Ok(domtblout)
}

/// Parse an hmmsearch `--domtblout` table from the dbCAN HMM search and apply the
/// canonical dbCAN cutoffs: HMM-profile coverage > 0.35, and domain i-E-value
/// < 1e-15 (tightened to < 1e-18 when the aligned HMM span is <= 80 residues).
/// Keeps the single best (lowest-E) hit per (gene, family). The CAZyme family is
/// the HMM profile name with any trailing ".hmm" removed (e.g. `GH5.hmm` → `GH5`).
pub fn parse_dbcan_domtblout(path: &Path) -> Result<HashMap<String, Vec<CazymeHit>>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);

    // Best hit per (gene_id, family), by lowest e-value.
    let mut best: HashMap<(String, String), CazymeHit> = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 22 {
            continue;
        }

        // hmmsearch --domtblout: col 0 = protein (target), col 3 = HMM (query),
        // col 5 = query length (HMM model length), col 12 = domain i-E-value,
        // col 15 = hmm coord from, col 16 = hmm coord to.
        let gene_id = f[0].to_string();
        let raw_name = f[3];
        let family = raw_name
            .strip_suffix(".hmm")
            .unwrap_or(raw_name)
            .to_string();
        let hmm_len: f64 = f[5].parse().unwrap_or(0.0);
        let evalue: f64 = f[12].parse().unwrap_or(f64::MAX);
        let hmm_from: f64 = f[15].parse().unwrap_or(0.0);
        let hmm_to: f64 = f[16].parse().unwrap_or(0.0);

        if hmm_len <= 0.0 {
            continue;
        }
        let aln_span = hmm_to - hmm_from; // residues aligned on the HMM model
        let coverage = aln_span / hmm_len;
        let e_cutoff = if aln_span > 80.0 { 1e-15 } else { 1e-18 };

        if evalue > e_cutoff || coverage <= 0.35 {
            continue;
        }

        let key = (gene_id.clone(), family.clone());
        let entry = best.entry(key).or_insert_with(|| CazymeHit {
            gene_id: gene_id.clone(),
            family: family.clone(),
            evalue,
            coverage,
            method: CazymeMethod::HMMER,
        });
        if evalue < entry.evalue {
            entry.evalue = evalue;
            entry.coverage = coverage;
        }
    }

    let mut map: HashMap<String, Vec<CazymeHit>> = HashMap::new();
    for ((gene, _family), hit) in best {
        map.entry(gene).or_default().push(hit);
    }
    Ok(map)
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    // Build one hmmsearch --domtblout data row (22 whitespace-delimited fields).
    // Only the columns the parser reads need to be meaningful:
    //   0 target(protein), 3 query(HMM name), 5 qlen(HMM len),
    //   12 dom i-Evalue, 15 hmm_from, 16 hmm_to.
    fn domtbl_row(
        protein: &str,
        hmm: &str,
        qlen: &str,
        ievalue: &str,
        hfrom: &str,
        hto: &str,
    ) -> String {
        let mut f = vec!["-"; 23];
        f[0] = protein;
        f[3] = hmm;
        f[5] = qlen;
        f[12] = ievalue;
        f[15] = hfrom;
        f[16] = hto;
        f[22] = "desc";
        f.join(" ")
    }

    #[test]
    fn dbcan_domtblout_applies_cutoffs_and_strips_hmm_suffix() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("dbcan_hmmer.domtblout");
        let mut fh = std::fs::File::create(&p).unwrap();
        writeln!(fh, "# comment line is skipped").unwrap();
        // PASS: cov = (300-10)/400 = 0.725 > 0.35, aln span 290 > 80 so E<1e-15; E=1e-30 ok.
        writeln!(
            fh,
            "{}",
            domtbl_row("GENE1", "GH5.hmm", "400", "1e-30", "10", "300")
        )
        .unwrap();
        // FAIL coverage: cov = (60-10)/400 = 0.125 <= 0.35.
        writeln!(
            fh,
            "{}",
            domtbl_row("GENE2", "GT2.hmm", "400", "1e-30", "10", "60")
        )
        .unwrap();
        // FAIL e-value: aln span 290 > 80 needs E<1e-15; 1e-12 fails.
        writeln!(
            fh,
            "{}",
            domtbl_row("GENE3", "AA9.hmm", "400", "1e-12", "10", "300")
        )
        .unwrap();
        // PASS short-domain rule: aln span 60 <= 80 needs E<1e-18; cov=(70-10)/70=0.857; E=1e-20 ok.
        writeln!(
            fh,
            "{}",
            domtbl_row("GENE4", "CBM1.hmm", "70", "1e-20", "10", "70")
        )
        .unwrap();

        let map = parse_dbcan_domtblout(&p).unwrap();
        // GENE1 and GENE4 pass; GENE2 and GENE3 are filtered out.
        assert_eq!(
            map.len(),
            2,
            "exactly two genes should pass the dbCAN cutoffs"
        );
        assert_eq!(map["GENE1"][0].family, "GH5", ".hmm suffix stripped");
        assert_eq!(map["GENE1"][0].method, CazymeMethod::HMMER);
        assert_eq!(map["GENE4"][0].family, "CBM1");
        assert!(!map.contains_key("GENE2"), "low-coverage hit rejected");
        assert!(!map.contains_key("GENE3"), "weak-e-value hit rejected");
    }

    #[test]
    fn dbcan_domtblout_keeps_best_evalue_per_gene_family() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("d.domtblout");
        let mut fh = std::fs::File::create(&p).unwrap();
        // Same gene+family twice; the lower e-value must win.
        writeln!(
            fh,
            "{}",
            domtbl_row("GENE1", "GH18.hmm", "400", "1e-20", "10", "300")
        )
        .unwrap();
        writeln!(
            fh,
            "{}",
            domtbl_row("GENE1", "GH18.hmm", "400", "1e-40", "10", "300")
        )
        .unwrap();
        let map = parse_dbcan_domtblout(&p).unwrap();
        assert_eq!(map["GENE1"].len(), 1, "one entry per (gene, family)");
        assert!(
            (map["GENE1"][0].evalue - 1e-40).abs() < 1e-50,
            "best (lowest) e-value kept"
        );
    }
}
