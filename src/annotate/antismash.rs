/// antiSMASH integration — secondary metabolite cluster annotation
///
/// antiSMASH (antibiotics & secondary metabolite analysis shell) identifies
/// biosynthetic gene clusters (BGCs) in fungal genomes.  This is essential
/// for:
///   - Identifying polyketide synthases (PKS), non-ribosomal peptide
///     synthetases (NRPS), terpene clusters, etc.
///   - Comparative genomics of secondary metabolite potential
///   - Connecting gene clusters to known metabolites via MIBiG database
///
/// Two integration modes:
///   1. **Run antiSMASH locally** — requires antismash ≥ 6 in PATH,
///      takes a GenBank file (GFF3 + FASTA) as input
///   2. **Parse existing antiSMASH results** — parse the JSON output from
///      a previous antiSMASH run (antismash output dir → *.json)
///
/// BGC annotations are merged back into the GFF3 as:
///   bgc_type    — cluster type (e.g. "T1PKS", "NRPS", "terpene")
///   bgc_id      — antiSMASH cluster ID
///   bgc_mibig   — closest MIBiG hit (if any)
///   bgc_product — predicted product (if available)
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Data types
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct BgcCluster {
    /// antiSMASH cluster ID (e.g. "ctg1_c1")
    pub cluster_id: String,
    /// Sequence / contig ID
    pub seq_id: String,
    /// Cluster start (1-based)
    pub start: u64,
    /// Cluster end (1-based, inclusive)
    pub end: u64,
    /// BGC type(s) — comma-separated if hybrid (e.g. "T1PKS-NRPS")
    pub bgc_type: String,
    /// Closest known cluster from MIBiG (if similarity > threshold)
    pub mibig_hit: Option<String>,
    /// Predicted product (if known)
    pub product: Option<String>,
    /// Genes within this cluster
    pub gene_ids: Vec<String>,
    /// Similarity to closest known cluster (0–100)
    pub similarity: f32,
}

// ─────────────────────────────────────────────────────────────────────────────
// Run antiSMASH
// ─────────────────────────────────────────────────────────────────────────────

pub fn antismash_available() -> bool {
    Command::new("which")
        .arg("antismash")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Run antiSMASH on a GenBank file and return the output directory.
pub fn run_antismash(
    genbank_file: &Path,
    out_dir: &Path,
    taxon: &str, // "fungi" | "bacteria" | "plants"
    threads: usize,
    extra_args: &[&str],
) -> Result<PathBuf> {
    if !antismash_available() {
        return Err(MycoNoteError::ExternalTool(
            "antiSMASH not found. Install: conda install -c bioconda antismash".to_string(),
        ));
    }

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let mut cmd = Command::new("antismash");
    cmd.arg(genbank_file)
        .arg("--taxon")
        .arg(taxon)
        .arg("--output-dir")
        .arg(out_dir)
        .arg("--cpus")
        .arg(threads.to_string())
        .arg("--genefinding-tool")
        .arg("none") // we have our own gene calls
        .arg("--output-basename")
        .arg("antismash");

    // Common useful flags
    cmd.arg("--cb-general") // ClusterBlast vs GenBank
        .arg("--cb-knownclusters") // ClusterBlast vs MIBiG
        .arg("--asf") // active site finder
        .arg("--pfam2go"); // Pfam → GO (useful for integration)

    for arg in extra_args {
        cmd.arg(arg);
    }

    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("antismash: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "antiSMASH exited with non-zero status".to_string(),
        ));
    }

    Ok(out_dir.to_path_buf())
}

// ─────────────────────────────────────────────────────────────────────────────
// Parse antiSMASH GFF output
// ─────────────────────────────────────────────────────────────────────────────

/// Parse antiSMASH cluster annotations from the GFF3 it writes to the output dir.
/// antiSMASH writes one GFF3 per contig; we scan them all.
pub fn parse_antismash_gff(out_dir: &Path) -> Result<Vec<BgcCluster>> {
    let mut clusters: Vec<BgcCluster> = Vec::new();

    let entries = std::fs::read_dir(out_dir).map_err(MycoNoteError::Io)?;
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) != Some("gff") {
            continue;
        }

        let file = match std::fs::File::open(&path) {
            Ok(f) => f,
            Err(_) => continue,
        };
        let reader = BufReader::new(file);

        for line_res in reader.lines() {
            let line = line_res.map_err(MycoNoteError::Io)?;
            let trimmed = line.trim();
            if trimmed.starts_with('#') || trimmed.is_empty() {
                continue;
            }

            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 9 {
                continue;
            }

            // antiSMASH uses "region" or "cluster" as feature type
            if fields[2] != "region" && fields[2] != "cluster" {
                continue;
            }

            let seq_id = fields[0].to_string();
            let start: u64 = fields[3].parse().unwrap_or(0);
            let end: u64 = fields[4].parse().unwrap_or(0);
            let attrs = parse_attrs(fields[8]);

            let bgc_type = attrs
                .get("product")
                .or(attrs.get("rules"))
                .cloned()
                .unwrap_or_else(|| "unknown".to_string());
            let cluster_id = attrs
                .get("ID")
                .cloned()
                .unwrap_or_else(|| format!("{}_{}-{}", seq_id, start, end));

            clusters.push(BgcCluster {
                cluster_id,
                seq_id,
                start,
                end,
                bgc_type,
                mibig_hit: attrs.get("knownclusterblast").cloned(),
                product: attrs.get("product").cloned(),
                gene_ids: vec![],
                similarity: 0.0,
            });
        }
    }

    Ok(clusters)
}

fn parse_attrs(attr_str: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for part in attr_str.split(';') {
        if let Some((k, v)) = part.split_once('=') {
            map.insert(k.trim().to_lowercase(), v.trim().to_string());
        }
    }
    map
}

// ─────────────────────────────────────────────────────────────────────────────
// Map genes to BGC clusters
// ─────────────────────────────────────────────────────────────────────────────

/// Given a list of gene coordinates (id, seqid, start, end), assign each
/// gene that falls within a BGC cluster to that cluster's gene_ids list.
pub fn assign_genes_to_clusters(
    clusters: &mut Vec<BgcCluster>,
    genes: &[(String, String, u64, u64)], // (id, seqid, start, end)
) {
    for (gene_id, seqid, g_start, g_end) in genes {
        for cluster in clusters.iter_mut() {
            if &cluster.seq_id == seqid && *g_start >= cluster.start && *g_end <= cluster.end {
                cluster.gene_ids.push(gene_id.clone());
            }
        }
    }
}

/// Return a map of gene_id → BgcCluster for all genes inside a cluster.
pub fn gene_to_cluster_map(clusters: &[BgcCluster]) -> HashMap<String, &BgcCluster> {
    let mut map = HashMap::new();
    for cluster in clusters {
        for gene_id in &cluster.gene_ids {
            map.insert(gene_id.clone(), cluster);
        }
    }
    map
}

// ─────────────────────────────────────────────────────────────────────────────
// Write results
// ─────────────────────────────────────────────────────────────────────────────

pub fn write_bgc_table(clusters: &[BgcCluster], output: &Path) -> Result<usize> {
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(
        out,
        "cluster_id\tseq_id\tstart\tend\tbgc_type\tgenes\tmibig_hit\tproduct\tsimilarity"
    )
    .map_err(MycoNoteError::Io)?;

    for c in clusters {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            c.cluster_id,
            c.seq_id,
            c.start,
            c.end,
            c.bgc_type,
            c.gene_ids.len(),
            c.mibig_hit.as_deref().unwrap_or("-"),
            c.product.as_deref().unwrap_or("-"),
            c.similarity,
        )
        .map_err(MycoNoteError::Io)?;
    }

    Ok(clusters.len())
}

pub fn print_bgc_summary(clusters: &[BgcCluster]) {
    let mut type_counts: HashMap<&str, usize> = HashMap::new();
    for c in clusters {
        *type_counts.entry(c.bgc_type.as_str()).or_insert(0) += 1;
    }

    println!("  BGC clusters found: {}", clusters.len());
    let mut sorted: Vec<(&&str, &usize)> = type_counts.iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(a.1));
    for (typ, count) in &sorted {
        println!("    {:3}  {}", count, bgc_type_description(typ));
    }
}

fn bgc_type_description(t: &str) -> String {
    let desc = match t {
        "T1PKS" => "Type I polyketide synthase",
        "T2PKS" => "Type II polyketide synthase",
        "T3PKS" => "Type III polyketide synthase",
        "NRPS" => "Non-ribosomal peptide synthetase",
        "terpene" => "Terpene cluster",
        "indole" => "Indole alkaloid",
        "RiPP" => "Ribosomally synthesised & post-translationally modified peptide",
        "siderophore" => "Siderophore",
        "betalactone" => "Beta-lactone",
        _ => t,
    };
    format!("[{}] {}", t, desc)
}
