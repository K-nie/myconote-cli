/// GlimmerHMM gene predictor wrapper
///
/// GlimmerHMM uses a generalised hidden Markov model trained on
/// known gene structures.  It is particularly good at predicting
/// genes in AT-rich fungal genomes and performs well alongside Augustus.
///
/// Pre-trained models are distributed with GlimmerHMM for several
/// organisms.  Myconote ships kingdom-aware defaults:
///   fungi    → saccharomyces_cerevisiae_S288C (closest general fungal model)
///   plant    → arabidopsis
///   animal   → human
///
/// GlimmerHMM outputs GFF-like format which we convert to GFF3.
use crate::utils::error::{MycoNoteError, Result};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Availability + training model detection
// ─────────────────────────────────────────────────────────────────────────────

pub fn glimmerhmm_available() -> bool {
    Command::new("which")
        .arg("glimmerhmm")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Find the GlimmerHMM training directory for a given species.
/// Searches common install locations.
pub fn find_training_dir(species: &str) -> Option<PathBuf> {
    let candidates = [
        format!("/usr/share/glimmerhmm/trained_dir/{}", species),
        format!("/opt/conda/share/glimmerhmm/trained_dir/{}", species),
        format!("/usr/local/share/glimmerhmm/trained_dir/{}", species),
        format!("{}/.myconote/glimmerhmm/{}", dirs_home_str(), species),
    ];

    for c in &candidates {
        let p = PathBuf::from(c);
        if p.exists() {
            return Some(p);
        }
    }

    // Try using `glimmerhmm` to find its own data dir
    if let Ok(out) = Command::new("glimmerhmm").arg("--help").output() {
        let help = String::from_utf8_lossy(&out.stderr);
        for line in help.lines() {
            if line.contains("trained_dir") {
                if let Some(dir) = line.split_whitespace().last() {
                    let p = PathBuf::from(dir).join(species);
                    if p.exists() {
                        return Some(p);
                    }
                }
            }
        }
    }

    None
}

fn dirs_home_str() -> String {
    std::env::var("HOME").unwrap_or_else(|_| "/tmp".to_string())
}

/// Map kingdom to a GlimmerHMM training species name
pub fn default_training_species(kingdom: &str) -> &'static str {
    match kingdom.to_lowercase().as_str() {
        "fungi" | "ascomycota" | "basidiomycota" => "saccharomyces_cerevisiae_S288C",
        "plant" | "plants" => "arabidopsis",
        "animal" | "animals" | "mammals" => "human",
        "insect" | "insects" => "drosophila",
        _ => "saccharomyces_cerevisiae_S288C",
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Run GlimmerHMM
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_glimmerhmm(
    genome_fasta: &Path,
    training_dir: &Path,
    output_gff3: &Path,
    threads: usize,
) -> Result<usize> {
    if !glimmerhmm_available() {
        return Err(MycoNoteError::ExternalTool(
            "glimmerhmm not found. Install: conda install -c bioconda glimmerhmm".to_string(),
        ));
    }

    let raw_out = output_gff3.with_extension("glimmer_raw.txt");

    // GlimmerHMM is single-threaded; for multi-thread we split by contig
    if threads > 1 {
        run_glimmerhmm_parallel(genome_fasta, training_dir, &raw_out, threads)?;
    } else {
        let out_f = std::fs::File::create(&raw_out).map_err(MycoNoteError::Io)?;
        let status = Command::new("glimmerhmm")
            .arg(genome_fasta)
            .arg(training_dir)
            .arg("-f") // output in GFF format
            .stdout(std::process::Stdio::from(out_f))
            .status()
            .map_err(|e| MycoNoteError::ExternalTool(format!("glimmerhmm: {}", e)))?;

        if !status.success() {
            return Err(MycoNoteError::ExternalTool(
                "GlimmerHMM exited with non-zero status".to_string(),
            ));
        }
    }

    // Convert GlimmerHMM native output → GFF3
    let count = convert_glimmer_to_gff3(&raw_out, output_gff3)?;
    let _ = std::fs::remove_file(&raw_out);

    Ok(count)
}

/// Run GlimmerHMM in parallel by splitting the genome into per-contig FASTAs.
#[allow(unused_assignments)]
fn run_glimmerhmm_parallel(
    genome_fasta: &Path,
    training_dir: &Path,
    raw_out: &Path,
    threads: usize,
) -> Result<()> {
    use std::io::{BufRead, BufReader};

    // Split FASTA into contigs
    let tmp_dir = raw_out
        .parent()
        .unwrap_or(Path::new("."))
        .join("glimmer_tmp");
    std::fs::create_dir_all(&tmp_dir).map_err(MycoNoteError::Io)?;

    let file = std::fs::File::open(genome_fasta).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut contig_files: Vec<PathBuf> = vec![];
    let mut current_id = String::new();
    let mut current_out: Option<std::fs::File> = None;

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        if line.starts_with('>') {
            current_id = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            let p = tmp_dir.join(format!("{}.fa", current_id));
            let mut f = std::fs::File::create(&p).map_err(MycoNoteError::Io)?;
            writeln!(f, "{}", line).map_err(MycoNoteError::Io)?;
            contig_files.push(p.clone());
            current_out = Some(f);
        } else if let Some(ref mut f) = current_out {
            writeln!(f, "{}", line).map_err(MycoNoteError::Io)?;
        }
    }
    drop(current_out);

    // Run GlimmerHMM on each contig (simple sequential batching)
    let all_output = std::fs::File::create(raw_out).map_err(MycoNoteError::Io)?;
    let batch_size = (contig_files.len() + threads - 1) / threads;

    for chunk in contig_files.chunks(batch_size.max(1)) {
        for fasta in chunk {
            let status = Command::new("glimmerhmm")
                .arg(fasta)
                .arg(training_dir)
                .arg("-f")
                .stdout(std::process::Stdio::from(
                    all_output.try_clone().map_err(MycoNoteError::Io)?,
                ))
                .status()
                .map_err(|e| MycoNoteError::ExternalTool(format!("glimmerhmm: {}", e)))?;
            if !status.success() {
                eprintln!("  ⚠  GlimmerHMM failed on {}", fasta.display());
            }
        }
    }

    let _ = std::fs::remove_dir_all(&tmp_dir);
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// GlimmerHMM output → GFF3 conversion
// ─────────────────────────────────────────────────────────────────────────────

/// Convert GlimmerHMM GFF2-like output to proper GFF3 with gene/mRNA/exon/CDS hierarchy.
pub fn convert_glimmer_to_gff3(input: &Path, output: &Path) -> Result<usize> {
    let file = std::fs::File::open(input).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;

    writeln!(out, "##gff-version 3").map_err(MycoNoteError::Io)?;
    writeln!(out, "# Source: GlimmerHMM").map_err(MycoNoteError::Io)?;

    // GlimmerHMM output groups features by gene (lines starting with "##" are gene separators)
    // Format: seqname source feature start end score strand frame attributes
    let mut gene_idx = 0usize;
    let mut exon_idx = 0usize;
    let mut gene_count = 0usize;

    // Buffer all records for a single gene
    let mut current_gene: Vec<(String, u64, u64, char, String)> = vec![]; // (seqid, start, end, strand, type)
    let mut current_seqid = String::new();

    let flush_gene = |out: &mut std::fs::File,
                      seqid: &str,
                      records: &[(String, u64, u64, char, String)],
                      gene_idx: &mut usize,
                      exon_idx: &mut usize,
                      gene_count: &mut usize|
     -> std::io::Result<()> {
        if records.is_empty() {
            return Ok(());
        }

        let start = records.iter().map(|r| r.1).min().unwrap_or(0);
        let end = records.iter().map(|r| r.2).max().unwrap_or(0);
        let strand = records[0].3;

        *gene_idx += 1;
        let gene_id = format!("GLIMMER_{:06}", gene_idx);
        let mrna_id = format!("{}.mRNA1", gene_id);

        writeln!(
            out,
            "{}\tGlimmerHMM\tgene\t{}\t{}\t.\t{}\t.\tID={}",
            seqid, start, end, strand, gene_id
        )?;
        writeln!(
            out,
            "{}\tGlimmerHMM\tmRNA\t{}\t{}\t.\t{}\t.\tID={};Parent={}",
            seqid, start, end, strand, mrna_id, gene_id
        )?;

        for (_, fstart, fend, fstrand, ftype) in records {
            *exon_idx += 1;
            let ftype_gff3 = if ftype == "CDS" { "CDS" } else { "exon" };
            let feat_id = format!("{}.{}{}", mrna_id, ftype_gff3, exon_idx);
            writeln!(
                out,
                "{}\tGlimmerHMM\t{}\t{}\t{}\t.\t{}\t.\tID={};Parent={}",
                seqid, ftype_gff3, fstart, fend, fstrand, feat_id, mrna_id
            )?;
        }

        *gene_count += 1;
        Ok(())
    };

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();

        if trimmed.starts_with("##") {
            // Gene boundary — flush previous gene
            flush_gene(
                &mut out,
                &current_seqid,
                &current_gene,
                &mut gene_idx,
                &mut exon_idx,
                &mut gene_count,
            )
            .map_err(MycoNoteError::Io)?;
            current_gene.clear();
            continue;
        }
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let seqid = fields[0].to_string();
        let feat = fields[2].to_string();
        let start: u64 = fields[3].parse().unwrap_or(0);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand: char = fields[6].chars().next().unwrap_or('+');

        if feat == "mRNA" {
            // New gene starts — flush previous
            flush_gene(
                &mut out,
                &current_seqid,
                &current_gene,
                &mut gene_idx,
                &mut exon_idx,
                &mut gene_count,
            )
            .map_err(MycoNoteError::Io)?;
            current_gene.clear();
            current_seqid = seqid;
        } else if feat == "CDS" || feat == "exon" {
            current_seqid = seqid;
            current_gene.push((current_seqid.clone(), start, end, strand, feat));
        }
    }

    // Flush last gene
    flush_gene(
        &mut out,
        &current_seqid,
        &current_gene,
        &mut gene_idx,
        &mut exon_idx,
        &mut gene_count,
    )
    .map_err(MycoNoteError::Io)?;

    Ok(gene_count)
}
