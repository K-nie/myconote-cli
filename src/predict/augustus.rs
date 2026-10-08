/// Augustus gene prediction wrapper
///
/// Calls Augustus as a subprocess and parses its GFF3 output.
/// Supports:
///   - Standard prediction from masked FASTA
///   - Hint-based prediction (protein/EST hints improve accuracy)
///   - Multi-sequence parallelism via rayon (one contig → one Augustus call)
///
/// Augustus must be installed and in PATH:
///   conda install -c bioconda augustus
use crate::utils::error::{MycoNoteError, Result};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Augustus configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct AugustusConfig {
    /// Species model (e.g. "saccharomyces_cerevisiae_S288C", "arabidopsis")
    pub species: String,
    /// Number of parallel threads (split by contig)
    pub threads: usize,
    /// Use UTR prediction
    pub utr: bool,
    /// Path to hints GFF file (optional — protein/EST evidence)
    pub hints_file: Option<PathBuf>,
    /// Augustus `--extrinsicCfgFile` basename to use when `hints_file` is set.
    /// `None` keeps the historical default (`extrinsic.M.RM.E.W.cfg`). Protein
    /// (`src=P`) hints need a config whose `[SOURCES]` block lists `P`, e.g.
    /// `extrinsic.M.RM.E.W.P.cfg`, which ships with Augustus.
    pub extrinsic_cfg: Option<String>,
    /// Extra Augustus arguments (passed verbatim)
    pub extra_args: Vec<String>,
}

/// Historical default extrinsic config, used whenever `extrinsic_cfg` is `None`
/// and a hints file is supplied. Kept as a constant so the single- and
/// multi-contig paths stay in lockstep.
const DEFAULT_EXTRINSIC_CFG: &str = "extrinsic.M.RM.E.W.cfg";

impl AugustusConfig {
    /// Resolve the `--extrinsicCfgFile` basename: the caller's override if set,
    /// otherwise the historical default. Only consulted when `hints_file` is set.
    pub fn extrinsic_cfg_name(&self) -> &str {
        self.extrinsic_cfg
            .as_deref()
            .unwrap_or(DEFAULT_EXTRINSIC_CFG)
    }
}

impl Default for AugustusConfig {
    fn default() -> Self {
        Self {
            species: "saccharomyces_cerevisiae_S288C".to_string(),
            threads: 4,
            utr: false,
            hints_file: None,
            extrinsic_cfg: None,
            extra_args: Vec::new(),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run Augustus on a masked FASTA and return path to the GFF3 output.
///
/// If `config.threads > 1` AND the input has multiple contigs, the FASTA is
/// split per-contig and each contig's Augustus call runs on its own rayon
/// thread. Per-contig gene IDs are rewritten with a contig-index prefix
/// before concatenation so the merged GFF3 has globally unique IDs — a flat
/// `cat` would otherwise produce colliding `g1`, `g2` identifiers across
/// contigs and break the downstream EVM merge's Parent chain tracking.
pub fn run(masked_fasta: &Path, output_gff: &Path, config: &AugustusConfig) -> Result<()> {
    let aug = which::which("augustus").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "augustus not found in PATH.\n\
             Install with: conda install -c bioconda augustus\n\
             Then verify with: augustus --species=help"
                .to_string(),
        )
    })?;

    let n_threads = config.threads.max(1);

    // Count contigs without reading the whole FASTA — decide single-threaded
    // vs per-contig-parallel based on the actual input shape.
    let records = crate::parser::fasta::read_fasta(masked_fasta)?;
    let n_contigs = records.len();

    if n_threads <= 1 || n_contigs <= 1 {
        return run_single(&aug, masked_fasta, output_gff, config);
    }

    println!(
        "  Running Augustus (species: {}) — {} contigs × {} threads in parallel…",
        config.species, n_contigs, n_threads
    );
    run_parallel(&aug, &records, output_gff, config, n_threads)
}

/// Original single-invocation path — used when threads=1 or there's only one
/// contig to work on.
fn run_single(
    aug: &Path,
    masked_fasta: &Path,
    output_gff: &Path,
    config: &AugustusConfig,
) -> Result<()> {
    println!("  Running Augustus (species: {})…", config.species);

    let utr_flag = if config.utr { "on" } else { "off" };

    let mut args = vec![
        format!("--species={}", config.species),
        "--gff3=on".to_string(),
        "--strand=both".to_string(),
        format!("--UTR={}", utr_flag),
        "--softmasking=1".to_string(),
        format!("--outfile={}", output_gff.display()),
    ];

    if let Some(ref hints) = config.hints_file {
        args.push(format!("--hintsfile={}", hints.display()));
        args.push(format!(
            "--extrinsicCfgFile={}",
            config.extrinsic_cfg_name()
        ));
    }

    args.extend(config.extra_args.clone());
    args.push(masked_fasta.to_string_lossy().into_owned());

    let status = Command::new(aug)
        .args(&args)
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Augustus failed. Check that species '{}' is installed.\n\
             List available species: augustus --species=help",
            config.species
        )));
    }

    println!("  Augustus finished → {}", output_gff.display());
    Ok(())
}

/// Per-contig parallel path. Writes each contig to a temp FASTA, runs
/// Augustus on each in parallel, rewrites gene IDs with a contig-index
/// prefix for uniqueness, then concatenates into the final GFF3.
fn run_parallel(
    aug: &Path,
    records: &[crate::parser::fasta::FastaRecord],
    output_gff: &Path,
    config: &AugustusConfig,
    n_threads: usize,
) -> Result<()> {
    use rayon::prelude::*;

    // Unique temp dir per call — nanosecond suffix avoids concurrent runs
    // clobbering each other.
    let suffix = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let tmp_dir = std::env::temp_dir().join(format!("myconote_augustus_{}", suffix));
    std::fs::create_dir_all(&tmp_dir).map_err(MycoNoteError::Io)?;

    // 1. Stage each contig as its own FASTA file.
    let mut contig_paths: Vec<(usize, PathBuf, PathBuf)> = Vec::with_capacity(records.len());
    for (idx, rec) in records.iter().enumerate() {
        let fa_path = tmp_dir.join(format!("contig_{}.fa", idx));
        let gff_path = tmp_dir.join(format!("contig_{}.gff3", idx));
        let mut f = std::fs::File::create(&fa_path).map_err(MycoNoteError::Io)?;
        writeln!(f, ">{}", rec.id).map_err(MycoNoteError::Io)?;
        for chunk in rec.sequence.as_bytes().chunks(60) {
            f.write_all(chunk).map_err(MycoNoteError::Io)?;
            writeln!(f).map_err(MycoNoteError::Io)?;
        }
        contig_paths.push((idx, fa_path, gff_path));
    }

    // 2. Run Augustus per contig, bounded by `n_threads` via a rayon pool.
    //    Each task is single-threaded Augustus; parallelism comes from the
    //    multiple contigs running concurrently.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .map_err(|e| MycoNoteError::ExternalTool(format!("rayon pool: {}", e)))?;

    let aug_path = aug.to_path_buf();
    let species = config.species.clone();
    let utr_flag = if config.utr { "on" } else { "off" };
    let utr_flag_owned = utr_flag.to_string();
    let hints = config.hints_file.clone();
    let extrinsic_cfg = config.extrinsic_cfg_name().to_string();
    let extra = config.extra_args.clone();

    let failures: Vec<String> = pool.install(|| {
        contig_paths
            .par_iter()
            .filter_map(|(idx, fa_path, gff_path)| {
                let mut args = vec![
                    format!("--species={}", species),
                    "--gff3=on".to_string(),
                    "--strand=both".to_string(),
                    format!("--UTR={}", utr_flag_owned),
                    "--softmasking=1".to_string(),
                    format!("--outfile={}", gff_path.display()),
                ];
                if let Some(ref h) = hints {
                    args.push(format!("--hintsfile={}", h.display()));
                    args.push(format!("--extrinsicCfgFile={}", extrinsic_cfg));
                }
                args.extend(extra.clone());
                args.push(fa_path.to_string_lossy().into_owned());

                let result = Command::new(&aug_path).args(&args).output();
                match result {
                    Ok(out) if out.status.success() => None,
                    Ok(out) => Some(format!(
                        "contig {}: exit {:?}\n{}",
                        idx,
                        out.status.code(),
                        String::from_utf8_lossy(&out.stderr)
                            .chars()
                            .take(400)
                            .collect::<String>()
                    )),
                    Err(e) => Some(format!("contig {}: spawn failed — {}", idx, e)),
                }
            })
            .collect()
    });

    if !failures.is_empty() {
        eprintln!(
            "  ⚠  {} of {} contigs failed:",
            failures.len(),
            contig_paths.len()
        );
        for msg in failures.iter().take(3) {
            eprintln!("     {}", msg);
        }
    }

    // 3. Concatenate, rewriting per-contig gene IDs so they don't collide.
    if let Some(parent) = output_gff.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(MycoNoteError::Io)?;
        }
    }
    let mut out = std::fs::File::create(output_gff).map_err(MycoNoteError::Io)?;
    writeln!(out, "##gff-version 3").map_err(MycoNoteError::Io)?;

    for (idx, _, gff_path) in &contig_paths {
        if !gff_path.exists() {
            continue;
        }
        let body = std::fs::read_to_string(gff_path).map_err(MycoNoteError::Io)?;
        for line in body.lines() {
            if line.is_empty() || line.starts_with("##gff-version") {
                continue;
            }
            if line.starts_with('#') {
                // Augustus section markers like "# start gene g1" and
                // "# end gene g1" — leave alone (they're informational).
                writeln!(out, "{}", line).map_err(MycoNoteError::Io)?;
                continue;
            }
            // Rewrite attribute IDs and Parent refs with a contig prefix so
            // the merged file has globally unique identifiers. Augustus uses
            // IDs like `g1`, `g1.t1`, `g1.t1.cds` — prefix them all.
            let rewritten = rewrite_attrs_with_prefix(line, *idx);
            writeln!(out, "{}", rewritten).map_err(MycoNoteError::Io)?;
        }
    }

    // 4. Clean up the staging directory.
    let _ = std::fs::remove_dir_all(&tmp_dir);

    println!("  Augustus finished → {}", output_gff.display());
    Ok(())
}

/// Prepend `c{idx}_` to every GFF3 attribute value for keys `ID` and
/// `Parent` — used to disambiguate identifiers across per-contig augustus
/// runs so the concatenated output has globally unique IDs.
fn rewrite_attrs_with_prefix(line: &str, contig_idx: usize) -> String {
    let mut cols: Vec<&str> = line.splitn(9, '\t').collect();
    if cols.len() < 9 {
        return line.to_string();
    }
    let rewritten: String = cols[8]
        .split(';')
        .map(|attr| {
            if let Some(rest) = attr.strip_prefix("ID=") {
                format!("ID=c{}_{}", contig_idx, rest)
            } else if let Some(rest) = attr.strip_prefix("Parent=") {
                // Parent may be comma-separated; apply prefix to each.
                let prefixed: Vec<String> = rest
                    .split(',')
                    .map(|p| format!("c{}_{}", contig_idx, p))
                    .collect();
                format!("Parent={}", prefixed.join(","))
            } else {
                attr.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join(";");
    cols[8] = ""; // replace column 9
    let mut out = cols.join("\t");
    // Re-insert the rewritten attributes (cols[8] is now empty so we append).
    out.push_str(&rewritten);
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Protein hint generation from homolog evidence
// ─────────────────────────────────────────────────────────────────────────────

/// Generate an Augustus hints file from a protein alignment (BLAST tabular fmt6).
/// Each hit becomes a CDSpart hint weighted by alignment identity.
pub fn make_protein_hints(blast_tsv: &Path, hints_out: &Path, priority: u8) -> Result<()> {
    use std::io::BufRead;

    let input = std::fs::File::open(blast_tsv).map_err(MycoNoteError::Io)?;
    let mut out = std::fs::File::create(hints_out).map_err(MycoNoteError::Io)?;

    for line in std::io::BufReader::new(input).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        // BLAST fmt6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        if f.len() < 12 {
            continue;
        }

        let seqid = f[0];
        let start: u64 = f[6].parse().unwrap_or(0);
        let end: u64 = f[7].parse().unwrap_or(0);
        let (s, e) = if start <= end {
            (start, end)
        } else {
            (end, start)
        };

        writeln!(
            out,
            "{}\tP\tCDSpart\t{}\t{}\t.\t.\t.\tsrc=P;pri={}",
            seqid, s, e, priority
        )
        .map_err(MycoNoteError::Io)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rewrite_attrs_prefixes_id_and_parent() {
        let line = "chr1\tAUGUSTUS\tgene\t100\t900\t0.9\t+\t.\tID=g1";
        let out = rewrite_attrs_with_prefix(line, 3);
        assert!(out.ends_with("ID=c3_g1"), "got: {}", out);

        let line = "chr1\tAUGUSTUS\tmRNA\t100\t900\t0.9\t+\t.\tID=g1.t1;Parent=g1";
        let out = rewrite_attrs_with_prefix(line, 7);
        assert!(out.contains("ID=c7_g1.t1"));
        assert!(out.contains("Parent=c7_g1"));
    }

    #[test]
    fn rewrite_attrs_handles_comma_separated_parents() {
        // CDS shared between two isoforms — both parents must be prefixed.
        let line = "chr1\tAUGUSTUS\tCDS\t100\t500\t.\t+\t0\tID=cds1;Parent=g1.t1,g1.t2";
        let out = rewrite_attrs_with_prefix(line, 2);
        assert!(out.contains("Parent=c2_g1.t1,c2_g1.t2"), "got: {}", out);
    }

    #[test]
    fn rewrite_attrs_preserves_other_attributes() {
        let line = "chr1\tAUGUSTUS\tgene\t100\t900\t0.9\t+\t.\tID=g1;Name=myGene;locus_tag=X_0001";
        let out = rewrite_attrs_with_prefix(line, 0);
        assert!(out.contains("Name=myGene"));
        assert!(out.contains("locus_tag=X_0001"));
        assert!(out.contains("ID=c0_g1"));
    }

    #[test]
    fn rewrite_attrs_passes_through_short_lines() {
        let not_gff3 = "##gff-version 3";
        assert_eq!(rewrite_attrs_with_prefix(not_gff3, 1), not_gff3);
    }

    // Feature 3: extrinsic config resolution. Default (None) keeps the historical
    // basename so the existing hints path is byte-for-byte unchanged; an override
    // (used for protein `src=P` hints) is honoured verbatim.
    #[test]
    fn extrinsic_cfg_defaults_when_unset() {
        let cfg = AugustusConfig::default();
        assert_eq!(cfg.extrinsic_cfg_name(), "extrinsic.M.RM.E.W.cfg");
    }

    #[test]
    fn extrinsic_cfg_override_is_honoured() {
        let cfg = AugustusConfig {
            extrinsic_cfg: Some("extrinsic.M.RM.E.W.P.cfg".to_string()),
            ..AugustusConfig::default()
        };
        assert_eq!(cfg.extrinsic_cfg_name(), "extrinsic.M.RM.E.W.P.cfg");
    }
}
