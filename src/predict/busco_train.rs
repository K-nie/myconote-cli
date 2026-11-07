//! BUSCO-seeded per-genome Augustus training (v0.9.0, Feature #1).
//!
//! This is funannotate's accuracy mechanism, and most likely its real edge over
//! a stock-species prediction: rather than training Augustus on the tool's own
//! first-pass calls (the circular `--augustus-training self` mode), it trains on
//! an INDEPENDENT high-quality gene set — the genome's own COMPLETE single-copy
//! BUSCO genes. Those are conserved, well-structured, single-transcript loci, so
//! they make honest training material for a genome-specific Augustus model.
//!
//! Pipeline:
//!   1. Run BUSCO in genome mode with the Augustus gene finder
//!      (`-m genome --augustus`) so BUSCO emits per-gene structured models.
//!   2. Read `full_table.tsv`, take the Complete single-copy BUSCOs.
//!   3. Convert BUSCO's Augustus gene models that overlap those complete hits
//!      into a GFF3 training set.
//!   4. Hand that GFF3 + the assembly to the existing Augustus-training plumbing
//!      (`predict::train::run_training`), which does
//!      `gff → mini-GenBank → new_species.pl → etraining → (optional optimize)`.
//!   5. Predict with the per-genome-trained model.
//!
//! BUSCO is already a MycoNote dependency (the annotate stage uses it), so this
//! adds no new hard dependency. The whole feature is opt-in
//! (`--augustus-training busco`) and falls back loudly to the stock species on
//! any failure, so a default run is unaffected.
//!
//! VALIDATION CAVEAT: the live BUSCO directory layout and the exact Augustus
//! GFF dialect vary across BUSCO versions and gene-finder backends. The parsing
//! and conversion logic here is unit-tested on synthetic fixtures, but the
//! end-to-end path (real `busco -m genome --augustus` output) must be confirmed
//! on the cluster before the trained model is trusted for a benchmark.

use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// 0. Pre-flight: BUSCO ↔ Augustus version compatibility
// ─────────────────────────────────────────────────────────────────────────────

/// First BUSCO release whose `--augustus` pipeline needs a newer Augustus.
/// BUSCO 5.7.x calls Augustus options that 3.3.3 does not have, so the genome
/// run "completes" without usable gene models and training falls back only
/// after the expensive BUSCO run.
const BUSCO_NEEDS_NEW_AUGUSTUS: (u32, u32, u32) = (5, 7, 0);
/// Minimum Augustus those BUSCO releases work with.
const MIN_AUGUSTUS_FOR_NEW_BUSCO: (u32, u32, u32) = (3, 4, 0);

/// Pull the first `X.Y[.Z]` version number out of a tool's `--version` text,
/// e.g. "BUSCO 5.7.1" or "AUGUSTUS (3.3.3) is a gene prediction tool".
fn parse_version(text: &str) -> Option<(u32, u32, u32)> {
    for token in text.split(|c: char| !(c.is_ascii_digit() || c == '.')) {
        let parts: Vec<&str> = token.split('.').filter(|p| !p.is_empty()).collect();
        if parts.len() < 2 {
            continue;
        }
        let nums: Vec<u32> = parts
            .iter()
            .take(3)
            .filter_map(|p| p.parse().ok())
            .collect();
        if nums.len() >= 2 {
            return Some((nums[0], nums[1], *nums.get(2).unwrap_or(&0)));
        }
    }
    None
}

fn fmt_version(v: (u32, u32, u32)) -> String {
    format!("{}.{}.{}", v.0, v.1, v.2)
}

/// Decide whether a BUSCO/Augustus pair can run BUSCO's `--augustus` mode.
/// Returns a human-readable reason when it cannot. Unknown versions are not
/// treated as a problem: the existing post-run fallback still covers them.
fn version_incompatibility(
    busco: Option<(u32, u32, u32)>,
    augustus: Option<(u32, u32, u32)>,
) -> Option<String> {
    let (b, a) = (busco?, augustus?);
    if b >= BUSCO_NEEDS_NEW_AUGUSTUS && a < MIN_AUGUSTUS_FOR_NEW_BUSCO {
        return Some(format!(
            "BUSCO {} needs Augustus >= {} for `--augustus` gene finding, but the \
             Augustus BUSCO would use is {}. Either install Augustus >= {} into \
             BUSCO's environment, or pin BUSCO below {} (e.g. `busco=5.6`)",
            fmt_version(b),
            fmt_version(MIN_AUGUSTUS_FOR_NEW_BUSCO),
            fmt_version(a),
            fmt_version(MIN_AUGUSTUS_FOR_NEW_BUSCO),
            fmt_version(BUSCO_NEEDS_NEW_AUGUSTUS),
        ));
    }
    None
}

/// Run a command and return its combined stdout + stderr (Augustus prints its
/// version banner to stderr on some builds and exits non-zero).
fn version_text(mut cmd: std::process::Command) -> Option<String> {
    let out = cmd.output().ok()?;
    let mut s = String::from_utf8_lossy(&out.stdout).into_owned();
    s.push('\n');
    s.push_str(&String::from_utf8_lossy(&out.stderr));
    Some(s)
}

/// The Augustus binary BUSCO itself will pick up: the one inside BUSCO's
/// conda env (or beside a direct `busco` binary), else `augustus` on PATH.
fn augustus_for_busco(busco: &crate::annotate::busco::BuscoCmd) -> Option<PathBuf> {
    use crate::annotate::busco::BuscoCmd;
    let beside = match busco {
        BuscoCmd::Direct(p) => p.parent().map(|d| d.join("augustus")),
        BuscoCmd::CondaRun {
            conda_path,
            env_name,
        } => conda_path
            .parent()
            .and_then(|bin| bin.parent())
            .map(|root| {
                root.join("envs")
                    .join(env_name)
                    .join("bin")
                    .join("augustus")
            }),
    };
    beside
        .filter(|p| p.is_file())
        .or_else(|| which::which("augustus").ok())
}

/// Pre-flight check run before the (expensive) BUSCO genome run. Returns
/// `Some(reason)` when BUSCO-seeded training is known not to work with the
/// installed BUSCO/Augustus pair.
pub fn busco_augustus_incompatibility() -> Option<String> {
    let busco = crate::annotate::busco::find_busco().ok()?;
    let busco_v = version_text(busco.command(&["--version"])).and_then(|t| parse_version(&t));
    let aug_v = augustus_for_busco(&busco).and_then(|p| {
        let mut cmd = std::process::Command::new(p);
        cmd.arg("--version");
        version_text(cmd).and_then(|t| parse_version(&t))
    });
    version_incompatibility(busco_v, aug_v)
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Run BUSCO in genome mode (Augustus gene finder)
// ─────────────────────────────────────────────────────────────────────────────

/// Run `busco -m genome --augustus` on an assembly and return the run directory
/// that holds `full_table.tsv` and the Augustus gene models. Tries `--offline`
/// first (the cluster convention) and retries once with download allowed if the
/// lineage dataset is missing.
pub fn run_busco_genome(
    assembly: &Path,
    lineage: &str,
    out_dir: &Path,
    threads: usize,
    seed_species: &str,
) -> Result<PathBuf> {
    let busco = crate::annotate::busco::find_busco().map_err(MycoNoteError::UnsupportedFormat)?;
    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    let run_name = "myconote_busco_genome";
    let abs_fa = std::fs::canonicalize(assembly).map_err(MycoNoteError::Io)?;
    let fa_str = abs_fa.to_string_lossy().into_owned();
    let threads_s = threads.to_string();

    println!(
        "  Running BUSCO genome mode (lineage: {}, seed species: {})…",
        lineage, seed_species
    );

    // --augustus forces the Augustus gene-finding pipeline (vs. the v5 metaeuk
    // default), which is what produces the structured per-gene models we train
    // from. --augustus_species seeds that first pass with a clade-appropriate
    // model so divergent fungi do not bootstrap from S288C.
    let base_args: Vec<String> = vec![
        "-m".into(),
        "genome".into(),
        "-i".into(),
        fa_str.clone(),
        "-o".into(),
        run_name.into(),
        "-l".into(),
        lineage.into(),
        "--cpu".into(),
        threads_s.clone(),
        "--augustus".into(),
        "--augustus_species".into(),
        seed_species.into(),
        "--force".into(),
    ];

    // First attempt: offline (cluster default — datasets pre-provisioned).
    let mut offline_args: Vec<&str> = base_args.iter().map(|s| s.as_str()).collect();
    offline_args.push("--offline");
    let _ = busco
        .command(&offline_args)
        .current_dir(out_dir)
        .status()
        .map_err(MycoNoteError::Io)?;

    if let Some(d) = find_run_dir(out_dir, run_name, lineage) {
        return Ok(d);
    }

    // Retry allowing dataset download.
    println!("  BUSCO offline run produced no results dir — retrying with dataset download…");
    let online_args: Vec<&str> = base_args.iter().map(|s| s.as_str()).collect();
    let status = busco
        .command(&online_args)
        .current_dir(out_dir)
        .status()
        .map_err(MycoNoteError::Io)?;
    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "BUSCO genome-mode run failed. Check that the lineage dataset is \
             available and that `--augustus` is supported by this BUSCO build."
                .to_string(),
        ));
    }

    find_run_dir(out_dir, run_name, lineage).ok_or_else(|| {
        MycoNoteError::InvalidFormat(format!(
            "BUSCO finished but no run_<lineage> directory was found under {}/{}",
            out_dir.display(),
            run_name
        ))
    })
}

/// Locate the `run_<lineage>` directory BUSCO 5 writes under `<out>/<run>/`.
/// Falls back to scanning for any child directory whose name starts with
/// `run_`, which covers minor version-to-version layout drift.
fn find_run_dir(out_dir: &Path, run_name: &str, lineage: &str) -> Option<PathBuf> {
    let base = out_dir.join(run_name);
    let direct = base.join(format!("run_{}", lineage));
    if direct.is_dir() {
        return Some(direct);
    }
    if let Ok(entries) = std::fs::read_dir(&base) {
        for e in entries.flatten() {
            let p = e.path();
            if p.is_dir() {
                if let Some(name) = p.file_name().and_then(|n| n.to_str()) {
                    if name.starts_with("run_") {
                        return Some(p);
                    }
                }
            }
        }
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. full_table.tsv — the Complete single-copy BUSCO hits
// ─────────────────────────────────────────────────────────────────────────────

/// One Complete single-copy BUSCO's genomic location, from `full_table.tsv`.
#[derive(Debug, Clone)]
struct BuscoHit {
    seqid: String,
    start: u64,
    end: u64,
}

/// Parse `full_table.tsv`, keeping only `Complete` (single-copy) BUSCOs with a
/// usable sequence + coordinate. `Duplicated` is deliberately excluded: a
/// duplicated BUSCO is complete but multi-copy, and training on paralogs muddies
/// the model (funannotate uses complete single-copy only for the same reason).
fn parse_full_table(path: &Path) -> Result<Vec<BuscoHit>> {
    let content = std::fs::read_to_string(path).map_err(MycoNoteError::Io)?;
    let mut hits = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        // Genome-mode layout: Busco id | Status | Sequence | Start | End | Strand | …
        if cols.len() < 5 {
            continue;
        }
        if cols[1] != "Complete" {
            continue;
        }
        let seqid = cols[2].to_string();
        let (start, end) = match (cols[3].parse::<u64>(), cols[4].parse::<u64>()) {
            (Ok(a), Ok(b)) => (a.min(b), a.max(b)),
            _ => continue,
        };
        if seqid.is_empty() || start == 0 {
            continue;
        }
        hits.push(BuscoHit { seqid, start, end });
    }
    Ok(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Augustus / GFF gene-model parsing + conversion to a training GFF3
// ─────────────────────────────────────────────────────────────────────────────

/// A transcript's CDS structure, keyed while parsing by a transcript id.
#[derive(Debug, Clone, Default)]
struct TxModel {
    seqid: String,
    strand: char,
    source: String,
    /// CDS segments as (start, end, phase); 1-based inclusive.
    cds: Vec<(u64, u64, u8)>,
}

/// Extract a transcript key from a GFF/GTF attribute column, tolerating both
/// dialects BUSCO's Augustus output can use:
///   * GTF:   `transcript_id "g1.t1"; gene_id "g1";`
///   * GFF3:  `ID=g1.t1;Parent=g1`  (for the CDS: `Parent=g1.t1`)
///   * bare:  Augustus gene/transcript lines whose column 9 is just `g1.t1`.
fn extract_tx_key(attrs: &str, feature: &str) -> Option<String> {
    // GTF transcript_id "…"
    if let Some(p) = attrs.find("transcript_id") {
        let rest = &attrs[p + "transcript_id".len()..];
        if let Some(q1) = rest.find('"') {
            if let Some(q2) = rest[q1 + 1..].find('"') {
                return Some(rest[q1 + 1..q1 + 1 + q2].to_string());
            }
        }
    }
    // GFF3: a CDS/exon points at its transcript via Parent; a transcript names
    // itself via ID.
    let want = if feature == "CDS" || feature == "exon" {
        "Parent="
    } else {
        "ID="
    };
    if let Some(p) = attrs.find(want) {
        let rest = &attrs[p + want.len()..];
        let end = rest.find(';').unwrap_or(rest.len());
        let v = rest[..end].trim();
        if !v.is_empty() {
            return Some(v.to_string());
        }
    }
    // Bare Augustus token (e.g. `g1.t1`), no key=value structure. The GFF3
    // empty-attribute placeholder "." is not a key.
    let bare = attrs.trim();
    if !bare.is_empty() && bare != "." && !bare.contains('=') && !bare.contains(' ') {
        return Some(bare.to_string());
    }
    None
}

/// Parse one GFF/GTF file into transcript models keyed by transcript id,
/// collecting only CDS segments (all the trainer needs). Non-9-column and
/// comment lines are skipped. Returns transcripts in first-seen order for
/// deterministic output.
fn parse_gff_transcripts(path: &Path) -> Vec<TxModel> {
    let content = match std::fs::read_to_string(path) {
        Ok(c) => c,
        Err(_) => return Vec::new(),
    };
    let mut order: Vec<String> = Vec::new();
    let mut map: HashMap<String, TxModel> = HashMap::new();

    for line in content.lines() {
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 9 {
            continue;
        }
        if f[2] != "CDS" {
            continue; // the trainer derives everything from CDS
        }
        let (start, end) = match (f[3].parse::<u64>(), f[4].parse::<u64>()) {
            (Ok(a), Ok(b)) if a <= b => (a, b),
            _ => continue,
        };
        let strand = f[6].chars().next().unwrap_or('.');
        let phase = f[7].parse::<u8>().unwrap_or(0);
        let key = match extract_tx_key(f[8], "CDS") {
            Some(k) => k,
            None => continue,
        };
        let entry = map.entry(key.clone()).or_insert_with(|| {
            order.push(key.clone());
            TxModel {
                seqid: f[0].to_string(),
                strand,
                source: f[1].to_string(),
                cds: Vec::new(),
            }
        });
        entry.cds.push((start, end, phase));
    }

    order
        .into_iter()
        .filter_map(|k| map.remove(&k))
        .filter(|t| !t.cds.is_empty())
        .collect()
}

/// Recursively collect every `.gff` / `.gff3` file under a directory.
fn collect_gff_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let entries = match std::fs::read_dir(dir) {
        Ok(e) => e,
        Err(_) => return,
    };
    for e in entries.flatten() {
        let p = e.path();
        if p.is_dir() {
            collect_gff_files(&p, out);
        } else if let Some(ext) = p.extension().and_then(|x| x.to_str()) {
            if ext == "gff" || ext == "gff3" {
                out.push(p);
            }
        }
    }
}

/// True when a transcript's span overlaps any Complete BUSCO hit on the same
/// sequence. Keeps the training set to complete single-copy loci.
fn overlaps_complete_hit(tx: &TxModel, hits: &[BuscoHit]) -> bool {
    let tx_start = tx.cds.iter().map(|c| c.0).min().unwrap_or(0);
    let tx_end = tx.cds.iter().map(|c| c.1).max().unwrap_or(0);
    hits.iter()
        .any(|h| h.seqid == tx.seqid && tx_start <= h.end && tx_end >= h.start)
}

/// Build the GFF3 training set from a BUSCO genome-mode run directory.
///
/// Reads `full_table.tsv` for the Complete single-copy BUSCOs, parses every
/// Augustus GFF under the run dir, keeps the transcripts overlapping a complete
/// hit, and writes one `gene → mRNA → CDS` GFF3 block per kept transcript (the
/// shape `predict::train::load_gene_models` expects). Returns the number of
/// training genes written.
pub fn collect_busco_training_gff(run_dir: &Path, out_gff: &Path) -> Result<usize> {
    let full_table = run_dir.join("full_table.tsv");
    if !full_table.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "BUSCO full_table.tsv not found under {} — cannot build a training set.",
            run_dir.display()
        )));
    }
    let hits = parse_full_table(&full_table)?;
    if hits.is_empty() {
        return Err(MycoNoteError::InvalidFormat(
            "BUSCO found no Complete single-copy genes — nothing to train on. \
             The assembly may be too fragmented, or the wrong lineage was used."
                .to_string(),
        ));
    }

    let mut gff_files = Vec::new();
    collect_gff_files(run_dir, &mut gff_files);
    gff_files.sort(); // deterministic order
    if gff_files.is_empty() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "BUSCO run dir {} holds no Augustus GFF files. Was BUSCO run with \
             `--augustus`? (The metaeuk backend does not emit trainable models.)",
            run_dir.display()
        )));
    }

    // Parse all transcripts, keep complete-hit overlappers, de-duplicate by
    // genomic span so the same BUSCO gene parsed from two files is not written
    // twice.
    let mut written = 0usize;
    let mut seen: std::collections::HashSet<(String, u64, u64, char)> =
        std::collections::HashSet::new();
    let mut buf = String::from("##gff-version 3\n");

    for gff in &gff_files {
        for mut tx in parse_gff_transcripts(gff) {
            if tx.cds.is_empty() || !overlaps_complete_hit(&tx, &hits) {
                continue;
            }
            tx.cds.sort_by_key(|c| c.0);
            let gstart = tx.cds.first().unwrap().0;
            let gend = tx.cds.last().unwrap().1;
            let key = (tx.seqid.clone(), gstart, gend, tx.strand);
            if !seen.insert(key) {
                continue;
            }
            written += 1;
            let gid = format!("busco_g{}", written);
            let source = if tx.source.is_empty() {
                "BUSCO"
            } else {
                tx.source.as_str()
            };
            let strand = if tx.strand == '+' || tx.strand == '-' {
                tx.strand
            } else {
                '+'
            };
            buf.push_str(&format!(
                "{seq}\t{src}\tgene\t{s}\t{e}\t.\t{st}\t.\tID={gid}\n",
                seq = tx.seqid,
                src = source,
                s = gstart,
                e = gend,
                st = strand,
                gid = gid
            ));
            buf.push_str(&format!(
                "{seq}\t{src}\tmRNA\t{s}\t{e}\t.\t{st}\t.\tID={gid}.t1;Parent={gid}\n",
                seq = tx.seqid,
                src = source,
                s = gstart,
                e = gend,
                st = strand,
                gid = gid
            ));
            for (i, (cs, ce, ph)) in tx.cds.iter().enumerate() {
                buf.push_str(&format!(
                    "{seq}\t{src}\tCDS\t{s}\t{e}\t.\t{st}\t{ph}\tID={gid}.cds{n};Parent={gid}.t1\n",
                    seq = tx.seqid,
                    src = source,
                    s = cs,
                    e = ce,
                    st = strand,
                    ph = ph,
                    gid = gid,
                    n = i + 1
                ));
            }
        }
    }

    std::fs::write(out_gff, buf).map_err(MycoNoteError::Io)?;
    Ok(written)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests (synthetic fixtures — the parsing/conversion logic, not the live run)
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_tx_key_handles_gtf_gff3_and_bare() {
        assert_eq!(
            extract_tx_key("transcript_id \"g1.t1\"; gene_id \"g1\";", "CDS").as_deref(),
            Some("g1.t1")
        );
        assert_eq!(
            extract_tx_key("ID=g2.cds1;Parent=g2.t1", "CDS").as_deref(),
            Some("g2.t1")
        );
        assert_eq!(
            extract_tx_key("ID=g3.t1;Parent=g3", "mRNA").as_deref(),
            Some("g3.t1")
        );
        assert_eq!(
            extract_tx_key("g4.t1", "transcript").as_deref(),
            Some("g4.t1")
        );
        assert_eq!(extract_tx_key(".", "CDS"), None);
    }

    #[test]
    fn parse_full_table_keeps_only_complete_single_copy() {
        let dir = tempfile::tempdir().unwrap();
        let ft = dir.path().join("full_table.tsv");
        std::fs::write(
            &ft,
            "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n\
             BUSCOa\tComplete\tchr1\t100\t900\t+\t500.0\t800\n\
             BUSCOb\tDuplicated\tchr1\t2000\t2500\t+\t400.0\t500\n\
             BUSCOc\tFragmented\tchr2\t10\t50\t-\t50.0\t40\n\
             BUSCOd\tMissing\t\t\t\t\t\t\n\
             BUSCOe\tComplete\tchr2\t3000\t3600\t-\t600.0\t600\n",
        )
        .unwrap();
        let hits = parse_full_table(&ft).unwrap();
        assert_eq!(hits.len(), 2, "only the two Complete entries are kept");
        assert_eq!(hits[0].seqid, "chr1");
        assert_eq!((hits[0].start, hits[0].end), (100, 900));
        assert_eq!(hits[1].seqid, "chr2");
    }

    #[test]
    fn parse_gff_transcripts_groups_cds_by_transcript() {
        let dir = tempfile::tempdir().unwrap();
        let gff = dir.path().join("BUSCOa.gff");
        // Augustus GTF-style output (two CDS exons, one transcript).
        std::fs::write(
            &gff,
            "# start gene g1\n\
             chr1\tAUGUSTUS\tgene\t100\t900\t.\t+\t.\tg1\n\
             chr1\tAUGUSTUS\ttranscript\t100\t900\t.\t+\t.\tg1.t1\n\
             chr1\tAUGUSTUS\tCDS\t100\t400\t.\t+\t0\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n\
             chr1\tAUGUSTUS\tCDS\t600\t900\t.\t+\t2\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n",
        )
        .unwrap();
        let txs = parse_gff_transcripts(&gff);
        assert_eq!(txs.len(), 1);
        assert_eq!(txs[0].seqid, "chr1");
        assert_eq!(txs[0].strand, '+');
        assert_eq!(txs[0].cds.len(), 2);
    }

    #[test]
    fn collect_training_gff_writes_complete_overlappers_only() {
        let dir = tempfile::tempdir().unwrap();
        let run = dir.path();
        std::fs::write(
            run.join("full_table.tsv"),
            "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\n\
             BUSCOa\tComplete\tchr1\t100\t900\t+\n",
        )
        .unwrap();
        // Two transcripts: one overlaps the complete hit, one does not.
        let gdir = run.join("augustus_output").join("gff");
        std::fs::create_dir_all(&gdir).unwrap();
        std::fs::write(
            gdir.join("models.gff"),
            "chr1\tAUGUSTUS\tCDS\t150\t400\t.\t+\t0\ttranscript_id \"g1.t1\";\n\
             chr1\tAUGUSTUS\tCDS\t600\t850\t.\t+\t0\ttranscript_id \"g1.t1\";\n\
             chr9\tAUGUSTUS\tCDS\t50\t200\t.\t-\t0\ttranscript_id \"gX.t1\";\n",
        )
        .unwrap();

        let out = run.join("busco_training.gff3");
        let n = collect_busco_training_gff(run, &out).unwrap();
        assert_eq!(n, 1, "only the complete-hit overlapper is written");
        let gff = std::fs::read_to_string(&out).unwrap();
        assert!(gff.contains("##gff-version 3"));
        assert_eq!(gff.matches("\tgene\t").count(), 1);
        assert_eq!(gff.matches("\tmRNA\t").count(), 1);
        assert_eq!(gff.matches("\tCDS\t").count(), 2);
        assert!(gff.contains("ID=busco_g1"));
        assert!(gff.contains("Parent=busco_g1.t1"));
    }

    #[test]
    fn collect_training_gff_errors_without_full_table() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("t.gff3");
        assert!(collect_busco_training_gff(dir.path(), &out).is_err());
    }

    #[test]
    fn parse_version_reads_busco_and_augustus_banners() {
        assert_eq!(parse_version("BUSCO 5.7.1\n"), Some((5, 7, 1)));
        assert_eq!(
            parse_version("AUGUSTUS (3.3.3) is a gene prediction tool\nwritten by M. Stanke"),
            Some((3, 3, 3))
        );
        assert_eq!(parse_version("AUGUSTUS (3.5.0) is a gene"), Some((3, 5, 0)));
        assert_eq!(parse_version("busco 5.4"), Some((5, 4, 0)));
        assert_eq!(parse_version("no version here"), None);
    }

    #[test]
    fn new_busco_with_old_augustus_is_rejected() {
        let why = version_incompatibility(Some((5, 7, 1)), Some((3, 3, 3))).unwrap();
        assert!(why.contains("5.7.1") && why.contains("3.3.3") && why.contains("3.4.0"));
    }

    #[test]
    fn compatible_or_unknown_versions_pass() {
        assert!(version_incompatibility(Some((5, 7, 1)), Some((3, 4, 0))).is_none());
        assert!(version_incompatibility(Some((5, 6, 1)), Some((3, 3, 3))).is_none());
        assert!(version_incompatibility(None, Some((3, 3, 3))).is_none());
        assert!(version_incompatibility(Some((5, 7, 1)), None).is_none());
    }
}
