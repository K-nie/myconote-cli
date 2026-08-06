/// Secretome prediction: signal peptide + TM-helix filtering
///
/// Identifies secreted proteins — those with a classical N-terminal
/// signal peptide but no transmembrane helices (i.e. genuinely secreted,
/// not membrane-anchored).  This "secretome" list is especially important
/// for fungal pathogenicity studies (effectors, cell-wall-degrading enzymes).
///
/// Pipeline:
///   1. **DeepSig** (free, pip) or SignalP 6/5/4 (academic licence) —
///      predict signal peptide cleavage site
///   2. **DeepTMHMM** via biolib (free, pip) or TMHMM 2 (academic licence) —
///      predict transmembrane helices
///   3. Retain proteins with signal peptide AND zero TM helices
///
/// If neither signal-peptide tool is available, Phobius (predicts both
/// signal peptide and TM topology in one run) is used as a fallback.
///
/// All tools are optional — we always gracefully degrade:
///   - Signal tool only  → flag "has_signal_peptide" but no TM filter
///   - TM tool only      → flag "transmembrane" proteins
///   - Neither           → skip secretome prediction entirely
use crate::utils::error::{MycoNoteError, Result};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Data types
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct SignalPResult {
    pub protein_id: String,
    /// Predicted signal peptide cleavage site (1-based)
    pub cleavage_site: Option<usize>,
    /// SignalP probability score
    pub sp_score: f64,
    /// Tool-level prediction ("SP", "LIPO", "TAT", "OTHER", "NO_SP")
    pub prediction: String,
    pub has_signal: bool,
}

#[derive(Debug, Clone)]
pub struct TmhmmResult {
    pub protein_id: String,
    /// Number of predicted transmembrane helices
    pub tm_count: usize,
    pub is_inside: bool, // N-terminus orientation
}

// ─────────────────────────────────────────────────────────────────────────────
// Run SignalP
// ─────────────────────────────────────────────────────────────────────────────

/// Detect the best available signal-peptide tool.
/// Priority: deepsig (free) > signalp6 > signalp (licensed)
pub fn signalp_available() -> Option<String> {
    for name in &["deepsig", "signalp6", "signalp"] {
        if Command::new("which")
            .arg(name)
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
        {
            return Some(name.to_string());
        }
    }
    None
}

/// Select the SignalP 6 `--mode` value that matches what's installed.
///
/// Resolves the `signalp6` binary through `which`, walks up to its
/// conda env prefix, and checks whether the distilled model weight
/// file (`site-packages/signalp/model_weights/distilled_model_signalp6.pt`)
/// exists. Fast mode requires that file; every other mode needs only
/// the base weights that ship with every DTU download variant.
///
/// Returns `"fast"` when the distilled weights are present,
/// `"slow_sequential"` otherwise.
fn signalp6_mode_for_install() -> &'static str {
    if let Ok(out) = Command::new("which").arg("signalp6").output() {
        let bin_path = String::from_utf8_lossy(&out.stdout).trim().to_string();
        if !bin_path.is_empty() {
            let mut p = PathBuf::from(&bin_path);
            // Resolve symlink so `~/.conda/envs/myconote/bin/signalp6 ->
            // ~/.conda/envs/myconote_signalp6/bin/signalp6` points at the
            // real env, not the umbrella myconote env.
            if let Ok(canon) = std::fs::canonicalize(&p) {
                p = canon;
            }
            // signalp6 lives at <env>/bin/signalp6, so climb two dirs
            // to reach <env> and probe every python3.* site-packages.
            if let Some(env_root) = p.parent().and_then(|bin| bin.parent()) {
                if let Ok(entries) = std::fs::read_dir(env_root.join("lib")) {
                    for entry in entries.flatten() {
                        let candidate = entry.path().join(
                            "site-packages/signalp/model_weights/distilled_model_signalp6.pt",
                        );
                        if candidate.exists() {
                            return "fast";
                        }
                    }
                }
            }
        }
    }
    "slow_sequential"
}

pub fn run_signalp(
    protein_fasta: &Path,
    out_dir: &Path,
    organism: &str, // "euk" for fungi / "gram+" / "gram-"
) -> Result<PathBuf> {
    let tool = signalp_available().ok_or_else(|| {
        MycoNoteError::ExternalTool(
            "No signal-peptide tool found. Run: myconote-cli install annotate".to_string(),
        )
    })?;

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    if tool == "deepsig" {
        // DeepSig CLI: deepsig -f <fasta> -o <output_file> -k euk|gram+|gram-
        // Output: TSV  ID \t Prediction \t Score
        let out_file = out_dir.join("deepsig_results.tsv");
        let status = Command::new("deepsig")
            .arg("-f")
            .arg(protein_fasta)
            .arg("-o")
            .arg(&out_file)
            .arg("-k")
            .arg(organism)
            .status()
            .map_err(|e| MycoNoteError::ExternalTool(format!("deepsig: {}", e)))?;

        if !status.success() {
            return Err(MycoNoteError::ExternalTool(
                "deepsig exited with non-zero status".to_string(),
            ));
        }
        return Ok(out_file);
    }

    // SignalP (licensed versions 4/5/6).
    //
    // SignalP 6's default `--mode fast` requires a distilled model
    // weights file (`distilled_model_signalp6.pt`) that DTU ships in
    // the "fast + distilled" tarball, but NOT in the plain "fast"
    // tarball that most first-time users download. Without it, the
    // tool crashes with `FileNotFoundError: Fast mode requires model
    // to be installed at .../distilled_model_signalp6.pt`. Fixes F10
    // by preflighting the model file: when the distilled weights are
    // present we pass `--mode fast`; otherwise we fall back to
    // `--mode slow_sequential`, which works with the base weights
    // that every SignalP 6 variant ships.
    let mut cmd = Command::new(&tool);
    cmd.arg("--fastafile")
        .arg(protein_fasta)
        .arg("--organism")
        .arg(organism)
        .arg("--output_dir")
        .arg(out_dir)
        .arg("--format")
        .arg("txt");
    if tool == "signalp6" {
        let mode = signalp6_mode_for_install();
        eprintln!(
            "  SignalP 6: using --mode {} (distilled-model detection)",
            mode
        );
        cmd.arg("--mode").arg(mode);
    }
    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("signalp: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "SignalP exited with non-zero status".to_string(),
        ));
    }

    // SignalP writes output_dir/prediction_results.txt (v6) or output.txt (v5/v4)
    for candidate in &[
        out_dir.join("prediction_results.txt"),
        out_dir.join("output.gff3"),
        out_dir.join("output.txt"),
    ] {
        if candidate.exists() {
            return Ok(candidate.clone());
        }
    }

    // Fallback: any *.txt in out_dir
    if let Ok(entries) = std::fs::read_dir(out_dir) {
        for entry in entries.flatten() {
            if entry.path().extension().and_then(|e| e.to_str()) == Some("txt") {
                return Ok(entry.path());
            }
        }
    }

    Err(MycoNoteError::ExternalTool(
        "Could not find SignalP output file".to_string(),
    ))
}

/// Parse signal-peptide output.
/// Handles DeepSig TSV, SignalP v4/v5/v6 formats (auto-detected).
pub fn parse_signalp_output(path: &Path) -> Result<HashMap<String, SignalPResult>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }

        let protein_id = fields[0].to_string();

        // SignalP 6: ID \t Prediction \t SP(Sec/SPI) \t LIPO(Sec/SPII) \t ... \t CS Position
        // SignalP 5: ID \t Prediction \t SP_prob \t ... \t CS
        // Detect format by column count
        let (prediction, sp_score, has_signal, cleavage_site) = if fields.len() >= 6 {
            let pred = fields[1].to_string();
            let score: f64 = fields[2].parse().unwrap_or(0.0);
            let has = pred.contains("SP") || pred == "SP" || pred.starts_with("Signal");
            // CS position: last field or specific column
            let cs_field = fields.last().unwrap_or(&"");
            let cs = cs_field
                .trim_start_matches("CS pos: ")
                .split('-')
                .next()
                .and_then(|s| s.trim().parse::<usize>().ok());
            (pred, score, has, cs)
        } else {
            let pred = fields.get(1).unwrap_or(&"OTHER").to_string();
            let has = pred == "Y" || pred.contains("SP");
            (pred, 0.0, has, None)
        };

        map.insert(
            protein_id.clone(),
            SignalPResult {
                protein_id,
                cleavage_site,
                sp_score,
                prediction,
                has_signal,
            },
        );
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// Run TMHMM
// ─────────────────────────────────────────────────────────────────────────────

/// Detect the best available TM-helix predictor.
/// Priority: biolib/DeepTMHMM (free) > deeptmhmm > tmhmm (licensed)
pub fn tmhmm_available() -> Option<String> {
    for name in &["biolib", "deeptmhmm", "tmhmm", "tmhmm2.0c"] {
        if Command::new("which")
            .arg(name)
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
        {
            return Some(name.to_string());
        }
    }
    None
}

pub fn run_tmhmm(protein_fasta: &Path, out_dir: &Path) -> Result<PathBuf> {
    let tool = tmhmm_available().ok_or_else(|| {
        MycoNoteError::ExternalTool(
            "No TM-helix tool found. Run: myconote-cli install annotate".to_string(),
        )
    })?;

    std::fs::create_dir_all(out_dir).map_err(MycoNoteError::Io)?;

    if tool == "biolib" {
        // DeepTMHMM via pybiolib:
        //   biolib run DTU/DeepTMHMM --fasta <file>
        // Writes results to ./biolib_results/ in the working directory,
        // so we run it from out_dir and return the topology file.
        let status = Command::new("biolib")
            .args(["run", "DTU/DeepTMHMM", "--fasta"])
            .arg(
                protein_fasta
                    .canonicalize()
                    .unwrap_or_else(|_| protein_fasta.to_path_buf()),
            )
            .current_dir(out_dir)
            .status()
            .map_err(|e| MycoNoteError::ExternalTool(format!("biolib/DeepTMHMM: {}", e)))?;

        if !status.success() {
            return Err(MycoNoteError::ExternalTool(
                "DeepTMHMM (biolib) exited with non-zero status".to_string(),
            ));
        }

        // DeepTMHMM writes: biolib_results/TMRs.gff3 and
        //                   biolib_results/predicted_topologies.3line
        // We convert it to the same short-format TSV that parse_tmhmm_output expects.
        let gff = out_dir.join("biolib_results").join("TMRs.gff3");
        let out_file = out_dir.join("tmhmm.txt");
        convert_deeptmhmm_to_short_fmt(&gff, &out_file)?;
        return Ok(out_file);
    }

    // TMHMM (licensed) short format
    let out_file = out_dir.join("tmhmm.txt");
    let out_f = std::fs::File::create(&out_file).map_err(MycoNoteError::Io)?;

    let status = Command::new(&tool)
        .arg(protein_fasta)
        .arg("--short")
        .stdout(std::process::Stdio::from(out_f))
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("tmhmm: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "TMHMM exited with non-zero status".to_string(),
        ));
    }

    Ok(out_file)
}

/// Convert DeepTMHMM GFF3 output → TMHMM short-format TSV so that the
/// existing `parse_tmhmm_output` function works unchanged.
///
/// DeepTMHMM GFF3 format (TMRs.gff3):
///   ##sequence-region  <protein_id>  1  <len>
///   <protein_id>  DeepTMHMM  TMhelix  <start>  <end>  .  .  .  .
///
/// Target short format (one line per protein):
///   <id>  len=N  ExpAA=0  First60=0  PredHel=K  Topology=...
fn convert_deeptmhmm_to_short_fmt(gff: &Path, out: &Path) -> Result<()> {
    use std::collections::HashMap;
    use std::io::BufRead;

    let reader = BufReader::new(std::fs::File::open(gff).map_err(MycoNoteError::Io)?);

    // protein_id → TM helix count
    let mut tm_counts: HashMap<String, usize> = HashMap::new();
    // protein_id → sequence length (from sequence-region directives)
    let mut lengths: HashMap<String, usize> = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let line = line.trim();

        if line.starts_with("##sequence-region") {
            // ##sequence-region  <id>  1  <len>
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                let id = parts[1].to_string();
                let len: usize = parts[3].parse().unwrap_or(0);
                lengths.entry(id.clone()).or_insert(len);
                tm_counts.entry(id).or_insert(0);
            }
        } else if !line.starts_with('#') && !line.is_empty() {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 && parts[2] == "TMhelix" {
                let id = parts[0].to_string();
                *tm_counts.entry(id).or_insert(0) += 1;
            }
        }
    }

    let mut outf = std::fs::File::create(out).map_err(MycoNoteError::Io)?;
    let mut ids: Vec<&String> = tm_counts.keys().collect();
    ids.sort();
    for id in ids {
        let k = tm_counts[id];
        let len = lengths.get(id).copied().unwrap_or(0);
        writeln!(
            outf,
            "{}\tlen={}\tExpAA=0\tFirst60=0\tPredHel={}\tTopology=",
            id, len, k
        )
        .map_err(MycoNoteError::Io)?;
    }

    Ok(())
}

pub fn parse_tmhmm_output(path: &Path) -> Result<HashMap<String, TmhmmResult>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(MycoNoteError::Io)?;
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }

        // Short format: ID \t len=NNN \t ExpAA=X \t First60=X \t PredHel=N \t Topology=...
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let protein_id = fields[0].to_string();
        let tm_count: usize = fields[4]
            .trim_start_matches("PredHel=")
            .parse()
            .unwrap_or(0);
        let is_inside = fields
            .get(5)
            .map(|t| t.trim_start_matches("Topology=").starts_with('i'))
            .unwrap_or(false);

        map.insert(
            protein_id.clone(),
            TmhmmResult {
                protein_id,
                tm_count,
                is_inside,
            },
        );
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// Secretome filter: signal peptide + no TM helices
// ─────────────────────────────────────────────────────────────────────────────

/// Return the set of protein IDs that are predicted secreted:
///   - Has a signal peptide (SignalP)
///   - Has 0 transmembrane helices (TMHMM)
pub fn filter_secretome(
    signalp: &HashMap<String, SignalPResult>,
    tmhmm: Option<&HashMap<String, TmhmmResult>>,
) -> HashSet<String> {
    signalp
        .iter()
        .filter(|(id, sp)| {
            if !sp.has_signal {
                return false;
            }
            if let Some(tm_map) = tmhmm {
                if let Some(tm) = tm_map.get(*id) {
                    return tm.tm_count == 0;
                }
            }
            true // No TMHMM data → pass through on signal alone
        })
        .map(|(id, _)| id.clone())
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Write secretome results
// ─────────────────────────────────────────────────────────────────────────────

pub fn write_secretome_table(
    secretome: &HashSet<String>,
    signalp: &HashMap<String, SignalPResult>,
    tmhmm: Option<&HashMap<String, TmhmmResult>>,
    output: &Path,
) -> Result<usize> {
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(
        out,
        "protein_id\thas_signal\tsp_score\tcleavage_site\ttm_helices\tsecreted"
    )
    .map_err(MycoNoteError::Io)?;

    let mut ids: Vec<&String> = signalp.keys().collect();
    ids.sort();

    for id in &ids {
        let sp = &signalp[*id];
        let tm = tmhmm.and_then(|m| m.get(*id));
        let tm_count = tm.map(|t| t.tm_count).unwrap_or(0);
        let cs = sp
            .cleavage_site
            .map(|n| n.to_string())
            .unwrap_or_else(|| "-".to_string());
        let is_secreted = secretome.contains(*id);

        writeln!(
            out,
            "{}\t{}\t{:.3}\t{}\t{}\t{}",
            id, sp.has_signal, sp.sp_score, cs, tm_count, is_secreted
        )
        .map_err(MycoNoteError::Io)?;
    }

    Ok(secretome.len())
}

pub fn phobius_available() -> bool {
    Command::new("which")
        .arg("phobius.pl")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}
