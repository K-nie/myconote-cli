/// Database download and management
///
/// Handles downloading and indexing all databases needed for annotation:
///   - Swiss-Prot (UniProt reviewed) — for MMseqs2 homology
///   - Pfam-A.hmm                    — for hmmscan domain search
///
/// All databases are stored in ~/.myconote/dbs/ by default.
/// Run once with: myconote annotate --download-dbs
use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Download configuration
// ─────────────────────────────────────────────────────────────────────────────

const SWISSPROT_URL: &str =
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";

const PFAM_URL: &str = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz";

// ─────────────────────────────────────────────────────────────────────────────
// Status report
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug)]
pub struct DbStatus {
    pub swissprot_ready: bool,
    pub pfam_ready: bool,
    pub db_dir: PathBuf,
}

impl DbStatus {
    pub fn print(&self) {
        println!("  Database directory: {}", self.db_dir.display());
        println!(
            "  Swiss-Prot (MMseqs2) : {}",
            status_str(self.swissprot_ready)
        );
        println!("  Pfam-A (hmmscan)     : {}", status_str(self.pfam_ready));
    }
}

fn status_str(ready: bool) -> &'static str {
    if ready {
        "✓ ready"
    } else {
        "✗ not found — run: myconote annotate --download-dbs"
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Check status
// ─────────────────────────────────────────────────────────────────────────────

pub fn check_status(db_dir: &Path) -> DbStatus {
    let sp_db = db_dir.join("swissprot").join("swissprot");
    let pfam = db_dir.join("pfam").join("Pfam-A.hmm");

    DbStatus {
        swissprot_ready: sp_db.exists(),
        pfam_ready: pfam.exists(),
        db_dir: db_dir.to_path_buf(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Download all databases
// ─────────────────────────────────────────────────────────────────────────────

/// Download and index all databases.  Safe to re-run (skips already-present files).
pub fn download_all(db_dir: &Path, threads: usize) -> Result<()> {
    std::fs::create_dir_all(db_dir).map_err(MycoNoteError::Io)?;

    println!("── Downloading annotation databases ─────────────────────────");
    println!("  Target directory: {}", db_dir.display());

    download_swissprot(db_dir, threads)?;
    download_pfam(db_dir)?;

    println!("── All databases ready ───────────────────────────────────────");
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Swiss-Prot download + MMseqs2 indexing
// ─────────────────────────────────────────────────────────────────────────────

fn download_swissprot(db_dir: &Path, threads: usize) -> Result<()> {
    let sp_dir = db_dir.join("swissprot");
    let sp_db = sp_dir.join("swissprot");

    if sp_db.exists() {
        println!("  Swiss-Prot MMseqs2 database already present — skipping download.");
        return Ok(());
    }

    std::fs::create_dir_all(&sp_dir).map_err(MycoNoteError::Io)?;

    let gz_path = sp_dir.join("uniprot_sprot.fasta.gz");
    let fa_path = sp_dir.join("uniprot_sprot.fasta");

    // Download
    println!("  Downloading Swiss-Prot FASTA from UniProt…");
    curl_download(SWISSPROT_URL, &gz_path)?;

    // Decompress
    println!("  Decompressing Swiss-Prot FASTA…");
    gunzip(&gz_path, &fa_path)?;
    let _ = std::fs::remove_file(&gz_path);

    // Build MMseqs2 database
    println!("  Building MMseqs2 database…");
    mmseqs_createdb(&fa_path, &sp_db, threads)?;
    let _ = std::fs::remove_file(&fa_path);

    println!("  Swiss-Prot ready → {}", sp_db.display());
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Pfam download
// ─────────────────────────────────────────────────────────────────────────────

fn download_pfam(db_dir: &Path) -> Result<()> {
    let pfam_dir = db_dir.join("pfam");
    let pfam_hmm = pfam_dir.join("Pfam-A.hmm");

    if pfam_hmm.exists() {
        println!("  Pfam-A.hmm already present — skipping download.");
        return Ok(());
    }

    std::fs::create_dir_all(&pfam_dir).map_err(MycoNoteError::Io)?;

    let gz_path = pfam_dir.join("Pfam-A.hmm.gz");

    println!("  Downloading Pfam-A HMM from EBI (this may take a few minutes)…");
    curl_download(PFAM_URL, &gz_path)?;

    println!("  Decompressing Pfam-A.hmm…");
    gunzip(&gz_path, &pfam_hmm)?;
    let _ = std::fs::remove_file(&gz_path);

    // Press the HMM database for hmmscan
    println!("  Pressing Pfam-A.hmm for hmmscan…");
    press_hmm(&pfam_hmm)?;

    println!("  Pfam-A ready → {}", pfam_hmm.display());
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Shell helpers
// ─────────────────────────────────────────────────────────────────────────────

fn curl_download(url: &str, dest: &Path) -> Result<()> {
    // Prefer wget, fall back to curl
    let downloader = if which::which("wget").is_ok() {
        let status = Command::new("wget")
            .args(["-q", "-O", dest.to_str().unwrap_or(""), url])
            .status()
            .map_err(MycoNoteError::Io)?;
        status
    } else {
        let status = Command::new("curl")
            .args(["-L", "-s", "-o", dest.to_str().unwrap_or(""), url])
            .status()
            .map_err(MycoNoteError::Io)?;
        status
    };

    if !downloader.success() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Failed to download {}. Check your internet connection.",
            url
        )));
    }
    Ok(())
}

fn gunzip(gz: &Path, dest: &Path) -> Result<()> {
    let status = Command::new("gunzip")
        .args(["-c", gz.to_str().unwrap_or("")])
        .stdout(std::fs::File::create(dest).map_err(MycoNoteError::Io)?)
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat("gunzip failed.".to_string()));
    }
    Ok(())
}

fn mmseqs_createdb(fasta: &Path, db_out: &Path, threads: usize) -> Result<()> {
    let mmseqs = which::which("mmseqs").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "mmseqs not found in PATH. Install with: conda install -c bioconda mmseqs2".to_string(),
        )
    })?;

    let tmp_dir = tempfile::TempDir::new().map_err(MycoNoteError::Io)?;

    // createdb
    let s1 = Command::new(&mmseqs)
        .args([
            "createdb",
            fasta.to_str().unwrap_or(""),
            db_out.to_str().unwrap_or(""),
        ])
        .status()
        .map_err(MycoNoteError::Io)?;
    if !s1.success() {
        return Err(MycoNoteError::InvalidFormat(
            "mmseqs createdb failed.".to_string(),
        ));
    }

    // createindex for faster search
    let s2 = Command::new(&mmseqs)
        .args([
            "createindex",
            db_out.to_str().unwrap_or(""),
            tmp_dir.path().to_str().unwrap_or("/tmp"),
            "--threads",
            &threads.to_string(),
        ])
        .status()
        .map_err(MycoNoteError::Io)?;
    if !s2.success() {
        eprintln!("  ⚠  mmseqs createindex failed (non-fatal — search will still work)");
    }

    Ok(())
}

fn press_hmm(hmm: &Path) -> Result<()> {
    let hmmpress = which::which("hmmpress").map_err(|_| {
        MycoNoteError::UnsupportedFormat(
            "hmmpress not found. Install with: conda install -c bioconda hmmer".to_string(),
        )
    })?;

    let status = Command::new(hmmpress)
        .args(["-f", hmm.to_str().unwrap_or("")])
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat("hmmpress failed.".to_string()));
    }
    Ok(())
}
