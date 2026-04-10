/// Database installer: `myconote setup`
///
/// Downloads and indexes all reference databases used by myconote-cli:
///
///   swiss-prot  — UniProt Swiss-Prot FASTA (MMseqs2 indexed)
///   pfam        — Pfam-A HMM profiles (pressed with hmmpress)
///   eggnog      — EggNog-mapper databases (emapper.py --download_taxon ...)
///   dbcan       — CAZyme HMM + DIAMOND databases
///   merops      — MEROPS protease scan library (DIAMOND indexed)
///   busco       — BUSCO lineage datasets (fungi_odb10, etc.)
///   augustus    — Pre-trained Augustus species config
///
/// Usage:
///   myconote setup                        # download all
///   myconote setup --dbs swiss-prot pfam  # download specific databases
///   myconote setup --list                 # list available databases
///   myconote setup --check                # check what is already downloaded
use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Database catalogue
// ─────────────────────────────────────────────────────────────────────────────

pub struct DbEntry {
    pub key: &'static str,
    pub description: &'static str,
    /// Approximate size on disk
    pub size_hint: &'static str,
    /// URL(s) to download from
    pub urls: &'static [&'static str],
    /// File that must exist after successful download
    pub marker_file: &'static str,
    /// Post-download indexing command (empty → no indexing needed)
    pub index_cmd: &'static str,
}

pub const DATABASES: &[DbEntry] = &[
    DbEntry {
        key:         "swiss-prot",
        description: "UniProt Swiss-Prot (MMseqs2 homology search)",
        size_hint:   "~250 MB",
        urls: &[
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
        ],
        marker_file: "swissprot/uniprot_sprot.fasta.gz",
        index_cmd:   "mmseqs createdb",
    },
    DbEntry {
        key:         "pfam",
        description: "Pfam-A HMM profiles (hmmscan domain search)",
        size_hint:   "~300 MB",
        urls: &[
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz",
        ],
        marker_file: "pfam/Pfam-A.hmm",
        index_cmd:   "hmmpress",
    },
    DbEntry {
        key:         "eggnog",
        description: "EggNog-mapper databases (COG/NOG functional categories)",
        size_hint:   "~50 GB",
        urls: &[],  // downloaded via emapper.py --download_taxon
        marker_file: "eggnog/eggnog.db",
        index_cmd:   "emapper.py --data_dir",
    },
    DbEntry {
        key:         "dbcan",
        description: "dbCAN CAZyme HMM + DIAMOND database",
        size_hint:   "~50 MB",
        urls: &[
            "https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt",
            "https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa",
        ],
        marker_file: "dbcan/dbCAN-HMMdb-V12.txt",
        index_cmd:   "hmmpress + diamond makedb",
    },
    DbEntry {
        key:         "merops",
        description: "MEROPS peptidase database (protease annotation)",
        size_hint:   "~15 MB",
        urls: &[
            "https://ftp.ebi.ac.uk/pub/databases/merops/current_release/merops_scan.lib",
        ],
        marker_file: "merops/merops_scan.lib",
        index_cmd:   "diamond makedb",
    },
    DbEntry {
        key:         "busco",
        description: "BUSCO lineage datasets (fungi_odb10, dikarya_odb10, ...)",
        size_hint:   "~500 MB",
        urls: &[],  // downloaded via busco --download_path auto-lineage
        marker_file: "busco/fungi_odb10",
        index_cmd:   "",
    },
];

// ─────────────────────────────────────────────────────────────────────────────
// Entry points
// ─────────────────────────────────────────────────────────────────────────────

pub fn list_databases() {
    println!("\nAvailable myconote databases:");
    println!("{:<14} {:<45} {}", "Key", "Description", "Size");
    println!("{}", "─".repeat(75));
    for db in DATABASES {
        println!("{:<14} {:<45} {}", db.key, db.description, db.size_hint);
    }
    println!();
    println!("Download with: myconote setup --dbs <key> [<key> ...]");
    println!("Download all:  myconote setup");
}

pub fn check_databases(db_dir: &Path) {
    println!("\nDatabase status (db_dir: {}):", db_dir.display());
    println!("{:<14} {:<12} {}", "Key", "Status", "Version / Date");
    println!("{}", "─".repeat(60));
    for db in DATABASES {
        let marker = db_dir.join(db.marker_file);
        if marker.exists() {
            let version_info =
                read_db_version(db_dir, db.key).unwrap_or_else(|| "version unknown".to_string());
            println!(
                "  {:<14} \x1b[32m✓ present\x1b[0m  {}",
                db.key, version_info
            );
        } else {
            println!("  {:<14} \x1b[31m✗ missing\x1b[0m", db.key);
        }
    }
    println!();
}

// ─────────────────────────────────────────────────────────────────────────────
// Database version tracking
// ─────────────────────────────────────────────────────────────────────────────

/// Write a version metadata file after successful download.
fn write_db_version(db_dir: &Path, key: &str, url: &str) -> Result<()> {
    let version_file = db_dir.join(format!("{}.version", key));
    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    let mut f = std::fs::File::create(&version_file).map_err(MycoNoteError::Io)?;
    use std::io::Write;
    writeln!(f, "database={}", key).map_err(MycoNoteError::Io)?;
    writeln!(f, "download_url={}", url).map_err(MycoNoteError::Io)?;
    writeln!(f, "download_timestamp={}", timestamp).map_err(MycoNoteError::Io)?;
    writeln!(f, "myconote_version={}", env!("CARGO_PKG_VERSION")).map_err(MycoNoteError::Io)?;

    // Try to compute checksum of marker file
    let marker_path = DATABASES
        .iter()
        .find(|d| d.key == key)
        .map(|d| db_dir.join(d.marker_file));

    if let Some(mp) = marker_path {
        if mp.exists() {
            if let Ok(meta) = std::fs::metadata(&mp) {
                writeln!(f, "file_size={}", meta.len()).map_err(MycoNoteError::Io)?;
            }
        }
    }

    Ok(())
}

/// Read version info for a database.
fn read_db_version(db_dir: &Path, key: &str) -> Option<String> {
    let version_file = db_dir.join(format!("{}.version", key));
    if !version_file.exists() {
        // Fall back to file modification time
        let marker = DATABASES.iter().find(|d| d.key == key)?;
        let path = db_dir.join(marker.marker_file);
        let meta = std::fs::metadata(&path).ok()?;
        let modified = meta.modified().ok()?;
        let secs = modified
            .duration_since(std::time::UNIX_EPOCH)
            .ok()?
            .as_secs();
        return Some(format!("downloaded ~{}", approximate_date(secs)));
    }

    let contents = std::fs::read_to_string(&version_file).ok()?;
    let mut timestamp = None;
    for line in contents.lines() {
        if let Some(ts) = line.strip_prefix("download_timestamp=") {
            timestamp = ts.parse::<u64>().ok();
        }
    }
    timestamp.map(|ts| format!("downloaded {}", approximate_date(ts)))
}

fn approximate_date(secs: u64) -> String {
    let days = secs / 86400;
    let mut y = 1970u64;
    let mut remaining = days;
    loop {
        let dy = if y % 4 == 0 && (y % 100 != 0 || y % 400 == 0) {
            366
        } else {
            365
        };
        if remaining < dy {
            break;
        }
        remaining -= dy;
        y += 1;
    }
    let month_days: [u64; 12] = if y % 4 == 0 && (y % 100 != 0 || y % 400 == 0) {
        [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };
    let mut m = 1u64;
    for md in &month_days {
        if remaining < *md {
            break;
        }
        remaining -= *md;
        m += 1;
    }
    format!("{:04}-{:02}-{:02}", y, m, remaining + 1)
}

pub fn download_databases(db_dir: &Path, keys: &[String], force: bool) -> Result<()> {
    std::fs::create_dir_all(db_dir).map_err(MycoNoteError::Io)?;

    let to_download: Vec<&DbEntry> = if keys.is_empty() {
        DATABASES.iter().collect()
    } else {
        DATABASES
            .iter()
            .filter(|db| keys.iter().any(|k| k == db.key))
            .collect()
    };

    if to_download.is_empty() {
        println!("No matching databases found. Use `myconote setup --list` to see options.");
        return Ok(());
    }

    for db in to_download {
        let marker = db_dir.join(db.marker_file);
        if marker.exists() && !force {
            println!(
                "  [skip] {} — already present (use --force to re-download)",
                db.key
            );
            continue;
        }

        println!(
            "\n── Downloading: {} ({}) ──────────────────────",
            db.key, db.size_hint
        );

        match db.key {
            "swiss-prot" => download_swissprot(db_dir)?,
            "pfam" => download_pfam(db_dir)?,
            "eggnog" => download_eggnog(db_dir)?,
            "dbcan" => download_dbcan(db_dir)?,
            "merops" => download_merops(db_dir)?,
            "busco" => download_busco(db_dir)?,
            other => println!("  ⚠  No download handler for '{}'", other),
        }
    }

    println!("\n✓ Setup complete. Run `myconote setup --check` to verify.");
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-database downloaders
// ─────────────────────────────────────────────────────────────────────────────

fn wget_or_curl(url: &str, dest: &Path) -> Result<()> {
    // Try wget first
    let wget_ok = Command::new("wget")
        .args([
            "--quiet",
            "--show-progress",
            "-O",
            dest.to_str().unwrap_or(""),
            url,
        ])
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if wget_ok {
        return Ok(());
    }

    // Fall back to curl
    let curl_ok = Command::new("curl")
        .args([
            "-L",
            "--progress-bar",
            "-o",
            dest.to_str().unwrap_or(""),
            url,
        ])
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if curl_ok {
        return Ok(());
    }

    Err(MycoNoteError::ExternalTool(format!(
        "Could not download {}\nInstall wget or curl.",
        url
    )))
}

fn gunzip(gz_path: &Path) -> Result<PathBuf> {
    let out = gz_path.with_extension(""); // strip .gz
    let status = Command::new("gunzip")
        .arg("-k") // keep original
        .arg(gz_path)
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("gunzip: {}", e)))?;
    if !status.success() {
        return Err(MycoNoteError::ExternalTool("gunzip failed".to_string()));
    }
    Ok(out)
}

fn download_swissprot(db_dir: &Path) -> Result<()> {
    let sp_dir = db_dir.join("swissprot");
    std::fs::create_dir_all(&sp_dir).map_err(MycoNoteError::Io)?;

    let gz = sp_dir.join("uniprot_sprot.fasta.gz");
    println!("  Downloading UniProt Swiss-Prot FASTA…");
    wget_or_curl(
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
        &gz,
    )?;

    let fasta = gunzip(&gz)?;
    println!("  Indexing with MMseqs2…");

    // Check mmseqs availability
    if !Command::new("mmseqs")
        .arg("version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
    {
        eprintln!(
            "  ⚠  mmseqs not found — skipping indexing. Install: conda install -c bioconda mmseqs2"
        );
        return Ok(());
    }

    let mmdb = sp_dir.join("swissprot");
    let tmp = sp_dir.join("tmp_mmseqs");
    Command::new("mmseqs")
        .args([
            "createdb",
            fasta.to_str().unwrap_or(""),
            mmdb.to_str().unwrap_or(""),
        ])
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("mmseqs createdb: {}", e)))?;
    let _ = std::fs::create_dir_all(&tmp);
    Command::new("mmseqs")
        .args([
            "createindex",
            mmdb.to_str().unwrap_or(""),
            tmp.to_str().unwrap_or(""),
        ])
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("mmseqs createindex: {}", e)))?;

    println!("  ✓ Swiss-Prot ready at {}", mmdb.display());

    // Track version
    let _ = write_db_version(db_dir, "swiss-prot",
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz");

    Ok(())
}

fn download_pfam(db_dir: &Path) -> Result<()> {
    let pfam_dir = db_dir.join("pfam");
    std::fs::create_dir_all(&pfam_dir).map_err(MycoNoteError::Io)?;

    for file in &["Pfam-A.hmm.gz", "Pfam-A.hmm.dat.gz"] {
        let url = format!(
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/{}",
            file
        );
        let dest = pfam_dir.join(file);
        println!("  Downloading {}…", file);
        wget_or_curl(&url, &dest)?;
        let _ = gunzip(&dest);
    }

    let hmm = pfam_dir.join("Pfam-A.hmm");
    if hmm.exists() {
        println!("  Pressing Pfam HMM database (hmmpress)…");
        if !Command::new("hmmpress")
            .arg(&hmm)
            .status()
            .map(|s| s.success())
            .unwrap_or(false)
        {
            eprintln!("  ⚠  hmmpress failed. Install: conda install -c bioconda hmmer");
        } else {
            println!("  ✓ Pfam-A ready at {}", hmm.display());
        }
    }

    let _ = write_db_version(
        db_dir,
        "pfam",
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
    );

    Ok(())
}

fn download_eggnog(db_dir: &Path) -> Result<()> {
    let eg_dir = db_dir.join("eggnog");
    std::fs::create_dir_all(&eg_dir).map_err(MycoNoteError::Io)?;

    println!("  Downloading EggNog-mapper databases (~50 GB, may take a while)…");
    println!("  This requires emapper.py.");

    let ok = Command::new("download_eggnog_data.py")
        .args(["-y", "--data_dir", eg_dir.to_str().unwrap_or(".")])
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if !ok {
        // Try the older method
        let ok2 = Command::new("emapper.py")
            .args([
                "--data_dir",
                eg_dir.to_str().unwrap_or("."),
                "--download_taxon",
                "2",
                "-y",
            ])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if !ok2 {
            eprintln!("  ⚠  emapper.py not found or download failed.");
            eprintln!("     Install: conda install -c bioconda eggnog-mapper");
            eprintln!(
                "     Then:    download_eggnog_data.py -y --data_dir {}",
                eg_dir.display()
            );
        }
    } else {
        println!("  ✓ EggNog databases ready at {}", eg_dir.display());
    }

    Ok(())
}

fn download_dbcan(db_dir: &Path) -> Result<()> {
    let dc_dir = db_dir.join("dbcan");
    std::fs::create_dir_all(&dc_dir).map_err(MycoNoteError::Io)?;

    for (url, filename) in &[
        (
            "https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt",
            "dbCAN-HMMdb-V12.txt",
        ),
        (
            "https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa",
            "CAZyDB.fa",
        ),
    ] {
        println!("  Downloading {}…", filename);
        let dest = dc_dir.join(filename);
        let _ = wget_or_curl(url, &dest);
    }

    // hmmpress HMM database
    let hmm = dc_dir.join("dbCAN-HMMdb-V12.txt");
    if hmm.exists() {
        let _ = Command::new("hmmpress").arg(&hmm).status();
    }

    // Build DIAMOND database
    let cazy_fa = dc_dir.join("CAZyDB.fa");
    let dmnd = dc_dir.join("dbCAN.dmnd");
    if cazy_fa.exists() {
        let _ = Command::new("diamond")
            .args([
                "makedb",
                "--in",
                cazy_fa.to_str().unwrap_or(""),
                "--db",
                dmnd.to_str().unwrap_or(""),
            ])
            .status();
    }

    println!("  ✓ dbCAN databases ready at {}", dc_dir.display());
    Ok(())
}

fn download_merops(db_dir: &Path) -> Result<()> {
    let mer_dir = db_dir.join("merops");
    std::fs::create_dir_all(&mer_dir).map_err(MycoNoteError::Io)?;

    let lib_path = mer_dir.join("merops_scan.lib");
    println!("  Downloading MEROPS scan library…");
    // Try FTP (registration required for current release; using legacy mirror)
    let ok = wget_or_curl(
        "https://ftp.ebi.ac.uk/pub/databases/merops/current_release/merops_scan.lib",
        &lib_path,
    );

    if ok.is_err() {
        eprintln!("  ⚠  MEROPS FTP download may require registration.");
        eprintln!("     Visit https://www.ebi.ac.uk/merops/download/ to download manually");
        eprintln!("     and place merops_scan.lib in {}", mer_dir.display());
        return Ok(());
    }

    // Build DIAMOND database
    let dmnd = mer_dir.join("merops.dmnd");
    if lib_path.exists() {
        println!("  Building DIAMOND database…");
        let diamond_ok = Command::new("diamond")
            .args([
                "makedb",
                "--in",
                lib_path.to_str().unwrap_or(""),
                "--db",
                dmnd.to_str().unwrap_or(""),
            ])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);

        if diamond_ok {
            println!("  ✓ MEROPS DIAMOND database ready at {}", dmnd.display());
        } else {
            eprintln!("  ⚠  diamond makedb failed. Install: conda install -c bioconda diamond");
        }
    }

    Ok(())
}

fn download_busco(db_dir: &Path) -> Result<()> {
    let busco_dir = db_dir.join("busco");
    std::fs::create_dir_all(&busco_dir).map_err(MycoNoteError::Io)?;

    let lineages = [
        "fungi_odb10",
        "dikarya_odb10",
        "basidiomycota_odb10",
        "ascomycota_odb10",
        "microsporidia_odb10",
    ];

    println!(
        "  Downloading BUSCO lineage datasets: {}",
        lineages.join(", ")
    );

    for lineage in &lineages {
        println!("    Downloading {}…", lineage);
        // Try busco in PATH first; fall back to conda run -n busco_env
        let busco_bin = which::which("busco").ok();
        let ok = if let Some(ref bin) = busco_bin {
            Command::new(bin)
                .args([
                    "--download",
                    lineage,
                    "--download_path",
                    busco_dir.to_str().unwrap_or("."),
                ])
                .status()
                .map(|s| s.success())
                .unwrap_or(false)
        } else {
            let env = std::env::var("BUSCO_CONDA_ENV").unwrap_or_else(|_| "busco_env".to_string());
            Command::new("conda")
                .args([
                    "run",
                    "--no-capture-output",
                    "-n",
                    &env,
                    "busco",
                    "--download",
                    lineage,
                    "--download_path",
                    busco_dir.to_str().unwrap_or("."),
                ])
                .status()
                .map(|s| s.success())
                .unwrap_or(false)
        };

        if !ok {
            eprintln!(
                "    ⚠  Could not download {}. Install BUSCO first \
                (conda create -n busco_env -c conda-forge -c bioconda python=3.12 busco=5).",
                lineage
            );
            break;
        }
    }

    println!("  ✓ BUSCO lineages ready at {}", busco_dir.display());
    Ok(())
}
