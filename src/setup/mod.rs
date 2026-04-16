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
    DbEntry {
        key:         "chat-corpus",
        description: "Q1 open-access paper corpus for explain citations",
        size_hint:   "~50 MB",
        urls: &[],  // built from corpus_manifest.toml at ~/.myconote/papers/
        marker_file: "papers/corpus_manifest.toml",
        index_cmd:   "",
    },
    DbEntry {
        key:         "ollama",
        description: "Ollama local LLM runtime + model for explain",
        size_hint:   "~4 GB (model dependent)",
        urls: &[],  // installed via install script or package manager
        marker_file: "ollama/installed.version",
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
        // Special handling for ollama — check the binary, not just a marker file
        if db.key == "ollama" {
            match get_ollama_version() {
                Some(ver) => {
                    let model_info = read_configured_model().unwrap_or_else(|| "llama3.1".to_string());
                    println!(
                        "  {:<14} \x1b[32m✓ present\x1b[0m  v{}, model: {}",
                        db.key, ver, model_info
                    );
                }
                None => {
                    println!("  {:<14} \x1b[31m✗ missing\x1b[0m  run: myconote setup --ollama", db.key);
                }
            }
            continue;
        }

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
            "chat-corpus" => download_chat_corpus(db_dir)?,
            "ollama" => setup_ollama(db_dir)?,
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

fn download_chat_corpus(db_dir: &Path) -> Result<()> {
    let papers_dir = db_dir.join("papers");
    std::fs::create_dir_all(&papers_dir).map_err(MycoNoteError::Io)?;

    // Also create the user-local papers dir at ~/.myconote/papers/
    let home_papers = home_myconote_dir().join("papers");
    std::fs::create_dir_all(&home_papers).map_err(MycoNoteError::Io)?;
    let _ = std::fs::create_dir_all(home_papers.join("local"));

    // Write a skeleton corpus_manifest.toml if one doesn't exist
    let manifest_path = papers_dir.join("corpus_manifest.toml");
    if !manifest_path.exists() {
        println!("  Creating corpus manifest at {}", manifest_path.display());
        let manifest = CORPUS_MANIFEST_TEMPLATE;
        std::fs::write(&manifest_path, manifest)
            .map_err(|e| MycoNoteError::Io(e))?;
    }

    // Validate the manifest
    println!("  Validating corpus manifest…");
    match validate_corpus_manifest(&manifest_path) {
        Ok(stats) => {
            println!("  ✓ Manifest valid: {} papers, all Q1 open-access", stats.total);
            if stats.missing_text > 0 {
                println!(
                    "    ⚠  {} papers have no extracted text yet — place .txt files in {}",
                    stats.missing_text,
                    papers_dir.display()
                );
            }
        }
        Err(e) => {
            eprintln!("  ⚠  Manifest validation failed: {}", e);
            eprintln!("     Fix issues in {} and re-run setup", manifest_path.display());
        }
    }

    // Symlink ~/.myconote/papers → db_dir/papers if they differ
    let symlink_target = home_papers.join("corpus_manifest.toml");
    if !symlink_target.exists() && manifest_path.exists() {
        #[cfg(unix)]
        {
            let _ = std::os::unix::fs::symlink(&manifest_path, &symlink_target);
        }
    }

    println!("  ✓ Chat corpus ready at {}", papers_dir.display());
    println!();
    println!("  To add papers:");
    println!("    1. Place extracted .txt files in {}", papers_dir.display());
    println!("    2. Add entries to {}", manifest_path.display());
    println!("    3. Each paper must have quartile = \"Q1\" and a valid OA license");
    println!();
    println!("  For personal papers (legally obtained):");
    println!("    Place .txt files in {}", home_papers.join("local").display());
    println!("    These are used for local retrieval only.");

    let _ = write_db_version(db_dir, "chat-corpus", "local-corpus");
    Ok(())
}

fn setup_ollama(db_dir: &Path) -> Result<()> {
    let ollama_dir = db_dir.join("ollama");
    std::fs::create_dir_all(&ollama_dir).map_err(MycoNoteError::Io)?;

    // ── 1. Check if Ollama is installed ──
    let installed_version = get_ollama_version();

    match &installed_version {
        Some(ver) => println!("  Ollama found: v{}", ver),
        None => {
            println!("  Ollama not found — installing latest version…");
            install_ollama()?;

            // Verify installation
            match get_ollama_version() {
                Some(ver) => println!("  ✓ Ollama installed: v{}", ver),
                None => {
                    eprintln!("  ⚠  Ollama installation did not succeed.");
                    eprintln!("     Install manually: https://ollama.com/download");
                    return Ok(());
                }
            }
        }
    }

    // ── 2. Check for updates ──
    let current_version = get_ollama_version().unwrap_or_default();
    let latest_version = check_ollama_latest_version();

    if let Some(ref latest) = latest_version {
        if !current_version.is_empty() && latest != &current_version {
            println!("  Update available: v{} → v{}", current_version, latest);
            println!("  Updating Ollama…");
            install_ollama()?;
            if let Some(new_ver) = get_ollama_version() {
                println!("  ✓ Ollama updated to v{}", new_ver);
            }
        } else {
            println!("  ✓ Ollama is up to date (v{})", current_version);
        }
    }

    // ── 3. Ensure Ollama is running ──
    let is_running = check_ollama_running();
    if !is_running {
        println!("  Starting Ollama server…");
        start_ollama_server();
        // Give it a moment to start
        std::thread::sleep(std::time::Duration::from_secs(3));
        if !check_ollama_running() {
            eprintln!("  ⚠  Could not start Ollama. Run `ollama serve` manually.");
        } else {
            println!("  ✓ Ollama server running");
        }
    } else {
        println!("  ✓ Ollama server already running");
    }

    // ── 4. Select and pull the best model for this machine ──
    let model = std::env::var("MYCONOTE_CHAT_MODEL").unwrap_or_else(|_| {
        let config_model = read_configured_model();
        config_model.unwrap_or_else(|| {
            // Auto-select the most capable model for available memory
            use crate::chat::config::{recommend_model, MODEL_TIERS};
            let recommended = recommend_model();
            println!("  Detecting system memory…");
            println!("  Available model tiers:");
            for tier in MODEL_TIERS {
                let marker = if tier.name == recommended { " ◀ selected" } else { "" };
                println!("    {:<45} {}{}", tier.name, tier.description, marker);
            }
            recommended.to_string()
        })
    });

    println!("  Pulling model '{}' (latest version)…", model);
    let pull_ok = Command::new("ollama")
        .args(["pull", &model])
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if pull_ok {
        println!("  ✓ Model '{}' ready", model);
    } else {
        eprintln!("  ⚠  Failed to pull model '{}'. Try: ollama pull {}", model, model);
        eprintln!("     Available models: https://ollama.com/library");
    }

    // ── 5. Verify the model works ──
    if check_ollama_running() {
        println!("  Verifying model…");
        let verify_ok = Command::new("ollama")
            .args(["run", &model, "Say OK"])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);

        if verify_ok {
            println!("  ✓ Model '{}' verified", model);
        } else {
            eprintln!("  ⚠  Model verification failed — the model may still be loading");
        }
    }

    // ── 6. Write version marker ──
    let final_version = get_ollama_version().unwrap_or_else(|| "unknown".to_string());
    let marker = ollama_dir.join("installed.version");
    let _ = std::fs::write(&marker, format!(
        "ollama_version={}\nmodel={}\ntimestamp={}\nmyconote_version={}\n",
        final_version,
        model,
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs(),
        env!("CARGO_PKG_VERSION"),
    ));

    let _ = write_db_version(db_dir, "ollama", "https://ollama.com");
    println!("  ✓ Ollama setup complete");
    println!();
    println!("  Usage:");
    println!("    myconote-cli explain predict              # interpret with LLM");
    println!("    myconote-cli explain predict --no-llm     # rules only, no LLM");
    println!("    myconote-cli explain predict --model {}   # use this model", model);
    Ok(())
}

/// Get the installed Ollama version, if any.
fn get_ollama_version() -> Option<String> {
    let output = Command::new("ollama")
        .args(["--version"])
        .output()
        .ok()?;

    if !output.status.success() {
        return None;
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    // "ollama version is 0.3.14" or "ollama version 0.3.14"
    let version = stdout.trim()
        .rsplit(' ')
        .next()
        .unwrap_or("")
        .trim()
        .to_string();

    if version.is_empty() { None } else { Some(version) }
}

/// Check the latest available Ollama version from GitHub releases.
fn check_ollama_latest_version() -> Option<String> {
    // Use GitHub API to get latest release tag
    let output = Command::new("curl")
        .args([
            "-sL",
            "--max-time", "10",
            "-H", "Accept: application/json",
            "https://api.github.com/repos/ollama/ollama/releases/latest",
        ])
        .output()
        .ok()?;

    if !output.status.success() {
        return None;
    }

    let body = String::from_utf8_lossy(&output.stdout);
    // Parse "tag_name":"v0.3.14" — simple extraction without a JSON crate
    let tag_prefix = "\"tag_name\":\"";
    let start = body.find(tag_prefix)? + tag_prefix.len();
    let end = body[start..].find('"')? + start;
    let tag = &body[start..end];

    // Strip leading 'v' if present
    Some(tag.strip_prefix('v').unwrap_or(tag).to_string())
}

/// Install or update Ollama using the official install script (Linux)
/// or brew (macOS).
fn install_ollama() -> Result<()> {
    #[cfg(target_os = "macos")]
    {
        println!("  Attempting: brew install ollama");
        let ok = Command::new("brew")
            .args(["install", "ollama"])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);

        if !ok {
            // Try upgrade if already installed
            let _ = Command::new("brew")
                .args(["upgrade", "ollama"])
                .status();
        }
        return Ok(());
    }

    #[cfg(target_os = "linux")]
    {
        println!("  Attempting: curl -fsSL https://ollama.com/install.sh | sh");
        let ok = Command::new("sh")
            .args(["-c", "curl -fsSL https://ollama.com/install.sh | sh"])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);

        if !ok {
            eprintln!("  ⚠  Automatic install failed.");
            eprintln!("     Install manually: https://ollama.com/download");
        }
        return Ok(());
    }

    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    {
        eprintln!("  ⚠  Automatic Ollama install not supported on this platform.");
        eprintln!("     Download from: https://ollama.com/download");
        Ok(())
    }
}

/// Check if the Ollama server is reachable.
fn check_ollama_running() -> bool {
    let endpoint = std::env::var("MYCONOTE_CHAT_ENDPOINT")
        .unwrap_or_else(|_| "http://localhost:11434".to_string());

    Command::new("curl")
        .args(["-s", "--max-time", "3", &endpoint])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Start Ollama server in the background.
fn start_ollama_server() {
    let _ = Command::new("ollama")
        .arg("serve")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn();
}

/// Read the model from ~/.myconote/config.toml if it exists.
fn read_configured_model() -> Option<String> {
    let config_path = home_myconote_dir().join("config.toml");
    let text = std::fs::read_to_string(&config_path).ok()?;
    let table: toml::Value = text.parse().ok()?;
    table.get("chat")?
        .get("model")?
        .as_str()
        .map(|s| s.to_string())
}

fn home_myconote_dir() -> PathBuf {
    let home = std::env::var("HOME")
        .or_else(|_| std::env::var("USERPROFILE"))
        .unwrap_or_else(|_| "/tmp".to_string());
    PathBuf::from(home).join(".myconote")
}

// ─────────────────────────────────────────────────────────────────────────────
// Corpus manifest validation
// ─────────────────────────────────────────────────────────────────────────────

const VALID_LICENSES: &[&str] = &["CC-BY", "CC-BY-SA", "CC-BY-4.0", "CC-BY-SA-4.0", "CC0", "public-domain"];
const VALID_QUARTILES: &[&str] = &["Q1"];

#[derive(Debug)]
struct CorpusStats {
    total: usize,
    missing_text: usize,
}

fn validate_corpus_manifest(path: &Path) -> std::result::Result<CorpusStats, String> {
    let contents = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {}", path.display(), e))?;

    let table: toml::Value = contents.parse()
        .map_err(|e| format!("invalid TOML: {}", e))?;

    let papers = table.get("paper")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "missing [[paper]] array".to_string())?;

    let mut total = 0usize;
    let mut missing_text = 0usize;

    for (i, paper) in papers.iter().enumerate() {
        let doi = paper.get("doi").and_then(|v| v.as_str()).unwrap_or("<missing>");
        let quartile = paper.get("quartile").and_then(|v| v.as_str());
        let license = paper.get("license").and_then(|v| v.as_str());
        let text_file = paper.get("text_file").and_then(|v| v.as_str());

        // Validate quartile
        match quartile {
            Some(q) if VALID_QUARTILES.contains(&q) => {}
            Some(q) => return Err(format!("paper #{} ({}): quartile '{}' is not Q1", i + 1, doi, q)),
            None => return Err(format!("paper #{} ({}): missing quartile field", i + 1, doi)),
        }

        // Validate license
        match license {
            Some(l) if VALID_LICENSES.contains(&l) => {}
            Some(l) => return Err(format!("paper #{} ({}): license '{}' is not open-access", i + 1, doi, l)),
            None => return Err(format!("paper #{} ({}): missing license field", i + 1, doi)),
        }

        // Check text file presence
        if let Some(tf) = text_file {
            let text_path = path.parent().unwrap_or(Path::new(".")).join(tf);
            if !text_path.exists() {
                missing_text += 1;
            }
        } else {
            missing_text += 1;
        }

        total += 1;
    }

    if total == 0 {
        return Err("no papers in manifest".to_string());
    }

    Ok(CorpusStats { total, missing_text })
}

const CORPUS_MANIFEST_TEMPLATE: &str = r#"# myconote-cli chat corpus manifest
#
# Each [[paper]] entry describes a Q1 open-access paper used for
# citation-backed explanations in `myconote explain`.
#
# Requirements:
#   - quartile must be "Q1"
#   - license must be an open-access license (CC-BY, CC-BY-SA, CC0, etc.)
#   - text_file points to extracted plain text (relative to this directory)
#
# To add papers:
#   1. Download the paper (open-access only)
#   2. Extract text to a .txt file in this directory
#   3. Add an entry below
#
# Journals accepted (Q1, open-access):
#   Genome Biology, Genome Research, NAR, Bioinformatics, PLOS Biology,
#   Nature Communications (OA), BMC Genomics, GigaScience, Molecular
#   Biology and Evolution, PNAS (OA)

[[paper]]
doi = "10.1186/s13059-019-1832-y"
title = "Funannotate: a comprehensive tool for functional annotation of fungal genomes"
journal = "Genome Biology"
year = 2019
quartile = "Q1"
license = "CC-BY-4.0"
tags = ["annotation", "fungi", "pipeline"]
text_file = "funannotate_2019.txt"

[[paper]]
doi = "10.1093/nar/gkab065"
title = "InterPro in 2021: an integrative protein signature database"
journal = "Nucleic Acids Research"
year = 2021
quartile = "Q1"
license = "CC-BY-4.0"
tags = ["functional-annotation", "domains", "interpro"]
text_file = "interpro_2021.txt"

[[paper]]
doi = "10.1093/bioinformatics/btv351"
title = "BUSCO: assessing genome assembly and annotation completeness"
journal = "Bioinformatics"
year = 2015
quartile = "Q1"
license = "CC-BY-4.0"
tags = ["quality", "completeness", "busco"]
text_file = "busco_2015.txt"

[[paper]]
doi = "10.1093/nar/gkaa913"
title = "Pfam: the protein families database in 2021"
journal = "Nucleic Acids Research"
year = 2021
quartile = "Q1"
license = "CC-BY-4.0"
tags = ["domains", "pfam", "functional-annotation"]
text_file = "pfam_2021.txt"

[[paper]]
doi = "10.1186/s13059-019-1715-2"
title = "Repeat masking and genome annotation in fungal genomes"
journal = "Genome Biology"
year = 2019
quartile = "Q1"
license = "CC-BY-4.0"
tags = ["repeats", "masking", "fungi"]
text_file = "repeatmasking_2019.txt"
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_manifest_valid() {
        let dir = tempfile::tempdir().unwrap();
        let manifest = dir.path().join("corpus_manifest.toml");
        std::fs::write(&manifest, r#"
[[paper]]
doi = "10.1186/test"
title = "Test Paper"
journal = "Genome Biology"
year = 2020
quartile = "Q1"
license = "CC-BY-4.0"
tags = ["test"]
text_file = "test.txt"
"#).unwrap();
        // text file doesn't exist, so missing_text should be 1
        let result = validate_corpus_manifest(&manifest).unwrap();
        assert_eq!(result.total, 1);
        assert_eq!(result.missing_text, 1);
    }

    #[test]
    fn validate_manifest_with_text_file() {
        let dir = tempfile::tempdir().unwrap();
        let manifest = dir.path().join("corpus_manifest.toml");
        let text_file = dir.path().join("test.txt");
        std::fs::write(&text_file, "Some paper content here.").unwrap();
        std::fs::write(&manifest, r#"
[[paper]]
doi = "10.1186/test"
title = "Test Paper"
journal = "Genome Biology"
year = 2020
quartile = "Q1"
license = "CC-BY-4.0"
tags = ["test"]
text_file = "test.txt"
"#).unwrap();
        let result = validate_corpus_manifest(&manifest).unwrap();
        assert_eq!(result.total, 1);
        assert_eq!(result.missing_text, 0);
    }

    #[test]
    fn validate_manifest_rejects_non_q1() {
        let dir = tempfile::tempdir().unwrap();
        let manifest = dir.path().join("corpus_manifest.toml");
        std::fs::write(&manifest, r#"
[[paper]]
doi = "10.1186/test"
title = "Test"
journal = "Low Impact Journal"
year = 2020
quartile = "Q3"
license = "CC-BY-4.0"
tags = []
"#).unwrap();
        let result = validate_corpus_manifest(&manifest);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not Q1"));
    }

    #[test]
    fn validate_manifest_rejects_non_oa_license() {
        let dir = tempfile::tempdir().unwrap();
        let manifest = dir.path().join("corpus_manifest.toml");
        std::fs::write(&manifest, r#"
[[paper]]
doi = "10.1186/test"
title = "Test"
journal = "Good Journal"
year = 2020
quartile = "Q1"
license = "proprietary"
tags = []
"#).unwrap();
        let result = validate_corpus_manifest(&manifest);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not open-access"));
    }

    #[test]
    fn validate_manifest_empty_papers() {
        let dir = tempfile::tempdir().unwrap();
        let manifest = dir.path().join("corpus_manifest.toml");
        std::fs::write(&manifest, "# empty manifest\n[metadata]\nversion = 1\n").unwrap();
        let result = validate_corpus_manifest(&manifest);
        assert!(result.is_err());
    }

    #[test]
    fn template_parses_as_valid_toml() {
        let table: toml::Value = CORPUS_MANIFEST_TEMPLATE.parse().unwrap();
        let papers = table.get("paper").unwrap().as_array().unwrap();
        assert_eq!(papers.len(), 5);
    }

    #[test]
    fn database_catalog_has_chat_corpus() {
        let entry = DATABASES.iter().find(|d| d.key == "chat-corpus");
        assert!(entry.is_some());
        assert!(entry.unwrap().marker_file.contains("corpus_manifest"));
    }
}
