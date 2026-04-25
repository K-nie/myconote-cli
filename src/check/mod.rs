/// Dependency checker: `myconote check`
///
/// Scans PATH for all external tools required by myconote-cli and reports
/// which are available, which are missing, and where to install them.
///
/// Output is a coloured table with three columns:
///   [✓] tool    version   required_by
///   [✗] tool    —         required_by   → install suggestion
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Tool catalogue
// ─────────────────────────────────────────────────────────────────────────────

struct Tool {
    /// Binary name in PATH
    name: &'static str,
    /// Which myconote sub-command uses it
    used_by: &'static str,
    /// Installation hint
    install_cmd: &'static str,
    /// Flag to get a version string (empty → just run `which`)
    version_arg: &'static str,
}

const TOOLS: &[Tool] = &[
    // ── Core annotation ──────────────────────────────────────────────────────
    Tool {
        name: "augustus",
        used_by: "predict",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda augustus",
    },
    Tool {
        name: "snap",
        used_by: "predict/train",
        version_arg: "",
        install_cmd: "conda install -c bioconda snap",
    },
    Tool {
        name: "glimmerhmm",
        used_by: "predict",
        version_arg: "--help",
        install_cmd: "conda install -c bioconda glimmerhmm",
    },
    Tool {
        name: "gmes_petap.pl",
        used_by: "predict",
        version_arg: "--version",
        install_cmd: "License required — http://topaz.gatech.edu/GeneMark/",
    },
    Tool {
        name: "prothint.py",
        used_by: "predict",
        version_arg: "--version",
        install_cmd: "Bundled with GeneMark-ES tarball or bioconda braker3 — needed for --genemark-mode ep|etp",
    },
    Tool {
        name: "braker.pl",
        used_by: "predict",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda braker3 — needed for --use-braker",
    },
    // ── Repeat masking ───────────────────────────────────────────────────────
    Tool {
        name: "RepeatMasker",
        used_by: "mask",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda repeatmasker",
    },
    Tool {
        name: "RepeatModeler",
        used_by: "mask",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda repeatmodeler",
    },
    Tool {
        name: "tantan",
        used_by: "mask",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda tantan",
    },
    // ── RNA-seq / training ───────────────────────────────────────────────────
    Tool {
        name: "Trinity",
        used_by: "train/update",
        version_arg: "--cite",
        install_cmd: "conda install -c bioconda trinity",
    },
    Tool {
        name: "minimap2",
        used_by: "train/update/clean",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda minimap2",
    },
    Tool {
        name: "samtools",
        used_by: "train/update",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda samtools",
    },
    Tool {
        name: "Launch_PASA_pipeline.pl",
        used_by: "train/update",
        version_arg: "",
        install_cmd: "conda install -c bioconda pasa",
    },
    Tool {
        name: "TransDecoder.LongOrfs",
        used_by: "train",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda transdecoder",
    },
    // ── Homology / domain search ─────────────────────────────────────────────
    Tool {
        name: "mmseqs",
        used_by: "annotate",
        version_arg: "version",
        install_cmd: "conda install -c bioconda mmseqs2",
    },
    Tool {
        name: "diamond",
        used_by: "annotate",
        version_arg: "version",
        install_cmd: "conda install -c bioconda diamond",
    },
    Tool {
        name: "hmmscan",
        used_by: "annotate",
        version_arg: "-h",
        install_cmd: "conda install -c bioconda hmmer",
    },
    // ── Functional annotation ─────────────────────────────────────────────────
    Tool {
        name: "emapper.py",
        used_by: "annotate",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda eggnog-mapper",
    },
    Tool {
        name: "run_dbcan",
        used_by: "annotate",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda dbcan",
    },
    Tool {
        name: "signalp6",
        used_by: "annotate",
        version_arg: "--version",
        install_cmd: "Free download: https://services.healthtech.dtu.dk/services/SignalP-6.0/",
    },
    Tool {
        name: "biolib",
        used_by: "annotate",
        version_arg: "--version",
        install_cmd: "pip install pybiolib  (provides DeepTMHMM)",
    },
    Tool {
        name: "antismash",
        used_by: "annotate",
        version_arg: "",
        install_cmd: "conda install -c bioconda antismash",
    },
    Tool {
        name: "busco",
        used_by: "annotate",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda busco",
    },
    // ── Secondary tools ───────────────────────────────────────────────────────
    Tool {
        name: "miniprot",
        used_by: "predict",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda miniprot",
    },
    Tool {
        name: "tRNAscan-SE",
        used_by: "annotate",
        version_arg: "-h",
        install_cmd: "conda install -c bioconda trnascan-se",
    },
    Tool {
        name: "table2asn",
        used_by: "submit",
        version_arg: "-version",
        install_cmd: "NCBI binary — download from https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/",
    },
    Tool {
        name: "orthofinder",
        used_by: "compare",
        version_arg: "-h",
        install_cmd: "conda install -c bioconda orthofinder",
    },
    Tool {
        name: "salmon",
        used_by: "quant",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda salmon",
    },
    Tool {
        name: "fastp",
        used_by: "quant",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda fastp",
    },
    Tool {
        name: "prefetch",
        used_by: "fetch-rna",
        version_arg: "--version",
        install_cmd: "conda install -c bioconda sra-tools  # optional — ENA backend is default",
    },
];

// ─────────────────────────────────────────────────────────────────────────────
// Check entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_check(filter_cmd: Option<&str>) {
    let divider = "─".repeat(72);
    println!("\nmyconote dependency check");
    println!("{}", divider);
    println!(
        "{:<22} {:<10} {:<16} {}",
        "Tool", "Status", "Used by", "Version / Note"
    );
    println!("{}", divider);

    let mut n_ok = 0usize;
    let mut n_miss = 0usize;

    for tool in TOOLS {
        // Filter by sub-command if requested
        if let Some(cmd) = filter_cmd {
            if !tool.used_by.contains(cmd) {
                continue;
            }
        }

        let (found, version) = check_tool(tool);

        if found {
            n_ok += 1;
            println!(
                "  \x1b[32m✓\x1b[0m {:<20} {:<10} {:<16} {}",
                tool.name,
                "OK",
                tool.used_by,
                version.as_deref().unwrap_or("")
            );
        } else {
            n_miss += 1;
            println!(
                "  \x1b[31m✗\x1b[0m {:<20} {:<10} {:<16} {}",
                tool.name, "MISSING", tool.used_by, tool.install_cmd
            );
        }
    }

    println!("{}", divider);
    println!("  {n_ok} tools available, {n_miss} missing.\n");

    if n_miss > 0 {
        println!("  Tip: run `myconote-cli install` to install all missing tools automatically.");
        println!("  All tools install via conda/mamba; pybiolib (DeepTMHMM) installs via pip.\n");
    }

    // ── Extra database-readiness checks ──────────────────────────────────────
    check_eggnog_db();
}

/// Check whether the eggNOG-mapper data directory has been populated.
/// The databases are separate from the binary and must be downloaded once.
fn check_eggnog_db() {
    // emapper stores data in ~/.eggnog_mapper/data  OR  a user-specified path
    let default_dir = dirs_home_str() + "/.eggnog_mapper/data";
    let marker = std::path::Path::new(&default_dir).join("eggnog.db");

    if !marker.exists() {
        println!(
            "  \x1b[33m⚠\x1b[0m  eggNOG-mapper databases not found at {}",
            default_dir
        );
        println!("     Run once to download (~50 GB) before using --eggnog:");
        println!(
            "       download_eggnog_data.py -y --data_dir {}",
            default_dir
        );
        println!();
    }
}

fn dirs_home_str() -> String {
    std::env::var("HOME").unwrap_or_else(|_| "/tmp".to_string())
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-tool check
// ─────────────────────────────────────────────────────────────────────────────

fn check_tool(tool: &Tool) -> (bool, Option<String>) {
    // First: confirm it's in PATH via `which`
    let in_path = Command::new("which")
        .arg(tool.name)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);

    if !in_path {
        return (false, None);
    }

    // Get version string if possible
    if tool.version_arg.is_empty() {
        return (true, None);
    }

    // Support multi-word version args (e.g. "version" vs "--version --flag")
    let args: Vec<&str> = tool.version_arg.split_whitespace().collect();

    let version = Command::new(tool.name)
        .args(&args)
        .output()
        .ok()
        .and_then(|o| {
            // Prefer stdout, fall back to stderr
            let stdout = String::from_utf8_lossy(&o.stdout).to_string();
            let stderr = String::from_utf8_lossy(&o.stderr).to_string();

            // Combine: try stdout lines first, then stderr lines
            let combined = format!("{}\n{}", stdout, stderr);

            // Find the first line that looks like a version / tool description,
            // skipping lines that are clearly errors, tracebacks, or warnings.
            let version_line = combined
                .lines()
                .filter(|l| {
                    let t = l.trim().to_lowercase();
                    !t.is_empty()
                        && !t.starts_with("error")
                        && !t.starts_with("warning")
                        && !t.starts_with("traceback")
                        && !t.starts_with("can't locate")
                        && !t.starts_with("there was an error")
                        && !t.contains("perhaps no data has been downloaded")
                        && !t.starts_with("  file ") // Python traceback frames
                })
                .next()
                .map(|l| {
                    let l = l.trim().to_string();
                    if l.len() > 40 {
                        l[..40].to_string() + "…"
                    } else {
                        l
                    }
                });

            version_line
        });

    (true, version)
}
