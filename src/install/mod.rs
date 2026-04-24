/// Dependency installer: `myconote install`
///
/// Checks which external tools required by myconote-cli are missing and
/// installs them automatically via conda/mamba or pip.  Tools that require a
/// manual licence download (GeneMark) are flagged clearly so the user knows
/// exactly what to do.  SignalP and TMHMM are replaced by their free
/// successors: DeepSig (pip) and DeepTMHMM via pybiolib (pip).
///
/// Usage:
///   myconote-cli install                 # install everything missing
///   myconote-cli install predict         # only tools needed for predict
///   myconote-cli install --yes           # skip confirmation prompt
///   myconote-cli install --mamba         # prefer mamba over conda
use std::io::{self, Write};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Tool catalogue
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Clone)]
#[allow(dead_code)]
struct Tool {
    /// Binary name searched in PATH
    name: &'static str,
    /// Which myconote sub-command(s) use it
    used_by: &'static str,
    /// Conda package name, or None if not available via conda
    conda_pkg: Option<&'static str>,
    /// Conda channel (ignored when conda_pkg is None)
    conda_chan: &'static str,
    /// PyPI package name for pip install, or None if not available via pip
    pip_pkg: Option<&'static str>,
    /// Manual install note shown when no automated install is possible
    manual_note: &'static str,
    /// Flag to pass for a version check (empty = skip version)
    version_arg: &'static str,
}

const TOOLS: &[Tool] = &[
    // ── Ab-initio gene predictors ─────────────────────────────────────────────
    Tool {
        name: "augustus",
        used_by: "predict/train",
        conda_pkg: Some("augustus"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/Gaius-Augustus/Augustus",
        version_arg: "--version",
    },
    Tool {
        name: "snap",
        used_by: "predict/train",
        conda_pkg: Some("snap"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/KorfLab/SNAP",
        version_arg: "",
    },
    Tool {
        name: "glimmerhmm",
        used_by: "predict",
        conda_pkg: Some("glimmerhmm"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://ccb.jhu.edu/software/glimmerhmm/",
        version_arg: "--help",
    },
    Tool {
        name: "gmes_petap.pl",
        used_by: "predict",
        conda_pkg: None,
        conda_chan: "",
        pip_pkg: None,
        manual_note: "Licence required — http://topaz.gatech.edu/GeneMark/",
        version_arg: "--version",
    },
    // ── Repeat masking ───────────────────────────────────────────────────────
    Tool {
        name: "RepeatMasker",
        used_by: "mask",
        conda_pkg: Some("repeatmasker"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://www.repeatmasker.org/",
        version_arg: "--version",
    },
    Tool {
        name: "RepeatModeler",
        used_by: "mask",
        conda_pkg: Some("repeatmodeler"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://www.repeatmasker.org/RepeatModeler/",
        version_arg: "--version",
    },
    Tool {
        name: "tantan",
        used_by: "mask",
        conda_pkg: Some("tantan"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://gitlab.com/mcfrith/tantan",
        version_arg: "--version",
    },
    // ── RNA-seq / training ───────────────────────────────────────────────────
    Tool {
        name: "Trinity",
        used_by: "train/update",
        conda_pkg: Some("trinity"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/trinityrnaseq/trinityrnaseq",
        version_arg: "--version",
    },
    Tool {
        name: "minimap2",
        used_by: "train/update/synteny",
        conda_pkg: Some("minimap2"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/lh3/minimap2",
        version_arg: "--version",
    },
    Tool {
        name: "samtools",
        used_by: "train/update",
        conda_pkg: Some("samtools"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://www.htslib.org/",
        version_arg: "--version",
    },
    Tool {
        name: "Launch_PASA_pipeline.pl",
        used_by: "train/update",
        conda_pkg: Some("pasa"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/PASApipeline/PASApipeline",
        version_arg: "--version",
    },
    Tool {
        name: "TransDecoder.LongOrfs",
        used_by: "train",
        conda_pkg: Some("transdecoder"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/TransDecoder/TransDecoder",
        version_arg: "--version",
    },
    // ── Homology / domain search ─────────────────────────────────────────────
    Tool {
        name: "mmseqs",
        used_by: "annotate",
        conda_pkg: Some("mmseqs2"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/soedinglab/MMseqs2",
        version_arg: "version",
    },
    Tool {
        name: "diamond",
        used_by: "annotate",
        conda_pkg: Some("diamond"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/bbuchfink/diamond",
        version_arg: "version",
    },
    Tool {
        name: "hmmscan",
        used_by: "annotate",
        conda_pkg: Some("hmmer"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://hmmer.org/",
        version_arg: "--version",
    },
    // ── Functional annotation ─────────────────────────────────────────────────
    Tool {
        name: "emapper.py",
        used_by: "annotate",
        conda_pkg: Some("eggnog-mapper"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/eggnogdb/eggnog-mapper",
        version_arg: "--version",
    },
    Tool {
        name: "run_dbcan",
        used_by: "annotate",
        conda_pkg: Some("dbcan"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://bcb.unl.edu/dbCAN2/",
        version_arg: "--version",
    },
    // Signal peptide prediction. DeepSig (bioconda) requires TF 2.2.0 which is
    // no longer available. Best free option: SignalP 6.0 from DTU (free after
    // registration). The secretome pipeline also accepts signalp or signalp6
    // if already installed. Without any of these, SP prediction is skipped.
    Tool {
        name: "signalp6",
        used_by: "annotate",
        conda_pkg: None,
        conda_chan: "",
        pip_pkg: None,
        manual_note:
            "Free download (registration): https://services.healthtech.dtu.dk/services/SignalP-6.0/",
        version_arg: "--version",
    },
    // DeepTMHMM replaces TMHMM (academic licence). pip install pybiolib.
    // Provides the `biolib` CLI: `biolib run DTU/DeepTMHMM --fasta <file>`.
    Tool {
        name: "biolib",
        used_by: "annotate",
        conda_pkg: None,
        conda_chan: "",
        pip_pkg: Some("pybiolib"),
        manual_note: "https://github.com/biolib/pybiolib",
        version_arg: "--version",
    },
    Tool {
        name: "antismash",
        used_by: "annotate",
        conda_pkg: Some("antismash"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://antismash.secondarymetabolites.org/",
        version_arg: "--version",
    },
    Tool {
        name: "busco",
        used_by: "annotate",
        conda_pkg: Some("busco"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://busco.ezlab.org/",
        version_arg: "--version",
    },
    Tool {
        name: "tRNAscan-SE",
        used_by: "annotate",
        conda_pkg: Some("trnascan-se"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "http://lowelab.ucsc.edu/tRNAscan-SE/",
        version_arg: "--version",
    },
    Tool {
        name: "miniprot",
        used_by: "predict",
        conda_pkg: Some("miniprot"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/lh3/miniprot",
        version_arg: "--version",
    },
    Tool {
        name: "table2asn",
        used_by: "submit",
        conda_pkg: None,
        conda_chan: "",
        pip_pkg: None,
        manual_note:
            "https://ftp.ncbi.nlm.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn/",
        version_arg: "--help",
    },
    // ── BLAST / alignment ────────────────────────────────────────────────────
    Tool {
        name: "blastp",
        used_by: "blast",
        conda_pkg: Some("blast"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://blast.ncbi.nlm.nih.gov/",
        version_arg: "-version",
    },
    // ── Comparative genomics ─────────────────────────────────────────────────
    Tool {
        name: "orthofinder",
        used_by: "compare",
        conda_pkg: Some("orthofinder"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/davidemms/OrthoFinder",
        version_arg: "-h",
    },
    // ── RNA-seq expression quantification ────────────────────────────────────
    Tool {
        name: "salmon",
        used_by: "quant",
        conda_pkg: Some("salmon"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://salmon.readthedocs.io/",
        version_arg: "--version",
    },
    Tool {
        name: "fastp",
        used_by: "quant",
        conda_pkg: Some("fastp"),
        conda_chan: "bioconda",
        pip_pkg: None,
        manual_note: "https://github.com/OpenGene/fastp",
        version_arg: "--version",
    },
];

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub struct InstallOptions {
    /// Only install tools required by this subcommand (None = all)
    pub filter_cmd: Option<String>,
    /// Skip the confirmation prompt
    pub yes: bool,
    /// Prefer mamba over conda
    pub mamba: bool,
}

pub fn run_install(opts: InstallOptions) {
    let divider = "─".repeat(70);

    println!("\nmyconote dependency installer");
    println!("{divider}");

    // ── Step 0: heal broken mamba/libarchive before doing anything else ───────
    let pkg_mgr = heal_and_detect_package_manager(opts.mamba);

    // ── Step 1: scan for missing tools ───────────────────────────────────────
    let mut conda_install: Vec<&Tool> = Vec::new();
    let mut pip_install: Vec<&Tool> = Vec::new();
    let mut manual: Vec<&Tool> = Vec::new();
    let mut already_ok: Vec<&Tool> = Vec::new();

    for tool in TOOLS {
        if let Some(ref cmd) = opts.filter_cmd {
            if !tool.used_by.contains(cmd.as_str()) {
                continue;
            }
        }
        if tool_in_path(tool.name) {
            already_ok.push(tool);
        } else if tool.conda_pkg.is_some() && pkg_mgr.is_some() {
            conda_install.push(tool);
        } else if tool.pip_pkg.is_some() {
            pip_install.push(tool);
        } else {
            manual.push(tool);
        }
    }

    // ── Step 2: report status ────────────────────────────────────────────────
    if !already_ok.is_empty() {
        println!("\n  Already installed ({}):", already_ok.len());
        for t in &already_ok {
            println!("    \x1b[32m✓\x1b[0m {:<28} ({})", t.name, t.used_by);
        }
    }

    if !conda_install.is_empty() {
        let mgr = pkg_mgr.as_deref().unwrap_or("conda");
        println!("\n  Will install via {} ({}):", mgr, conda_install.len());
        for t in &conda_install {
            println!(
                "    \x1b[33m○\x1b[0m {:<28} {} install -c {} {}",
                t.name,
                mgr,
                t.conda_chan,
                t.conda_pkg.unwrap_or("")
            );
        }
    }

    if !pip_install.is_empty() {
        println!("\n  Will install via pip ({}):", pip_install.len());
        for t in &pip_install {
            println!(
                "    \x1b[33m○\x1b[0m {:<28} pip install {}",
                t.name,
                t.pip_pkg.unwrap_or("")
            );
        }
    }

    if !manual.is_empty() {
        println!("\n  Require manual download ({}):", manual.len());
        for t in &manual {
            println!("    \x1b[31m✗\x1b[0m {:<28} {}", t.name, t.manual_note);
        }
    }

    if conda_install.is_empty() && pip_install.is_empty() && manual.is_empty() {
        println!("\n  \x1b[32mAll tools are already installed. You're good to go!\x1b[0m\n");
        return;
    }

    // ── Step 3: confirm ──────────────────────────────────────────────────────
    let has_auto = !conda_install.is_empty() || !pip_install.is_empty();
    if has_auto {
        if !opts.yes {
            println!("\n{divider}");
            print!("  Proceed with installation? [y/N] ");
            io::stdout().flush().ok();
            let mut input = String::new();
            io::stdin().read_line(&mut input).ok();
            if !matches!(input.trim().to_lowercase().as_str(), "y" | "yes") {
                println!("  Aborted — no changes made.\n");
                return;
            }
        }

        // ── Step 4a: install conda packages one at a time ────────────────────
        //
        // Bundling packages lets conda find a globally consistent solution but
        // the classic solver is exponentially slower with each additional pkg.
        // Trinity alone has 200+ transitive deps and causes 20+ min hangs when
        // combined with others.  Installing one-by-one keeps each solve tiny
        // (30–90 s).  Conflicts are rare with bioconda since all packages share
        // a curated common base.
        if !conda_install.is_empty() {
            println!();
            let mgr = pkg_mgr.as_deref().unwrap_or("conda");

            // Deduplicate package names, heavy packages last so lighter ones
            // succeed even if a heavy one fails.
            let heavy = ["trinity", "antismash", "busco", "eggnog-mapper", "dbcan"];
            let mut light_pkgs: Vec<&str> = Vec::new();
            let mut heavy_pkgs: Vec<&str> = Vec::new();
            for tool in &conda_install {
                let pkg = tool.conda_pkg.unwrap();
                if heavy.contains(&pkg) {
                    if !heavy_pkgs.contains(&pkg) {
                        heavy_pkgs.push(pkg);
                    }
                } else {
                    if !light_pkgs.contains(&pkg) {
                        light_pkgs.push(pkg);
                    }
                }
            }
            let pkgs: Vec<&str> = light_pkgs.into_iter().chain(heavy_pkgs).collect();

            println!("  Installing {} package(s) one at a time…\n", pkgs.len());
            println!("  (light packages first, heavy packages last)\n");

            let mut ok_count = 0usize;
            let mut fail_count = 0usize;

            for (i, pkg) in pkgs.iter().enumerate() {
                print!("  [{:>2}/{}] {:<28} … ", i + 1, pkgs.len(), pkg);
                io::stdout().flush().ok();

                // mamba/micromamba use libsolv natively — --solver flag not supported.
                // Only pass --solver classic for plain conda.
                let mut cmd = Command::new(mgr);
                cmd.arg("install").arg("-y");
                if mgr == "conda" {
                    cmd.arg("--solver").arg("classic");
                }
                cmd.arg("-c")
                    .arg("bioconda")
                    .arg("-c")
                    .arg("conda-forge")
                    .arg(pkg);
                let status = cmd.status();

                match status {
                    Ok(s) if s.success() => {
                        println!("\x1b[32mdone\x1b[0m");
                        ok_count += 1;
                    }
                    _ => {
                        println!("\x1b[31mFAILED\x1b[0m");
                        println!("      → Retry: {} install -y --solver classic -c bioconda -c conda-forge {}", mgr, pkg);
                        fail_count += 1;
                    }
                }
            }

            println!("\n{divider}");
            if fail_count == 0 {
                println!(
                    "  \x1b[32mAll {} conda package(s) installed.\x1b[0m",
                    ok_count
                );
            } else {
                println!(
                    "  {} installed, \x1b[31m{} failed\x1b[0m.",
                    ok_count, fail_count
                );
                println!(
                    "  Re-running `myconote-cli install --yes` will retry only what is missing."
                );
            }
        }

        // ── Step 4b: install pip-only tools one at a time ─────────────────────
        if !pip_install.is_empty() {
            println!();
            let mut pip_ok = 0usize;
            let mut pip_fail = 0usize;

            for tool in &pip_install {
                let pkg = tool.pip_pkg.unwrap();
                print!("  pip install {:<24} ", pkg);
                io::stdout().flush().ok();

                let status = Command::new("pip")
                    .args(["install", "--quiet", pkg])
                    .status();

                match status {
                    Ok(s) if s.success() => {
                        println!("\x1b[32mdone\x1b[0m");
                        pip_ok += 1;
                    }
                    _ => {
                        println!("\x1b[31mFAILED\x1b[0m");
                        println!("    → Try manually: pip install {}", pkg);
                        pip_fail += 1;
                    }
                }
            }

            println!("\n{divider}");
            if pip_fail == 0 {
                println!("  \x1b[32m{pip_ok} pip package(s) installed.\x1b[0m");
            } else {
                println!("  {pip_ok} pip installed, \x1b[31m{pip_fail} failed\x1b[0m.");
            }
        }
    }

    // ── Step 6: manual tools reminder ────────────────────────────────────────
    if !manual.is_empty() {
        println!("\n  The following tools require a manual download:");
        for t in &manual {
            println!("    \x1b[31m✗\x1b[0m {}  →  {}", t.name, t.manual_note);
        }
    }

    println!();
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Check whether a binary is present in PATH.
fn tool_in_path(name: &str) -> bool {
    Command::new("which")
        .arg(name)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Returns true when mamba is in PATH but its libarchive dylib is missing,
/// which happens when RepeatModeler (or another package) downgrades libarchive
/// below the version mamba was compiled against.
fn mamba_libarchive_broken() -> bool {
    if !tool_in_path("mamba") {
        return false;
    }
    // Run `mamba --version`; if it fails with a dylib error, mamba is broken.
    let out = Command::new("mamba").arg("--version").output();
    match out {
        Ok(o) => !o.status.success(),
        Err(_) => true,
    }
}

/// Attempt to restore the correct libarchive version so that mamba can load.
/// Uses conda with `--solver classic` (which does NOT need libmamba).
/// Returns true if the fix succeeded or was not needed.
fn try_fix_libarchive(conda: &str) -> bool {
    println!("  \x1b[33m⚠\x1b[0m  mamba's libarchive is out of date — attempting auto-fix…");
    println!(
        "       Running: {} install -y --solver classic -c conda-forge 'libarchive>=3.7'",
        conda
    );

    // Only pass --solver classic for conda (mamba/micromamba don't support it)
    let is_conda = conda == "conda";
    let mut cmd = Command::new(conda);
    cmd.arg("install").arg("-y");
    if is_conda {
        cmd.arg("--solver").arg("classic");
    }
    cmd.arg("-c").arg("conda-forge").arg("libarchive>=3.7");
    let status = cmd.status();

    match status {
        Ok(s) if s.success() => {
            println!("  \x1b[32m✓\x1b[0m  libarchive restored — mamba is now usable.\n");
            true
        }
        _ => {
            println!("  \x1b[31m✗\x1b[0m  auto-fix failed. Falling back to conda classic solver.");
            println!("       To fix manually:  conda install -c conda-forge 'libarchive>=3.7' --solver classic\n");
            false
        }
    }
}

/// Detect an available package manager, healing broken mamba first if needed.
/// Preference order: mamba (fast) → conda (reliable fallback).
fn heal_and_detect_package_manager(prefer_mamba: bool) -> Option<String> {
    // If the caller explicitly wants mamba, or mamba is present, check health.
    if tool_in_path("mamba") || prefer_mamba {
        if mamba_libarchive_broken() {
            // Try to fix using conda (which doesn't need libmamba).
            if tool_in_path("conda") {
                let fixed = try_fix_libarchive("conda");
                if fixed && !mamba_libarchive_broken() {
                    return Some("mamba".to_string());
                }
            }
            // Fix failed or conda not present — fall through to conda.
        } else if tool_in_path("mamba") {
            // mamba is healthy
            return Some("mamba".to_string());
        }
    }
    if tool_in_path("conda") {
        return Some("conda".to_string());
    }
    None
}
