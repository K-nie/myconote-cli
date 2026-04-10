/// List available trained Augustus species: `myconote species`
///
/// Scans the Augustus configuration directory for trained species and
/// prints a formatted table showing:
///   - Species name (directory name used with --species=...)
///   - Organism / display name (from species.cfg if present)
///   - Completeness (whether all required HMM files are present)
///   - Whether it was trained by the user vs. shipped with Augustus
///
/// Augustus config directory is found via:
///   1. AUGUSTUS_CONFIG_PATH environment variable
///   2. `augustus --species=help` output
///   3. Common install locations (/opt/conda, /usr/share, etc.)
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Species entry
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct AugustusSpecies {
    pub name: String,
    pub display: String,    // from meta/*.cfg
    pub complete: bool,     // has parameters.cfg + all HMM files
    pub user_trained: bool, // found in user config dir, not system
}

// ─────────────────────────────────────────────────────────────────────────────
// Required files for a valid species
// ─────────────────────────────────────────────────────────────────────────────

const REQUIRED_FILES: &[&str] = &[
    "parameters.cfg",
    "exon_probs.pbl",
    "intron_probs.pbl",
    "igenic_probs.pbl",
];

// ─────────────────────────────────────────────────────────────────────────────
// Entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn list_species(filter: Option<&str>) {
    // Find config directory
    let config_dirs = find_augustus_config_dirs();

    if config_dirs.is_empty() {
        eprintln!("Could not find Augustus configuration directory.");
        eprintln!("Set AUGUSTUS_CONFIG_PATH or install Augustus:");
        eprintln!("  conda install -c bioconda augustus");
        return;
    }

    let mut all_species: Vec<AugustusSpecies> = Vec::new();
    let mut seen: HashMap<String, ()> = HashMap::new();

    for (dir, is_user) in &config_dirs {
        let species_dir = dir.join("species");
        if !species_dir.exists() {
            continue;
        }

        if let Ok(entries) = std::fs::read_dir(&species_dir) {
            for entry in entries.flatten() {
                let name = entry.file_name().to_string_lossy().to_string();
                if seen.contains_key(&name) {
                    continue;
                }
                if name.starts_with('.') {
                    continue;
                }

                let sp_path = entry.path();
                if !sp_path.is_dir() {
                    continue;
                }

                let complete = check_species_complete(&sp_path);
                let display = read_display_name(&sp_path, &name);

                all_species.push(AugustusSpecies {
                    name: name.clone(),
                    display,
                    complete,
                    user_trained: *is_user,
                });
                seen.insert(name, ());
            }
        }
    }

    // Sort: user species first, then alphabetical
    all_species.sort_by(|a, b| {
        b.user_trained
            .cmp(&a.user_trained)
            .then(a.name.cmp(&b.name))
    });

    // Apply filter
    let filtered: Vec<&AugustusSpecies> = if let Some(f) = filter {
        let f_lower = f.to_lowercase();
        all_species
            .iter()
            .filter(|s| {
                s.name.to_lowercase().contains(&f_lower)
                    || s.display.to_lowercase().contains(&f_lower)
            })
            .collect()
    } else {
        all_species.iter().collect()
    };

    if filtered.is_empty() {
        println!(
            "No Augustus species found matching {:?}",
            filter.unwrap_or("")
        );
        return;
    }

    // Print table
    println!("\nAvailable Augustus species ({} total):", filtered.len());
    println!("{}", "─".repeat(78));
    println!(
        "{:<35} {:<8} {:<8} {}",
        "Species name", "Status", "Source", "Description"
    );
    println!("{}", "─".repeat(78));

    for sp in &filtered {
        let status = if sp.complete { "✓" } else { "⚠ incomplete" };
        let source = if sp.user_trained { "user" } else { "system" };
        println!(
            "  {:<33} {:<8} {:<8} {}",
            sp.name,
            status,
            source,
            if sp.display == sp.name {
                ""
            } else {
                &sp.display
            },
        );
    }

    println!("{}", "─".repeat(78));
    println!("  Use with: myconote predict --species <name>");
    println!("  Train new: myconote train --rna-r1 reads_R1.fastq.gz ...\n");
}

// ─────────────────────────────────────────────────────────────────────────────
// Well-known grouped species for quick reference
// ─────────────────────────────────────────────────────────────────────────────

pub fn list_species_grouped() {
    let groups: &[(&str, &[&str])] = &[
        (
            "Fungi (recommended for myconote)",
            &[
                "aspergillus_fumigatus",
                "aspergillus_nidulans",
                "aspergillus_oryzae",
                "aspergillus_terreus",
                "botrytis_cinerea",
                "candida_albicans",
                "candida_guilliermondii",
                "candida_tropicalis",
                "chaetomium_globosum",
                "coccidioides_immitis",
                "coprinus_cinereus",
                "cryptococcus_neoformans_gattii",
                "cryptococcus_neoformans_neoformans_B",
                "encephalitozoon_cuniculi_GB",
                "fusarium_graminearum",
                "histoplasma_capsulatum",
                "neurospora_crassa",
                "pneumocystis_jirovecii",
                "rhizopus_oryzae",
                "saccharomyces_cerevisiae_S288C",
                "schizosaccharomyces_pombe",
                "ustilago_maydis",
                "yarrowia_lipolytica",
            ],
        ),
        (
            "Plants",
            &[
                "arabidopsis",
                "maize",
                "maize5",
                "rice",
                "tomato",
                "wheat",
                "Solanaceae",
            ],
        ),
        (
            "Animals / Metazoa",
            &[
                "fly",
                "honeybee1",
                "human",
                "mouse",
                "zebrafish",
                "nematode",
                "caenorhabditis",
                "Drosophila",
            ],
        ),
        ("Oomycetes / Protists", &["phytophthora", "Chlamydomonas"]),
    ];

    println!("\nAugustus species — grouped by kingdom");
    println!("{}", "═".repeat(60));
    for (group, species) in groups {
        println!("\n  {}:", group);
        for sp in *species {
            println!("    • {}", sp);
        }
    }
    println!("\n  Use: myconote predict --species <name>");
    println!("  Full list: myconote species --list\n");
}

// ─────────────────────────────────────────────────────────────────────────────
// Config directory finders
// ─────────────────────────────────────────────────────────────────────────────

fn find_augustus_config_dirs() -> Vec<(PathBuf, bool)> {
    let mut dirs: Vec<(PathBuf, bool)> = Vec::new();

    // 1. Environment variable
    if let Ok(v) = std::env::var("AUGUSTUS_CONFIG_PATH") {
        let p = PathBuf::from(&v);
        if p.exists() {
            dirs.push((p, false));
        }
    }

    // 2. User-trained species (~/.myconote/augustus_config)
    if let Ok(home) = std::env::var("HOME") {
        let user_cfg = PathBuf::from(home)
            .join(".myconote")
            .join("augustus_config");
        if user_cfg.exists() {
            dirs.push((user_cfg, true));
        }
    }

    // 3. Ask augustus
    if let Ok(out) = Command::new("augustus").arg("--species=help").output() {
        let text = String::from_utf8_lossy(&out.stderr).to_string()
            + &String::from_utf8_lossy(&out.stdout);
        for line in text.lines() {
            if line.contains("AUGUSTUS_CONFIG_PATH") || line.contains("config") {
                if let Some(p) = line.split_whitespace().last() {
                    let pb = PathBuf::from(p);
                    if pb.exists() && !dirs.iter().any(|(d, _)| d == &pb) {
                        dirs.push((pb, false));
                    }
                }
            }
        }
    }

    // 4. Common locations
    let common = [
        "/opt/conda/config",
        "/usr/share/augustus/config",
        "/usr/local/share/augustus/config",
        "/opt/augustus/config",
    ];
    for c in &common {
        let p = PathBuf::from(c);
        if p.exists() && !dirs.iter().any(|(d, _)| d == &p) {
            dirs.push((p, false));
        }
    }

    dirs
}

fn check_species_complete(sp_path: &Path) -> bool {
    REQUIRED_FILES.iter().all(|f| sp_path.join(f).exists())
}

fn read_display_name(sp_path: &Path, fallback: &str) -> String {
    // Look for species.cfg or meta/*.cfg
    let cfg = sp_path.join("species.cfg");
    if cfg.exists() {
        if let Ok(text) = std::fs::read_to_string(&cfg) {
            for line in text.lines() {
                if line.trim_start().starts_with("Species") || line.contains("organism") {
                    if let Some(val) = line.split('=').nth(1) {
                        let v = val.trim().to_string();
                        if !v.is_empty() {
                            return v;
                        }
                    }
                }
            }
        }
    }

    // Try meta/ directory
    if let Ok(entries) = std::fs::read_dir(sp_path.join("meta")) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.extension().map(|e| e == "cfg").unwrap_or(false) {
                if let Ok(text) = std::fs::read_to_string(&path) {
                    for line in text.lines() {
                        if line.to_lowercase().contains("name") && line.contains('=') {
                            if let Some(val) = line.split('=').nth(1) {
                                let v = val.trim().trim_matches('"').to_string();
                                if !v.is_empty() {
                                    return v;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fallback.to_string()
}
