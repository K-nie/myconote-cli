use crate::utils::error::Result;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

pub mod alignment;
pub mod phylogeny;
pub mod synteny;

#[derive(Debug, Clone)]
pub struct CompareConfig {
    pub output_dir: PathBuf,
    pub method: String, // "blast", "mmseqs", "mummer"
    pub min_identity: f64,
    pub min_length: usize,
    pub threads: usize,
    pub generate_tree: bool,
    pub generate_synteny: bool,
}

impl Default for CompareConfig {
    fn default() -> Self {
        Self {
            output_dir: PathBuf::from("compare_results"),
            method: "mmseqs".to_string(),
            min_identity: 30.0,
            min_length: 100,
            threads: 4,
            generate_tree: false,
            generate_synteny: true,
        }
    }
}

impl CompareConfig {
    pub fn new(output_dir: &str) -> Self {
        Self {
            output_dir: PathBuf::from(output_dir),
            ..Default::default()
        }
    }

    pub fn with_method(mut self, method: &str) -> Self {
        self.method = method.to_string();
        self
    }

    pub fn with_min_identity(mut self, min_identity: f64) -> Self {
        self.min_identity = min_identity;
        self
    }

    pub fn with_min_length(mut self, min_length: usize) -> Self {
        self.min_length = min_length;
        self
    }

    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    pub fn with_tree(mut self) -> Self {
        self.generate_tree = true;
        self
    }

    pub fn with_synteny(mut self) -> Self {
        self.generate_synteny = true;
        self
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct GenomeInfo {
    pub path: PathBuf,
    pub name: String,
    pub gene_count: usize,
    pub total_length: u64,
}

pub fn compare_genomes<P: AsRef<Path>>(genome_paths: &[P], config: &CompareConfig) -> Result<()> {
    println!("\n🔬 Comparing {} genomes", genome_paths.len());
    println!("Method: {}", config.method);
    println!("Output directory: {}", config.output_dir.display());

    // Create output directory
    std::fs::create_dir_all(&config.output_dir)?;

    // Load genome information
    let mut genomes = Vec::new();
    for path in genome_paths {
        let path = path.as_ref();
        let name = path
            .file_stem()
            .unwrap_or_else(|| path.as_os_str())
            .to_string_lossy()
            .to_string();

        println!("Loading: {} ({})", name, path.display());

        // Count genes (simplified - just count gene features)
        let mut gene_count = 0;
        let mut total_length = 0;
        let reader = crate::parser::gff::GFFReader::from_path(path)?;

        for result in reader {
            if let Ok(record) = result {
                if record.feature_type == "gene" {
                    gene_count += 1;
                    total_length += record.length();
                }
            }
        }

        genomes.push(GenomeInfo {
            path: path.to_path_buf(),
            name,
            gene_count,
            total_length,
        });
    }

    // Save genome info
    let info_path = config.output_dir.join("genomes.json");
    let info_json = serde_json::to_string_pretty(&genomes)?;
    std::fs::write(info_path, info_json)?;

    // Run alignment based on method
    match config.method.as_str() {
        "blast" => alignment::run_blast(&genomes, config)?,
        "mmseqs" => alignment::run_mmseqs(&genomes, config)?,
        "mummer" => alignment::run_mummer(&genomes, config)?,
        _ => println!("Unknown alignment method: {}", config.method),
    }

    // Generate synteny plot if requested
    if config.generate_synteny {
        synteny::generate_synteny_plot(&genomes, config)?;
    }

    // Generate phylogenetic tree if requested
    if config.generate_tree {
        phylogeny::build_tree(&genomes, config)?;
    }

    println!(
        "\n✅ Comparison complete! Results in: {}",
        config.output_dir.display()
    );

    Ok(())
}
