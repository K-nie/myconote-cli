use crate::utils::error::Result;
use std::collections::HashMap;
use std::path::PathBuf;

pub mod jbrowse;
pub mod synteny;
pub mod ucsc;

/// Which genome browser backend to use for visualisation
#[derive(Debug, Clone, PartialEq)]
pub enum BrowserType {
    JBrowse2,
    UCSC,
    NCBI,
}

impl BrowserType {
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "ucsc" => BrowserType::UCSC,
            "ncbi" => BrowserType::NCBI,
            _ => BrowserType::JBrowse2,
        }
    }
}

/// Configuration for the `view` command
#[derive(Debug, Clone)]
pub struct ViewConfig {
    /// Path to the GFF3 annotation file
    pub gff_path: PathBuf,
    /// Optional FASTA reference (enables sequence track)
    pub fasta_path: Option<PathBuf>,
    /// Output file path (HTML)
    pub output: PathBuf,
    /// Which browser backend to use
    pub browser: BrowserType,
    /// Optional genomic region to focus on (e.g. "contig1:1000-50000")
    pub region: Option<String>,
    /// Title shown in the viewer header
    pub title: Option<String>,
    /// Species / assembly name (defaults to filename stem)
    pub assembly_name: Option<String>,
    /// If true, print a UCSC custom-track URL in addition to the HTML output
    pub also_ucsc: bool,
    /// Optional gene ID → display name mapping (from --names file or API fetch).
    /// When present, gene IDs in GFF3 are replaced with these display labels
    /// in the browser view.
    pub names: HashMap<String, String>,
}

impl Default for ViewConfig {
    fn default() -> Self {
        Self {
            gff_path: PathBuf::new(),
            fasta_path: None,
            output: PathBuf::from("genome_view.html"),
            browser: BrowserType::JBrowse2,
            region: None,
            title: None,
            assembly_name: None,
            also_ucsc: false,
            names: HashMap::new(),
        }
    }
}

/// Top-level entry point called from main
pub fn generate_view(config: &ViewConfig) -> Result<()> {
    match &config.browser {
        BrowserType::JBrowse2 => {
            jbrowse::generate_jbrowse_html(config)?;
            if config.also_ucsc {
                ucsc::print_ucsc_url(config)?;
            }
        }
        BrowserType::UCSC => {
            ucsc::generate_ucsc_track_file(config)?;
            ucsc::print_ucsc_url(config)?;
        }
        BrowserType::NCBI => {
            println!("ℹ  NCBI Sequence Viewer integration is planned for a future release.");
            println!("   For now, use --browser jbrowse2 (default) or --browser ucsc.");
        }
    }
    Ok(())
}
