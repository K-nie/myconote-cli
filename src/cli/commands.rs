use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about = "MycoNote: Blazing-fast genome annotation analysis", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Calculate statistics from annotation files
    Stats(StatsArgs),

    /// Convert between genomic formats
    Convert(ConvertArgs),

    /// Clean and validate annotation files
    Clean(CleanArgs),

    /// BLAST search against NCBI databases
    Blast(BlastArgs),

    /// Align sequences (BLAST, MMseqs2, MUMmer, minimap2)
    Align(AlignArgs),

    /// Build phylogenetic trees with IQ-TREE
    Phylogeny(PhylogenyArgs),

    /// Compare multiple genomes (alignment + visualization)
    Compare(CompareArgs),
}

#[derive(Args)]
pub struct StatsArgs {
    /// Input annotation file (GFF, GBK, FASTA)
    pub input: PathBuf,

    /// Output format (json, csv, human)
    #[arg(short, long, default_value = "human")]
    pub format: String,

    /// Taxonomic group for benchmarking (fungi, plants, animals, etc.)
    #[arg(long)]
    pub taxon: Option<String>,

    /// Use primary transcripts only (collapse isoforms)
    #[arg(long)]
    pub primary_only: bool,

    /// Specific chromosomes to analyze
    #[arg(long, value_delimiter = ',')]
    pub chromosomes: Option<Vec<String>>,

    /// Region in format "chr:start-end"
    #[arg(long)]
    pub region: Option<String>,

    /// Generate statistical plots
    #[arg(long)]
    pub plot: bool,

    /// Output directory for plots
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Args)]
pub struct BlastArgs {
    /// Input FASTA file with protein sequences
    pub input: PathBuf,

    /// Number of top genes to BLAST (max 50, warns above)
    #[arg(short, long, default_value_t = 10)]
    pub top_n: usize,

    /// BLAST database (nr, swissprot, etc.)
    #[arg(long, default_value = "nr")]
    pub database: String,

    /// Maximum hits per query
    #[arg(long, default_value_t = 5)]
    pub max_hits: usize,

    /// E-value threshold
    #[arg(long, default_value_t = 1e-5)]
    pub evalue: f64,

    /// Output directory
    #[arg(short, long, default_value = "blast_results")]
    pub output: PathBuf,

    /// Your email (required by NCBI)
    #[arg(long, env = "NCBI_EMAIL")]
    pub email: String,

    /// Force run even with >50 genes
    #[arg(long)]
    pub force: bool,
}

#[derive(Args)]
pub struct AlignArgs {
    /// Input files (FASTA, GFF, GBK)
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    /// Alignment method (blast, mmseqs, mummer, minimap2)
    #[arg(short, long, default_value = "mmseqs")]
    pub method: String,

    /// Sequence type (protein, nucleotide)
    #[arg(long, default_value = "protein")]
    pub seqtype: String,

    /// Output file
    #[arg(short, long)]
    pub output: PathBuf,

    /// Minimum identity threshold (%)
    #[arg(long, default_value_t = 30.0)]
    pub min_identity: f64,

    /// Minimum alignment length
    #[arg(long, default_value_t = 50)]
    pub min_length: usize,

    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Args)]
pub struct PhylogenyArgs {
    /// Input alignment file (FASTA, PHYLIP)
    pub input: PathBuf,

    /// Model specification (MFP for ModelFinder Plus)
    #[arg(short, long, default_value = "MFP")]
    pub model: String,

    /// Number of bootstrap replicates
    #[arg(short, long, default_value_t = 1000)]
    pub bootstrap: usize,

    /// Output file (Newick, PDF, SVG)
    #[arg(short, long)]
    pub output: PathBuf,

    /// Partition file for multi-gene analysis
    #[arg(long)]
    pub partition: Option<PathBuf>,

    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Args)]
pub struct CompareArgs {
    /// Input genome files
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    /// Specific chromosomes (format: file:chr)
    #[arg(long)]
    pub chromosomes: Option<Vec<String>>,

    /// Regions (format: file:chr:start-end)
    #[arg(long)]
    pub regions: Option<Vec<String>>,

    /// Alignment method
    #[arg(long, default_value = "mmseqs")]
    pub method: String,

    /// Build phylogenetic tree
    #[arg(long)]
    pub build_tree: bool,

    /// Output directory
    #[arg(short, long)]
    pub output: PathBuf,

    /// Plot type (linear, circular, matrix)
    #[arg(long, default_value = "linear")]
    pub plot_type: String,

    /// Generate HTML report
    #[arg(long)]
    pub html: bool,

    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Args)]
pub struct ConvertArgs {
    /// Input file
    pub input: PathBuf,

    /// Input format (auto-detect if not specified)
    #[arg(long)]
    pub from: Option<String>,

    /// Output format
    #[arg(short, long)]
    pub to: String,

    /// Output file (stdout if not specified)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Validate after conversion
    #[arg(long)]
    pub validate: bool,
}

#[derive(Args)]
pub struct CleanArgs {
    /// Input annotation file
    pub input: PathBuf,

    /// Output file (modifies in place if not specified)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Fix common errors automatically
    #[arg(long)]
    pub fix: bool,

    /// Remove orphaned features
    #[arg(long)]
    pub remove_orphans: bool,

    /// Use primary transcripts only
    #[arg(long)]
    pub primary_only: bool,
}

