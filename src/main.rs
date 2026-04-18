use anyhow::Result;
use std::env;
use std::path::{Path, PathBuf};

pub mod align;
pub mod annotate;
pub mod batch;
pub mod chat;
pub mod check;
pub mod cli;
pub mod compare;
pub mod convert;
pub mod fix;
pub mod install;
pub mod learn;
pub mod mask;
pub mod names;
pub mod parser;
pub mod phylogeny;
pub mod predict;
pub mod progress;
pub mod remote;
pub mod setup;
pub mod sort;
pub mod species;
pub mod stats;
pub mod submit;
pub mod train;
pub mod update;
pub mod utils;
pub mod y1000plus;

use parser::region::RegionSelector;

/// Returns true if `--help` or `-h` appears anywhere in the slice.
fn has_help_flag(args: &[String]) -> bool {
    args.iter().any(|a| a == "--help" || a == "-h")
}

fn print_banner() {
    // Colour / style codes — gracefully degrade on non-ANSI terminals
    let g = "\x1b[32m"; // green  (MYCONOTE)
    let c = "\x1b[36m"; // cyan   (_CLI)
    let b = "\x1b[1m"; // bold
    let y = "\x1b[33m"; // yellow (version)
    let d = "\x1b[2m"; // dim
    let r = "\x1b[0m"; // reset

    // ── MYCONOTE (green) ── and ── _CLI (cyan) ── side-by-side on one band ──
    let mc1 = "███╗   ███╗██╗   ██╗ ██████╗  ██████╗ ███╗   ██╗ ██████╗ ████████╗███████╗";
    let mc2 = "████╗ ████║╚██╗ ██╔╝██╔════╝ ██╔═══██╗████╗  ██║██╔═══██╗╚══██╔══╝██╔════╝";
    let mc3 = "██╔████╔██║ ╚████╔╝ ██║      ██║   ██║██╔██╗ ██║██║   ██║   ██║   █████╗  ";
    let mc4 = "██║╚██╔╝██║  ╚██╔╝  ██║      ██║   ██║██║╚██╗██║██║   ██║   ██║   ██╔══╝  ";
    let mc5 = "██║ ╚═╝ ██║   ██║   ╚██████╗ ╚██████╔╝██║ ╚████║╚██████╔╝   ██║   ███████╗";
    let mc6 = "╚═╝     ╚═╝   ╚═╝    ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝    ╚═╝   ╚══════╝";

    let cl1 = "        ██████╗ ██╗     ██╗";
    let cl2 = "       ██╔════╝ ██║     ██║";
    let cl3 = "       ██║      ██║     ██║";
    let cl4 = "       ██║      ██║     ██║";
    let cl5 = "───────╚██████╗ ███████╗██║";
    let cl6 = "        ╚═════╝ ╚══════╝╚═╝";

    println!();
    println!("{b}{g}  {mc1}{r}{b}{c}{cl1}{r}");
    println!("{b}{g}  {mc2}{r}{b}{c}{cl2}{r}");
    println!("{b}{g}  {mc3}{r}{b}{c}{cl3}{r}");
    println!("{b}{g}  {mc4}{r}{b}{c}{cl4}{r}");
    println!("{b}{g}  {mc5}{r}{b}{c}{cl5}{r}");
    println!("{b}{g}  {mc6}{r}{b}{c}{cl6}{r}");
    println!();
    println!("  {b}Genome Annotation Pipeline{r}");
    println!(
        "  {y}v{}{r}  ·  Hittinger Lab  ·  Laboratory of Genetics  ·  UW–Madison",
        env!("CARGO_PKG_VERSION")
    );
    println!("  {d}Benjamin Narh-Madey  ·  narhmadey@wisc.edu{r}");
    println!();
}

fn print_version() {
    print_banner();
}

fn print_main_help() {
    print_banner();
    println!("Usage: myconote-cli <command> [options]");
    println!("\nPipeline commands (run in order):");
    println!("  sort     Sort + rename genome contigs by length (pre-processing)");
    println!("  mask     Identify and soft-mask repeats in a genome FASTA");
    println!("  train    RNA-seq mediated training of Augustus/SNAP (Trinity + PASA)");
    println!("  predict  Predict genes (Augustus + SNAP + GlimmerHMM + GeneMark + EVM)");
    println!("  update   Refine gene models with RNA-seq evidence (PASA UTR extension)");
    println!("  annotate Functionally annotate genes (MMseqs2 + Pfam + EggNog + CAZyme + MEROPS + tRNAscan + ...)");
    println!("  submit   Prepare NCBI GenBank submission (validation + table2asn)");
    println!("  remote   Submit proteins to remote annotation servers (Phobius, InterProScan)");
    println!("  batch    Annotate multiple genomes (directory or sample sheet, HTCondor support)");
    println!("\nAnalysis commands:");
    println!("  stats    Calculate statistics from annotation files");
    println!("  phylogeny Build a maximum-likelihood tree with IQ-TREE (from an alignment)");
    println!("  compare  N-genome ortholog inference + pan-genome summary (OrthoFinder)");
    println!("  place    Place your genome in the Y1000+ 1,154-yeast reference (functional)");
    println!("  convert  Convert between genome annotation and sequence formats");
    println!("  clean    Validate and fix a GFF3 annotation file");
    println!("  fix      Repair errors in GenBank (.gbk) files");
    println!("\nUtility commands:");
    println!("  install  Install missing external tools via conda/mamba");
    println!("  check    Check which external tools are installed");
    println!("  setup    Download and index reference databases");
    println!("  species  List available trained Augustus species");
    println!("  learn    Interactive tutorial — learn myconote step by step (like R swirl)");
    println!("\nOptions for stats:");
    println!("  --format <json|csv|human>     Output format (default: human)");
    println!("  --taxon <group>               Taxonomic group for benchmarking");
    println!("  --chromosome <chr>            Focus on specific chromosome(s)");
    println!("  --region <chr:start-end>      Focus on specific region");
    println!("  --exclude <chr>               Exclude a chromosome");
    println!("  --primary-only                Collapse isoforms");
    println!("\nTaxonomic groups: fungi, ascomycota, basidiomycota, plants, animals, mammals");
    println!("\nExamples:");
    println!("  myconote-cli stats genome.gff3");
    println!("  myconote-cli annotate genes.gff3 --fasta genome.fa --kingdom fungi");
    println!("\nFor visualisation, pipe myconote outputs into external tools:");
    println!("  • Proksee (web)  — upload the GenBank from `convert --to genbank`");
    println!("  • IGV (desktop) — drop GFF3 + FASTA for interactive browsing");
    println!("  • clinker (pip) — cross-species gene-cluster synteny from .gbk files");
    println!("\nRun 'myconote-cli <command> --help' for detailed help on any command.");
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || matches!(args[1].as_str(), "--help" | "-h" | "help") {
        print_main_help();
        return Ok(());
    }

    if matches!(args[1].as_str(), "--version" | "-V" | "version") {
        print_version();
        return Ok(());
    }

    let command = &args[1];

    match command.as_str() {
        "sort" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli sort <genome.fa> [options]");
                println!("\nSorts contigs by length (longest first) and renames headers to");
                println!("clean sequential IDs (scaffold_1, scaffold_2, ...).  Run this");
                println!("before masking to ensure consistent identifiers throughout the pipeline.");
                println!("\nOptions:");
                println!("  --output <file>      Output FASTA (default: <input>_sorted.fa)");
                println!("  --prefix <str>       Contig ID prefix (default: scaffold)");
                println!("  --min-length <bp>    Discard contigs shorter than this (default: 0)");
                println!("  --rename-table <f>   Write old→new ID table to TSV file");
                println!("  --keep-desc          Keep original description after ID (default: strip)");
                println!("  --sort-by-name       Sort alphabetically instead of by length");
                println!("\nExamples:");
                println!("  myconote-cli sort assembly.fa --min-length 500");
                println!("  myconote-cli sort assembly.fa --prefix chr --rename-table rename.tsv");
                return Ok(());
            }
            let path = &args[2];
            handle_sort(path, &args[3..])?;
        }
        "train" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli train <masked_genome.fa> [options]");
                println!("\nRNA-seq mediated training of Augustus and SNAP gene predictors.");
                println!("Requires Trinity and PASA. Run after masking, before predicting.");
                println!("\nOptions (RNA-seq input — at least one required):");
                println!("  --left <file(s)>     R1 reads (comma-separated for multiple libs)");
                println!("  --right <file(s)>    R2 reads (paired, same order as --left)");
                println!("  --single <file(s)>   Single-end reads");
                println!("  --trinity <file>     Pre-assembled Trinity FASTA (skip Trinity step)");
                println!("\nOptions:");
                println!("  --species <name>     Augustus species name (default: myconote_trained)");
                println!("  --output <dir>       Output directory (default: train_out)");
                println!("  --max-intron <bp>    Maximum intron size (default: 3000)");
                println!("  --min-models <n>     Minimum training models required (default: 200)");
                println!("  --strand <RF|FR>     RNA-seq strand specificity (default: unstranded)");
                println!("  --memory <str>       Trinity memory (default: 50G)");
                println!("  --no-snap            Skip SNAP training");
                println!("  --threads <n>        Threads (default: 4)");
                println!("\nExamples:");
                println!("  myconote-cli train genome_masked.fa --left R1.fastq --right R2.fastq --species myorg");
                println!("  myconote-cli train genome_masked.fa --trinity trinity.fasta --species myorg");
                return Ok(());
            }
            let path = &args[2];
            handle_train(path, &args[3..])?;
        }
        "mask" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli mask <genome.fa> [options]");
                println!("\nOptions:");
                println!("  --output <masked.fa>        Output masked FASTA (default: <input>_masked.fa)");
                println!("  --engine <engine>           Masking engine (default: self)");
                println!("  --species <name>            RepeatMasker species (e.g. fungi, arabidopsis)");
                println!("  --repeat-lib <file>         Custom RepeatMasker library FASTA");
                println!("  --hard-mask                 Use N instead of lowercase (hard-masking)");
                println!("  --min-length <bp>           Minimum repeat length to mask (default: 200)");
                println!("  --threads <n>               Threads (default: 4)");
                println!("\nEngines:");
                println!("  self         minimap2 self-alignment + native TRF (no database needed)");
                println!("  repeatmasker RepeatMasker with RepBase species library");
                println!("  both         RepeatMasker + self-alignment merged");
                println!("  repeatmodeler  De novo library (RepeatModeler) → RepeatMasker [RECOMMENDED]");
                println!("  full         RepeatModeler + RepeatMasker + self-align (most thorough)");
                println!("\nExamples:");
                println!("  myconote-cli mask genome.fa");
                println!("  myconote-cli mask genome.fa --engine repeatmasker --species fungi");
                println!("  myconote-cli mask genome.fa --engine both --threads 8");
                return Ok(());
            }
            let path = &args[2];
            handle_mask(path, &args[3..])?;
        }
        "predict" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli predict <masked.fa> [options]");
                println!("\nOptions:");
                println!("  --output <dir>              Output directory (default: predict_out)");
                println!("  --kingdom <kingdom>         fungi|plant|animal|insect|protist (default: fungi)");
                println!("  --species <augustus-species> Override Augustus species model");
                println!("  --no-snap                   Disable SNAP (use only Augustus)");
                println!("  --snap-hmm <hmm>            SNAP HMM model (default: auto from kingdom)");
                println!("  --protein-evidence <blast>  BLAST tabular fmt6 protein hits for hints");
                println!("  --locus-prefix <prefix>     Gene ID prefix (default: GENE)");
                println!("  --train                     Auto-train Augustus from first-pass prediction");
                println!("  --train-species <name>      Species name for trained model (default: <prefix>_trained)");
                println!("  --glimmerhmm                Also run GlimmerHMM (adds third ab initio predictor)");
                println!("  --glimmer-dir <dir>         GlimmerHMM training directory (default: auto-detect)");
                println!("  --genemark                  Also run GeneMark-ES (self-training, no model needed)");
                println!("  --genemark-hints <gff>      Use GeneMark-ET with RNA-seq intron hints");
                println!("  --threads <n>               Threads (default: 4)");
                println!("\nKingdoms and their default Augustus species:");
                println!("  fungi    saccharomyces_cerevisiae_S288C");
                println!("  plant    arabidopsis");
                println!("  animal   human");
                println!("  insect   fly");
                println!("  protist  toxoplasma");
                println!("\nExamples:");
                println!("  myconote-cli predict genome_masked.fa --kingdom fungi --locus-prefix AFUB");
                println!("  myconote-cli predict genome_masked.fa --kingdom plant --species arabidopsis");
                println!("  myconote-cli predict genome_masked.fa --protein-evidence blastx.tsv");
                return Ok(());
            }
            let path = &args[2];
            handle_predict(path, &args[3..])?;
        }
        "annotate" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli annotate <genes.gff3> --fasta <genome.fa> [options]");
                println!("\nOptions:");
                println!("  --fasta <file>              Genome FASTA (required)");
                println!("  --output <dir>              Output directory (default: annotate_out)");
                println!("  --kingdom <kingdom>         fungi|plant|animal|insect|protist");
                println!("  --organism <name>           Organism name for report");
                println!("  --locus-prefix <prefix>     Must match prefix used in predict");
                println!("  --swissprot-db <dir>        Path to Swiss-Prot MMseqs2 database");
                println!("  --pfam-db <file>            Path to Pfam-A.hmm");
                println!("  --db-dir <dir>              Database directory (default: ~/.myconote/dbs)");
                println!("  --evalue <float>            E-value cutoff (default: 1e-5)");
                println!("  --min-identity <float>      Min sequence identity for MMseqs2 (default: 0.3)");
                println!("  --no-pfam                   Skip Pfam domain search");
                println!("  --no-busco                  Skip BUSCO completeness check");
                println!("  --interproscan              Run InterProScan via EBI REST API");
                println!("  --email <address>           Email for InterProScan (required by EBI)");
                println!("  --eggnog                    Run EggNog-mapper (COG/NOG categories)");
                println!("  --eggnog-db <dir>           EggNog-mapper database directory");
                println!("  --eggnog-results <file>     Pre-computed emapper.annotations file");
                println!("  --cazyme                    Run CAZyme annotation (dbCAN / DIAMOND)");
                println!("  --cazyme-db <file>          dbCAN DIAMOND database (.dmnd)");
                println!("  --secretome                 Run secretome prediction (DeepSig + DeepTMHMM)");
                println!("  --signalp-organism <str>    Signal-peptide organism: euk|gram+|gram- (default: euk)");
                println!("  --antismash                 Run antiSMASH BGC cluster prediction");
                println!("  --antismash-dir <dir>       Pre-computed antiSMASH output directory");
                println!("  --antismash-taxon <str>     antiSMASH taxon: fungi|bacteria|plants (default: fungi)");
                println!("  --merops                    Run MEROPS protease annotation (DIAMOND vs merops.dmnd)");
                println!("  --merops-db <file>          MEROPS DIAMOND database (default: auto-detect from db-dir)");
                println!("  --trnascan                  Run tRNAscan-SE for tRNA gene prediction");
                println!("  --trnascan-mode <mode>      tRNAscan mode: eukaryotic|mitochondrial|general (default: eukaryotic)");
                println!("  --genetic-code <n>          Translation table (1=standard, 12=Candida CTG, etc.)");
                println!("  --threads <n>               Threads (default: 4)");
                println!("  --download-dbs              Download Swiss-Prot and Pfam databases");
                println!("\nOutputs:");
                println!("  annotated.gff3              GFF3 with all annotations in attributes");
                println!("  annotations.tsv             Full annotation table (all sources)");
                println!("  proteins.fa                 Extracted protein sequences");
                println!("  annotation_report.txt       Human-readable summary");
                println!("  eggnog_hits.tsv             EggNog COG annotations (if --eggnog)");
                println!("  cazyme_hits.tsv             CAZyme family assignments (if --cazyme)");
                println!("  secretome.tsv               Secreted protein predictions (if --secretome)");
                println!("  bgc_clusters.tsv            BGC cluster table (if --antismash)");
                println!("  merops_hits.tsv             MEROPS protease families (if --merops)");
                println!("  busco/                      BUSCO completeness results");
                println!("\nExamples:");
                println!("  myconote-cli annotate --download-dbs");
                println!("  myconote-cli annotate predict_out/consensus.gff3 --fasta genome.fa");
                println!("  myconote-cli annotate genes.gff3 --fasta genome.fa --eggnog --cazyme --secretome --antismash --merops");
                return Ok(());
            }
            if args[2] == "--download-dbs" {
                handle_annotate_download_dbs(&args[3..])?;
            } else {
                let path = &args[2];
                handle_annotate(path, &args[3..])?;
            }
        }
        "stats" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli stats <annotation.gff3> [options]");
                println!("\nCalculates comprehensive statistics from a GFF3 annotation file.");
                println!("\nOptions:");
                println!("  --format <json|csv|human>     Output format (default: human)");
                println!("  --taxon <group>               Taxonomic group for benchmarking context");
                println!("  --chromosome <chr>            Focus on one or more specific sequences");
                println!("  --region <chr:start-end>      Focus on a specific genomic region");
                println!("  --exclude <chr>               Exclude a sequence from stats");
                println!("  --primary-only                Collapse isoforms to primary transcript only");
                println!("  --benchmark y1000plus         Append the Y1000+ reference distributions");
                println!("                                (requires `setup --y1000plus --preset starter`)");
                println!("  --annotated <path.tsv>        Combined with --benchmark: place the user's");
                println!("                                KEGG-KO count on the 1,154-yeast percentile.");
                println!("                                Accepts the TSV emitted by `myconote-cli annotate`.");
                println!("\nTaxonomic groups: fungi, ascomycota, basidiomycota, plants, animals, mammals");
                println!("\nOutputs (human format):");
                println!("  Total features, genes, transcripts, CDS, exons");
                println!("  Mean/median/min/max gene length, N50");
                println!("  Per-chromosome breakdown");
                println!("\nExamples:");
                println!("  myconote-cli stats genes.gff3");
                println!("  myconote-cli stats genes.gff3 --taxon fungi");
                println!("  myconote-cli stats genes.gff3 --format json > stats.json");
                println!("  myconote-cli stats genes.gff3 --chromosome scaffold_1");
                println!("  myconote-cli stats genes.gff3 --primary-only --format csv");
                return Ok(());
            }
            let path = &args[2];
            handle_stats(path, &args[3..])?;
        }
        "phylogeny" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli phylogeny <alignment.fa> [options]");
                println!("\nInfers a maximum-likelihood phylogenetic tree using IQ-TREE 2.");
                println!("\nOptions:");
                println!("  --model <model>       Substitution model, or MFP for auto-select (default: MFP)");
                println!("  --bootstrap <n>       Ultrafast bootstrap replicates (default: 1000, 0 to skip)");
                println!("  --partition <file>    Partition file for multi-gene/multi-locus analysis");
                println!("  --threads <n>         Threads (default: 4)");
                println!("  --prefix <name>       Output prefix (default: alignment filename stem)");
                println!("  --force               Overwrite an existing <prefix>.treefile (default: skip)");
                println!("\nOutput:");
                println!("  <prefix>.treefile     Best-fit ML tree in Newick format");
                println!("  <prefix>.iqtree       Full IQ-TREE run report");
                println!("  <prefix>.log          Console log");
                println!("\nExamples:");
                println!("  myconote-cli phylogeny aligned_cds.fa");
                println!("  myconote-cli phylogeny aligned_cds.fa --model GTR+G --bootstrap 1000");
                println!("  myconote-cli phylogeny multilocus.fa --partition parts.txt --threads 8");
                return Ok(());
            }
            let input = &args[2];
            let mut config = phylogeny::PhylogenyConfig::default();
            let mut i = 3usize;
            while i < args.len() {
                match args[i].as_str() {
                    "--model" | "-m" if i + 1 < args.len() => {
                        config.model = args[i + 1].clone(); i += 2;
                    }
                    "--bootstrap" | "-b" if i + 1 < args.len() => {
                        if let Ok(n) = args[i + 1].parse::<usize>() { config.bootstrap = n; }
                        i += 2;
                    }
                    "--partition" if i + 1 < args.len() => {
                        config.partition = Some(std::path::PathBuf::from(&args[i + 1])); i += 2;
                    }
                    "--threads" | "-T" if i + 1 < args.len() => {
                        if let Ok(n) = args[i + 1].parse::<usize>() { config.threads = n; }
                        i += 2;
                    }
                    "--prefix" if i + 1 < args.len() => {
                        config.prefix = Some(args[i + 1].clone()); i += 2;
                    }
                    "--force" => {
                        config.force = true; i += 1;
                    }
                    _ => { i += 1; }
                }
            }
            phylogeny::build_tree(input, &config)
                .map_err(|e| anyhow::anyhow!("{}", e))?;
        }
        "compare" => {
            if args.len() < 4 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli compare <g1.gff3> <g1.fa> <g2.gff3> <g2.fa> [...] [options]");
                println!("\nN-genome ortholog inference via OrthoFinder (Emms & Kelly 2019).");
                println!("Produces an ortholog table, pan-genome summary (core / soft-core / shell / cloud),");
                println!("rooted species tree, and a set of single-copy orthogroups ready for `phylogeny`.");
                println!("\nOptions:");
                println!("  --output <dir>         Output directory (default: compare_out)");
                println!("  --threads <n>          Threads for OrthoFinder (default: all cores)");
                println!("  --sensitive            Use diamond_ultra_sens search (default: on)");
                println!("  --fast                 Use default diamond search (faster, less accurate)");
                println!("  --msa                  MSA-based tree refinement (2–3× slower)");
                println!("  --species-tree         Also build a concatenated supermatrix alignment");
                println!("                         for downstream `phylogeny` (MAFFT + NEXUS partitions)");
                println!("  --genetic-code <n>     NCBI translation table (default: 1)");
                println!("  --soft-core <frac>     Soft-core threshold fraction (default: 0.95)");
                println!("  --cloud <frac>         Cloud upper bound fraction (default: 0.15)");
                println!("\nGenome-count caps (auto-detected from protein count):");
                println!("  Small  (≤15 000 proteins/genome, e.g. fungi): cap = 5");
                println!("  Medium (15–30 k, e.g. small plants):          cap = 3");
                println!("  Large  (>30 k, e.g. crops):                    cap = 2");
                println!("\nRequires `orthofinder` on PATH:");
                println!("  conda install -c bioconda orthofinder");
                println!("\nExamples:");
                println!("  myconote-cli compare s1.gff3 s1.fa s2.gff3 s2.fa s3.gff3 s3.fa");
                println!("  myconote-cli compare a.gff3 a.fa b.gff3 b.fa --threads 8 --msa");
                return Ok(());
            }
            handle_compare(&args[2..])?;
        }
        "place" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli place --annotated <annotated.tsv> [options]");
                println!("\nPlaces your genome against the Y1000+ bundle of 1,154 yeasts by");
                println!("comparing functional profiles (KEGG-KO Jaccard similarity).");
                println!("\nOptions:");
                println!("  --annotated <file.tsv>   myconote `annotate` output (required)");
                println!("  --top <n>                How many closest species to report (default: 10)");
                println!("  --format <human|tsv>     Output format (default: human)");
                println!("\nDependencies:");
                println!("  Requires the `kegg` subset of the Y1000+ bundle:");
                println!("    myconote-cli setup --y1000plus --include kegg");
                println!("\nExamples:");
                println!("  myconote-cli place --annotated annotated.tsv");
                println!("  myconote-cli place --annotated annotated.tsv --top 20 --format tsv");
                println!("\nCitation: Opulente DA et al. (2024). Science 384(6694): eadj4503.");
                return Ok(());
            }
            handle_place(&args[2..])?;
        }
        "convert" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli convert <input> --to <format> [options]");
                println!("\nGFF3 conversions (input: .gff3):");
                println!("  --to gtf          → GTF (GENCODE/Ensembl gene transfer format)");
                println!("  --to bed          → BED6 (simple browser track)");
                println!("  --to bed12        → BED12 (exon-block structure per transcript)");
                println!("  --to bedgraph     → BEDGraph (per-feature coverage depth)");
                println!("  --to table        → TSV feature table (attributes expanded)");
                println!("  --to protein      → Protein FASTA (.faa)  [requires --fasta]");
                println!("  --to genbank      → GenBank (.gbk)        [requires --fasta]");
                println!("\nSequence conversions (input: .fasta/.fa/.fastq/.fq/.aln/.phy/.nex):");
                println!("  --to fastq        → FASTQ  (FASTA → FASTQ, dummy quality scores)");
                println!("  --to fasta        → FASTA  (FASTQ → FASTA, strip quality)");
                println!("  --to table        → TSV    (id, length, sequence)");
                println!("  --to phylip       → PHYLIP alignment format");
                println!("  --to nexus        → NEXUS  alignment format");
                println!("  --to clustal      → CLUSTAL alignment format");
                println!("\nVCF conversions (input: .vcf):");
                println!("  --to bed          → BED (variant positions)");
                println!("  --to table        → TSV (all fields expanded)");
                println!("  --to consensus    → Consensus FASTA  [requires --fasta]");
                println!("  --to annovar      → ANNOVAR input format");
                println!("  --to maf          → MAF (Mutation Annotation Format)");
                println!("\nOptions:");
                println!("  --to <format>           Target format (required)");
                println!("  --output <file>  -o     Output file (default: auto-named)");
                println!("  --fasta <file>          Reference FASTA (for protein/genbank/consensus)");
                println!("  --organism <name>       Organism name (for genbank output)");
                println!("  --sample <name>         Tumor sample name (for maf output)");
                println!("  --qual <char>           Default quality character for FASTQ (default: I)");
                println!("\nExamples:");
                println!("  myconote-cli convert genes.gff3 --to gtf");
                println!("  myconote-cli convert genes.gff3 --to protein --fasta genome.fa");
                println!("  myconote-cli convert reads.fastq --to fasta -o reads.fa");
                println!("  myconote-cli convert align.aln --to phylip");
                println!("  myconote-cli convert variants.vcf --to table");
                return Ok(());
            }
            let path = &args[2];
            handle_convert(path, &args[3..])?;
        }
        "clean" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli clean <annotations.gff3> [options]");
                println!("\nOptions:");
                println!("  --output <file.gff3>    Output file (default: <input>_clean.gff3)");
                println!("  --fix-coords            Fix off-by-one coordinate errors");
                println!("  --remove-orphans        Remove features with missing parents");
                println!("  --min-length <bp>       Remove features shorter than this (default: 1)");
                println!("\nExamples:");
                println!("  myconote-cli clean genes.gff3");
                println!("  myconote-cli clean genes.gff3 --remove-orphans --output clean.gff3");
                return Ok(());
            }
            let path = &args[2];
            handle_clean(path, &args[3..])?;
        }
        "update" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli update <genes.gff3> --fasta <genome.fa> [options]");
                println!("\nRefines initial gene models by adding UTRs and correcting boundaries");
                println!("using RNA-seq transcript evidence (PASA or lightweight fallback).");
                println!("\nOptions:");
                println!("  --fasta <file>            Genome FASTA (required)");
                println!("  --output <dir>            Output directory (default: update_out)");
                println!("  --transcripts <file>      Pre-assembled transcript FASTA (Trinity output)");
                println!("  --rna-r1 <file>           R1 reads (Trinity will assemble them)");
                println!("  --rna-r2 <file>           R2 reads (paired)");
                println!("  --rna-bam <file>          Pre-aligned BAM (skip alignment step)");
                println!("  --locus-prefix <str>      Must match prefix used in predict (default: GENE)");
                println!("  --organism <name>         Organism name for PASA config");
                println!("  --max-utr-ext <bp>        Maximum UTR extension (default: 2000)");
                println!("  --min-identity <float>    Min transcript alignment identity (default: 0.95)");
                println!("  --threads <n>             Threads (default: 4)");
                println!("\nOutputs:");
                println!("  updated.gff3              Refined gene models with UTR features");
                println!("\nExamples:");
                println!("  myconote-cli update predict_out/consensus.gff3 --fasta genome.fa --rna-bam rnaseq.bam");
                println!("  myconote-cli update genes.gff3 --fasta genome.fa --transcripts trinity.fasta");
                return Ok(());
            }
            let path = &args[2];
            handle_update(path, &args[3..])?;
        }
        "install" => {
            if has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli install [command] [options]");
                println!("\nChecks which external tools are missing and installs them");
                println!("automatically via conda or mamba.");
                println!("\nOptions:");
                println!("  [command]    Only install tools needed by this pipeline step");
                println!("               e.g. predict, mask, train, annotate");
                println!("  --yes  -y    Skip the confirmation prompt");
                println!("  --mamba      Force mamba instead of conda");
                println!("\nExamples:");
                println!("  myconote-cli install              # install all missing tools");
                println!("  myconote-cli install predict      # only tools for predict");
                println!("  myconote-cli install --yes        # non-interactive install");
                return Ok(());
            }

            let mut opts = install::InstallOptions {
                filter_cmd: None,
                yes:        false,
                mamba:      false,
            };
            let mut i = 2usize;
            while i < args.len() {
                match args[i].as_str() {
                    "--yes" | "-y"  => { opts.yes   = true; i += 1; }
                    "--mamba"       => { opts.mamba  = true; i += 1; }
                    other if !other.starts_with('-') => {
                        opts.filter_cmd = Some(other.to_string()); i += 1;
                    }
                    _ => { i += 1; }
                }
            }
            install::run_install(opts);
        }
        "check" => {
            if has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli check [command_name]");
                println!("\nChecks whether each external tool required by myconote-cli is");
                println!("installed and available in PATH, and reports its version.");
                println!("\nOptions:");
                println!("  [command_name]   Filter to tools used by a specific command");
                println!("                   e.g. 'annotate', 'predict', 'mask'");
                println!("\nTools checked: augustus, snap, glimmerhmm, genemark, evm, trinity,");
                println!("  minimap2, samtools, diamond, hmmscan, mmseqs, pasa, busco,");
                println!("  deepsig, biolib (DeepTMHMM), deeploc, repeatmasker, repeatmodeler,");
                println!("  antismash, eggnog-mapper, dbcan, interproscan, phobius, wget, curl");
                println!("\nExamples:");
                println!("  myconote-cli check");
                println!("  myconote-cli check predict");
                println!("  myconote-cli check annotate");
                return Ok(());
            }
            let filter = args.iter().skip(2)
                .find(|a| !a.starts_with('-'))
                .map(|s| s.as_str());
            check::run_check(filter);
        }
        "setup" => {
            if has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli setup [options]");
                println!("\nDownloads and indexes the reference databases used by myconote-cli.");
                println!("Databases are stored in ~/.myconote/dbs by default.");
                println!("\nOptions:");
                println!("  --list                   Show all databases and their download status");
                println!("  --check                  Verify existing databases are intact");
                println!("  --db <name> [name ...]   Download only specific databases");
                println!("  --db-dir <dir>           Custom database directory (default: ~/.myconote/dbs)");
                println!("  --force                  Re-download even if the database already exists");
                println!("\nAvailable databases:");
                println!("  swiss-prot   UniProt/Swiss-Prot (MMseqs2 indexed) — used by annotate");
                println!("  pfam         Pfam-A HMM profiles (hmmpress indexed) — used by annotate");
                println!("  eggnog       EggNog-mapper database (COG/NOG) — used by annotate --eggnog");
                println!("  dbcan        CAZyme DIAMOND database — used by annotate --cazyme");
                println!("  merops       MEROPS protease DIAMOND database — used by annotate --merops");
                println!("  busco        BUSCO fungi lineage data — used by annotate");
                println!("  chat-corpus  Q1 open-access paper corpus — used by explain");
                println!("  ollama       Ollama LLM runtime + model — used by explain");
                println!("\nY1000+ yeast reference bundle (Opulente et al. 2024, Science):");
                println!("  --y1000plus                  Enter the Y1000+ subsystem");
                println!("  --y1000plus --list           Show all subsets + what's installed");
                println!("  --y1000plus --dry-run --preset starter   Preview the ~175 MB starter bundle");
                println!("  --y1000plus --include kegg,busco         Install just these subsets");
                println!("  --y1000plus --preset phylogeny           Enables phylogenetic placement");
                println!("  --y1000plus --uninstall <subsets>        Remove specific subsets");
                println!("\nExamples:");
                println!("  myconote-cli setup --list");
                println!("  myconote-cli setup                     # download everything");
                println!("  myconote-cli setup --db swiss-prot pfam");
                println!("  myconote-cli setup --check");
                println!("  myconote-cli setup --y1000plus --list");
                return Ok(());
            }
            handle_setup(&args[2..])?;
        }
        "fix" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli fix <input.gbk> [options]");
                println!("\nDetects and repairs common errors in GenBank (.gbk/.gb) files.");
                println!("\nOptions:");
                println!("  --output <file>  -o    Output file (default: <input>_fixed.gbk)");
                println!("  --report <file>        Write a repair report to this file");
                println!("  --dry-run              Report problems without writing output");
                println!("  --no-fix-tags          Don't renumber duplicate locus tags");
                println!("  --no-fix-product       Don't add missing /product qualifiers");
                println!("  --no-fix-stops         Don't remove internal stop codons from /translation");
                println!("\nProblems fixed:");
                println!("  • Duplicate /locus_tag values (renumbered with suffix _2, _3, ...)");
                println!("  • Missing /product qualifier (filled with 'hypothetical protein')");
                println!("  • Internal stop codons in /translation (removed)");
                println!("  • Invalid characters in /locus_tag and /gene (replaced with _)");
                println!("  • Invalid /codon_start values (reset to 1)");
                println!("  • Unclosed records missing // terminator (added)");
                println!("\nExamples:");
                println!("  myconote-cli fix annotation.gbk -o annotation_fixed.gbk");
                println!("  myconote-cli fix annotation.gbk --dry-run --report fix_report.txt");
                return Ok(());
            }
            let path = &args[2];
            handle_fix(path, &args[3..])?;
        }
        "remote" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli remote <proteins.fa> [options]");
                println!("\nSubmits protein sequences to remote annotation servers.");
                println!("\nOptions:");
                println!("  --output <dir>          Output directory (default: remote_out)");
                println!("  --phobius               Run Phobius (signal peptide + TM topology)");
                println!("  --interproscan          Run InterProScan (domain/family search)");
                println!("  --deeploc               Run DeepLoc 2 (subcellular localisation)");
                println!("  --email <address>       Email for InterProScan (required by EBI)");
                println!("  --batch-size <n>        Proteins per batch (default: 100)");
                println!("  --poll-interval <secs>  Seconds between status polls (default: 60)");
                println!("  --max-retries <n>       Max poll attempts per batch (default: 60)");
                println!("  --threads <n>           Threads (default: 4)");
                println!("\nOutputs:");
                println!("  remote_annotations.tsv  Merged results from all services");
                println!("  remote_cache/           Cached results (re-run is fast)");
                println!("\nNote: Remote services require an internet connection.");
                println!("      Results are cached — interrupted runs can be resumed.");
                println!("\nExamples:");
                println!("  myconote-cli remote annotate_out/proteins.fa --phobius");
                println!("  myconote-cli remote proteins.fa --interproscan --email me@uni.edu");
                return Ok(());
            }
            let path = &args[2];
            handle_remote(path, &args[3..])?;
        }
        "species" => {
            if has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli species [filter] [options]");
                println!("\nLists Augustus species models available on this system.");
                println!("Scans AUGUSTUS_CONFIG_PATH, ~/.myconote/augustus_config, and");
                println!("common install locations for trained species profiles.");
                println!("\nOptions:");
                println!("  [filter]           Show only species whose name contains this string");
                println!("  --filter <str>     Same as positional filter");
                println!("  --grouped  -g      Group species by kingdom (fungi/plant/animal/other)");
                println!("  --list     -l      List all species (same as no filter)");
                println!("\nColumns shown:");
                println!("  Name         Augustus species identifier (used with --species flag)");
                println!("  Complete     Whether all required HMM files are present");
                println!("  Trained by   'user' for custom-trained, 'reference' for built-in");
                println!("\nExamples:");
                println!("  myconote-cli species");
                println!("  myconote-cli species fungi");
                println!("  myconote-cli species --grouped");
                println!("  myconote-cli species --filter aspergillus");
                return Ok(());
            }

            let mut _list_all   = false;
            let mut grouped     = false;
            let mut filter: Option<String> = None;

            let mut i = 2usize;
            while i < args.len() {
                match args[i].as_str() {
                    "--list" | "-l"   => { _list_all = true; i += 1; }
                    "--grouped" | "-g"=> { grouped  = true; i += 1; }
                    "--filter" | "--search" if i + 1 < args.len() => {
                        filter = Some(args[i+1].clone()); i += 2;
                    }
                    other if !other.starts_with('-') => {
                        filter = Some(other.to_string()); i += 1;
                    }
                    _ => { i += 1; }
                }
            }

            if grouped {
                species::list_species_grouped();
            } else {
                species::list_species(filter.as_deref());
            }
        }
        "learn" | "tutorial" | "swirl" => {
            learn::run_learn(&args[2..]).map_err(|e| anyhow::anyhow!("{}", e))?;
        }
        "explain" => {
            chat::run_explain(&args[2..]).map_err(|e| anyhow::anyhow!("{}", e))?;
        }
        "batch" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli batch <genomes_dir|sample_sheet.tsv> [options]");
                println!("\nAnnotate multiple genomes in one command. Accepts a directory of");
                println!("FASTA files or a sample sheet (TSV) with per-genome settings.");
                println!("\nSample sheet columns (TSV, header required):");
                println!("  name    fasta    kingdom    species    genetic_code    locus_prefix");
                println!("  (only 'fasta' is required; others use defaults)");
                println!("\nPipeline options:");
                println!("  --output <dir>          Batch output directory (default: batch_out)");
                println!("  --stages <list>         Comma-separated stages (default: sort,mask,predict,annotate,submit)");
                println!("  --kingdom <k>           Default kingdom (default: fungi)");
                println!("  --threads <n>           Threads per genome (default: 4)");
                println!("  --parallel <n>          Max genomes in parallel (default: 2)");
                println!("  --min-length <bp>       Min contig length for sort (default: 500)");
                println!("  --mask-engine <engine>  Masking engine (default: repeatmodeler)");
                println!("  --genetic-code <n>      Default translation table (default: 1)");
                println!("  --locus-prefix <str>    Default locus prefix (default: GENE)");
                println!("  --resume <dir>          Resume a previous batch run");
                println!("\nHTCondor options:");
                println!("  --condor                Generate HTCondor submit files (don't run locally)");
                println!("  --condor-cpus <n>       CPUs per job (default: 8)");
                println!("  --condor-mem <size>     Memory per job (default: 32G)");
                println!("  --condor-disk <size>    Disk per job (default: 50G)");
                println!("  --condor-queue <name>   HTCondor accounting group");
                println!("  --condor-extra <file>   Extra submit directives to append");
                println!("\nExamples:");
                println!("  myconote-cli batch genomes/");
                println!("  myconote-cli batch samples.tsv --threads 8 --parallel 4");
                println!("  myconote-cli batch genomes/ --condor --condor-mem 64G");
                println!("  myconote-cli batch genomes/ --stages sort,mask,predict");
                println!("  myconote-cli batch --resume batch_out/");
                return Ok(());
            }
            handle_batch(&args[2..])?;
        }
        "submit" => {
            if args.len() < 3 || has_help_flag(&args[2..]) {
                println!("Usage: myconote-cli submit <annotated.gff3> --fasta <genome.fa> [options]");
                println!("\nPrepares genome annotations for NCBI GenBank submission.");
                println!("\nOptions:");
                println!("  --fasta <file>              Genome FASTA (required)");
                println!("  --output <dir>              Output directory (default: submit_out)");
                println!("  --organism <name>           Organism name (required)");
                println!("  --strain <name>             Strain name");
                println!("  --bioproject <acc>          BioProject accession");
                println!("  --biosample <acc>           BioSample accession");
                println!("  --locus-prefix <str>        Locus tag prefix (default: MYCO)");
                println!("  --genetic-code <n>          Translation table (default: 1)");
                println!("  --email <address>           Contact email");
                println!("  --validate-only             Only validate, do not generate files");
                println!("\nOutputs:");
                println!("  annotation.tbl              NCBI feature table");
                println!("  annotation.fsa              Genome FASTA copy");
                println!("  annotation.sqn              Sequin file (if table2asn available)");
                println!("  template.sbt                Submission template");
                println!("\nExamples:");
                println!("  myconote-cli submit genes.gff3 --fasta genome.fa --organism 'Aspergillus niger'");
                println!("  myconote-cli submit genes.gff3 --fasta genome.fa --organism 'Candida albicans' --genetic-code 12");
                return Ok(());
            }
            let path = &args[2];
            handle_submit(path, &args[3..])?;
        }
        _ => println!("Unknown command: {}. Try: sort | mask | train | predict | update | annotate | submit | batch | explain | remote | stats | phylogeny | compare | place | convert | clean | fix | install | check | setup | species | learn", command),
    }

    Ok(())
}

fn handle_stats(path: &str, args: &[String]) -> Result<()> {
    // Parse options
    let mut format = "human".to_string();
    let mut taxon = None;
    let mut chromosomes = Vec::new();
    let mut regions = Vec::new();
    let mut exclude = Vec::new();
    let mut primary_only = false;
    let mut benchmark_y1000plus = false;
    let mut annotated_path: Option<PathBuf> = None;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--format" if i + 1 < args.len() => {
                format = args[i + 1].clone();
                i += 2;
            }
            "--taxon" if i + 1 < args.len() => {
                taxon = Some(args[i + 1].clone());
                i += 2;
            }
            "--chromosome" if i + 1 < args.len() => {
                chromosomes.push(args[i + 1].clone());
                i += 2;
            }
            "--region" if i + 1 < args.len() => {
                regions.push(args[i + 1].clone());
                i += 2;
            }
            "--exclude" if i + 1 < args.len() => {
                exclude.push(args[i + 1].clone());
                i += 2;
            }
            "--primary-only" => {
                primary_only = true;
                i += 1;
            }
            "--benchmark" if i + 1 < args.len() => {
                if args[i + 1].eq_ignore_ascii_case("y1000plus") {
                    benchmark_y1000plus = true;
                } else {
                    eprintln!(
                        "Unknown --benchmark target '{}'. Only 'y1000plus' is supported.",
                        args[i + 1]
                    );
                }
                i += 2;
            }
            "--annotated" if i + 1 < args.len() => {
                annotated_path = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            _ => i += 1,
        }
    }

    // Build region selector
    let mut selector = if !regions.is_empty() {
        RegionSelector::with_regions(&regions)?
    } else if !chromosomes.is_empty() {
        RegionSelector::with_chromosomes(&chromosomes)
    } else {
        RegionSelector::new()
    };

    for chr in &exclude {
        selector.exclude_chromosome(chr);
    }

    // Calculate stats
    let stats = stats::GenomeStatistics::from_gff_with_selector(path, &selector, primary_only)?;

    match format.as_str() {
        "json" => {
            let json = serde_json::to_string_pretty(&stats)?;
            println!("{}", json);
        }
        "csv" => {
            println!("Category,Value");
            println!("Total Features,{}", stats.total_features);
            println!("Genes,{}", stats.total_genes);
            println!("Transcripts,{}", stats.total_transcripts);
            println!("CDS,{}", stats.total_cds);
            println!("Exons,{}", stats.total_exons);
            println!("Mean Gene Length,{:.2}", stats.mean_gene_length());
            println!("Median Gene Length,{:.2}", stats.median_gene_length());
            println!("Min Gene Length,{}", stats.min_gene_length());
            println!("Max Gene Length,{}", stats.max_gene_length());
            println!("N50,{}", stats.n50());

            for (chr, chr_stats) in &stats.chromosome_stats {
                println!("Chromosome {}-Genes,{}", chr, chr_stats.gene_count);
                println!(
                    "Chromosome {}-Transcripts,{}",
                    chr, chr_stats.transcript_count
                );
                println!("Chromosome {}-CDS,{}", chr, chr_stats.cds_count);
                println!("Chromosome {}-Exons,{}", chr, chr_stats.exon_count);
            }
        }
        _ => {
            if !chromosomes.is_empty() {
                println!("\n🔍 Focusing on chromosomes: {}", chromosomes.join(", "));
            }
            if !regions.is_empty() {
                println!("\n🔍 Focusing on regions: {}", regions.join(", "));
            }
            if !exclude.is_empty() {
                println!("\n🔍 Excluding chromosomes: {}", exclude.join(", "));
            }
            if primary_only {
                println!("\n🔍 Showing PRIMARY transcripts only");
            }

            stats.print_summary();

            if let Some(t) = taxon {
                if let Some(warning) = stats.get_taxon_warning(&t) {
                    println!("\n{}", warning);
                }
            }

            if benchmark_y1000plus {
                print_y1000plus_benchmark(&stats, annotated_path.as_deref())?;
            }
        }
    }

    Ok(())
}

/// Append a Y1000+ reference-distribution block to the human-format stats
/// output. Pulls from the installed bundle; silently skips distributions
/// whose subset isn't on disk. If `annotated_tsv` is provided, reads the
/// user's distinct KEGG KO count from it and places them on the percentile.
fn print_y1000plus_benchmark(
    stats: &stats::GenomeStatistics,
    annotated_tsv: Option<&Path>,
) -> Result<()> {
    use y1000plus::benchmark::{count_user_kos, load_reference, DistStats};

    let reference = match load_reference() {
        Ok(r) => r,
        Err(e) => {
            eprintln!("\n⚠  --benchmark y1000plus skipped: {e}");
            return Ok(());
        }
    };

    // Optional: count user KOs from their annotated TSV for percentile rank.
    let user_ko_count: Option<f64> = annotated_tsv.and_then(|p| match count_user_kos(p) {
        Ok((kos, rows)) => {
            println!(
                "   ↪ scanned {} gene rows from {}; {} distinct KEGG KOs",
                rows,
                p.display(),
                kos
            );
            Some(kos as f64)
        }
        Err(e) => {
            eprintln!("   ⚠  couldn't parse --annotated: {e}");
            None
        }
    });

    println!("\n══════════════════════════════════════════════════════════════════");
    println!(
        "📊 Y1000+ benchmark  ({} species · Opulente et al. 2024, Science)",
        reference.species_count
    );
    println!("══════════════════════════════════════════════════════════════════");

    let row = |label: &str, d: &DistStats, user: Option<f64>, fmt: &str| {
        let user_col = match user {
            Some(u) => {
                let pct = d.percentile_of(u);
                let marker = if pct < 25.0 {
                    "below p25 — lean"
                } else if pct < 75.0 {
                    "typical range"
                } else {
                    "above p75 — rich"
                };
                if fmt == "int" {
                    format!("{:>7} ({:>5.1} pct · {})", u as u64, pct, marker)
                } else {
                    format!("{:>7.1}% ({:>5.1} pct · {})", u, pct, marker)
                }
            }
            None => "(not measured)".to_string(),
        };
        let pretty = if fmt == "int" {
            format!(
                "{:<22} min {:>6.0}   p25 {:>6.0}   median {:>6.0}   p75 {:>6.0}   max {:>6.0}",
                label, d.min, d.p25, d.median, d.p75, d.max
            )
        } else {
            format!(
                "{:<22} min {:>5.1}%  p25 {:>5.1}%  median {:>5.1}%  p75 {:>5.1}%  max {:>5.1}%",
                label, d.min, d.p25, d.median, d.p75, d.max
            )
        };
        println!("{pretty}\n  your genome: {user_col}");
    };

    let _ = stats; // gene-count benchmarking arrives with the `annotations` subset

    if let Some(ref d) = reference.busco_completeness {
        row("BUSCO completeness", d, None, "pct");
        println!(
            "  (run BUSCO on your genome to see where you land; we surface \
             the Y1000+ distribution so you can interpret it in context.)"
        );
    }
    if let Some(ref d) = reference.kegg_ko_count {
        row("Distinct KEGG KOs", d, user_ko_count, "int");
        if user_ko_count.is_none() {
            println!(
                "  (pass --annotated <path.tsv> pointing at the output of \
                 `myconote-cli annotate` to see your percentile.)"
            );
        }
    }

    if reference.busco_completeness.is_none() && reference.kegg_ko_count.is_none() {
        println!(
            "No usable reference subsets installed. Install with:\n  \
             myconote-cli setup --y1000plus --include busco,kegg"
        );
    }

    println!();
    println!(
        "Reference cache: {}\nCitation: Opulente DA et al. (2024). Science 384(6694): eadj4503.",
        reference.source_root.display()
    );
    Ok(())
}

fn handle_compare(args: &[String]) -> Result<()> {
    use compare::{parse_positional_inputs, run_compare, CompareConfig};

    let mut config = CompareConfig::default();

    // Split positional args (before the first `--flag`) from option args.
    let mut positional: Vec<String> = Vec::new();
    let mut i = 0usize;
    while i < args.len() {
        if args[i].starts_with("--") {
            break;
        }
        positional.push(args[i].clone());
        i += 1;
    }

    // Parse option args starting from where positional parsing stopped.
    while i < args.len() {
        match args[i].as_str() {
            "--output" | "-o" if i + 1 < args.len() => {
                config.output_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--threads" | "-t" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.threads = n.max(1);
                }
                i += 2;
            }
            "--sensitive" => {
                config.sensitive = true;
                i += 1;
            }
            "--fast" => {
                config.sensitive = false;
                i += 1;
            }
            "--msa" => {
                config.msa = true;
                i += 1;
            }
            "--species-tree" => {
                config.species_tree = true;
                i += 1;
            }
            "--no-primary-only" => {
                config.primary_only = false;
                i += 1;
            }
            "--primary-only" => {
                config.primary_only = true;
                i += 1;
            }
            "--genetic-code" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<u8>() {
                    config.genetic_code = n;
                }
                i += 2;
            }
            "--soft-core" if i + 1 < args.len() => {
                if let Ok(f) = args[i + 1].parse::<f64>() {
                    config.soft_core_frac = f;
                }
                i += 2;
            }
            "--cloud" if i + 1 < args.len() => {
                if let Ok(f) = args[i + 1].parse::<f64>() {
                    config.cloud_frac = f;
                }
                i += 2;
            }
            "--force-cap" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.force_cap = Some(n);
                }
                i += 2;
            }
            _ => {
                eprintln!("Unknown compare option: {}", args[i]);
                i += 1;
            }
        }
    }

    config.inputs = parse_positional_inputs(&positional)?;
    run_compare(&config)?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// place command: species placement against the Y1000+ bundle
// ─────────────────────────────────────────────────────────────────────────────

fn handle_place(args: &[String]) -> Result<()> {
    use y1000plus::place::{place_functional, PlaceOptions};

    let mut annotated: Option<PathBuf> = None;
    let mut top_n: usize = 10;
    let mut format = "human".to_string();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--annotated" if i + 1 < args.len() => {
                annotated = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--top" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    top_n = n;
                }
                i += 2;
            }
            "--format" if i + 1 < args.len() => {
                format = args[i + 1].clone();
                i += 2;
            }
            _ => i += 1,
        }
    }

    let annotated = annotated.ok_or_else(|| {
        anyhow::anyhow!(
            "--annotated <annotated.tsv> is required. Run `myconote-cli annotate` first."
        )
    })?;

    println!(
        "🧬 Placing your genome against Y1000+ (1,154 yeasts, Opulente et al. 2024, Science)…"
    );
    let report = place_functional(&PlaceOptions {
        annotated_tsv: &annotated,
        top_n: top_n.max(1),
    })?;

    match format.as_str() {
        "tsv" => {
            println!("species\tjaccard\tshared_kos\tref_kos\tuser_kos");
            for p in &report.top {
                println!(
                    "{}\t{:.4}\t{}\t{}\t{}",
                    p.species, p.jaccard, p.shared, p.ref_kos, report.user_kos,
                );
            }
            if let Some(ref c) = report.top_codon_table {
                println!(
                    "# top_codon_table\tspecies={}\tctg_is_ser={}\tdeviations={}",
                    c.top_species,
                    c.ctg_is_ser,
                    c.deviations.len()
                );
            }
            if let Some(ref m) = report.metabolic_prediction {
                println!(
                    "# metabolic_prediction\tcarbon={}\tcarbon_conf={:.3}\tnitrogen={}\tnitrogen_conf={:.3}\tn_neighbours={}",
                    m.carbon_vote.label,
                    m.carbon_vote.confidence,
                    m.nitrogen_vote.label,
                    m.nitrogen_vote.confidence,
                    m.neighbours.len(),
                );
            }
            if let Some(ref t) = report.thermotolerance {
                println!(
                    "# thermotolerance_37C\tlabel={}\tconfidence={:.3}\tn_neighbours={}",
                    t.vote.label,
                    t.vote.confidence,
                    t.neighbours.len(),
                );
            }
        }
        _ => {
            println!();
            println!(
                "  user KEGG-KO set : {} distinct KOs across {} gene rows",
                report.user_kos, report.user_gene_rows
            );
            println!(
                "  reference space  : {} species with KEGG annotations",
                report.species_total
            );
            println!();
            println!(
                "{:<45}  {:>8}  {:>8}  {:>8}",
                "Closest Y1000+ species", "Jaccard", "shared", "ref KOs"
            );
            println!("{}", "─".repeat(80));
            for (i, p) in report.top.iter().enumerate() {
                let badge = if i == 0 { " ⭐" } else { "   " };
                println!(
                    "{badge}{:<42}  {:>8.4}  {:>8}  {:>8}",
                    prettify_species(&p.species),
                    p.jaccard,
                    p.shared,
                    p.ref_kos,
                );
            }
            // Metabolic lifestyle prediction — rendered when the `metabolism`
            // subset is installed and at least one neighbour was classified.
            if let Some(ref m) = report.metabolic_prediction {
                println!();
                println!("── Predicted metabolic lifestyle ──");
                println!(
                    "  carbon  : {:<12}  confidence {:.0}%  (weighted vote across {} neighbours)",
                    m.carbon_vote.label,
                    m.carbon_vote.confidence * 100.0,
                    m.neighbours.len()
                );
                println!(
                    "  nitrogen: {:<12}  confidence {:.0}%",
                    m.nitrogen_vote.label,
                    m.nitrogen_vote.confidence * 100.0
                );
                if m.carbon_vote.tallies.len() > 1 || m.nitrogen_vote.tallies.len() > 1 {
                    println!("  neighbour classifications:");
                    for n in &m.neighbours {
                        println!(
                            "    {:<36}  J={:.4}  C={:<10} N={}",
                            n.species, n.jaccard, n.carbon_class, n.nitrogen_class
                        );
                    }
                }
            }

            // Thermotolerance prediction — rendered when the `phenotypes`
            // subset is installed and at least one neighbour had a growth-at-37 label.
            if let Some(ref t) = report.thermotolerance {
                use y1000plus::place::thermo_label_description;
                println!();
                println!("── Predicted thermotolerance (growth at 37 °C) ──");
                println!(
                    "  prediction : {}  ({})  confidence {:.0}%",
                    t.vote.label,
                    thermo_label_description(&t.vote.label),
                    t.vote.confidence * 100.0
                );
                if t.vote.tallies.len() > 1 {
                    println!("  vote tallies:");
                    for (label, w) in &t.vote.tallies {
                        println!(
                            "    {:<2} ({}): weight {:.4}",
                            label,
                            thermo_label_description(label),
                            w
                        );
                    }
                }
                println!(
                    "  neighbours consulted: {}/{} of top-N had phenotype data",
                    t.neighbours.len(),
                    report.top.len()
                );
            }

            // Codon-table advice — only rendered when the `codontable` subset
            // is installed. Warns about CTG-Ser and other non-standard codes.
            if let Some(ref c) = report.top_codon_table {
                println!();
                println!("── Inferred genetic code (codetta, top hit) ──");
                println!(
                    "  {:<42}  {} chars",
                    prettify_species(&c.top_species),
                    c.code.len()
                );
                if c.deviations.is_empty() {
                    println!("  Standard NCBI genetic code — no reassignments detected.");
                } else {
                    println!(
                        "  {} codon reassignment(s) vs NCBI standard:",
                        c.deviations.len()
                    );
                    for (codon, std, inferred) in c.deviations.iter().take(8) {
                        println!("    {codon}: {std} → {inferred}");
                    }
                    if c.deviations.len() > 8 {
                        println!("    … and {} more", c.deviations.len() - 8);
                    }
                }
                if c.ctg_is_ser {
                    println!(
                        "  ⚠  CTG-Ser clade detected — your closest yeast translates CTG as Serine,"
                    );
                    println!(
                        "     not Leucine. Re-run `predict` / `convert --to protein` with the"
                    );
                    println!("     alternative yeast mitochondrial code (NCBI transl_table=12).");
                }
            }
            println!();
            println!(
                "Note: functional placement via KEGG-KO Jaccard. True phylogenetic placement\n\
                 (1,403 marker genes + EPA-ng) will land once the `phylogeny-place` subset\n\
                 and EPA-ng integration are wired."
            );
            println!("Cite: Opulente DA et al. (2024). Science 384(6694): eadj4503.");
        }
    }

    Ok(())
}

fn prettify_species(raw: &str) -> String {
    // Per-species filenames in the Y1000+ kegg subset look like
    // `candida_tropicalis.txt` after our walk yields `candida_tropicalis`.
    // Render `Candida tropicalis` for humans, but leave underscores in the
    // stem alone for readability (some names contain strain suffixes).
    let mut parts = raw.splitn(2, '_');
    match (parts.next(), parts.next()) {
        (Some(genus), Some(rest)) if !genus.is_empty() => {
            let mut g = genus.chars();
            let upper = g.next().map(|c| c.to_ascii_uppercase()).unwrap_or_default();
            let tail: String = g.collect();
            format!("{upper}{tail} {rest}")
        }
        _ => raw.to_string(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// convert command: multi-format file conversion
// ─────────────────────────────────────────────────────────────────────────────

fn handle_convert(input: &str, args: &[String]) -> Result<()> {
    use std::path::Path;

    // ── Parse options ──────────────────────────────────────────────────────
    let mut to_format: Option<String> = None;
    let mut output_path: Option<String> = None;
    let mut fasta_path: Option<String> = None;
    let mut organism: Option<String> = None;
    let mut sample: Option<String> = None;
    let mut qual_char: char = 'I';

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--to" if i + 1 < args.len() => {
                to_format = Some(args[i + 1].to_lowercase());
                i += 2;
            }
            "--output" | "-o" if i + 1 < args.len() => {
                output_path = Some(args[i + 1].clone());
                i += 2;
            }
            "--fasta" if i + 1 < args.len() => {
                fasta_path = Some(args[i + 1].clone());
                i += 2;
            }
            "--organism" if i + 1 < args.len() => {
                organism = Some(args[i + 1].clone());
                i += 2;
            }
            "--sample" if i + 1 < args.len() => {
                sample = Some(args[i + 1].clone());
                i += 2;
            }
            "--qual" if i + 1 < args.len() => {
                qual_char = args[i + 1].chars().next().unwrap_or('I');
                i += 2;
            }
            _ => i += 1,
        }
    }

    let to = match to_format {
        Some(ref f) => f.as_str(),
        None => {
            eprintln!("Error: --to <format> is required.");
            eprintln!("Run 'myconote-cli convert' with no arguments to see all formats.");
            return Ok(());
        }
    };

    let input_path = Path::new(input);
    let ext = input_path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

    // ── Auto-derive output path if not specified ───────────────────────────
    let out_path: PathBuf = match output_path {
        Some(ref p) => PathBuf::from(p),
        None => {
            let new_ext = match to {
                "gtf" => "gtf",
                "bed" | "bed6" | "bed12" | "bedgraph" => "bed",
                "table" => "tsv",
                "protein" => "faa",
                "genbank" => "gbk",
                "fastq" => "fastq",
                "fasta" => "fasta",
                "phylip" => "phy",
                "nexus" => "nex",
                "clustal" => "aln",
                "consensus" => "fasta",
                "annovar" => "avinput",
                "maf" => "maf",
                other => other,
            };
            input_path.with_extension(new_ext)
        }
    };

    println!("Converting: {} → {} ({})", input, to, out_path.display());

    // ── Dispatch by input type + target format ─────────────────────────────

    // GFF3 conversions
    if matches!(ext.as_str(), "gff3" | "gff") {
        let n = match to {
            "gtf" => convert::gff3_to_gtf(input_path, &out_path)?,
            "bed" | "bed6" => convert::gff3_to_bed(input_path, &out_path, &[])?,
            "bed12" => convert::gff3_to_bed12(input_path, &out_path)?,
            "bedgraph" => convert::gff3_to_bedgraph(input_path, &out_path, "gene")?,
            "table" | "tsv" => convert::gff3_to_table(input_path, &out_path)?,
            "protein" | "faa" => {
                let fa = require_fasta(&fasta_path, to)?;
                convert::gff3_to_protein(input_path, Path::new(&fa), &out_path)?
            }
            "genbank" | "gbk" => {
                handle_convert_genbank(
                    input,
                    fasta_path.as_deref(),
                    &out_path,
                    organism.as_deref(),
                )?;
                return Ok(());
            }
            other => {
                eprintln!(
                    "Error: unsupported target format '{}' for GFF3 input.",
                    other
                );
                eprintln!("Supported: gtf, bed, bed12, bedgraph, table, protein, genbank");
                return Ok(());
            }
        };
        println!("✓ {} records written → {}", n, out_path.display());
        return Ok(());
    }

    // VCF conversions
    if ext == "vcf" {
        let n = match to {
            "bed" => convert::vcf_to_bed(input_path, &out_path)?,
            "table" | "tsv" => convert::vcf_to_table(input_path, &out_path)?,
            "consensus" => {
                let fa = require_fasta(&fasta_path, to)?;
                convert::vcf_to_consensus(input_path, Path::new(&fa), &out_path)?
            }
            "annovar" => convert::vcf_to_annovar(input_path, &out_path)?,
            "maf" => {
                let samp = sample.as_deref().unwrap_or("TUMOR");
                convert::vcf_to_maf(input_path, &out_path, samp)?
            }
            other => {
                eprintln!(
                    "Error: unsupported target format '{}' for VCF input.",
                    other
                );
                eprintln!("Supported: bed, table, consensus, annovar, maf");
                return Ok(());
            }
        };
        println!("✓ {} records written → {}", n, out_path.display());
        return Ok(());
    }

    // Sequence / alignment conversions
    if matches!(ext.as_str(), "fasta" | "fa" | "fna" | "faa" | "fas") {
        let n = match to {
            "fastq" | "fq" => convert::fasta_to_fastq(input_path, &out_path, qual_char)?,
            "table" | "tsv" => convert::fasta_to_table(input_path, &out_path)?,
            other => {
                eprintln!(
                    "Error: unsupported target format '{}' for FASTA input.",
                    other
                );
                eprintln!("Supported: fastq, table");
                return Ok(());
            }
        };
        println!("✓ {} sequences written → {}", n, out_path.display());
        return Ok(());
    }

    if matches!(ext.as_str(), "fastq" | "fq") {
        let n = match to {
            "fasta" | "fa" => convert::fastq_to_fasta(input_path, &out_path)?,
            other => {
                eprintln!(
                    "Error: unsupported target format '{}' for FASTQ input.",
                    other
                );
                eprintln!("Supported: fasta");
                return Ok(());
            }
        };
        println!("✓ {} sequences written → {}", n, out_path.display());
        return Ok(());
    }

    // Alignment format conversions
    if matches!(
        ext.as_str(),
        "aln" | "phy" | "nex" | "nxs" | "nexus" | "phylip" | "clustal"
    ) {
        if !matches!(
            to,
            "phylip" | "phy" | "nexus" | "nex" | "clustal" | "aln" | "fasta" | "fa"
        ) {
            eprintln!(
                "Error: unsupported target format '{}' for alignment input.",
                to
            );
            eprintln!("Supported: phylip, nexus, clustal, fasta");
            return Ok(());
        }
        let n = convert::convert_alignment(input_path, &out_path, &ext, to)?;
        println!("✓ {} sequences written → {}", n, out_path.display());
        return Ok(());
    }

    eprintln!("Error: unrecognised input file extension '.{}'.", ext);
    eprintln!(
        "Supported input types: .gff3, .gff, .vcf, .fasta, .fa, .fastq, .fq, .aln, .phy, .nex"
    );
    Ok(())
}

/// Require --fasta flag; emit a helpful error if missing.
fn require_fasta(fasta_path: &Option<String>, to: &str) -> Result<String> {
    match fasta_path {
        Some(p) => Ok(p.clone()),
        None => {
            eprintln!(
                "Error: --fasta <reference.fa> is required for '{}' output.",
                to
            );
            Err(anyhow::anyhow!("missing --fasta argument"))
        }
    }
}

/// Inner helper for GFF3 → GenBank (keeps the detailed seqid-mismatch reporting).
fn handle_convert_genbank(
    gff_path: &str,
    fasta_path: Option<&str>,
    out_path: &std::path::Path,
    organism: Option<&str>,
) -> Result<()> {
    use parser::genbank::write_genbank;
    use parser::{read_fasta_index, GFFReader};
    use std::fs::File;
    use std::io::BufWriter;

    let fasta_path = match fasta_path {
        Some(p) => p,
        None => {
            eprintln!("Error: --fasta <genome.fa> is required for GenBank output.");
            return Ok(());
        }
    };

    println!("  GFF3:     {}", gff_path);
    println!("  FASTA:    {}", fasta_path);
    println!("  Output:   {}", out_path.display());
    if let Some(org) = organism {
        println!("  Organism: {}", org);
    }
    println!();

    print!("Reading GFF3 annotations... ");
    let records: Vec<_> = GFFReader::from_path(gff_path)?
        .filter_map(|r| r.ok())
        .collect();
    println!("{} records", records.len());

    print!("Indexing FASTA sequences...  ");
    let fasta_index = read_fasta_index(fasta_path)?;
    println!("{} sequences", fasta_index.len());

    let gff_seqids: std::collections::BTreeSet<&str> =
        records.iter().map(|r| r.seqid.as_str()).collect();
    let mut missing = 0usize;
    for seqid in &gff_seqids {
        if !fasta_index.contains_key(*seqid) {
            eprintln!("  ⚠  GFF3 seqid '{}' has no matching FASTA sequence", seqid);
            missing += 1;
        }
    }
    if missing > 0 {
        eprintln!("\n  {} seqid(s) will be skipped.", missing);
        eprintln!("  Tip: seqid in GFF3 col 1 must match the first token of the FASTA '>' header.");
        eprintln!("  E.g. '>NODE_1 length=12345 ...' → seqid must be 'NODE_1'\n");
    }

    print!("Writing GenBank file...      ");
    let out_file = File::create(out_path)?;
    let mut writer = BufWriter::new(out_file);
    write_genbank(&mut writer, &records, &fasta_index, organism)?;
    println!("done");

    println!("\n✓ GenBank output: {}", out_path.display());
    println!("  Open with: Geneious, Benchling, SnapGene, or BioPython SeqIO");
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// clean command: validate and fix GFF3
// ─────────────────────────────────────────────────────────────────────────────

fn handle_clean(gff_path: &str, args: &[String]) -> Result<()> {
    use parser::gff::GFFRecord;
    use std::collections::HashSet;
    use std::io::{BufRead, BufReader, Write as IoWrite};

    let mut output_path = String::new();
    let mut fix_coords = false;
    let mut remove_orphans = false;
    let mut min_length: u64 = 1;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--output" | "-o" if i + 1 < args.len() => {
                output_path = args[i + 1].clone();
                i += 2;
            }
            "--fix-coords" => {
                fix_coords = true;
                i += 1;
            }
            "--remove-orphans" => {
                remove_orphans = true;
                i += 1;
            }
            "--min-length" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<u64>() {
                    min_length = n;
                }
                i += 2;
            }
            _ => i += 1,
        }
    }

    if output_path.is_empty() {
        let p = PathBuf::from(gff_path);
        let stem = p.file_stem().unwrap_or_default().to_string_lossy();
        let ext = p.extension().unwrap_or_default().to_string_lossy();
        let dir = p.parent().unwrap_or(std::path::Path::new("."));
        output_path = format!(
            "{}_clean.{}",
            dir.join(&*stem).display(),
            if ext.is_empty() {
                "gff3".to_string()
            } else {
                ext.to_string()
            }
        );
    }

    println!("Cleaning GFF3: {}", gff_path);

    // ── Pass 1: lenient read to collect all IDs ───────────────────────────
    // We do a lenient line-by-line parse that tries each data line, fixing
    // swapped coordinates before rejecting the record.
    let raw_file = std::fs::File::open(gff_path)?;
    let raw_reader = BufReader::new(raw_file);

    let mut all_records: Vec<GFFRecord> = Vec::new();
    let mut parse_errors = 0usize;

    for (line_num, line_res) in raw_reader.lines().enumerate() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // Attempt normal parse first
        match GFFRecord::from_line(trimmed, line_num + 1) {
            Ok(rec) => all_records.push(rec),
            Err(_) if fix_coords => {
                // Try swapping columns 3 and 4 (start/end) and re-parse
                let fields: Vec<&str> = trimmed.splitn(9, '\t').collect();
                if fields.len() == 9 {
                    let swapped = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        fields[0],
                        fields[1],
                        fields[2],
                        fields[4],
                        fields[3], // swap start/end
                        fields[5],
                        fields[6],
                        fields[7],
                        fields[8]
                    );
                    match GFFRecord::from_line(&swapped, line_num + 1) {
                        Ok(rec) => {
                            eprintln!(
                                "  🔧 Fixed swapped coords at line {} ({}..{} → {}..{})",
                                line_num + 1,
                                fields[3],
                                fields[4],
                                fields[4],
                                fields[3]
                            );
                            all_records.push(rec);
                        }
                        Err(_) => {
                            eprintln!("  ⚠  Skipping unparseable line {}", line_num + 1);
                            parse_errors += 1;
                        }
                    }
                }
            }
            Err(_) => {
                eprintln!("  ⚠  Skipping unparseable line {}", line_num + 1);
                parse_errors += 1;
            }
        }
    }

    let known_ids: HashSet<String> = all_records
        .iter()
        .filter_map(|r| r.id().map(|s| s.clone()))
        .collect();

    // ── Pass 2: filter and write ──────────────────────────────────────────
    let mut kept = 0usize;
    let mut removed = 0usize;
    let mut fixed = 0usize;

    let out_file = std::fs::File::create(&output_path)?;
    let mut writer = std::io::BufWriter::new(out_file);
    writeln!(writer, "##gff-version 3")?;

    for mut rec in all_records {
        // Min-length filter
        if rec.length() < min_length {
            removed += 1;
            continue;
        }

        // Orphan filter
        if remove_orphans {
            if let Some(parent) = rec.parent() {
                if !known_ids.contains(parent) {
                    eprintln!(
                        "  ✂  Removing orphan {} (parent '{}' not found)",
                        rec.id().map(|s| s.as_str()).unwrap_or("?"),
                        parent
                    );
                    removed += 1;
                    continue;
                }
            }
        }

        // Clamp start to 1 (GFF3 is 1-based)
        if fix_coords && rec.start == 0 {
            rec.start = 1;
            fixed += 1;
        }

        writeln!(writer, "{}", rec.to_gff3_line())?;
        kept += 1;
    }

    println!("\n  Features kept:    {}", kept);
    println!("  Features removed: {}", removed);
    if parse_errors > 0 {
        println!("  Unparseable lines skipped: {}", parse_errors);
    }
    if fix_coords {
        println!("  Coords fixed:     {}", fixed);
    }
    println!("\n✓ Clean GFF3: {}", output_path);

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// mask command: repeat identification and soft-masking
// ─────────────────────────────────────────────────────────────────────────────

fn handle_mask(fasta_path: &str, args: &[String]) -> Result<()> {
    use mask::{run_masking, MaskConfig, MaskEngine};

    let mut config = MaskConfig {
        input: std::path::PathBuf::from(fasta_path),
        ..MaskConfig::default()
    };

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--output" | "-o" if i + 1 < args.len() => {
                config.output = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--engine" if i + 1 < args.len() => {
                config.engine = MaskEngine::from_str(&args[i + 1]);
                i += 2;
            }
            "--species" if i + 1 < args.len() => {
                config.species = Some(args[i + 1].clone());
                i += 2;
            }
            "--repeat-lib" if i + 1 < args.len() => {
                config.repeat_lib = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--hard-mask" => {
                config.hard_mask = true;
                i += 1;
            }
            "--min-length" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.min_length = n;
                }
                i += 2;
            }
            "--threads" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.threads = n;
                }
                i += 2;
            }
            _ => i += 1,
        }
    }

    // Default output path
    if config.output == PathBuf::new() {
        let p = PathBuf::from(fasta_path);
        let stem = p.file_stem().unwrap_or_default().to_string_lossy();
        let ext = p
            .extension()
            .map(|e| format!(".{}", e.to_string_lossy()))
            .unwrap_or_default();
        let dir = p.parent().unwrap_or(std::path::Path::new("."));
        config.output = dir.join(format!("{}_masked{}", stem, ext));
    }

    let stats = run_masking(&config).map_err(|e| anyhow::anyhow!("{}", e))?;

    println!("\n── Masking summary ──────────────────────────────────────────");
    println!("  Total bases   : {}", stats.total_bases);
    println!(
        "  Masked bases  : {} ({:.1}%)",
        stats.masked_bases,
        stats.percent_masked()
    );
    println!("  Repeat regions: {}", stats.repeat_regions);
    println!("  Output FASTA  : {}", config.output.display());

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// predict command: gene prediction pipeline
// ─────────────────────────────────────────────────────────────────────────────

fn handle_predict(fasta_path: &str, args: &[String]) -> Result<()> {
    use predict::kingdom::Kingdom;
    use predict::{run_prediction, PredictConfig};

    let mut config = PredictConfig {
        masked_fasta: PathBuf::from(fasta_path),
        ..PredictConfig::default()
    };

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--output" | "-o" if i + 1 < args.len() => {
                config.out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--kingdom" if i + 1 < args.len() => {
                config.kingdom = Kingdom::from_str(&args[i + 1]);
                i += 2;
            }
            "--species" if i + 1 < args.len() => {
                config.augustus_species = Some(args[i + 1].clone());
                i += 2;
            }
            "--no-snap" => {
                config.use_snap = false;
                i += 1;
            }
            "--snap-hmm" if i + 1 < args.len() => {
                config.snap_hmm = Some(args[i + 1].clone());
                i += 2;
            }
            "--protein-evidence" if i + 1 < args.len() => {
                config.protein_evidence = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--locus-prefix" if i + 1 < args.len() => {
                config.locus_prefix = args[i + 1].clone();
                i += 2;
            }
            "--train" => {
                config.self_train = true;
                i += 1;
            }
            "--train-species" if i + 1 < args.len() => {
                config.train_species = Some(args[i + 1].clone());
                i += 2;
            }
            "--glimmerhmm" | "--glimmer" => {
                config.use_glimmerhmm = true;
                i += 1;
            }
            "--glimmer-dir" if i + 1 < args.len() => {
                config.glimmer_dir = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--genemark" => {
                config.use_genemark = true;
                i += 1;
            }
            "--genemark-hints" if i + 1 < args.len() => {
                config.genemark_hints = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--protein-fasta" if i + 1 < args.len() => {
                config.protein_fasta = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--max-intron" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.max_intron = n;
                }
                i += 2;
            }
            "--ploidy" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<u8>() {
                    config.ploidy = Some(n);
                }
                i += 2;
            }
            "--weights" if i + 1 < args.len() => {
                config.weights_file = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--threads" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.threads = n;
                }
                i += 2;
            }
            _ => i += 1,
        }
    }

    let (_gff_path, gene_count) = run_prediction(&config).map_err(|e| anyhow::anyhow!("{}", e))?;

    println!("\n── Prediction complete ──────────────────────────────────────");
    println!(
        "  {} genes in {}/consensus.gff3",
        gene_count,
        config.out_dir.display()
    );
    println!(
        "\n  Next step: myconote-cli annotate {}/consensus.gff3 \\",
        config.out_dir.display()
    );
    println!("               --fasta {} \\", fasta_path);
    println!(
        "               --kingdom {} \\",
        format!("{:?}", config.kingdom).to_lowercase()
    );
    println!("               --locus-prefix {} \\", config.locus_prefix);
    println!("               --interproscan --email <your@email.com>  # optional, very thorough");

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// annotate command: functional annotation pipeline
// ─────────────────────────────────────────────────────────────────────────────

fn handle_annotate(gff_path: &str, args: &[String]) -> Result<()> {
    use annotate::{run_annotation, AnnotateConfig};
    use predict::kingdom::Kingdom;

    let mut config = AnnotateConfig {
        gff: PathBuf::from(gff_path),
        ..AnnotateConfig::default()
    };

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--fasta" if i + 1 < args.len() => {
                config.fasta = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--output" | "-o" if i + 1 < args.len() => {
                config.out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--kingdom" if i + 1 < args.len() => {
                config.kingdom = Kingdom::from_str(&args[i + 1]);
                i += 2;
            }
            "--organism" if i + 1 < args.len() => {
                config.organism = Some(args[i + 1].clone());
                i += 2;
            }
            "--locus-prefix" if i + 1 < args.len() => {
                config.locus_prefix = args[i + 1].clone();
                i += 2;
            }
            "--swissprot-db" if i + 1 < args.len() => {
                config.swissprot_db = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--pfam-db" if i + 1 < args.len() => {
                config.pfam_db = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--db-dir" if i + 1 < args.len() => {
                config.db_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--evalue" if i + 1 < args.len() => {
                if let Ok(e) = args[i + 1].parse::<f64>() {
                    config.evalue = e;
                }
                i += 2;
            }
            "--min-identity" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse::<f64>() {
                    config.min_identity = v;
                }
                i += 2;
            }
            "--no-pfam" => {
                config.run_pfam = false;
                i += 1;
            }
            "--no-busco" => {
                config.run_busco = false;
                i += 1;
            }
            "--interproscan" => {
                config.run_interproscan = true;
                i += 1;
            }
            "--email" if i + 1 < args.len() => {
                config.interproscan_email = args[i + 1].clone();
                i += 2;
            }
            // ── EggNog ────────────────────────────────────────────────────
            "--eggnog" | "--eggnog-mapper" => {
                config.run_eggnog = true;
                i += 1;
            }
            "--eggnog-db" if i + 1 < args.len() => {
                config.eggnog_db = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--eggnog-results" if i + 1 < args.len() => {
                config.eggnog_results = Some(PathBuf::from(&args[i + 1]));
                config.run_eggnog = true;
                i += 2;
            }
            // ── CAZymes ───────────────────────────────────────────────────
            "--cazyme" | "--cazymes" => {
                config.run_cazyme = true;
                i += 1;
            }
            "--cazyme-db" if i + 1 < args.len() => {
                config.cazyme_db = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            // ── Secretome ─────────────────────────────────────────────────
            "--secretome" | "--signalp" => {
                config.run_secretome = true;
                i += 1;
            }
            "--signalp-organism" if i + 1 < args.len() => {
                config.signalp_organism = args[i + 1].clone();
                i += 2;
            }
            // ── antiSMASH ─────────────────────────────────────────────────
            "--antismash" => {
                config.run_antismash = true;
                i += 1;
            }
            "--antismash-dir" if i + 1 < args.len() => {
                config.antismash_dir = Some(PathBuf::from(&args[i + 1]));
                config.run_antismash = true;
                i += 2;
            }
            "--antismash-taxon" if i + 1 < args.len() => {
                config.antismash_taxon = args[i + 1].clone();
                i += 2;
            }
            "--merops" => {
                config.run_merops = true;
                i += 1;
            }
            "--merops-db" if i + 1 < args.len() => {
                config.merops_db = Some(PathBuf::from(&args[i + 1]));
                config.run_merops = true;
                i += 2;
            }
            // ── tRNAscan-SE ───────────────────────────────────────────────
            "--trnascan" | "--trna" => {
                config.run_trnascan = true;
                i += 1;
            }
            "--trnascan-mode" if i + 1 < args.len() => {
                config.trnascan_mode = args[i + 1].clone();
                i += 2;
            }
            // ── Genetic code ──────────────────────────────────────────────
            "--genetic-code" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<u8>() {
                    config.genetic_code = n;
                }
                i += 2;
            }
            "--threads" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.threads = n;
                }
                i += 2;
            }
            _ => i += 1,
        }
    }

    if config.fasta == PathBuf::new() {
        eprintln!("Error: --fasta <genome.fa> is required");
        eprintln!("Usage: myconote-cli annotate <genes.gff3> --fasta <genome.fa>");
        return Ok(());
    }

    // Show database status before starting
    let status = annotate::db::check_status(&config.db_dir);
    status.print();

    let results = run_annotation(&config).map_err(|e| anyhow::anyhow!("{}", e))?;

    println!("\n── Annotation complete ──────────────────────────────────────");
    println!(
        "  {:.1}% of genes functionally annotated",
        results.annotated_fraction() * 100.0
    );
    if let Some(ref b) = results.busco_summary {
        println!(
            "  BUSCO completeness: {:.1}% ({})",
            b.percent_complete(),
            b.lineage
        );
    }
    println!("  Output directory: {}", config.out_dir.display());

    Ok(())
}

fn handle_annotate_download_dbs(args: &[String]) -> Result<()> {
    use annotate::db::download_all;

    let home = std::env::var("HOME").unwrap_or_else(|_| "/tmp".to_string());
    let mut db_dir = PathBuf::from(home).join(".myconote").join("dbs");
    let mut threads = 4usize;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--db-dir" if i + 1 < args.len() => {
                db_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--threads" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    threads = n;
                }
                i += 2;
            }
            _ => i += 1,
        }
    }

    download_all(&db_dir, threads).map_err(|e| anyhow::anyhow!("{}", e))?;

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// sort command: sort + rename genome contigs
// ─────────────────────────────────────────────────────────────────────────────

fn handle_sort(input: &str, args: &[String]) -> Result<()> {
    use sort::{run_sort, SortConfig};

    let mut config = SortConfig {
        input: PathBuf::from(input),
        ..SortConfig::default()
    };

    // Default output: <stem>_sorted.<ext>
    let p = PathBuf::from(input);
    let stem = p.file_stem().unwrap_or_default().to_string_lossy();
    let ext = p
        .extension()
        .map(|e| format!(".{}", e.to_str().unwrap_or("")))
        .unwrap_or_default();
    config.output = p
        .parent()
        .unwrap_or(std::path::Path::new("."))
        .join(format!("{}_sorted{}", stem, ext));

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--output" | "-o" if i + 1 < args.len() => {
                config.output = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--prefix" if i + 1 < args.len() => {
                config.prefix = args[i + 1].clone();
                i += 2;
            }
            "--min-length" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.min_length = n;
                }
                i += 2;
            }
            "--rename-table" if i + 1 < args.len() => {
                config.rename_table = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--keep-desc" => {
                config.strip_desc = false;
                i += 1;
            }
            "--sort-by-name" => {
                config.sort_by_name = true;
                i += 1;
            }
            _ => i += 1,
        }
    }

    println!("Sorting genome: {}", input);
    run_sort(&config).map_err(|e| anyhow::anyhow!("{}", e))?;
    println!("\n✓ Sorted FASTA: {}", config.output.display());
    println!(
        "  Next step: myconote-cli mask {} --engine repeatmodeler",
        config.output.display()
    );
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// train command: RNA-seq mediated training
// ─────────────────────────────────────────────────────────────────────────────

fn handle_train(masked_fasta: &str, args: &[String]) -> Result<()> {
    use train::{run_training, TrainConfig};

    let mut config = TrainConfig {
        masked_fasta: PathBuf::from(masked_fasta),
        ..TrainConfig::default()
    };

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--left" if i + 1 < args.len() => {
                config.left_reads = args[i + 1]
                    .split(',')
                    .map(|p| PathBuf::from(p.trim()))
                    .collect();
                i += 2;
            }
            "--right" if i + 1 < args.len() => {
                config.right_reads = args[i + 1]
                    .split(',')
                    .map(|p| PathBuf::from(p.trim()))
                    .collect();
                i += 2;
            }
            "--single" if i + 1 < args.len() => {
                config.single_reads = args[i + 1]
                    .split(',')
                    .map(|p| PathBuf::from(p.trim()))
                    .collect();
                i += 2;
            }
            "--trinity" if i + 1 < args.len() => {
                config.trinity_fasta = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--species" if i + 1 < args.len() => {
                config.species = args[i + 1].clone();
                i += 2;
            }
            "--output" | "-o" if i + 1 < args.len() => {
                config.out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--max-intron" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.max_intron = n;
                }
                i += 2;
            }
            "--min-models" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.min_models = n;
                }
                i += 2;
            }
            "--strand" if i + 1 < args.len() => {
                config.strand = args[i + 1].clone();
                i += 2;
            }
            "--memory" if i + 1 < args.len() => {
                config.trinity_memory = args[i + 1].clone();
                i += 2;
            }
            "--no-snap" => {
                config.train_snap = false;
                i += 1;
            }
            "--genemark" => {
                config.train_genemark = true;
                i += 1;
            }
            "--threads" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.threads = n;
                }
                i += 2;
            }
            _ => i += 1,
        }
    }

    if config.left_reads.is_empty()
        && config.single_reads.is_empty()
        && config.trinity_fasta.is_none()
    {
        eprintln!("Error: RNA-seq reads or a pre-assembled Trinity FASTA are required.");
        eprintln!(
            "Usage: myconote-cli train <genome.fa> --left R1.fq --right R2.fq --species myorg"
        );
        eprintln!("       myconote-cli train <genome.fa> --trinity trinity.fasta --species myorg");
        return Ok(());
    }

    run_training(&config).map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// update command: refine gene models with RNA-seq evidence
// ─────────────────────────────────────────────────────────────────────────────

fn handle_update(gff_path: &str, args: &[String]) -> Result<()> {
    use update::{run_update, UpdateConfig};

    let mut config = UpdateConfig {
        gff: PathBuf::from(gff_path),
        ..UpdateConfig::default()
    };

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--fasta" if i + 1 < args.len() => {
                config.fasta = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--output" | "-o" if i + 1 < args.len() => {
                config.out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--transcripts" if i + 1 < args.len() => {
                config.transcripts = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--rna-r1" if i + 1 < args.len() => {
                config.rna_r1 = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--rna-r2" if i + 1 < args.len() => {
                config.rna_r2 = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--rna-bam" if i + 1 < args.len() => {
                config.rna_bam = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--locus-prefix" if i + 1 < args.len() => {
                config.locus_prefix = args[i + 1].clone();
                i += 2;
            }
            "--organism" if i + 1 < args.len() => {
                config.organism = Some(args[i + 1].clone());
                i += 2;
            }
            "--max-utr-ext" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<u64>() {
                    config.max_utr_extension = n;
                }
                i += 2;
            }
            "--min-identity" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse::<f64>() {
                    config.min_identity = v;
                }
                i += 2;
            }
            "--threads" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.threads = n;
                }
                i += 2;
            }
            _ => i += 1,
        }
    }

    if config.fasta.as_os_str().is_empty() {
        eprintln!("Error: --fasta <genome.fa> is required.");
        return Ok(());
    }

    run_update(&config).map_err(|e| anyhow::anyhow!("{}", e))?;
    println!(
        "\n  Next step: myconote-cli annotate update_out/updated.gff3 --fasta {}",
        config.fasta.display()
    );
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// setup command: database downloader
// ─────────────────────────────────────────────────────────────────────────────

fn handle_setup(args: &[String]) -> Result<()> {
    use setup::{check_databases, download_databases, list_databases};
    use y1000plus::commands as y1000_cmds;

    // ── Y1000+ branch ─────────────────────────────────────────────────────
    // Any `--y1000plus` flag (alone, or with --list/--dry-run/--include/...)
    // routes to the Y1000+ bundle subsystem and skips the regular DB flow.
    if args.iter().any(|a| a == "--y1000plus") {
        let filtered: Vec<String> = args
            .iter()
            .filter(|a| a.as_str() != "--y1000plus")
            .cloned()
            .collect();
        let (y_args, _unused) = y1000_cmds::parse_args(&filtered)?;
        // Default to --list when nothing else is asked for.
        if !y_args.list
            && !y_args.dry_run
            && y_args.preset.is_none()
            && y_args.include.is_empty()
            && y_args.uninstall.is_empty()
        {
            y1000_cmds::run_list()?;
            return Ok(());
        }
        y1000_cmds::dispatch(&y_args)?;
        return Ok(());
    }

    let db_dir_default = {
        let home = std::env::var("HOME").unwrap_or_else(|_| "/tmp".to_string());
        PathBuf::from(home).join(".myconote").join("dbs")
    };

    let mut db_dir = db_dir_default;
    let mut keys: Vec<String> = Vec::new();
    let mut force = false;
    let mut do_list = false;
    let mut do_check = false;

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--db-dir" if i + 1 < args.len() => {
                db_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--dbs" | "--databases" => {
                i += 1;
                while i < args.len() && !args[i].starts_with('-') {
                    keys.push(args[i].clone());
                    i += 1;
                }
            }
            "--force" => {
                force = true;
                i += 1;
            }
            "--list" => {
                do_list = true;
                i += 1;
            }
            "--check" => {
                do_check = true;
                i += 1;
            }
            // Convenience: myconote setup --chat-corpus
            "--chat-corpus" => {
                keys.push("chat-corpus".to_string());
                i += 1;
            }
            // Convenience: myconote setup --ollama
            "--ollama" => {
                keys.push("ollama".to_string());
                i += 1;
            }
            // Legacy: myconote annotate --download-dbs
            "--download-dbs" | "--download" => {
                i += 1;
            }
            _ => i += 1,
        }
    }

    if do_list {
        list_databases();
        return Ok(());
    }
    if do_check {
        check_databases(&db_dir);
        return Ok(());
    }

    download_databases(&db_dir, &keys, force).map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// fix command: repair GenBank files
// ─────────────────────────────────────────────────────────────────────────────

fn handle_fix(gbk_path: &str, args: &[String]) -> Result<()> {
    use fix::{run_fix, FixConfig};

    let input = PathBuf::from(gbk_path);
    let default_output = {
        let p = PathBuf::from(gbk_path);
        let stem = p.file_stem().unwrap_or_default().to_string_lossy();
        let ext = p
            .extension()
            .map(|e| format!(".{}", e.to_str().unwrap_or("")))
            .unwrap_or_else(|| ".gbk".to_string());
        p.parent()
            .unwrap_or(std::path::Path::new("."))
            .join(format!("{}_fixed{}", stem, ext))
    };

    let mut config = FixConfig {
        input,
        output: default_output,
        ..FixConfig::default()
    };

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--output" | "-o" if i + 1 < args.len() => {
                config.output = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--report" if i + 1 < args.len() => {
                config.report = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--dry-run" => {
                config.dry_run = true;
                i += 1;
            }
            "--no-fix-tags" => {
                config.fix_dup_tags = false;
                i += 1;
            }
            "--no-fix-product" => {
                config.fix_product = false;
                i += 1;
            }
            "--no-fix-stops" => {
                config.fix_stops = false;
                i += 1;
            }
            _ => i += 1,
        }
    }

    run_fix(&config).map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// remote command: remote annotation services
// ─────────────────────────────────────────────────────────────────────────────

fn handle_remote(proteins_fa: &str, args: &[String]) -> Result<()> {
    use remote::{run_remote, RemoteConfig};

    let mut config = RemoteConfig {
        proteins_fa: PathBuf::from(proteins_fa),
        ..RemoteConfig::default()
    };

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--output" | "-o" if i + 1 < args.len() => {
                config.out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--phobius" => {
                config.run_phobius = true;
                i += 1;
            }
            "--interproscan" => {
                config.run_interpro = true;
                i += 1;
            }
            "--deeploc" => {
                config.run_deeploc = true;
                i += 1;
            }
            "--email" if i + 1 < args.len() => {
                config.email = args[i + 1].clone();
                i += 2;
            }
            "--batch-size" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.batch_size = n;
                }
                i += 2;
            }
            "--poll-interval" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<u64>() {
                    config.poll_interval = n;
                }
                i += 2;
            }
            "--max-retries" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.max_retries = n;
                }
                i += 2;
            }
            _ => i += 1,
        }
    }

    // If no service specified, default to Phobius (doesn't require sign-up)
    if !config.run_phobius && !config.run_interpro && !config.run_deeploc {
        config.run_phobius = true;
    }

    run_remote(&config).map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}

fn handle_submit(gff_path: &str, args: &[String]) -> Result<()> {
    let mut config = submit::SubmitConfig {
        gff: PathBuf::from(gff_path),
        ..submit::SubmitConfig::default()
    };

    let mut validate_only = false;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--fasta" if i + 1 < args.len() => {
                config.fasta = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--output" | "-o" if i + 1 < args.len() => {
                config.out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--organism" if i + 1 < args.len() => {
                config.organism = args[i + 1].clone();
                i += 2;
            }
            "--strain" if i + 1 < args.len() => {
                config.strain = Some(args[i + 1].clone());
                i += 2;
            }
            "--bioproject" if i + 1 < args.len() => {
                config.bioproject = Some(args[i + 1].clone());
                i += 2;
            }
            "--biosample" if i + 1 < args.len() => {
                config.biosample = Some(args[i + 1].clone());
                i += 2;
            }
            "--locus-prefix" if i + 1 < args.len() => {
                config.locus_tag_prefix = args[i + 1].clone();
                i += 2;
            }
            "--genetic-code" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<u8>() {
                    config.genetic_code = n;
                }
                i += 2;
            }
            "--email" if i + 1 < args.len() => {
                config.email = args[i + 1].clone();
                i += 2;
            }
            "--validate-only" => {
                validate_only = true;
                i += 1;
            }
            _ => i += 1,
        }
    }

    if config.fasta.as_os_str().is_empty() {
        return Err(anyhow::anyhow!(
            "--fasta is required. Run 'myconote-cli submit --help' for usage."
        ));
    }

    if validate_only {
        println!("── NCBI Submission Validation ───────────────────────────────");
        let result = submit::validate_for_ncbi(&config.gff, &config.fasta)
            .map_err(|e| anyhow::anyhow!("{}", e))?;
        result.print_summary();
    } else {
        println!("── NCBI Submission Preparation ──────────────────────────────");
        submit::run_table2asn(&config).map_err(|e| anyhow::anyhow!("{}", e))?;
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// batch command: multi-genome annotation
// ─────────────────────────────────────────────────────────────────────────────

fn handle_batch(args: &[String]) -> Result<()> {
    use batch::{run_batch, BatchConfig};

    let mut config = BatchConfig::default();

    // First positional arg is input (directory or sample sheet)
    // Handle --resume specially
    let mut i = 0;
    let mut input_set = false;

    while i < args.len() {
        match args[i].as_str() {
            "--resume" if i + 1 < args.len() => {
                config.resume = Some(PathBuf::from(&args[i + 1]));
                config.output_dir = PathBuf::from(&args[i + 1]);
                // Input is not required for resume
                input_set = true;
                i += 2;
            }
            "--output" | "-o" if i + 1 < args.len() => {
                config.output_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--stages" if i + 1 < args.len() => {
                config.stages = args[i + 1]
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .collect();
                i += 2;
            }
            "--kingdom" if i + 1 < args.len() => {
                config.kingdom = args[i + 1].clone();
                i += 2;
            }
            "--threads" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.threads = n;
                }
                i += 2;
            }
            "--parallel" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.max_parallel = n;
                }
                i += 2;
            }
            "--min-length" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.min_length = n;
                }
                i += 2;
            }
            "--mask-engine" if i + 1 < args.len() => {
                config.mask_engine = args[i + 1].clone();
                i += 2;
            }
            "--genetic-code" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<u8>() {
                    config.genetic_code = n;
                }
                i += 2;
            }
            "--locus-prefix" if i + 1 < args.len() => {
                config.locus_prefix = args[i + 1].clone();
                i += 2;
            }
            "--condor" => {
                config.condor = true;
                i += 1;
            }
            "--condor-cpus" if i + 1 < args.len() => {
                if let Ok(n) = args[i + 1].parse::<usize>() {
                    config.condor_cpus = n;
                }
                i += 2;
            }
            "--condor-mem" if i + 1 < args.len() => {
                config.condor_mem = args[i + 1].clone();
                i += 2;
            }
            "--condor-disk" if i + 1 < args.len() => {
                config.condor_disk = args[i + 1].clone();
                i += 2;
            }
            "--condor-queue" if i + 1 < args.len() => {
                config.condor_queue = Some(args[i + 1].clone());
                i += 2;
            }
            "--condor-extra" if i + 1 < args.len() => {
                config.condor_extra = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            other if !other.starts_with('-') && !input_set => {
                config.input = PathBuf::from(other);
                input_set = true;
                i += 1;
            }
            _ => i += 1,
        }
    }

    if !input_set {
        return Err(anyhow::anyhow!(
            "No input provided. Pass a directory of FASTAs or a sample sheet TSV.\n\
             Run 'myconote-cli batch --help' for usage."
        ));
    }

    println!("── Batch Annotation ─────────────────────────────────────────");
    run_batch(&config).map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}
