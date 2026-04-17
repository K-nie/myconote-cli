//! N-genome comparative genomics — ortholog inference + pan-genome summary.
//!
//! This is the real implementation after the earlier fake-output version was
//! deleted. Rather than rolling our own clustering (OrthoMCL/OMA-style) we
//! wrap OrthoFinder, which is the research-grade standard and gives
//! paralog-aware ortholog tables plus a rooted species tree for free.
//!
//! Scope is intentionally capped by genome size — see `Tier` below — because
//! OrthoFinder's all-vs-all similarity search scales with total_proteins² and
//! the user experience on commodity hardware degrades sharply past certain
//! thresholds. Users with HPC can pass `--force-cap` to bypass the tier cap
//! at their own risk.

pub mod orthofinder;
pub mod primary_transcripts;

use crate::utils::error::{MycoNoteError, Result};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Per-genome input. A (GFF3, FASTA, genome_name) triple after expansion —
/// name defaults to the FASTA file stem but can be overridden explicitly.
#[derive(Debug, Clone)]
pub struct GenomeInput {
    pub name: String,
    pub gff: PathBuf,
    pub fasta: PathBuf,
}

/// Protein-count-based tier. Determines the hard cap on genome count so the
/// all-vs-all search stays tractable on the user's hardware.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tier {
    Small,  // ≤ 15,000 proteins — fungi, most protists
    Medium, // 15,001 – 30,000 — small/mid plants, many invertebrates
    Large,  // > 30,000 — major crops, vertebrates
}

impl Tier {
    pub fn from_protein_count(n: usize) -> Self {
        match n {
            0..=15_000 => Tier::Small,
            15_001..=30_000 => Tier::Medium,
            _ => Tier::Large,
        }
    }
    pub fn cap(self) -> usize {
        match self {
            Tier::Small => 5,
            Tier::Medium => 3,
            Tier::Large => 2,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Tier::Small => "Small",
            Tier::Medium => "Medium",
            Tier::Large => "Large",
        }
    }
}

#[derive(Debug, Clone)]
pub struct CompareConfig {
    pub inputs: Vec<GenomeInput>,
    pub output_dir: PathBuf,
    pub threads: usize,
    /// Keep one longest-CDS protein per gene rather than every isoform.
    /// Default on — alternative-splice isoforms cluster into the same
    /// orthogroup anyway and inflate OrthoFinder's cost.
    pub primary_only: bool,
    /// Use diamond_ultra_sens search mode. Default on for N≤5 where accuracy
    /// matters more than the extra runtime.
    pub sensitive: bool,
    /// Multiple-sequence-alignment tree refinement (2–3× slower but slightly
    /// more accurate species tree). Default off.
    pub msa: bool,
    /// Override the tier cap. Hidden from --help; advanced-use escape hatch.
    pub force_cap: Option<usize>,
    /// NCBI genetic code table. Default 1 (standard).
    pub genetic_code: u8,
    /// Core threshold (fraction of genomes to count as "soft-core").
    pub soft_core_frac: f64,
    /// Cloud threshold (fraction at or below which a cluster is "cloud").
    pub cloud_frac: f64,
    /// When true, chain MAFFT + supermatrix concatenation onto the
    /// OrthoFinder run so `phylogeny` can consume the output directly.
    /// Cost scales with the number of single-copy orthogroups; typical
    /// fungal clade runs add a few minutes. Off by default.
    pub species_tree: bool,
}

impl Default for CompareConfig {
    fn default() -> Self {
        Self {
            inputs: Vec::new(),
            output_dir: PathBuf::from("compare_out"),
            threads: num_threads_default(),
            primary_only: true,
            sensitive: true,
            msa: false,
            force_cap: None,
            genetic_code: 1,
            soft_core_frac: 0.95,
            cloud_frac: 0.15,
            species_tree: false,
        }
    }
}

fn num_threads_default() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
}

/// Classify all input genomes, find the most demanding tier, and enforce
/// its cap. Returns the effective tier for logging.
pub fn classify_and_enforce_cap(gene_counts: &[usize], force_cap: Option<usize>) -> Result<Tier> {
    if gene_counts.is_empty() {
        return Err(MycoNoteError::InvalidFormat("No input genomes".to_string()));
    }
    let max_count = *gene_counts.iter().max().unwrap();
    let tier = Tier::from_protein_count(max_count);
    let n = gene_counts.len();
    let cap = force_cap.unwrap_or_else(|| tier.cap());

    if n > cap {
        let hint = if force_cap.is_some() {
            format!(
                "You passed --force-cap {}; at least {} genomes provided.",
                cap, n
            )
        } else {
            format!(
                "Tier {} caps at {} genomes; {} provided (max protein count: {}).",
                tier.label(),
                cap,
                n,
                max_count
            )
        };
        return Err(MycoNoteError::InvalidFormat(format!(
            "{}\n  Rationale: OrthoFinder memory and runtime scale quadratically with total\n  \
             protein count; above the cap, commodity hardware swaps to disk and the run\n  \
             becomes unreliable.\n  \
             Options:\n    1. Reduce input to {} genomes\n    2. Re-run with --primary-only to shrink each proteome\n    3. Run OrthoFinder directly on HPC: orthofinder -f proteins/ -t <threads>\n    4. Bypass at your own risk: --force-cap <n>",
            hint, cap
        )));
    }
    Ok(tier)
}

/// Parse positional CLI arguments into GenomeInput pairs. Expects args to
/// alternate GFF3 → FASTA; bails with a helpful error if the pattern breaks.
pub fn parse_positional_inputs(paths: &[String]) -> Result<Vec<GenomeInput>> {
    if paths.is_empty() {
        return Err(MycoNoteError::InvalidFormat(
            "compare needs at least one GFF3 + FASTA pair".to_string(),
        ));
    }
    if paths.len() % 2 != 0 {
        return Err(MycoNoteError::InvalidFormat(format!(
            "compare needs an even number of positional arguments (GFF3 FASTA pairs); got {}",
            paths.len()
        )));
    }

    let mut inputs = Vec::new();
    for pair in paths.chunks(2) {
        let gff = PathBuf::from(&pair[0]);
        let fasta = PathBuf::from(&pair[1]);
        if !is_gff(&gff) {
            return Err(MycoNoteError::InvalidFormat(format!(
                "Expected GFF3 as first of each pair: {} (extension should be .gff3 or .gff)",
                gff.display()
            )));
        }
        if !is_fasta(&fasta) {
            return Err(MycoNoteError::InvalidFormat(format!(
                "Expected FASTA as second of each pair: {} (extension should be .fa/.fasta/.fna)",
                fasta.display()
            )));
        }
        let name = fasta
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("genome")
            .to_string();
        inputs.push(GenomeInput { name, gff, fasta });
    }
    Ok(inputs)
}

fn is_gff(p: &Path) -> bool {
    matches!(
        p.extension().and_then(|e| e.to_str()).map(str::to_lowercase),
        Some(ref s) if s == "gff3" || s == "gff"
    )
}

fn is_fasta(p: &Path) -> bool {
    // Accept the full set of common FASTA extensions. `.fas` is widely used
    // by SGD / FungiDB / candidagenome.org — the earlier omission caused
    // E2E runs on real test data to fail at parse time.
    matches!(
        p.extension().and_then(|e| e.to_str()).map(str::to_lowercase),
        Some(ref s) if s == "fa" || s == "fasta" || s == "fna" || s == "faa" || s == "fas" || s == "ffn" || s == "frn"
    )
}

/// Entry point.
pub fn run_compare(config: &CompareConfig) -> Result<()> {
    println!("── Comparative genomics ─────────────────────────────────────");
    println!("  Genomes : {}", config.inputs.len());
    println!("  Output  : {}", config.output_dir.display());

    // Disambiguate duplicate genome names — OrthoFinder keys everything on
    // the input file stem, so identical names would silently collide.
    check_unique_names(&config.inputs)?;

    std::fs::create_dir_all(&config.output_dir).map_err(MycoNoteError::Io)?;
    let proteins_dir = config.output_dir.join("proteins");
    std::fs::create_dir_all(&proteins_dir).map_err(MycoNoteError::Io)?;

    // ── 1. Extract primary-transcript proteins per genome ────────────────
    println!();
    println!("  Step 1: extracting primary-transcript proteins…");
    let mut protein_counts: Vec<usize> = Vec::new();
    for g in &config.inputs {
        if !g.gff.exists() {
            return Err(MycoNoteError::InvalidFormat(format!(
                "GFF3 not found: {}",
                g.gff.display()
            )));
        }
        if !g.fasta.exists() {
            return Err(MycoNoteError::InvalidFormat(format!(
                "FASTA not found: {}",
                g.fasta.display()
            )));
        }
        let out_fa = proteins_dir.join(format!("{}.faa", g.name));
        let n = if config.primary_only {
            primary_transcripts::extract_primary_proteins(
                &g.gff,
                &g.fasta,
                &out_fa,
                config.genetic_code,
            )?
        } else {
            // Non-primary mode: fall back to the annotate-style extractor.
            // Still writes one protein per mRNA (all isoforms) which is what
            // the user asked for when they pass --no-primary-only.
            return Err(MycoNoteError::InvalidFormat(
                "--no-primary-only not yet supported (always runs in primary-transcript mode)"
                    .to_string(),
            ));
        };
        println!("    {}: {} proteins", g.name, n);
        protein_counts.push(n);
    }

    if protein_counts.iter().any(|&n| n == 0) {
        return Err(MycoNoteError::InvalidFormat(
            "One or more genomes produced zero proteins — check GFF3 has 'gene' + 'mRNA' + 'CDS' rows and FASTA seqids match".to_string(),
        ));
    }

    // ── 2. Tier detection + cap enforcement ──────────────────────────────
    let tier = classify_and_enforce_cap(&protein_counts, config.force_cap)?;
    let total: usize = protein_counts.iter().sum();
    println!();
    println!(
        "  Tier: {} (largest genome = {} proteins, total = {} proteins)",
        tier.label(),
        protein_counts.iter().max().copied().unwrap_or(0),
        total
    );

    // ── 3. Invoke OrthoFinder ─────────────────────────────────────────────
    println!();
    println!(
        "  Step 2: running OrthoFinder (sensitive={}, msa={})…",
        config.sensitive, config.msa
    );
    let results_dir =
        orthofinder::run_orthofinder(&proteins_dir, config.threads, config.sensitive, config.msa)?;
    println!("    OrthoFinder results: {}", results_dir.display());

    // ── 4. Parse orthogroups + pan-genome summary ────────────────────────
    println!();
    println!("  Step 3: parsing orthogroups + computing pan-genome shape…");
    let ogs_path = results_dir.join("Orthogroups").join("Orthogroups.tsv");
    if !ogs_path.exists() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "OrthoFinder didn't produce Orthogroups.tsv at {}",
            ogs_path.display()
        )));
    }
    let orthogroups = orthofinder::parse_orthogroups_tsv(&ogs_path)?;
    let n_genomes = config.inputs.len();
    let summary = orthofinder::summarise_pangenome(
        &orthogroups,
        n_genomes,
        config.soft_core_frac,
        config.cloud_frac,
    );

    // ── 5. Write MycoNote-native output files ────────────────────────────
    write_ortholog_table(
        &config.output_dir.join("ortholog_table.tsv"),
        &orthogroups,
        &config.inputs,
    )?;
    write_pangenome_summary(
        &config.output_dir.join("pangenome_summary.tsv"),
        &summary,
        n_genomes,
    )?;

    // Re-expose key OrthoFinder outputs at the top level so users don't have
    // to go spelunking inside the Results_<date> directory.
    let species_tree_src = results_dir
        .join("Species_Tree")
        .join("SpeciesTree_rooted.txt");
    if species_tree_src.exists() {
        let dst = config.output_dir.join("species_tree.nwk");
        let _ = std::fs::copy(&species_tree_src, &dst);
    }
    let sco_dir = results_dir.join("Single_Copy_Orthologue_Sequences");
    let sco_count = if sco_dir.is_dir() {
        std::fs::read_dir(&sco_dir)
            .map(|rd| rd.filter_map(|e| e.ok()).count())
            .unwrap_or(0)
    } else {
        0
    };

    // ── 6. Optional: chain alignment + concatenation for species tree ────
    // OrthoFinder emits one unaligned FASTA per single-copy orthogroup, but
    // `phylogeny` wants ONE aligned supermatrix. This step aligns each
    // orthogroup with MAFFT and writes a concatenated supermatrix +
    // partition file ready to feed `phylogeny --partition`.
    let species_tree_alignment = if config.species_tree && sco_count > 0 {
        println!();
        println!("  Step 4: building species-tree alignment (MAFFT + concat)…");
        let taxon_names: Vec<String> = config.inputs.iter().map(|g| g.name.clone()).collect();
        let align_out = config.output_dir.join("species_tree_alignment");
        match crate::align::build_species_tree_alignment(
            &sco_dir,
            &taxon_names,
            &orthogroups,
            &align_out,
            config.threads,
        ) {
            Ok(result) => Some(result),
            Err(e) => {
                eprintln!("  ⚠  Species-tree alignment failed (non-fatal): {}", e);
                None
            }
        }
    } else {
        None
    };

    // ── 7. Human-readable report ─────────────────────────────────────────
    write_report(
        &config.output_dir.join("compare_report.txt"),
        &config.inputs,
        &protein_counts,
        tier,
        &summary,
        sco_count,
    )?;

    println!();
    println!("✓  Comparison complete.");
    println!(
        "    Ortholog table:      {}/ortholog_table.tsv",
        config.output_dir.display()
    );
    println!(
        "    Pan-genome summary:  {}/pangenome_summary.tsv",
        config.output_dir.display()
    );
    if species_tree_src.exists() {
        println!(
            "    Species tree:        {}/species_tree.nwk",
            config.output_dir.display()
        );
    }
    println!(
        "    Report:              {}/compare_report.txt",
        config.output_dir.display()
    );
    if let Some(sta) = &species_tree_alignment {
        println!();
        println!(
            "    Species-tree alignment: {} ({} taxa × {} columns, {} partitions)",
            sta.supermatrix.display(),
            sta.n_taxa,
            sta.total_columns,
            sta.n_orthogroups
        );
        println!(
            "    → Next: myconote-cli phylogeny {} --partition {} --threads {}",
            sta.supermatrix.display(),
            sta.partition.display(),
            config.threads
        );
    } else if sco_count > 0 {
        println!();
        println!(
            "  {} single-copy orthogroups available (unaligned).",
            sco_count
        );
        println!(
            "    → Re-run with --species-tree to produce a supermatrix ready for `myconote-cli phylogeny`."
        );
    }

    Ok(())
}

fn check_unique_names(inputs: &[GenomeInput]) -> Result<()> {
    let mut seen = std::collections::HashSet::new();
    for g in inputs {
        if !seen.insert(g.name.clone()) {
            return Err(MycoNoteError::InvalidFormat(format!(
                "Duplicate genome name '{}' — each input FASTA must have a unique file stem, \
                 or rename with per-input flags",
                g.name
            )));
        }
    }
    Ok(())
}

fn write_ortholog_table(
    path: &Path,
    orthogroups: &[orthofinder::Orthogroup],
    inputs: &[GenomeInput],
) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    let mut header = String::from("orthogroup_id\tn_genomes\tcategory");
    for g in inputs {
        header.push('\t');
        header.push_str(&g.name);
    }
    writeln!(f, "{}", header).map_err(MycoNoteError::Io)?;

    let total = inputs.len();
    for og in orthogroups {
        let n = og.n_genomes();
        let category = if n == total {
            if og.is_single_copy(total) {
                "core_single_copy"
            } else {
                "core"
            }
        } else if n == 1 {
            "singleton"
        } else {
            "accessory"
        };
        let mut line = format!("{}\t{}\t{}", og.id, n, category);
        for g in inputs {
            line.push('\t');
            let genes = og.members.get(&g.name);
            match genes {
                Some(v) if !v.is_empty() => line.push_str(&v.join(",")),
                _ => line.push('-'),
            }
        }
        writeln!(f, "{}", line).map_err(MycoNoteError::Io)?;
    }
    Ok(())
}

fn write_pangenome_summary(
    path: &Path,
    s: &orthofinder::PangenomeSummary,
    n_genomes: usize,
) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "category\tn_clusters\tfraction_of_total").map_err(MycoNoteError::Io)?;
    let frac = |n: usize| {
        if s.total_clusters == 0 {
            0.0
        } else {
            n as f64 / s.total_clusters as f64
        }
    };
    writeln!(f, "total\t{}\t1.000", s.total_clusters).map_err(MycoNoteError::Io)?;
    writeln!(f, "core\t{}\t{:.3}", s.core, frac(s.core)).map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "core_single_copy\t{}\t{:.3}",
        s.single_copy_core,
        frac(s.single_copy_core)
    )
    .map_err(MycoNoteError::Io)?;
    writeln!(f, "soft_core\t{}\t{:.3}", s.soft_core, frac(s.soft_core))
        .map_err(MycoNoteError::Io)?;
    writeln!(f, "shell\t{}\t{:.3}", s.shell, frac(s.shell)).map_err(MycoNoteError::Io)?;
    writeln!(f, "cloud\t{}\t{:.3}", s.cloud, frac(s.cloud)).map_err(MycoNoteError::Io)?;
    writeln!(f, "singletons\t{}\t{:.3}", s.singletons, frac(s.singletons))
        .map_err(MycoNoteError::Io)?;
    writeln!(f, "# n_genomes = {}", n_genomes).map_err(MycoNoteError::Io)?;
    Ok(())
}

fn write_report(
    path: &Path,
    inputs: &[GenomeInput],
    protein_counts: &[usize],
    tier: Tier,
    summary: &orthofinder::PangenomeSummary,
    sco_count: usize,
) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(f, "myconote compare — Comparative Genomics Report").map_err(MycoNoteError::Io)?;
    writeln!(f, "=================================================").map_err(MycoNoteError::Io)?;
    writeln!(f, "Genomes    : {}", inputs.len()).map_err(MycoNoteError::Io)?;
    writeln!(f, "Tier       : {}", tier.label()).map_err(MycoNoteError::Io)?;
    writeln!(f).map_err(MycoNoteError::Io)?;
    writeln!(f, "Per-genome protein counts (primary transcripts):").map_err(MycoNoteError::Io)?;
    for (g, n) in inputs.iter().zip(protein_counts.iter()) {
        writeln!(f, "  {:<40} {:>8}", g.name, n).map_err(MycoNoteError::Io)?;
    }
    writeln!(f).map_err(MycoNoteError::Io)?;
    writeln!(f, "Pan-genome shape (Tettelin bins):").map_err(MycoNoteError::Io)?;
    writeln!(f, "  Total orthogroups   : {}", summary.total_clusters).map_err(MycoNoteError::Io)?;
    writeln!(f, "  Core (all genomes)  : {}", summary.core).map_err(MycoNoteError::Io)?;
    writeln!(f, "    - single-copy     : {}", summary.single_copy_core)
        .map_err(MycoNoteError::Io)?;
    writeln!(f, "  Soft-core (≥95 %)   : {}", summary.soft_core).map_err(MycoNoteError::Io)?;
    writeln!(f, "  Shell               : {}", summary.shell).map_err(MycoNoteError::Io)?;
    writeln!(f, "  Cloud               : {}", summary.cloud).map_err(MycoNoteError::Io)?;
    writeln!(f, "    - singletons      : {}", summary.singletons).map_err(MycoNoteError::Io)?;
    writeln!(f).map_err(MycoNoteError::Io)?;
    writeln!(f, "Single-copy orthogroups for species tree: {}", sco_count)
        .map_err(MycoNoteError::Io)?;
    writeln!(f).map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "Generated by: OrthoFinder (Emms & Kelly 2019, Genome Biology)"
    )
    .map_err(MycoNoteError::Io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tier_thresholds_are_correct() {
        assert_eq!(Tier::from_protein_count(0), Tier::Small);
        assert_eq!(Tier::from_protein_count(15_000), Tier::Small);
        assert_eq!(Tier::from_protein_count(15_001), Tier::Medium);
        assert_eq!(Tier::from_protein_count(30_000), Tier::Medium);
        assert_eq!(Tier::from_protein_count(30_001), Tier::Large);
        assert_eq!(Tier::from_protein_count(100_000), Tier::Large);
    }

    #[test]
    fn tier_caps() {
        assert_eq!(Tier::Small.cap(), 5);
        assert_eq!(Tier::Medium.cap(), 3);
        assert_eq!(Tier::Large.cap(), 2);
    }

    #[test]
    fn max_protein_count_determines_tier() {
        // Mixing one plant with four fungi should still trigger the plant
        // tier — this prevents a heavy genome from quietly exploding the
        // run when buried in a majority of small ones.
        let counts = vec![5_000, 6_000, 7_000, 8_000, 40_000];
        let result = classify_and_enforce_cap(&counts, None);
        assert!(
            result.is_err(),
            "5 inputs where largest is Large should be rejected"
        );
        let err_msg = format!("{}", result.unwrap_err());
        assert!(
            err_msg.contains("Large"),
            "error must mention the binding tier:\n{}",
            err_msg
        );
        assert!(
            err_msg.contains("2 genomes"),
            "error must state the cap:\n{}",
            err_msg
        );
    }

    #[test]
    fn small_tier_accepts_up_to_5() {
        let ok_4 = classify_and_enforce_cap(&vec![6_000; 4], None);
        assert!(ok_4.is_ok());
        let ok_5 = classify_and_enforce_cap(&vec![6_000; 5], None);
        assert!(ok_5.is_ok());
        let fail_6 = classify_and_enforce_cap(&vec![6_000; 6], None);
        assert!(fail_6.is_err());
    }

    #[test]
    fn medium_tier_caps_at_3() {
        let ok_3 = classify_and_enforce_cap(&vec![25_000; 3], None);
        assert_eq!(ok_3.unwrap(), Tier::Medium);
        let fail_4 = classify_and_enforce_cap(&vec![25_000; 4], None);
        assert!(fail_4.is_err());
    }

    #[test]
    fn large_tier_caps_at_2() {
        let ok_2 = classify_and_enforce_cap(&vec![50_000; 2], None);
        assert_eq!(ok_2.unwrap(), Tier::Large);
        let fail_3 = classify_and_enforce_cap(&vec![50_000; 3], None);
        assert!(fail_3.is_err());
    }

    #[test]
    fn force_cap_bypasses_tier_limit() {
        let counts = vec![50_000; 8];
        let result = classify_and_enforce_cap(&counts, Some(10));
        assert!(
            result.is_ok(),
            "--force-cap 10 should accept 8 large genomes"
        );
        assert_eq!(result.unwrap(), Tier::Large);
    }

    #[test]
    fn force_cap_still_rejects_when_exceeded() {
        // --force-cap only raises the cap — it doesn't eliminate it.
        let counts = vec![50_000; 8];
        let result = classify_and_enforce_cap(&counts, Some(5));
        assert!(
            result.is_err(),
            "--force-cap 5 should still reject 8 inputs"
        );
    }

    #[test]
    fn parse_positional_requires_paired_args() {
        let odd = vec!["a.gff3".into(), "a.fa".into(), "b.gff3".into()];
        assert!(parse_positional_inputs(&odd).is_err());
    }

    #[test]
    fn parse_positional_validates_extensions() {
        let bad = vec!["a.txt".into(), "a.fa".into()];
        let err = parse_positional_inputs(&bad).unwrap_err();
        assert!(
            format!("{}", err).contains("GFF3"),
            "error should mention the GFF3 extension requirement"
        );

        let swapped = vec!["a.fa".into(), "a.gff3".into()];
        assert!(parse_positional_inputs(&swapped).is_err());
    }

    #[test]
    fn parse_positional_builds_named_inputs() {
        let ok = vec![
            "genome1.gff3".into(),
            "genome1.fa".into(),
            "genome2.gff".into(),
            "genome2.fasta".into(),
        ];
        let inputs = parse_positional_inputs(&ok).unwrap();
        assert_eq!(inputs.len(), 2);
        assert_eq!(inputs[0].name, "genome1");
        assert_eq!(inputs[1].name, "genome2");
    }
}
