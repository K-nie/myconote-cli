//! Allele-specific expression (ASE) for heterozygous / hybrid /
//! polyploid fungal genomes.
//!
//! Given a phased VCF + a reference CDS FASTA + the annotated GFF3
//! + a sample sheet, build personalized transcriptomes per haplotype
//! and quantify expression separately against each. Produces
//! allele-count matrices suitable for downstream cis/trans
//! regression or binomial ASE tests in R.
//!
//! See `scratch/ase_spec.md` for the full design; all 10 open
//! decisions were locked on 2026-04-24 before code started. The GFF3
//! flag (not in the initial spec) was added during implementation
//! because personalize.rs needs it to map genome → CDS coordinates.

pub mod bundle;
pub mod merge;
pub mod personalize;
pub mod quant;
pub mod vcf;

use crate::ase::merge::HapQuantInput;
use crate::ase::personalize::{
    personalize, read_cds_fasta, write_cds_fasta, TranscriptIndex, VariantApplication,
};
use crate::ase::quant::{build_haplotype_indices, quantify_sample, AseQuantSpec};
use crate::ase::vcf::{
    parse_phased_vcf, summarize_variants, PhasedVariant, SkipReason, VcfOptions,
};
use crate::quant::bundle::{extract_fastp_summary, sha256_file, ExternalTools, HashedPath};
use crate::quant::index::{resolve_cache_root, CacheEnv};
use crate::quant::sample_sheet::{self, Sample};
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct AseConfig {
    pub cds_fa: PathBuf,
    pub vcf: PathBuf,
    pub gff3: PathBuf,
    pub samples_sheet: PathBuf,
    pub genome_fa: PathBuf,
    pub output_dir: PathBuf,
    pub haplotype_names: [String; 2],
    pub k: usize,
    pub threads: usize,
    pub tmpdir: Option<PathBuf>,
    pub index_cache: Option<PathBuf>,
    pub keep_trimmed: Option<PathBuf>,
    pub max_indel_size: usize,
    pub asymmetry_threshold: f64,
    pub seed: u64,
    pub fastp_bin: Option<String>,
    pub salmon_bin: Option<String>,
    pub raw_cmdline: String,
}

impl AseConfig {
    fn parse(args: &[String]) -> Result<Self> {
        let cds_fa: Option<PathBuf>;
        let mut vcf: Option<PathBuf> = None;
        let mut gff3: Option<PathBuf> = None;
        let mut samples_sheet: Option<PathBuf> = None;
        let mut genome_fa: Option<PathBuf> = None;
        let mut output_dir = PathBuf::from("ase_out");
        let mut haplotype_names = ["hap0".to_string(), "hap1".to_string()];
        let mut k: usize = 31;
        let mut threads = threads_default();
        let mut tmpdir: Option<PathBuf> = None;
        let mut index_cache: Option<PathBuf> = None;
        let mut keep_trimmed: Option<PathBuf> = None;
        let mut max_indel_size: usize = 50;
        let mut asymmetry_threshold: f64 = 0.05;
        let mut seed: u64 = 42;
        let mut fastp_bin: Option<String> = None;
        let mut salmon_bin: Option<String> = None;
        let mut positional: Vec<String> = Vec::new();

        let mut i = 0usize;
        while i < args.len() {
            match args[i].as_str() {
                "--vcf" if i + 1 < args.len() => {
                    vcf = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--gff3" if i + 1 < args.len() => {
                    gff3 = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--samples" if i + 1 < args.len() => {
                    samples_sheet = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--genome" if i + 1 < args.len() => {
                    genome_fa = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--output" | "-o" if i + 1 < args.len() => {
                    output_dir = PathBuf::from(&args[i + 1]);
                    i += 2;
                }
                "--haplotype-names" if i + 1 < args.len() => {
                    let parts: Vec<&str> = args[i + 1].split(',').collect();
                    if parts.len() != 2 {
                        return Err(MycoNoteError::QuantSheet(format!(
                            "ase: --haplotype-names expects 'NAME1,NAME2', got '{}'",
                            args[i + 1]
                        )));
                    }
                    haplotype_names = [parts[0].trim().to_string(), parts[1].trim().to_string()];
                    i += 2;
                }
                "-k" if i + 1 < args.len() => {
                    k = args[i + 1].parse().unwrap_or(k);
                    i += 2;
                }
                "--threads" | "-t" if i + 1 < args.len() => {
                    threads = args[i + 1].parse().unwrap_or(threads).max(1);
                    i += 2;
                }
                "--tmpdir" if i + 1 < args.len() => {
                    tmpdir = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--index-cache" if i + 1 < args.len() => {
                    index_cache = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--keep-trimmed" if i + 1 < args.len() => {
                    keep_trimmed = Some(PathBuf::from(&args[i + 1]));
                    i += 2;
                }
                "--max-indel-size" if i + 1 < args.len() => {
                    max_indel_size = args[i + 1].parse().unwrap_or(max_indel_size);
                    i += 2;
                }
                "--asymmetry-threshold" if i + 1 < args.len() => {
                    asymmetry_threshold = args[i + 1].parse().unwrap_or(asymmetry_threshold);
                    i += 2;
                }
                "--seed" if i + 1 < args.len() => {
                    seed = args[i + 1].parse().unwrap_or(seed);
                    i += 2;
                }
                "--fastp" if i + 1 < args.len() => {
                    fastp_bin = Some(args[i + 1].clone());
                    i += 2;
                }
                "--salmon" if i + 1 < args.len() => {
                    salmon_bin = Some(args[i + 1].clone());
                    i += 2;
                }
                other if other.starts_with("--") => {
                    return Err(MycoNoteError::QuantSheet(format!(
                        "ase: unknown flag '{other}'"
                    )));
                }
                _ => {
                    positional.push(args[i].clone());
                    i += 1;
                }
            }
        }

        cds_fa = positional.into_iter().next().map(PathBuf::from);
        let cds_fa = cds_fa.ok_or_else(|| {
            MycoNoteError::QuantSheet(
                "ase: missing CDS FASTA. Usage:\n  \
                 myconote-cli ase <cds.fa> --vcf <phased.vcf.gz> --gff3 <annotated.gff3> \
                 --samples <sheet.tsv> --genome <genome.fa>"
                    .to_string(),
            )
        })?;
        let vcf = vcf.ok_or_else(|| {
            MycoNoteError::QuantSheet("ase: --vcf <phased.vcf.gz> is required".to_string())
        })?;
        let gff3 = gff3.ok_or_else(|| {
            MycoNoteError::QuantSheet(
                "ase: --gff3 <annotated.gff3> is required. Use the same GFF3 that \
                 `convert --to cds` consumed to produce <cds.fa>."
                    .to_string(),
            )
        })?;
        let samples_sheet = samples_sheet.ok_or_else(|| {
            MycoNoteError::QuantSheet("ase: --samples <sheet.tsv> is required".to_string())
        })?;
        let genome_fa = genome_fa.ok_or_else(|| {
            MycoNoteError::QuantSheet(
                "ase: --genome <genome.fa> is required (used as the salmon decoy)".to_string(),
            )
        })?;

        if haplotype_names[0] == haplotype_names[1] {
            return Err(MycoNoteError::QuantSheet(format!(
                "ase: haplotype names must differ; got '{}' twice",
                haplotype_names[0]
            )));
        }

        let raw_cmdline = std::iter::once("myconote-cli ase".to_string())
            .chain(args.iter().cloned())
            .collect::<Vec<_>>()
            .join(" ");

        Ok(Self {
            cds_fa,
            vcf,
            gff3,
            samples_sheet,
            genome_fa,
            output_dir,
            haplotype_names,
            k,
            threads,
            tmpdir,
            index_cache,
            keep_trimmed,
            max_indel_size,
            asymmetry_threshold,
            seed,
            fastp_bin,
            salmon_bin,
            raw_cmdline,
        })
    }
}

fn threads_default() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
}

// ─────────────────────────────────────────────────────────────────────────────
// Dispatcher
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_ase(args: &[String]) -> Result<()> {
    let cfg = AseConfig::parse(args)?;
    dispatch(&cfg)
}

fn dispatch(cfg: &AseConfig) -> Result<()> {
    // 1. Validate inputs.
    require_file(&cfg.cds_fa, "<cds.fa>")?;
    require_file(&cfg.vcf, "--vcf")?;
    require_file(&cfg.gff3, "--gff3")?;
    require_file(&cfg.samples_sheet, "--samples")?;
    require_file(&cfg.genome_fa, "--genome")?;
    std::fs::create_dir_all(&cfg.output_dir)?;

    // 2. Parse the sample sheet (reuses src/quant/sample_sheet.rs).
    let sheet = sample_sheet::parse_sheet(&cfg.samples_sheet)?;
    eprintln!(
        "  ✓ {} samples from {}",
        sheet.samples.len(),
        cfg.samples_sheet.display()
    );
    require_all_fastqs(&sheet.samples)?;

    // 3. Parse the phased VCF. Unphased heterozygous sites fail the
    //    parse with a concrete line number; here we just surface any
    //    error as-is.
    let vopts = VcfOptions {
        max_indel_size: cfg.max_indel_size,
    };
    let vcf_parsed = parse_phased_vcf(&cfg.vcf, &vopts)?;
    eprintln!(
        "  ✓ {} phased variants from {} ({} skipped during parse)",
        vcf_parsed.variants.len(),
        cfg.vcf.display(),
        vcf_parsed.skipped.len()
    );

    // 4. Build the transcript index from the GFF3. This walks the
    //    file once and keeps per-transcript CDS exons.
    let tx_index = TranscriptIndex::from_gff3(&cfg.gff3)?;
    let n_transcripts: usize = tx_index.transcripts().count();
    eprintln!(
        "  ✓ {} transcripts indexed from {}",
        n_transcripts,
        cfg.gff3.display()
    );

    // 5. Read reference CDS, then build one personalized CDS FASTA
    //    per haplotype by applying phased variants.
    let reference_cds = read_cds_fasta(&cfg.cds_fa)?;
    let personalized = personalize(
        &reference_cds,
        &tx_index,
        &vcf_parsed.variants,
        &cfg.haplotype_names,
    )?;
    let mut hap_fastas: Vec<(String, PathBuf)> = Vec::with_capacity(2);
    let mut hap_sha: Vec<String> = Vec::with_capacity(2);
    let mut hap_n_modified: Vec<usize> = Vec::with_capacity(2);
    let mut hap_n_applied: Vec<usize> = Vec::with_capacity(2);
    for (hap_idx, h) in personalized.haplotypes.iter().enumerate() {
        let out_path = cfg.output_dir.join(format!("cds_{}.fa", h.name));
        write_cds_fasta(&out_path, &h.sequences)?;
        let sha = sha256_file(&out_path)?;
        let n_applied = personalized
            .applied
            .iter()
            .filter(|a| a.haplotype == h.name)
            .count();
        let n_modified = personalized
            .applied
            .iter()
            .filter(|a| a.haplotype == h.name)
            .map(|a| a.transcript_id.clone())
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        eprintln!(
            "  ✓ {} personalized transcriptome → {} ({} variants applied across {} transcripts)",
            h.name,
            out_path.display(),
            n_applied,
            n_modified
        );
        let _ = hap_idx;
        hap_fastas.push((h.name.clone(), out_path));
        hap_sha.push(sha);
        hap_n_applied.push(n_applied);
        hap_n_modified.push(n_modified);
    }

    // 6. Write the audit trails (variants_applied.tsv, variants_skipped.tsv).
    write_applied_tsv(
        &cfg.output_dir.join("variants_applied.tsv"),
        &personalized.applied,
    )?;
    write_skipped_tsv(
        &cfg.output_dir.join("variants_skipped.tsv"),
        &personalized.skipped,
        &vcf_parsed.skipped,
    )?;

    // 7. Build salmon indices per haplotype. Reuses the existing
    //    SHA256-keyed cache.
    let cache_env = CacheEnv::from_process(cfg.index_cache.clone());
    let cache_root = resolve_cache_root(&cache_env)?;
    let tmpdir = cfg.tmpdir.clone().unwrap_or_else(std::env::temp_dir);

    let ase_spec = AseQuantSpec {
        haplotype_cds: &hap_fastas,
        genome_fa: &cfg.genome_fa,
        output_dir: &cfg.output_dir,
        index_cache_root: &cache_root,
        k: cfg.k,
        threads: cfg.threads,
        tmpdir: &tmpdir,
        fastp_bin: cfg.fastp_bin.as_deref(),
        salmon_bin: cfg.salmon_bin.as_deref(),
        asymmetry_threshold: cfg.asymmetry_threshold,
    };
    let indices = build_haplotype_indices(&ase_spec)?;

    // 8. Per-sample quant (fastp once + salmon twice).
    let fastp_subdir = cfg.output_dir.join("fastp");
    std::fs::create_dir_all(&fastp_subdir)?;
    let mut per_sample_results = Vec::with_capacity(sheet.samples.len());
    let mut merge_inputs: Vec<HapQuantInput> = Vec::new();

    for sample in &sheet.samples {
        eprintln!("  → sample {}", sample.sample_id);
        let r = quantify_sample(sample, &indices, &ase_spec)?;
        // Persist fastp JSON report at ase_out/fastp/<sample_id>.json.
        let persisted_json = fastp_subdir.join(format!("{}.json", sample.sample_id));
        std::fs::copy(&r.fastp.json_report, &persisted_json)?;
        let _ = std::fs::copy(
            &r.fastp.html_report,
            fastp_subdir.join(format!("{}.html", sample.sample_id)),
        );
        for hr in &r.hap_results {
            merge_inputs.push(HapQuantInput {
                sample_id: sample.sample_id.clone(),
                haplotype: hr.haplotype.clone(),
                quant_sf: hr.quant.quant_sf.clone(),
            });
        }
        per_sample_results.push((sample.clone(), r, persisted_json));
    }

    // 9. Merge → allele matrices + summary.
    let matrix = merge::merge_allele_quant(&merge_inputs)?;
    let counts_path = cfg.output_dir.join("ase_counts.tsv");
    let tpm_path = cfg.output_dir.join("ase_tpm.tsv");
    matrix.write_counts(&counts_path)?;
    matrix.write_tpm(&tpm_path)?;

    let (variants_per_tx_hap0, variants_per_tx_hap1) =
        count_variants_per_transcript(&personalized.applied, &cfg.haplotype_names);
    let summary_rows = merge::build_summary(
        &matrix,
        &variants_per_tx_hap0,
        &variants_per_tx_hap1,
        &cfg.haplotype_names,
    );
    merge::write_summary(&cfg.output_dir.join("ase_summary.tsv"), &summary_rows)?;
    let n_informative = summary_rows.iter().filter(|r| r.informative).count();
    eprintln!(
        "  ✓ merged {} transcripts × {} sample-hap columns ({} informative for ASE)",
        matrix.transcripts.len(),
        matrix.columns.len(),
        n_informative
    );

    // 10. Copy the sample sheet (matches src/quant/ pattern).
    std::fs::copy(&cfg.samples_sheet, cfg.output_dir.join("sample_sheet.tsv"))?;

    // 11. Reproducibility bundle.
    let bundle = build_bundle(
        cfg,
        &hap_fastas,
        &hap_sha,
        &hap_n_applied,
        &hap_n_modified,
        &vcf_parsed.variants,
        &vcf_parsed.skipped,
        &personalized.applied,
        &personalized.skipped,
        &per_sample_results,
        &indices,
    )?;
    bundle.write(&cfg.output_dir.join("ase_bundle.json"))?;

    eprintln!("\nNext step — run the binomial ASE test in R:");
    eprintln!(
        "  myconote-cli ase-template --ase-dir {} -o ase_test.R",
        cfg.output_dir.display()
    );
    eprintln!("  Rscript ase_test.R");

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn require_file(p: &Path, role: &str) -> Result<()> {
    if !p.is_file() {
        return Err(MycoNoteError::QuantSheet(format!(
            "{role}: file not found at {}",
            p.display()
        )));
    }
    Ok(())
}

fn require_all_fastqs(samples: &[Sample]) -> Result<()> {
    let mut missing: Vec<String> = Vec::new();
    for s in samples {
        if !s.fastq_r1.is_file() {
            missing.push(format!("{}: r1 {}", s.sample_id, s.fastq_r1.display()));
        }
        if let Some(ref r2) = s.fastq_r2 {
            if !r2.is_file() {
                missing.push(format!("{}: r2 {}", s.sample_id, r2.display()));
            }
        }
    }
    if !missing.is_empty() {
        return Err(MycoNoteError::QuantSheet(format!(
            "{} missing FASTQ file(s):\n  - {}",
            missing.len(),
            missing.join("\n  - ")
        )));
    }
    Ok(())
}

fn write_applied_tsv(path: &Path, applied: &[VariantApplication]) -> Result<()> {
    use std::io::Write as _;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "transcript\thaplotype\tchrom\tgenome_pos\tcds_pos\tref\talt\tkind\tstrand"
    )?;
    for a in applied {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{}",
            a.transcript_id,
            a.haplotype,
            a.chrom,
            a.genome_pos,
            a.cds_pos,
            a.ref_allele,
            a.alt_allele,
            a.kind,
            a.strand
        )?;
    }
    Ok(())
}

fn write_skipped_tsv(
    path: &Path,
    per_tx_skipped: &[personalize::VariantSkip],
    vcf_skipped: &[crate::ase::vcf::SkippedVariant],
) -> Result<()> {
    use std::io::Write as _;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "transcript\thaplotype\tchrom\tgenome_pos\treason\tdetail"
    )?;
    // VCF-level skips (no transcript context) — record once with
    // empty transcript/haplotype.
    for s in vcf_skipped {
        let (reason, detail) = match &s.reason {
            SkipReason::StructuralVariant => ("structural_variant".to_string(), String::new()),
            SkipReason::IndelTooLarge {
                max_allowed,
                observed,
            } => (
                "indel_too_large".to_string(),
                format!("max={max_allowed},observed={observed}"),
            ),
            SkipReason::MultiAllelic => ("multi_allelic".to_string(), String::new()),
            SkipReason::MissingGenotype => ("missing_genotype".to_string(), String::new()),
            SkipReason::NotDiploid => ("not_diploid".to_string(), String::new()),
        };
        writeln!(w, "\t\t{}\t{}\t{}\t{}", s.chrom, s.pos, reason, detail)?;
    }
    // Per-transcript skips.
    for s in per_tx_skipped {
        let (reason, detail) = match &s.reason {
            personalize::ApplicationSkipReason::OutsideCds => {
                ("outside_cds".to_string(), String::new())
            }
            personalize::ApplicationSkipReason::SpansExonBoundary => {
                ("spans_exon_boundary".to_string(), String::new())
            }
            personalize::ApplicationSkipReason::RefMismatch { observed, expected } => (
                "ref_mismatch".to_string(),
                format!("observed={observed},expected={expected}"),
            ),
            personalize::ApplicationSkipReason::InCisOverlap => {
                ("in_cis_overlap".to_string(), String::new())
            }
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            s.transcript_id, s.haplotype, s.chrom, s.genome_pos, reason, detail
        )?;
    }
    Ok(())
}

fn count_variants_per_transcript(
    applied: &[VariantApplication],
    hap_names: &[String; 2],
) -> (HashMap<String, usize>, HashMap<String, usize>) {
    let mut h0: HashMap<String, usize> = HashMap::new();
    let mut h1: HashMap<String, usize> = HashMap::new();
    for a in applied {
        let target = if a.haplotype == hap_names[0] {
            &mut h0
        } else if a.haplotype == hap_names[1] {
            &mut h1
        } else {
            continue;
        };
        *target.entry(a.transcript_id.clone()).or_insert(0) += 1;
    }
    (h0, h1)
}

#[allow(clippy::too_many_arguments)]
fn build_bundle(
    cfg: &AseConfig,
    hap_fastas: &[(String, PathBuf)],
    hap_sha: &[String],
    hap_n_applied: &[usize],
    hap_n_modified: &[usize],
    vcf_variants: &[PhasedVariant],
    vcf_skipped: &[crate::ase::vcf::SkippedVariant],
    applied: &[VariantApplication],
    per_tx_skipped: &[personalize::VariantSkip],
    per_sample_results: &[(Sample, quant::SampleAseResult, PathBuf)],
    indices: &[(String, crate::quant::index::IndexResult)],
) -> Result<bundle::AseBundle> {
    let vcf_sha = sha256_file(&cfg.vcf)?;
    let cds_sha = sha256_file(&cfg.cds_fa)?;
    let genome_sha = sha256_file(&cfg.genome_fa)?;
    let sheet_sha = sha256_file(&cfg.samples_sheet)?;

    let salmon_version = indices
        .first()
        .map(|(_, r)| r.salmon_version.clone())
        .unwrap_or_else(|| "unknown".to_string());
    let fastp_version = detect_fastp_version(cfg.fastp_bin.as_deref().unwrap_or("fastp"))
        .unwrap_or_else(|_| "unknown".to_string());

    let personalized_transcriptomes: Vec<bundle::HaplotypeTranscriptome> = hap_fastas
        .iter()
        .enumerate()
        .map(|(i, (name, path))| bundle::HaplotypeTranscriptome {
            name: name.clone(),
            path: path.to_string_lossy().into_owned(),
            sha256: hap_sha[i].clone(),
            n_variants_applied: hap_n_applied[i],
            n_transcripts_modified: hap_n_modified[i],
        })
        .collect();

    // Count VCF-level skips by reason.
    let mut s_unphased = 0usize;
    let mut s_structural = 0usize;
    let mut s_multi = 0usize;
    let mut s_missing_gt = 0usize;
    let mut s_indel_too_large = 0usize;
    for s in vcf_skipped {
        match s.reason {
            SkipReason::StructuralVariant => s_structural += 1,
            SkipReason::MultiAllelic => s_multi += 1,
            SkipReason::MissingGenotype => s_missing_gt += 1,
            SkipReason::IndelTooLarge { .. } => s_indel_too_large += 1,
            SkipReason::NotDiploid => s_unphased += 1,
        }
    }
    // Count per-transcript skip reasons.
    let mut s_outside_cds = 0usize;
    let mut s_spans_exon = 0usize;
    let mut s_ref_mismatch = 0usize;
    let mut s_in_cis = 0usize;
    for s in per_tx_skipped {
        match s.reason {
            personalize::ApplicationSkipReason::OutsideCds => s_outside_cds += 1,
            personalize::ApplicationSkipReason::SpansExonBoundary => s_spans_exon += 1,
            personalize::ApplicationSkipReason::RefMismatch { .. } => s_ref_mismatch += 1,
            personalize::ApplicationSkipReason::InCisOverlap => s_in_cis += 1,
        }
    }
    let applied_hap0 = applied
        .iter()
        .filter(|a| a.haplotype == cfg.haplotype_names[0])
        .count();
    let applied_hap1 = applied
        .iter()
        .filter(|a| a.haplotype == cfg.haplotype_names[1])
        .count();

    let total_in_vcf = vcf_variants.len() + vcf_skipped.len();
    let _ = summarize_variants(vcf_variants);

    let mut samples: Vec<bundle::AseSampleSummary> = Vec::new();
    for (_sample, r, persisted_json) in per_sample_results {
        let fsum = extract_fastp_summary(persisted_json)?;
        let h0 = r
            .hap_results
            .iter()
            .find(|h| h.haplotype == cfg.haplotype_names[0]);
        let h1 = r
            .hap_results
            .iter()
            .find(|h| h.haplotype == cfg.haplotype_names[1]);
        samples.push(bundle::AseSampleSummary {
            sample_id: r.sample_id.clone(),
            mapping_rate_hap0: h0.map(|h| h.quant.mapping_rate).unwrap_or(0.0),
            mapping_rate_hap1: h1.map(|h| h.quant.mapping_rate).unwrap_or(0.0),
            library_size_hap0: h0.map(|h| h.quant.library_size).unwrap_or(0),
            library_size_hap1: h1.map(|h| h.quant.library_size).unwrap_or(0),
            asymmetry_flag: r.asymmetry_flag,
            fastp: fsum,
        });
    }

    Ok(bundle::AseBundle {
        version: env!("CARGO_PKG_VERSION").to_string(),
        run_timestamp: timestamp_rfc3339(),
        command: cfg.raw_cmdline.clone(),
        external_tools: ExternalTools {
            fastp: fastp_version,
            salmon: salmon_version,
        },
        inputs: bundle::AseBundleInputs {
            cds_fa: HashedPath {
                path: cfg.cds_fa.to_string_lossy().into_owned(),
                sha256: cds_sha,
            },
            genome_fa: HashedPath {
                path: cfg.genome_fa.to_string_lossy().into_owned(),
                sha256: genome_sha,
            },
            vcf: HashedPath {
                path: cfg.vcf.to_string_lossy().into_owned(),
                sha256: vcf_sha,
            },
            sample_sheet: HashedPath {
                path: cfg.samples_sheet.to_string_lossy().into_owned(),
                sha256: sheet_sha,
            },
        },
        haplotypes: cfg.haplotype_names.to_vec(),
        personalized_transcriptomes,
        variants: bundle::VariantApplicationSummary {
            total_in_vcf,
            applied_hap0,
            applied_hap1,
            skipped_unphased: s_unphased,
            skipped_structural: s_structural,
            skipped_multiallelic: s_multi,
            skipped_missing_genotype: s_missing_gt,
            skipped_indel_too_large: s_indel_too_large,
            skipped_outside_cds: s_outside_cds,
            skipped_spans_exon_boundary: s_spans_exon,
            skipped_ref_mismatch: s_ref_mismatch,
            skipped_in_cis_overlap: s_in_cis,
        },
        samples,
        seed: cfg.seed,
    })
}

fn detect_fastp_version(fastp_bin: &str) -> Result<String> {
    let output = duct::cmd!(fastp_bin, "--version")
        .stderr_to_stdout()
        .read()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "fastp".to_string(),
            message: format!("running '{fastp_bin} --version': {e}"),
        })?;
    for line in output.lines() {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.first().map(|t| t.eq_ignore_ascii_case("fastp")) == Some(true) && toks.len() >= 2 {
            return Ok(toks[1].trim_start_matches('v').to_string());
        }
    }
    Err(MycoNoteError::QuantTool {
        tool: "fastp".to_string(),
        message: format!("unexpected '--version' output: {}", output.trim()),
    })
}

fn timestamp_rfc3339() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    format_rfc3339(secs)
}

fn format_rfc3339(secs: i64) -> String {
    let days = secs.div_euclid(86400);
    let sod = secs.rem_euclid(86400);
    let hour = sod / 3600;
    let minute = (sod / 60) % 60;
    let second = sod % 60;
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        y, m, d, hour, minute, second
    )
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(s: &str) -> Vec<String> {
        s.split('|').map(|p| p.trim().to_string()).collect()
    }

    #[test]
    fn parse_minimal_ok() {
        let a = argv("cds.fa|--vcf|v.vcf.gz|--gff3|annotated.gff3|--samples|s.tsv|--genome|g.fa");
        let cfg = AseConfig::parse(&a).unwrap();
        assert_eq!(cfg.cds_fa, PathBuf::from("cds.fa"));
        assert_eq!(cfg.vcf, PathBuf::from("v.vcf.gz"));
        assert_eq!(cfg.gff3, PathBuf::from("annotated.gff3"));
        assert_eq!(cfg.samples_sheet, PathBuf::from("s.tsv"));
        assert_eq!(cfg.genome_fa, PathBuf::from("g.fa"));
        assert_eq!(
            cfg.haplotype_names,
            ["hap0".to_string(), "hap1".to_string()]
        );
        assert_eq!(cfg.k, 31);
        assert_eq!(cfg.max_indel_size, 50);
        assert!((cfg.asymmetry_threshold - 0.05).abs() < 1e-9);
    }

    #[test]
    fn parse_custom_haplotype_names() {
        let a = argv("cds.fa|--vcf|v|--gff3|g3|--samples|s|--genome|g|--haplotype-names|cer,par");
        let cfg = AseConfig::parse(&a).unwrap();
        assert_eq!(cfg.haplotype_names, ["cer".to_string(), "par".to_string()]);
    }

    #[test]
    fn parse_rejects_duplicate_haplotype_names() {
        let a = argv("cds.fa|--vcf|v|--gff3|g3|--samples|s|--genome|g|--haplotype-names|same,same");
        let err = AseConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("must differ"));
    }

    #[test]
    fn parse_rejects_missing_vcf() {
        let a = argv("cds.fa|--gff3|g3|--samples|s|--genome|g");
        let err = AseConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--vcf"));
    }

    #[test]
    fn parse_rejects_missing_gff3() {
        let a = argv("cds.fa|--vcf|v|--samples|s|--genome|g");
        let err = AseConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("--gff3"));
    }

    #[test]
    fn parse_rejects_missing_cds() {
        let a = argv("--vcf|v|--gff3|g3|--samples|s|--genome|g");
        let err = AseConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("missing CDS FASTA"));
    }

    #[test]
    fn parse_rejects_bad_haplotype_names_format() {
        let a = argv("cds.fa|--vcf|v|--gff3|g3|--samples|s|--genome|g|--haplotype-names|justone");
        let err = AseConfig::parse(&a).unwrap_err();
        assert!(format!("{err}").contains("'NAME1,NAME2'"));
    }

    #[test]
    fn parse_all_flags() {
        let a = argv(
            "cds.fa|--vcf|v|--gff3|g3|--samples|s|--genome|g|\
             --output|out|--haplotype-names|A,B|-k|25|-t|8|\
             --tmpdir|/tmp/x|--index-cache|/c|--keep-trimmed|/kt|\
             --max-indel-size|100|--asymmetry-threshold|0.1|\
             --seed|7|--fastp|/bin/fastp|--salmon|/bin/salmon",
        );
        let cfg = AseConfig::parse(&a).unwrap();
        assert_eq!(cfg.output_dir, PathBuf::from("out"));
        assert_eq!(cfg.haplotype_names, ["A".to_string(), "B".to_string()]);
        assert_eq!(cfg.k, 25);
        assert_eq!(cfg.threads, 8);
        assert_eq!(cfg.max_indel_size, 100);
        assert!((cfg.asymmetry_threshold - 0.1).abs() < 1e-9);
        assert_eq!(cfg.seed, 7);
        assert_eq!(cfg.fastp_bin.as_deref(), Some("/bin/fastp"));
        assert_eq!(cfg.salmon_bin.as_deref(), Some("/bin/salmon"));
    }

    #[test]
    fn rfc3339_epoch_and_known() {
        assert_eq!(format_rfc3339(0), "1970-01-01T00:00:00Z");
        // 2026-04-24T12:00:00Z
        assert_eq!(format_rfc3339(1777032000), "2026-04-24T12:00:00Z");
    }
}
