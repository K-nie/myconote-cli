/// myconote-cli integration tests
///
/// Run with:  cargo test
/// Run single: cargo test test_stats_basic
///
/// All tests that invoke the CLI binary use `assert_cmd` so they work
/// regardless of whether the binary is in PATH.
use assert_cmd::Command;
use predicates::prelude::*;
use std::path::PathBuf;
use tempfile::TempDir;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn bin() -> Command {
    Command::cargo_bin("myconote-cli").expect("binary not found — run `cargo build` first")
}

fn gff3() -> PathBuf {
    PathBuf::from("tests/data/candida_tropicalis.final.gff3")
}

fn fasta() -> PathBuf {
    PathBuf::from("tests/data/candida_tropicalis.fas")
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Top-level help
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_help_no_args() {
    bin()
        .assert()
        .success()
        .stdout(predicate::str::contains("myconote-cli"))
        .stdout(predicate::str::contains("predict"))
        .stdout(predicate::str::contains("annotate"))
        .stdout(predicate::str::contains("check"))
        .stdout(predicate::str::contains("setup"))
        .stdout(predicate::str::contains("species"));
}

#[test]
fn test_unknown_command() {
    bin()
        .arg("notacommand")
        .assert()
        .stdout(predicate::str::contains("Unknown command"));
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. stats command
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_stats_basic() {
    bin()
        .arg("stats")
        .arg(gff3())
        .assert()
        .success()
        .stdout(predicate::str::contains("Gene"))
        .stdout(predicate::str::contains("6290").or(predicate::str::contains("6,290")));
}

#[test]
#[ignore = "requires stats JSON output fix"]
fn test_stats_json_format() {
    bin()
        .args(["stats", "--format", "json"])
        .arg(gff3())
        .assert()
        .success()
        .stdout(predicate::str::contains("{"))
        .stdout(predicate::str::contains("gene_count"));
}

#[test]
#[ignore = "requires stats CSV output fix"]
fn test_stats_csv_format() {
    bin()
        .args(["stats", "--format", "csv"])
        .arg(gff3())
        .assert()
        .success()
        .stdout(predicate::str::contains(","));
}

#[test]
fn test_stats_fungi_taxon() {
    bin()
        .arg("stats")
        .arg(gff3())
        .args(["--taxon", "fungi"])
        .assert()
        .success();
}

#[test]
fn test_stats_missing_file() {
    bin().args(["stats", "nonexistent.gff3"]).assert().failure();
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. convert command — GFF3 outputs (pure Rust, no external tools)
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_convert_gff3_to_gtf() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.gtf");

    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "gtf", "-o"])
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists(), "output GTF not created");
    let content = std::fs::read_to_string(&out).unwrap();
    assert!(
        content.contains("transcript_id"),
        "GTF missing transcript_id"
    );
}

#[test]
fn test_convert_gff3_to_bed() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.bed");

    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "bed", "-o"])
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists());
    let content = std::fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = content.lines().filter(|l| !l.is_empty()).collect();
    // BED should have at least as many lines as genes
    assert!(lines.len() >= 100, "BED has too few lines: {}", lines.len());
    // Each BED line has 6 tab-separated fields
    let first = lines[0];
    assert_eq!(
        first.split('\t').count(),
        6,
        "BED line doesn't have 6 fields: {}",
        first
    );
}

#[test]
fn test_convert_gff3_to_bed12() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.bed12");

    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "bed12", "-o"])
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists());
}

#[test]
fn test_convert_gff3_to_table() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.tsv");

    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "table", "-o"])
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists());
    let content = std::fs::read_to_string(&out).unwrap();
    // TSV table should have a header and data rows
    assert!(content.contains('\t'), "TSV has no tabs");
}

#[test]
fn test_convert_gff3_to_genbank() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.gbk");

    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "genbank", "--fasta"])
        .arg(fasta())
        .arg("-o")
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists(), "GenBank file not created");
    let content = std::fs::read_to_string(&out).unwrap();
    assert!(content.contains("LOCUS"), "GenBank missing LOCUS line");
    assert!(content.contains("CDS"), "GenBank missing CDS features");
}

#[test]
fn test_convert_gff3_to_protein() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("proteins.faa");

    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "protein", "--fasta"])
        .arg(fasta())
        .arg("-o")
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists(), "Protein FASTA not created");
    let content = std::fs::read_to_string(&out).unwrap();
    assert!(content.contains('>'), "Protein FASTA has no sequences");
    // Should have a reasonable number of proteins
    let n = content.lines().filter(|l| l.starts_with('>')).count();
    assert!(n > 100, "Too few proteins extracted: {}", n);
}

#[test]
fn test_convert_gff3_to_cds() {
    // Extracts spliced CDS nucleotide sequences.
    // Asserts: (a) at least one record written;
    //          (b) g001568.m1 length == sum of its two CDS segment lengths
    //              (133 bp + 500 bp = 633 bp, phase 0 on transcript-first CDS);
    //          (c) g001559.m1 (minus strand, single CDS 31–1398, phase 0)
    //              is the reverse-complement of the raw forward-strand slice.
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("cds.fna");

    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "cds", "--fasta"])
        .arg(fasta())
        .arg("-o")
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists(), "CDS FASTA not created");
    let content = std::fs::read_to_string(&out).unwrap();
    let n_records = content.lines().filter(|l| l.starts_with('>')).count();
    assert!(n_records > 100, "Too few CDS extracted: {}", n_records);

    // Parse the emitted FASTA into (id -> sequence). Simple inline parse —
    // the test has no need to depend on the internal FastaReader.
    let mut records: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut cur_id = String::new();
    let mut cur_seq = String::new();
    for line in content.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !cur_id.is_empty() {
                records.insert(cur_id.clone(), cur_seq.clone());
            }
            cur_id = rest.split_whitespace().next().unwrap_or("").to_string();
            cur_seq.clear();
        } else {
            cur_seq.push_str(line.trim());
        }
    }
    if !cur_id.is_empty() {
        records.insert(cur_id, cur_seq);
    }

    // (b) Multi-CDS positive-strand length check.
    // g001568.m1: CDS 11031–11163 (133 bp) + CDS 11246–11745 (500 bp) = 633 bp.
    let g568 = records
        .get("g001568.m1")
        .expect("g001568.m1 should be present");
    assert_eq!(
        g568.len(),
        633,
        "g001568.m1 spliced CDS length should be 633 bp (133 + 500), got {}",
        g568.len()
    );

    // (c) Negative-strand reverse-complement check.
    // g001559.m1: single CDS on NW_003020049.1, 31–1398, strand '-', phase 0.
    // The emitted sequence must equal RC of the forward-strand slice.
    let g559 = records
        .get("g001559.m1")
        .expect("g001559.m1 should be present");
    assert_eq!(
        g559.len(),
        1398 - 31 + 1,
        "g001559.m1 length should equal CDS span (1368 bp)"
    );

    // Pull the forward-strand slice directly from the fixture FASTA
    // and reverse-complement it; must match what the tool emitted.
    let genome_text = std::fs::read_to_string(fasta()).unwrap();
    let mut target_seq = String::new();
    let mut in_target = false;
    for line in genome_text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            let id = rest.split_whitespace().next().unwrap_or("");
            in_target = id == "NW_003020049.1";
            continue;
        }
        if in_target {
            target_seq.push_str(line.trim());
        }
    }
    // 1-based inclusive [31, 1398] → 0-based Rust slice [30, 1398]
    let forward = target_seq[30..1398].to_ascii_uppercase();
    let rc: String = forward
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'N' => 'N',
            other => other,
        })
        .collect();
    assert_eq!(
        g559, &rc,
        "g001559.m1 emitted sequence is not the reverse-complement of the forward slice"
    );
}

#[test]
#[ignore = "error message format differs across platforms"]
fn test_convert_missing_fasta_for_protein() {
    // Should fail gracefully with a message, not panic
    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "protein"])
        .assert()
        .stdout(predicate::str::contains("--fasta").or(predicate::str::contains("required")));
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. clean command
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_clean_basic() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("clean.gff3");

    bin()
        .arg("clean")
        .arg(gff3())
        .arg("-o")
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists());
    let n_genes = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|l| l.contains("\tgene\t"))
        .count();
    assert!(n_genes > 6000, "clean removed too many genes: {}", n_genes);
}

#[test]
fn test_clean_remove_orphans() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("clean_no_orphans.gff3");

    bin()
        .arg("clean")
        .arg(gff3())
        .args(["--remove-orphans", "-o"])
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists());
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. sort command (pure Rust, no external tools)
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_sort_basic() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("sorted.fa");

    bin()
        .arg("sort")
        .arg(fasta())
        .arg("-o")
        .arg(&out)
        .assert()
        .success();

    assert!(out.exists(), "Sorted FASTA not created");
    let content = std::fs::read_to_string(&out).unwrap();
    // Should still have all 24 sequences
    let n = content.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(n, 24, "Expected 24 sequences, got {}", n);
    // New default (v0.7.9): original seqids are PRESERVED — no scaffold_N
    // renaming happens unless --ncbi-clean is passed.
    assert!(
        !content.lines().any(|l| l.starts_with(">scaffold_")),
        "sort must preserve original seqids by default (no scaffold_ renaming)"
    );
}

#[test]
fn test_sort_custom_prefix() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("sorted_chr.fa");

    // --ncbi-clean opts into the scaffold_N-style renaming; --prefix sets the
    // prefix. Without --ncbi-clean the prefix is inert (names are preserved).
    bin()
        .arg("sort")
        .arg(fasta())
        .args(["--ncbi-clean", "--prefix", "chr", "-o"])
        .arg(&out)
        .assert()
        .success();

    let content = std::fs::read_to_string(&out).unwrap();
    assert!(content.contains(">chr_"), "Expected chr_ prefix");
}

#[test]
fn test_sort_min_length() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("sorted_long.fa");

    bin()
        .arg("sort")
        .arg(fasta())
        .args(["--min-length", "500000", "-o"])
        .arg(&out)
        .assert()
        .success();

    // With a high min-length, some short scaffolds should be filtered
    let content = std::fs::read_to_string(&out).unwrap();
    let n = content.lines().filter(|l| l.starts_with('>')).count();
    assert!(
        n < 24,
        "Expected some scaffolds filtered out, got {} (all 24)",
        n
    );
}

#[test]
#[ignore = "rename table count differs with duplicate scaffold names"]
fn test_sort_rename_table() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("sorted.fa");
    let table = tmp.path().join("rename.tsv");

    bin()
        .arg("sort")
        .arg(fasta())
        .args(["--ncbi-clean", "--rename-table"])
        .arg(&table)
        .arg("-o")
        .arg(&out)
        .assert()
        .success();

    assert!(table.exists(), "Rename table not created");
    let content = std::fs::read_to_string(&table).unwrap();
    // TSV: old_id <tab> new_id
    assert!(content.contains('\t'), "Rename table has no tabs");
    let n = content.lines().filter(|l| !l.is_empty()).count();
    assert_eq!(n, 24, "Expected 24 rename entries, got {}", n);
}

// ─────────────────────────────────────────────────────────────────────────────
// 8. fix command — GenBank repair (requires convert output first)
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_fix_genbank_roundtrip() {
    let tmp = TempDir::new().unwrap();
    let gbk = tmp.path().join("annotation.gbk");
    let fixed = tmp.path().join("annotation_fixed.gbk");
    let report = tmp.path().join("fix_report.txt");

    // First create a GenBank file from our test data
    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "genbank", "--fasta"])
        .arg(fasta())
        .arg("-o")
        .arg(&gbk)
        .assert()
        .success();

    assert!(gbk.exists());

    // Now run fix on it
    bin()
        .arg("fix")
        .arg(&gbk)
        .arg("-o")
        .arg(&fixed)
        .args(["--report"])
        .arg(&report)
        .assert()
        .success();

    assert!(fixed.exists(), "Fixed GenBank not created");
    assert!(report.exists(), "Fix report not created");

    let report_text = std::fs::read_to_string(&report).unwrap();
    assert!(
        report_text.contains("Records processed"),
        "Report missing stats"
    );
}

#[test]
fn test_fix_dry_run() {
    let tmp = TempDir::new().unwrap();
    let gbk = tmp.path().join("annotation.gbk");
    let fixed = tmp.path().join("should_not_exist.gbk");

    bin()
        .arg("convert")
        .arg(gff3())
        .args(["--to", "genbank", "--fasta"])
        .arg(fasta())
        .arg("-o")
        .arg(&gbk)
        .assert()
        .success();

    bin()
        .arg("fix")
        .arg(&gbk)
        .arg("-o")
        .arg(&fixed)
        .arg("--dry-run")
        .assert()
        .success();

    assert!(!fixed.exists(), "dry-run should not produce output file");
}

// ─────────────────────────────────────────────────────────────────────────────
// 9. check command (no external tools needed)
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_check_runs() {
    bin()
        .arg("check")
        .assert()
        .success()
        .stdout(predicate::str::contains("augustus").or(predicate::str::contains("Tool")));
}

#[test]
fn test_check_filter_by_command() {
    bin()
        .args(["check", "annotate"])
        .assert()
        .success()
        .stdout(predicate::str::contains("diamond").or(predicate::str::contains("hmmscan")));
}

// ─────────────────────────────────────────────────────────────────────────────
// 10. setup command (no download, just list/check)
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_setup_list() {
    bin()
        .args(["setup", "--list"])
        .assert()
        .success()
        .stdout(predicate::str::contains("swiss-prot"))
        .stdout(predicate::str::contains("pfam"))
        .stdout(predicate::str::contains("merops"));
}

#[test]
fn test_setup_check() {
    let tmp = TempDir::new().unwrap();

    bin()
        .args(["setup", "--check", "--db-dir"])
        .arg(tmp.path())
        .assert()
        .success()
        .stdout(predicate::str::contains("swiss-prot"));
}

// ─────────────────────────────────────────────────────────────────────────────
// 11. species command (no external tools needed unless Augustus is installed)
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_species_runs() {
    // May find 0 species if Augustus not installed — should still succeed
    bin().arg("species").assert().success();
}

#[test]
fn test_species_grouped() {
    bin()
        .args(["species", "--grouped"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Fungi").or(predicate::str::contains("fungi")));
}

// ─────────────────────────────────────────────────────────────────────────────
// 12. Subcommand help text
// ─────────────────────────────────────────────────────────────────────────────

macro_rules! help_test {
    ($name:ident, $cmd:expr, $expected:expr) => {
        #[test]
        fn $name() {
            bin()
                .arg($cmd)
                .assert()
                .stdout(predicate::str::contains($expected));
        }
    };
}

help_test!(test_help_sort, "sort", "Usage");
help_test!(test_help_mask, "mask", "Usage");
help_test!(test_help_train, "train", "Usage");
help_test!(test_help_predict, "predict", "Usage");
help_test!(test_help_update, "update", "Usage");
help_test!(test_help_annotate, "annotate", "Usage");
help_test!(test_help_remote, "remote", "Usage");
help_test!(test_help_fix, "fix", "Usage");
help_test!(test_help_convert, "convert", "Usage");
help_test!(test_help_clean, "clean", "Usage");

// ─────────────────────────────────────────────────────────────────────────────
// Explain command integration tests
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_explain_help() {
    bin()
        .args(["explain", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("myconote explain"))
        .stdout(predicate::str::contains("Stages"))
        .stdout(predicate::str::contains("predict"))
        .stdout(predicate::str::contains("--no-llm"))
        .stdout(predicate::str::contains("--paste"));
}

#[test]
fn test_explain_unknown_stage() {
    bin()
        .args(["explain", "bogus"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("unknown stage"));
}

#[test]
fn test_explain_predict_no_llm() {
    // Run explain in rules-only mode on a temp directory
    let dir = TempDir::new().unwrap();

    // Create a minimal predict directory structure
    let predict_dir = dir.path().join("02_predict");
    std::fs::create_dir_all(&predict_dir).unwrap();
    std::fs::write(
        predict_dir.join("predict_summary.txt"),
        "Gene predictions summary\n  Augustus: 5000 genes\n  Consensus: 4800 genes\n",
    )
    .unwrap();

    bin()
        .args([
            "explain",
            "predict",
            "--no-llm",
            "--dir",
            predict_dir.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Rules-only mode"))
        .stdout(predicate::str::contains("suggestive, not definitive"));
}

#[test]
fn test_explain_help_shows_disclaimer() {
    bin()
        .args(["explain", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("suggestive, not definitive"));
}

#[test]
fn test_explain_predict_no_llm_with_gff3() {
    // Test with an actual GFF3 file
    if !gff3().exists() {
        return; // skip if test data not available
    }

    let dir = TempDir::new().unwrap();
    let predict_dir = dir.path().join("02_predict");
    std::fs::create_dir_all(&predict_dir).unwrap();

    // Copy the test GFF3
    std::fs::copy(gff3(), predict_dir.join("consensus.gff3")).unwrap();

    bin()
        .args([
            "explain",
            "predict",
            "--no-llm",
            "--dir",
            predict_dir.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Findings"));
}

#[test]
fn test_explain_submit_no_llm() {
    let dir = TempDir::new().unwrap();
    let submit_dir = dir.path().join("06_submit");
    std::fs::create_dir_all(&submit_dir).unwrap();

    // Create errorsummary.val with some errors
    std::fs::write(
        submit_dir.join("errorsummary.val"),
        "ERROR: valid [SEQ_FEAT.NoStop] No stop codon found\nWARNING: valid [SEQ_FEAT.NotSpliceConsensus]\n",
    ).unwrap();

    bin()
        .args([
            "explain",
            "submit",
            "--no-llm",
            "--dir",
            submit_dir.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("NCBI validation"));
}

#[test]
fn test_explain_missing_dir() {
    bin()
        .args(["explain", "predict", "--dir", "/nonexistent/path/really"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("directory not found"));
}

#[test]
fn test_explain_dry_run() {
    let dir = TempDir::new().unwrap();

    bin()
        .args([
            "explain",
            "predict",
            "--dry-run",
            "--dir",
            dir.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("--dry-run"))
        .stdout(predicate::str::contains("Prompt length"));
}

#[test]
fn test_explain_mask_no_llm() {
    let dir = TempDir::new().unwrap();
    let mask_dir = dir.path().join("01_mask");
    std::fs::create_dir_all(&mask_dir).unwrap();

    // Create a minimal masked FASTA with some lowercase
    std::fs::write(
        mask_dir.join("genome_masked.fas"),
        ">scaffold_1\nACGTacgtACGTACGTACGTACGT\n",
    )
    .unwrap();

    bin()
        .args([
            "explain",
            "mask",
            "--no-llm",
            "--dir",
            mask_dir.to_str().unwrap(),
        ])
        .assert()
        .success();
}

#[test]
fn test_explain_annotate_no_llm() {
    let dir = TempDir::new().unwrap();
    let ann_dir = dir.path().join("05_annotate");
    std::fs::create_dir_all(&ann_dir).unwrap();

    bin()
        .args([
            "explain",
            "annotate",
            "--no-llm",
            "--dir",
            ann_dir.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Findings"));
}

#[test]
fn test_setup_list_includes_chat_corpus() {
    bin()
        .args(["setup", "--list"])
        .assert()
        .success()
        .stdout(predicate::str::contains("chat-corpus"));
}

// ─────────────────────────────────────────────────────────────────────────────
// Batch subcommand
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_batch_help() {
    bin()
        .args(["batch", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Annotate multiple genomes"))
        .stdout(predicate::str::contains("--condor"))
        .stdout(predicate::str::contains("HTCondor"))
        .stdout(predicate::str::contains("sample sheet"));
}

#[test]
fn test_batch_no_input_shows_help() {
    bin()
        .args(["batch"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Annotate multiple genomes"));
}

#[test]
fn test_batch_empty_directory() {
    let tmp = TempDir::new().unwrap();
    bin()
        .args(["batch", tmp.path().to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("no FASTA files found"));
}

#[test]
fn test_batch_condor_generates_submit_files() {
    let tmp = TempDir::new().unwrap();
    let genomes_dir = tmp.path().join("genomes");
    std::fs::create_dir_all(&genomes_dir).unwrap();
    std::fs::write(genomes_dir.join("test_a.fa"), ">seq1\nACGT\n").unwrap();
    std::fs::write(genomes_dir.join("test_b.fa"), ">seq1\nACGT\n").unwrap();

    let out_dir = tmp.path().join("batch_out");
    bin()
        .args([
            "batch",
            genomes_dir.to_str().unwrap(),
            "--condor",
            "--output",
            out_dir.to_str().unwrap(),
            "--condor-mem",
            "64G",
            "--condor-cpus",
            "16",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("HTCondor submit files generated"))
        .stdout(predicate::str::contains("condor_submit"));

    // Verify files were created
    assert!(out_dir.join("condor.sub").exists());
    assert!(out_dir.join("condor_genomes.txt").exists());
    assert!(out_dir.join("run_genome.sh").exists());
    assert!(out_dir.join("condor_logs").exists());

    // Check submit file content
    let sub = std::fs::read_to_string(out_dir.join("condor.sub")).unwrap();
    assert!(sub.contains("request_cpus = 16"));
    assert!(sub.contains("request_memory = 64G"));
    assert!(sub.contains("queue 2"));
}

#[test]
fn test_batch_sample_sheet_condor() {
    let tmp = TempDir::new().unwrap();

    // Create FASTA files
    let fasta_a = tmp.path().join("genome_a.fa");
    let fasta_b = tmp.path().join("genome_b.fa");
    std::fs::write(&fasta_a, ">seq1\nACGT\n").unwrap();
    std::fs::write(&fasta_b, ">seq1\nACGT\n").unwrap();

    // Create sample sheet
    let sheet = tmp.path().join("samples.tsv");
    let sheet_content = format!(
        "name\tfasta\tkingdom\tspecies\tgenetic_code\tlocus_prefix\n\
         isolate_A\t{}\tfungi\tsaccharomyces\t1\tISOA\n\
         isolate_B\t{}\tfungi\tauto\t12\tISOB\n",
        fasta_a.display(),
        fasta_b.display()
    );
    std::fs::write(&sheet, &sheet_content).unwrap();

    let out_dir = tmp.path().join("batch_out");
    bin()
        .args([
            "batch",
            sheet.to_str().unwrap(),
            "--condor",
            "--output",
            out_dir.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("HTCondor submit files generated"));

    // Verify genome list has both entries
    let list = std::fs::read_to_string(out_dir.join("condor_genomes.txt")).unwrap();
    assert!(list.contains("isolate_A"));
    assert!(list.contains("isolate_B"));
}

#[test]
fn test_batch_nonexistent_input_errors() {
    bin()
        .args(["batch", "/nonexistent/path/to/genomes"])
        .assert()
        .failure();
}

#[test]
fn test_batch_in_main_help() {
    bin()
        .assert()
        .success()
        .stdout(predicate::str::contains("batch"));
}
