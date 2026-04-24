//! End-to-end integration test for `myconote-cli quant`.
//!
//! All tests here are `#[ignore]` by default because they require
//! `salmon` and `fastp` on PATH. Run them with:
//!
//!   PATH=~/miniconda3/envs/fast_myco_env/bin:$PATH \
//!     env -u CC cargo test --release --test quant_tests -- --ignored
//!
//! The tiny synthetic fixture runs in under 5 seconds end-to-end.

use assert_cmd::Command;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

fn bin() -> Command {
    Command::cargo_bin("myconote-cli").expect("binary not found — run `cargo build` first")
}

fn tools_available() -> bool {
    which::which("salmon").is_ok() && which::which("fastp").is_ok()
}

/// Write a tiny gzipped FASTQ whose reads are contiguous 30-bp slices
/// of the supplied reference string. Produces reads that map cleanly
/// against a salmon index built on the same reference.
fn write_reads_derived_from(path: &Path, reference: &str, stride: usize, read_len: usize) {
    let f = fs::File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::fast());
    let mut i = 0;
    while i + read_len <= reference.len() {
        writeln!(enc, "@read_{i}").unwrap();
        writeln!(enc, "{}", &reference[i..i + read_len]).unwrap();
        writeln!(enc, "+").unwrap();
        writeln!(enc, "{}", "I".repeat(read_len)).unwrap();
        i += stride;
    }
    enc.finish().unwrap();
}

fn setup_quant_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    // Two synthetic transcripts that read as distinct sequences, each
    // longer than salmon's k=21. Long enough to permit ~10 30-bp
    // reads per transcript after striding by 3.
    let t1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let t2 = "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG";
    let cds = dir.join("cds.fa");
    fs::write(&cds, format!(">t1\n{t1}\n>t2\n{t2}\n")).unwrap();

    // Decoy genome: unrelated to the CDS, so mapped reads come only
    // from the transcriptome and we exercise the decoy path.
    let genome = dir.join("genome.fa");
    fs::write(
        &genome,
        ">chr1\n\
         GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC\n\
         >chr2\n\
         TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT\n",
    )
    .unwrap();

    // One single-end and one paired-end sample. Reads are exact
    // substrings of t1 (for s1) and t2 (for s2).
    let s1_r1 = dir.join("s1.fq.gz");
    write_reads_derived_from(&s1_r1, t1, 3, 30);
    let s2_r1 = dir.join("s2_R1.fq.gz");
    let s2_r2 = dir.join("s2_R2.fq.gz");
    write_reads_derived_from(&s2_r1, t2, 3, 30);
    write_reads_derived_from(&s2_r2, t2, 3, 30);

    let sheet = dir.join("samples.tsv");
    fs::write(
        &sheet,
        format!(
            "sample_id\tfastq_r1\tfastq_r2\tcondition\tstrandedness\n\
             s1\t{s1_r1}\t\tctrl\tunstranded\n\
             s2\t{s2_r1}\t{s2_r2}\ttrt\tunstranded\n",
            s1_r1 = s1_r1.display(),
            s2_r1 = s2_r1.display(),
            s2_r2 = s2_r2.display(),
        ),
    )
    .unwrap();

    (cds, genome, sheet)
}

#[test]
#[ignore]
fn quant_end_to_end_two_samples() {
    if !tools_available() {
        eprintln!("salmon or fastp missing from PATH — skipping");
        return;
    }

    let tmp = TempDir::new().unwrap();
    let (cds, genome, sheet) = setup_quant_fixture(tmp.path());

    let out = tmp.path().join("quant_out");
    let cache = tmp.path().join("salmon_cache");
    let fastp_tmp = tmp.path().join("fastp_tmp");

    bin()
        .arg("quant")
        .arg(&cds)
        .args(["--samples"])
        .arg(&sheet)
        .args(["--genome"])
        .arg(&genome)
        .args(["--output"])
        .arg(&out)
        .args(["--index-cache"])
        .arg(&cache)
        .args(["--tmpdir"])
        .arg(&fastp_tmp)
        .args(["-k", "21", "--threads", "2"])
        .assert()
        .success();

    // Required output files.
    let counts = out.join("counts.tsv");
    let tpm = out.join("tpm.tsv");
    let bundle = out.join("quant_bundle.json");
    let sheet_copy = out.join("sample_sheet.tsv");
    assert!(counts.is_file(), "counts.tsv missing");
    assert!(tpm.is_file(), "tpm.tsv missing");
    assert!(bundle.is_file(), "quant_bundle.json missing");
    assert!(sheet_copy.is_file(), "sample_sheet.tsv missing");

    // counts.tsv shape: header + at least one data row; two samples.
    let counts_text = fs::read_to_string(&counts).unwrap();
    let lines: Vec<&str> = counts_text.lines().collect();
    assert!(lines.len() >= 2, "counts.tsv too short: {}", counts_text);
    let header = lines[0];
    assert_eq!(
        header.split('\t').count(),
        3,
        "header should be transcript + 2 samples, got: {header}"
    );
    assert!(header.contains("s1") && header.contains("s2"));

    // Per-sample quant.sf present (tximport needs these directly).
    assert!(
        out.join("salmon").join("s1").join("quant.sf").is_file(),
        "s1 quant.sf missing"
    );
    assert!(
        out.join("salmon").join("s2").join("quant.sf").is_file(),
        "s2 quant.sf missing"
    );

    // fastp JSON persisted per sample.
    assert!(out.join("fastp").join("s1.json").is_file());
    assert!(out.join("fastp").join("s2.json").is_file());

    // Bundle sanity: parses as JSON, carries both samples with
    // nonzero library size (reads were drawn from the indexed
    // transcripts, so at least some must have mapped).
    let bundle_text = fs::read_to_string(&bundle).unwrap();
    let bundle_json: serde_json::Value = serde_json::from_str(&bundle_text).unwrap();
    assert_eq!(
        bundle_json["version"].as_str().unwrap_or(""),
        env!("CARGO_PKG_VERSION")
    );
    let samples = bundle_json["samples"].as_array().unwrap();
    assert_eq!(samples.len(), 2);
    for s in samples {
        let lib = s["library_size"].as_u64().unwrap_or(0);
        assert!(
            lib > 0,
            "sample {} has library_size 0 — reads didn't map",
            s["sample_id"]
        );
    }
    assert!(bundle_json["external_tools"]["salmon"]
        .as_str()
        .unwrap()
        .contains("1."));
    assert!(bundle_json["external_tools"]["fastp"]
        .as_str()
        .unwrap()
        .contains("1."));

    // Second run on same inputs → index cache hit (stdout should
    // mention [cached]).
    let second = bin()
        .arg("quant")
        .arg(&cds)
        .args(["--samples"])
        .arg(&sheet)
        .args(["--genome"])
        .arg(&genome)
        .args(["--output"])
        .arg(tmp.path().join("quant_out_rerun"))
        .args(["--index-cache"])
        .arg(&cache)
        .args(["--tmpdir"])
        .arg(&fastp_tmp)
        .args(["-k", "21", "--threads", "2"])
        .output()
        .unwrap();
    assert!(second.status.success(), "rerun failed");
    let stderr = String::from_utf8_lossy(&second.stderr);
    assert!(
        stderr.contains("[cached]"),
        "expected index cache hit on rerun, stderr was:\n{stderr}"
    );
}

#[test]
fn quant_help_shows_usage() {
    // Non-ignored — runs in CI without external tools.
    bin()
        .args(["quant", "--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains("Usage: myconote-cli quant"))
        .stdout(predicates::str::contains("--samples"))
        .stdout(predicates::str::contains("--genome"));
}

#[test]
fn quant_missing_samples_errors() {
    bin()
        .args(["quant", "cds.fa", "--genome", "g.fa"])
        .assert()
        .failure()
        .stderr(predicates::str::contains("--samples"));
}
