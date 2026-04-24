//! Integration tests for `myconote-cli fetch-rna`.
//!
//! CI-safe tests (no network) exercise CLI parsing + error paths.
//! The `#[ignore]`d live tests hit ENA's public REST API and (for the
//! full smoke) download a tiny single-end FASTQ; run with:
//!
//!   env -u CC cargo test --release --test fetch_rna_tests -- --ignored

use assert_cmd::Command;
use predicates::prelude::*;

fn bin() -> Command {
    Command::cargo_bin("myconote-cli").expect("binary not found — run `cargo build` first")
}

#[test]
fn fetch_rna_help_shows_usage() {
    bin()
        .args(["fetch-rna", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Usage: myconote-cli fetch-rna"))
        .stdout(predicate::str::contains("--backend"))
        .stdout(predicate::str::contains("samples.tsv"));
}

#[test]
fn fetch_rna_no_accession_errors() {
    bin()
        .args(["fetch-rna", "--output", "/tmp/rna"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("no accessions"));
}

#[test]
fn fetch_rna_unknown_backend_errors() {
    bin()
        .args(["fetch-rna", "SRR1", "--backend", "narwhal"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("unknown --backend"));
}

#[test]
#[ignore]
fn fetch_rna_dry_run_against_live_ena() {
    // Dry-run path hits ENA for metadata but does NOT download FASTQs.
    // Keeps the test fast (~1s) and independent of download bandwidth.
    bin()
        .args(["fetch-rna", "SRR453566", "--dry-run", "--output"])
        .arg(std::env::temp_dir().join("myconote_fetch_dryrun"))
        .assert()
        .success()
        .stderr(predicate::str::contains("dry-run"))
        .stderr(predicate::str::contains("https://"));
}
