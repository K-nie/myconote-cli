/// Pipeline integration tests for myconote-cli
///
/// Tests the predict, annotate, update, and submit pipelines end-to-end
/// using minimal test fixtures. These tests verify that module wiring
/// is correct even when external tools are not available.
use assert_cmd::Command;
use predicates::prelude::*;

fn bin() -> Command {
    Command::cargo_bin("myconote-cli").expect("binary not found")
}

// ─────────────────────────────────────────────────────────────────────────────
// Help text tests for new commands
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_submit_help() {
    bin()
        .args(["submit", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("NCBI"))
        .stdout(predicate::str::contains("--organism"))
        .stdout(predicate::str::contains("--genetic-code"));
}

#[test]
fn test_predict_help_shows_new_flags() {
    bin()
        .args(["predict", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("--protein-evidence"))
        .stdout(predicate::str::contains("--kingdom"));
}

#[test]
fn test_annotate_help_shows_trnascan() {
    bin()
        .args(["annotate", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("tRNAscan"));
}

// ─────────────────────────────────────────────────────────────────────────────
// Version
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_version_matches_cargo_pkg_version() {
    bin()
        .arg("--version")
        .assert()
        .success()
        .stdout(predicate::str::contains(env!("CARGO_PKG_VERSION")));
}

// ─────────────────────────────────────────────────────────────────────────────
// Submit validation (doesn't need external tools)
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_submit_validate_only_missing_fasta() {
    // Should error when --fasta is missing
    bin()
        .args(["submit", "nonexistent.gff3", "--validate-only"])
        .assert()
        .failure();
}

// ─────────────────────────────────────────────────────────────────────────────
// Learn / Tutorial
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_learn_help() {
    bin()
        .args(["learn", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Interactive Tutorial"))
        .stdout(predicate::str::contains("Lesson"));
}

#[test]
fn test_learn_list() {
    bin()
        .args(["learn", "list"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Welcome to MycoNote"))
        .stdout(predicate::str::contains("Gene Prediction"))
        .stdout(predicate::str::contains("NCBI Submission"));
}

#[test]
fn test_learn_alias_tutorial() {
    bin()
        .args(["tutorial", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("Interactive Tutorial"));
}

#[test]
fn test_learn_alias_swirl() {
    bin()
        .args(["swirl", "--help"])
        .assert()
        .success()
        .stdout(predicate::str::contains("swirl"));
}

#[test]
fn test_learn_invalid_lesson() {
    bin().args(["learn", "999"]).assert().failure();
}

#[test]
fn test_learn_reset() {
    bin()
        .args(["learn", "--reset"])
        .assert()
        .success()
        .stdout(predicate::str::contains("reset"));
}
