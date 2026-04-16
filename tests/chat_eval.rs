/// Eval harness for `myconote explain` — verifies rule recall, ethics precision,
/// format detection, and command validity against hand-crafted fixtures.
///
/// Each fixture is a JSON file describing an input scenario and expected outputs.
/// The harness deserializes the fixture, runs the relevant component (rules, ethics,
/// paste detection), and checks the output against expectations.
///
/// Metrics tracked:
/// - Rule recall: did expected findings fire?
/// - Rule precision: did unexpected findings NOT fire?
/// - Ethics precision: correct refuse vs. correct pass
/// - Format detection accuracy: did paste format match?
/// - Command validity: are recommended commands within expected bounds?
use myconote_cli::chat::{
    commands,
    context::{ArtifactSummary, StageContext},
    ethics,
    paste,
    rules,
};
use serde::Deserialize;
use std::path::PathBuf;

// ─────────────────────────────────────────────────────────────────────────────
// Fixture types
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Deserialize)]
#[allow(dead_code)]
struct RuleFixture {
    description: String,
    stage: String,
    context: FixtureContext,
    expected_findings: ExpectedFindings,
    expected_commands: Option<ExpectedCommands>,
    ethics: Option<String>,
}

#[derive(Deserialize)]
struct FixtureContext {
    stage: String,
    dir: String,
    artifacts: Vec<FixtureArtifact>,
    stats: Option<serde_json::Value>,
    notes: Vec<String>,
}

#[derive(Deserialize)]
struct FixtureArtifact {
    name: String,
    size_bytes: u64,
    line_count: Option<usize>,
    preview: Option<String>,
}

#[derive(Deserialize)]
struct ExpectedFindings {
    must_fire: Vec<String>,
    must_not_fire: Vec<String>,
}

#[derive(Deserialize)]
struct ExpectedCommands {
    min_count: usize,
    max_count: usize,
    must_contain_substring: Vec<String>,
}

#[derive(Deserialize)]
struct EthicsFixture {
    description: String,
    input: String,
    expected_verdict: String,
    #[serde(default)]
    expected_rule_id_prefix: Option<String>,
}

#[derive(Deserialize)]
struct PasteFixture {
    description: String,
    input: String,
    expected_format: String,
    #[allow(dead_code)]
    expected_stage: String,
}

// ─────────────────────────────────────────────────────────────────────────────
// Eval results
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
struct EvalScore {
    fixture: String,
    rule_recall_ok: usize,
    rule_recall_total: usize,
    rule_precision_ok: usize,
    rule_precision_total: usize,
    command_count_ok: bool,
    command_substring_ok: usize,
    command_substring_total: usize,
    passed: bool,
}

impl EvalScore {
    fn rule_recall_pct(&self) -> f64 {
        if self.rule_recall_total == 0 { 100.0 }
        else { self.rule_recall_ok as f64 / self.rule_recall_total as f64 * 100.0 }
    }

    fn rule_precision_pct(&self) -> f64 {
        if self.rule_precision_total == 0 { 100.0 }
        else { self.rule_precision_ok as f64 / self.rule_precision_total as f64 * 100.0 }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/chat/fixtures")
}

fn fixture_to_stage_context(fc: &FixtureContext) -> StageContext {
    StageContext {
        stage: fc.stage.clone(),
        dir: fc.dir.clone(),
        artifacts: fc.artifacts.iter().map(|a| ArtifactSummary {
            name: a.name.clone(),
            size_bytes: a.size_bytes,
            line_count: a.line_count,
            preview: a.preview.clone(),
        }).collect(),
        stats: fc.stats.clone(),
        notes: fc.notes.clone(),
    }
}

fn eval_rule_fixture(fixture_path: &str) -> EvalScore {
    let path = fixtures_dir().join(fixture_path);
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("cannot read fixture: {}", path.display()));
    let fixture: RuleFixture = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("cannot parse fixture {}: {}", path.display(), e));

    let ctx = fixture_to_stage_context(&fixture.context);
    let findings = rules::evaluate(&ctx);
    let finding_ids: Vec<&str> = findings.iter().map(|f| f.rule_id.as_str()).collect();

    let mut score = EvalScore {
        fixture: fixture.description.clone(),
        ..Default::default()
    };

    // Check recall: expected rules must fire
    score.rule_recall_total = fixture.expected_findings.must_fire.len();
    for expected_id in &fixture.expected_findings.must_fire {
        if finding_ids.iter().any(|id| id == expected_id) {
            score.rule_recall_ok += 1;
        } else {
            eprintln!(
                "  [MISS] {} — expected rule '{}' did not fire. Fired: {:?}",
                fixture.description, expected_id, finding_ids
            );
        }
    }

    // Check precision: unexpected rules must NOT fire
    score.rule_precision_total = fixture.expected_findings.must_not_fire.len();
    for unexpected_id in &fixture.expected_findings.must_not_fire {
        if !finding_ids.iter().any(|id| id == unexpected_id) {
            score.rule_precision_ok += 1;
        } else {
            eprintln!(
                "  [SPURIOUS] {} — rule '{}' fired but should not have",
                fixture.description, unexpected_id
            );
        }
    }

    // Check commands
    if let Some(ref ec) = fixture.expected_commands {
        let recs = commands::recommend(&ctx, &findings);
        score.command_count_ok = recs.len() >= ec.min_count && recs.len() <= ec.max_count;

        score.command_substring_total = ec.must_contain_substring.len();
        for sub in &ec.must_contain_substring {
            if recs.iter().any(|r| r.command.contains(sub) || r.rationale.contains(sub)) {
                score.command_substring_ok += 1;
            } else {
                eprintln!(
                    "  [CMD MISS] {} — no command/rationale contains '{}'",
                    fixture.description, sub
                );
            }
        }
    } else {
        score.command_count_ok = true;
    }

    score.passed = score.rule_recall_ok == score.rule_recall_total
        && score.rule_precision_ok == score.rule_precision_total
        && score.command_count_ok;

    score
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn eval_predict_fungal_normal() {
    let score = eval_rule_fixture("predict_fungal_normal.json");
    assert!(score.passed, "predict_fungal_normal failed: {:?}", score);
    assert_eq!(score.rule_recall_pct(), 100.0);
    assert_eq!(score.rule_precision_pct(), 100.0);
}

#[test]
fn eval_predict_fungal_overpredicted() {
    let score = eval_rule_fixture("predict_fungal_overpredicted.json");
    assert!(score.passed, "predict_fungal_overpredicted failed: {:?}", score);
    assert_eq!(score.rule_recall_pct(), 100.0);
}

#[test]
fn eval_annotate_low_coverage() {
    let score = eval_rule_fixture("annotate_low_coverage.json");
    assert!(score.passed, "annotate_low_coverage failed: {:?}", score);
}

#[test]
fn eval_submit_validation_errors() {
    let score = eval_rule_fixture("submit_validation_errors.json");
    assert!(score.passed, "submit_validation_errors failed: {:?}", score);
}

#[test]
fn eval_ethics_refuse_biosecurity() {
    let path = fixtures_dir().join("ethics_refuse_biosecurity.json");
    let content = std::fs::read_to_string(&path).unwrap();
    let fixture: EthicsFixture = serde_json::from_str(&content).unwrap();

    let verdict = ethics::classify(&fixture.input);
    match verdict {
        ethics::EthicsVerdict::Refuse { rule_id, message: _ } => {
            assert_eq!(fixture.expected_verdict, "refuse");
            if let Some(prefix) = &fixture.expected_rule_id_prefix {
                assert!(
                    rule_id.contains(prefix),
                    "rule_id '{}' does not contain prefix '{}'", rule_id, prefix
                );
            }
        }
        ethics::EthicsVerdict::Pass => {
            panic!("Expected refuse but got pass for: {}", fixture.description);
        }
    }
}

#[test]
fn eval_ethics_pass_amr_genes() {
    let path = fixtures_dir().join("ethics_pass_amr_genes.json");
    let content = std::fs::read_to_string(&path).unwrap();
    let fixture: EthicsFixture = serde_json::from_str(&content).unwrap();

    let verdict = ethics::classify(&fixture.input);
    match verdict {
        ethics::EthicsVerdict::Pass => {
            assert_eq!(fixture.expected_verdict, "pass");
        }
        ethics::EthicsVerdict::Refuse { rule_id, message } => {
            panic!(
                "Expected pass but got refuse for: {} (rule: {}, msg: {})",
                fixture.description, rule_id, message
            );
        }
    }
}

#[test]
fn eval_paste_gff_line() {
    let path = fixtures_dir().join("paste_gff_line.json");
    let content = std::fs::read_to_string(&path).unwrap();
    let fixture: PasteFixture = serde_json::from_str(&content).unwrap();

    let detected = paste::detect_format(&fixture.input);
    let format_name = format!("{:?}", detected);
    assert_eq!(
        format_name, fixture.expected_format,
        "Format mismatch for: {}", fixture.description
    );
}

#[test]
fn eval_paste_ncbi_error() {
    let path = fixtures_dir().join("paste_ncbi_error.json");
    let content = std::fs::read_to_string(&path).unwrap();
    let fixture: PasteFixture = serde_json::from_str(&content).unwrap();

    let detected = paste::detect_format(&fixture.input);
    let format_name = format!("{:?}", detected);
    assert_eq!(
        format_name, fixture.expected_format,
        "Format mismatch for: {}", fixture.description
    );
}

/// Run all rule fixtures and produce a summary scoreboard.
#[test]
fn eval_scoreboard_summary() {
    let rule_fixtures = [
        "predict_fungal_normal.json",
        "predict_fungal_overpredicted.json",
        "annotate_low_coverage.json",
        "submit_validation_errors.json",
    ];

    let mut scores: Vec<EvalScore> = Vec::new();
    let mut all_passed = true;

    for fixture_name in &rule_fixtures {
        let score = eval_rule_fixture(fixture_name);
        if !score.passed {
            all_passed = false;
        }
        scores.push(score);
    }

    // Print scoreboard
    println!("\n━━━ Eval Scoreboard ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("{:<40} {:>8} {:>8} {:>6}", "Fixture", "Recall%", "Prec%", "Pass");
    println!("{}", "─".repeat(66));
    for s in &scores {
        println!(
            "{:<40} {:>7.0}% {:>7.0}% {:>6}",
            s.fixture.chars().take(40).collect::<String>(),
            s.rule_recall_pct(),
            s.rule_precision_pct(),
            if s.passed { "yes" } else { "NO" },
        );
    }
    println!("{}", "─".repeat(66));

    assert!(all_passed, "One or more eval fixtures failed");
}
