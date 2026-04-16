use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use super::context::StageContext;
use super::rules::Finding;

// ─────────────────────────────────────────────────────────────────────────────
// Command recommendation
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Serialize, Debug, Clone)]
pub struct CommandRecommendation {
    pub command: String,
    pub rationale: String,
    pub finding_id: Option<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Tool catalog (TOML)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Deserialize, Debug)]
struct ToolCatalog {
    commands: Vec<ToolEntry>,
}

#[derive(Deserialize, Debug, Clone)]
struct ToolEntry {
    name: String,
    #[allow(dead_code)]
    usage: String,
    #[allow(dead_code)]
    when: String,
    #[allow(dead_code)]
    flags: Vec<String>,
    #[allow(dead_code)]
    inputs: Option<Vec<String>>,
    #[allow(dead_code)]
    outputs: Option<Vec<String>>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Generate command recommendations based on the stage, context, and findings.
pub fn recommend(ctx: &StageContext, findings: &[Finding]) -> Vec<CommandRecommendation> {
    let mut recs = Vec::new();

    // Stage-specific recommendations
    match ctx.stage.as_str() {
        "sort" => recommend_after_sort(ctx, findings, &mut recs),
        "mask" => recommend_after_mask(ctx, findings, &mut recs),
        "train" => recommend_after_train(ctx, findings, &mut recs),
        "predict" => recommend_after_predict(ctx, findings, &mut recs),
        "update" => recommend_after_update(ctx, findings, &mut recs),
        "annotate" => recommend_after_annotate(ctx, findings, &mut recs),
        "submit" => recommend_after_submit(ctx, findings, &mut recs),
        _ => {}
    }

    // Validate against tool catalog
    let catalog = load_tool_catalog();
    validate_recommendations(&mut recs, &catalog);

    // Cap at 3
    recs.truncate(3);
    recs
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-stage recommenders
// ─────────────────────────────────────────────────────────────────────────────

fn recommend_after_sort(
    _ctx: &StageContext,
    _findings: &[Finding],
    recs: &mut Vec<CommandRecommendation>,
) {
    recs.push(CommandRecommendation {
        command: "myconote-cli mask <sorted_genome.fa>".to_string(),
        rationale: "Next pipeline step: mask repeats before gene prediction.".to_string(),
        finding_id: None,
    });
}

fn recommend_after_mask(
    _ctx: &StageContext,
    findings: &[Finding],
    recs: &mut Vec<CommandRecommendation>,
) {
    let has_low_masking = findings
        .iter()
        .any(|f| f.rule_id.contains("low_masking") || f.rule_id.contains("very_low_masking"));

    if has_low_masking {
        recs.push(CommandRecommendation {
            command: "myconote-cli mask <genome.fa> --engine full".to_string(),
            rationale: "Low masking detected. Re-run with full engine (RepeatModeler + RepeatMasker) for better results.".to_string(),
            finding_id: Some("mask.low_masking".to_string()),
        });
    }

    recs.push(CommandRecommendation {
        command: "myconote-cli predict <masked_genome.fa> --kingdom fungi".to_string(),
        rationale: "Next pipeline step: predict genes on the masked genome.".to_string(),
        finding_id: None,
    });
}

fn recommend_after_train(
    _ctx: &StageContext,
    _findings: &[Finding],
    recs: &mut Vec<CommandRecommendation>,
) {
    recs.push(CommandRecommendation {
        command: "myconote-cli predict <masked_genome.fa> --species <trained_species>".to_string(),
        rationale: "Use your trained Augustus species model for gene prediction.".to_string(),
        finding_id: None,
    });
}

fn recommend_after_predict(
    _ctx: &StageContext,
    findings: &[Finding],
    recs: &mut Vec<CommandRecommendation>,
) {
    // Always recommend stats
    recs.push(CommandRecommendation {
        command: "myconote-cli stats <consensus.gff3> --format human --taxon fungi".to_string(),
        rationale: "Review gene statistics and compare against kingdom norms.".to_string(),
        finding_id: None,
    });

    let has_high_count = findings
        .iter()
        .any(|f| f.rule_id == "predict.gene_count_high");
    if has_high_count {
        recs.push(CommandRecommendation {
            command: "myconote-cli explain mask --dir <mask_dir>".to_string(),
            rationale: "High gene count may indicate masking failure. Check the mask stage."
                .to_string(),
            finding_id: Some("predict.gene_count_high".to_string()),
        });
    } else {
        recs.push(CommandRecommendation {
            command: "myconote-cli update <predicted.gff3> --fasta <genome.fa>".to_string(),
            rationale: "Next pipeline step: refine gene models with PASA UTR extension."
                .to_string(),
            finding_id: None,
        });
    }
}

fn recommend_after_update(
    _ctx: &StageContext,
    _findings: &[Finding],
    recs: &mut Vec<CommandRecommendation>,
) {
    recs.push(CommandRecommendation {
        command: "myconote-cli annotate <updated.gff3> --fasta <genome.fa>".to_string(),
        rationale: "Next pipeline step: functional annotation with all available databases."
            .to_string(),
        finding_id: None,
    });

    recs.push(CommandRecommendation {
        command: "myconote-cli stats <updated.gff3> --format json".to_string(),
        rationale: "Compare gene counts before and after PASA update.".to_string(),
        finding_id: None,
    });
}

fn recommend_after_annotate(
    _ctx: &StageContext,
    findings: &[Finding],
    recs: &mut Vec<CommandRecommendation>,
) {
    let missing_sources = findings
        .iter()
        .find(|f| f.rule_id == "annotate.sources_missing");
    if let Some(f) = missing_sources {
        recs.push(CommandRecommendation {
            command: "myconote-cli annotate <gff3> --fasta <genome.fa>".to_string(),
            rationale: format!(
                "Some annotation sources are missing. Re-run to add: {}",
                f.evidence
            ),
            finding_id: Some("annotate.sources_missing".to_string()),
        });
    }

    recs.push(CommandRecommendation {
        command: "myconote-cli submit <annotated.gff3> --fasta <genome.fa> --organism '<name>'"
            .to_string(),
        rationale: "Next pipeline step: prepare NCBI GenBank submission.".to_string(),
        finding_id: None,
    });
}

fn recommend_after_submit(
    _ctx: &StageContext,
    findings: &[Finding],
    recs: &mut Vec<CommandRecommendation>,
) {
    let has_errors = findings
        .iter()
        .any(|f| f.rule_id == "submit.validation_errors");

    if has_errors {
        recs.push(CommandRecommendation {
            command: "myconote-cli fix <annotated.gff3> --fasta <genome.fa>".to_string(),
            rationale: "NCBI validation errors found. Run fix to auto-repair common issues."
                .to_string(),
            finding_id: Some("submit.validation_errors".to_string()),
        });
        recs.push(CommandRecommendation {
            command: "myconote-cli submit <fixed.gff3> --fasta <genome.fa> --validate-only"
                .to_string(),
            rationale: "Re-validate after fixing to confirm errors are resolved.".to_string(),
            finding_id: Some("submit.validation_errors".to_string()),
        });
    } else {
        recs.push(CommandRecommendation {
            command: "myconote-cli submit <annotated.gff3> --fasta <genome.fa> --organism '<name>' --bioproject <acc>".to_string(),
            rationale: "Validation passed. Proceed with full submission including BioProject accession.".to_string(),
            finding_id: None,
        });
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Validation
// ─────────────────────────────────────────────────────────────────────────────

fn validate_recommendations(recs: &mut Vec<CommandRecommendation>, catalog: &[ToolEntry]) {
    for rec in recs.iter_mut() {
        // Extract the subcommand from the recommendation
        let parts: Vec<&str> = rec.command.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        let subcmd = parts[1]; // e.g., "mask", "predict", etc.

        // Check subcommand exists in catalog
        let entry = catalog.iter().find(|e| e.name == subcmd);
        if entry.is_none() {
            // Replace with a safe fallback — but our built-in recs should always be valid
            rec.rationale.push_str(&format!(
                " (Run `myconote-cli {} --help` for options.)",
                subcmd
            ));
        }
    }
}

fn load_tool_catalog() -> Vec<ToolEntry> {
    let path = find_tool_catalog();
    match path {
        Some(p) => {
            let text = std::fs::read_to_string(&p).unwrap_or_default();
            let catalog: ToolCatalog = toml::from_str(&text).unwrap_or(ToolCatalog {
                commands: Vec::new(),
            });
            catalog.commands
        }
        None => Vec::new(),
    }
}

fn find_tool_catalog() -> Option<PathBuf> {
    let candidates = [
        std::env::current_exe()
            .ok()
            .and_then(|p| p.parent().map(|d| d.join("assets/chat/tool_catalog.toml"))),
        std::env::var("CARGO_MANIFEST_DIR")
            .ok()
            .map(|d| PathBuf::from(d).join("assets/chat/tool_catalog.toml")),
        Some(PathBuf::from("assets/chat/tool_catalog.toml")),
    ];

    for candidate in &candidates {
        if let Some(ref path) = candidate {
            if path.exists() {
                return Some(path.clone());
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chat::rules::Severity;

    #[test]
    fn predict_recommendations() {
        let ctx = StageContext {
            stage: "predict".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![],
            stats: Some(serde_json::json!({"total_genes": 8000})),
            notes: vec![],
        };
        let findings = vec![Finding {
            rule_id: "predict.gene_count_normal".to_string(),
            severity: Severity::Info,
            message: "Normal count.".to_string(),
            evidence: "8000".to_string(),
            citation: None,
        }];
        let recs = recommend(&ctx, &findings);
        assert!(!recs.is_empty());
        assert!(recs.iter().any(|r| r.command.contains("stats")));
    }

    #[test]
    fn predict_high_count_recommends_mask_check() {
        let ctx = StageContext {
            stage: "predict".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![],
            stats: Some(serde_json::json!({"total_genes": 25000})),
            notes: vec![],
        };
        let findings = vec![Finding {
            rule_id: "predict.gene_count_high".to_string(),
            severity: Severity::Warning,
            message: "High.".to_string(),
            evidence: "25000".to_string(),
            citation: None,
        }];
        let recs = recommend(&ctx, &findings);
        assert!(recs.iter().any(|r| r.command.contains("explain mask")));
    }

    #[test]
    fn submit_errors_recommends_fix() {
        let ctx = StageContext {
            stage: "submit".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![],
            stats: None,
            notes: vec![],
        };
        let findings = vec![Finding {
            rule_id: "submit.validation_errors".to_string(),
            severity: Severity::Critical,
            message: "Errors found.".to_string(),
            evidence: "2 ERRORs".to_string(),
            citation: None,
        }];
        let recs = recommend(&ctx, &findings);
        assert!(recs.iter().any(|r| r.command.contains("fix")));
    }

    #[test]
    fn recommendations_capped_at_three() {
        let ctx = StageContext {
            stage: "annotate".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![],
            stats: None,
            notes: vec![],
        };
        let findings = vec![Finding {
            rule_id: "annotate.sources_missing".to_string(),
            severity: Severity::Warning,
            message: "Missing sources.".to_string(),
            evidence: "Missing: eggnog, cazyme".to_string(),
            citation: None,
        }];
        let recs = recommend(&ctx, &findings);
        assert!(recs.len() <= 3);
    }
}
