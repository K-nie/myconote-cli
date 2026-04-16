use serde::{Deserialize, Serialize};

use super::context::StageContext;

// ─────────────────────────────────────────────────────────────────────────────
// Finding — a single deterministic result from the rule engine
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Serialize, Clone, Debug)]
pub struct Finding {
    pub rule_id: String,
    pub severity: Severity,
    pub message: String,
    pub evidence: String,
    pub citation: Option<String>,
}

#[derive(Serialize, Clone, Debug, PartialEq)]
pub enum Severity {
    Info,
    Warning,
    Critical,
}

impl std::fmt::Display for Severity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Severity::Info => write!(f, "INFO"),
            Severity::Warning => write!(f, "WARNING"),
            Severity::Critical => write!(f, "CRITICAL"),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Rule definitions loaded from TOML
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Deserialize, Debug)]
pub struct RuleSet {
    pub rules: Vec<RuleDef>,
}

#[derive(Deserialize, Debug)]
pub struct RuleDef {
    pub id: String,
    pub stage: String,
    pub field: String,
    pub condition: String,
    pub threshold: f64,
    pub severity: String,
    pub message: String,
    pub citation: Option<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Built-in rules (compiled in; TOML rules extend these)
// ─────────────────────────────────────────────────────────────────────────────

/// Evaluate all applicable rules against a StageContext. Returns findings.
pub fn evaluate(ctx: &StageContext) -> Vec<Finding> {
    let mut findings = Vec::new();

    match ctx.stage.as_str() {
        "sort" => evaluate_sort(ctx, &mut findings),
        "mask" => evaluate_mask(ctx, &mut findings),
        "train" => evaluate_train(ctx, &mut findings),
        "predict" => evaluate_predict(ctx, &mut findings),
        "update" => evaluate_update(ctx, &mut findings),
        "annotate" => evaluate_annotate(ctx, &mut findings),
        "submit" => evaluate_submit(ctx, &mut findings),
        _ => {}
    }

    // Extend with TOML-defined rules if available
    if let Some(toml_findings) = evaluate_toml_rules(ctx) {
        findings.extend(toml_findings);
    }

    findings
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-stage rule evaluators
// ─────────────────────────────────────────────────────────────────────────────

fn evaluate_sort(ctx: &StageContext, findings: &mut Vec<Finding>) {
    // Check for rename table
    let has_rename = ctx
        .artifacts
        .iter()
        .any(|a| a.name.contains("name_map") || a.name.contains("rename"));
    if !has_rename {
        findings.push(Finding {
            rule_id: "sort.no_rename_table".to_string(),
            severity: Severity::Info,
            message: "No rename table found. If contigs were renamed, the mapping is not recorded."
                .to_string(),
            evidence: "No file matching 'name_map' or 'rename' found in directory.".to_string(),
            citation: None,
        });
    }

    // Check for very small genome
    for art in &ctx.artifacts {
        if (art.name.ends_with(".fa") || art.name.ends_with(".fas") || art.name.ends_with(".fasta"))
            && art.size_bytes < 1_000_000
        {
            findings.push(Finding {
                rule_id: "sort.small_genome".to_string(),
                severity: Severity::Warning,
                message: format!("Genome file '{}' is very small ({} bytes). Verify this is the correct assembly.", art.name, art.size_bytes),
                evidence: format!("{}: {} bytes", art.name, art.size_bytes),
                citation: None,
            });
        }
    }
}

fn evaluate_mask(ctx: &StageContext, findings: &mut Vec<Finding>) {
    // Parse masking percentage from notes
    for note in &ctx.notes {
        if let Some(pct_str) = note.strip_prefix("Soft-masked content: ") {
            if let Ok(pct) = pct_str.trim_end_matches('%').parse::<f64>() {
                if pct < 1.0 {
                    findings.push(Finding {
                        rule_id: "mask.very_low_masking".to_string(),
                        severity: Severity::Critical,
                        message: format!("Only {:.1}% of the genome is masked. This is unusually low and may indicate masking failed. Gene prediction on unmasked genomes produces severe over-prediction.", pct),
                        evidence: format!("Masked content: {:.1}%", pct),
                        citation: Some("knowledge:mask.repeat_content_norms".to_string()),
                    });
                } else if pct < 3.0 {
                    findings.push(Finding {
                        rule_id: "mask.low_masking".to_string(),
                        severity: Severity::Warning,
                        message: format!("Masked content is {:.1}%, which is low even for compact fungal genomes (typical: 3-10%). Verify RepeatModeler/RepeatMasker ran correctly.", pct),
                        evidence: format!("Masked content: {:.1}%", pct),
                        citation: Some("knowledge:mask.repeat_content_norms".to_string()),
                    });
                } else if pct > 85.0 {
                    findings.push(Finding {
                        rule_id: "mask.very_high_masking".to_string(),
                        severity: Severity::Warning,
                        message: format!("Masked content is {:.1}%, which is very high. This is expected for large plant genomes but unusual for fungi or most animals. Verify the organism.", pct),
                        evidence: format!("Masked content: {:.1}%", pct),
                        citation: Some("knowledge:mask.repeat_content_norms".to_string()),
                    });
                } else {
                    findings.push(Finding {
                        rule_id: "mask.normal_masking".to_string(),
                        severity: Severity::Info,
                        message: format!(
                            "Masked content: {:.1}%. This is within typical range.",
                            pct
                        ),
                        evidence: format!("Masked content: {:.1}%", pct),
                        citation: None,
                    });
                }
            }
        }
    }
}

fn evaluate_train(_ctx: &StageContext, findings: &mut Vec<Finding>) {
    // Training stage — mostly informational for now
    findings.push(Finding {
        rule_id: "train.check_species".to_string(),
        severity: Severity::Info,
        message: "Verify that the Augustus species used for training matches your organism or a close relative.".to_string(),
        evidence: "Training stage context".to_string(),
        citation: None,
    });
}

fn evaluate_predict(ctx: &StageContext, findings: &mut Vec<Finding>) {
    if let Some(ref stats_val) = ctx.stats {
        let gene_count = stats_val
            .get("total_genes")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);

        if gene_count == 0 {
            findings.push(Finding {
                rule_id: "predict.no_genes".to_string(),
                severity: Severity::Critical,
                message: "No genes found in the GFF3. Gene prediction may have failed entirely."
                    .to_string(),
                evidence: "total_genes: 0".to_string(),
                citation: None,
            });
        } else if gene_count > 20_000 {
            findings.push(Finding {
                rule_id: "predict.gene_count_high".to_string(),
                severity: Severity::Warning,
                message: format!("Gene count ({}) is high. For fungal genomes, >15,000 genes usually indicates over-prediction from unmasked repeats. For plants, 20,000-50,000 may be normal.", gene_count),
                evidence: format!("total_genes: {}", gene_count),
                citation: Some("knowledge:predict.gene_count_ranges".to_string()),
            });
        } else if gene_count < 3_000 {
            findings.push(Finding {
                rule_id: "predict.gene_count_low".to_string(),
                severity: Severity::Warning,
                message: format!("Gene count ({}) is unusually low. Most eukaryotic genomes have >5,000 genes. Check if the assembly is complete.", gene_count),
                evidence: format!("total_genes: {}", gene_count),
                citation: Some("knowledge:predict.gene_count_ranges".to_string()),
            });
        } else {
            findings.push(Finding {
                rule_id: "predict.gene_count_normal".to_string(),
                severity: Severity::Info,
                message: format!(
                    "Gene count: {}. Within typical range for eukaryotic genomes.",
                    gene_count
                ),
                evidence: format!("total_genes: {}", gene_count),
                citation: None,
            });
        }

        // Check transcript/gene ratio (isoform complexity)
        let transcript_count = stats_val
            .get("total_transcripts")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);
        if gene_count > 0 && transcript_count > 0 {
            let ratio = transcript_count as f64 / gene_count as f64;
            if ratio > 2.0 {
                findings.push(Finding {
                    rule_id: "predict.high_isoform_ratio".to_string(),
                    severity: Severity::Info,
                    message: format!("Transcript/gene ratio is {:.1}. Multiple isoforms per gene — expected if RNA-seq evidence was used.", ratio),
                    evidence: format!("transcripts: {}, genes: {}, ratio: {:.1}", transcript_count, gene_count, ratio),
                    citation: None,
                });
            }
        }

        // Check for summary file
        let has_summary = ctx.artifacts.iter().any(|a| a.name.contains("summary"));
        if has_summary {
            findings.push(Finding {
                rule_id: "predict.summary_present".to_string(),
                severity: Severity::Info,
                message: "Prediction summary file found.".to_string(),
                evidence: "predict_summary.txt present".to_string(),
                citation: None,
            });
        }
    } else if !ctx.artifacts.is_empty() {
        findings.push(Finding {
            rule_id: "predict.no_stats".to_string(),
            severity: Severity::Warning,
            message: "GFF3 files found but GenomeStatistics could not be computed. The GFF3 may be malformed.".to_string(),
            evidence: format!("Artifacts: {}", ctx.artifacts.iter().map(|a| a.name.as_str()).collect::<Vec<_>>().join(", ")),
            citation: None,
        });
    }
}

fn evaluate_update(ctx: &StageContext, findings: &mut Vec<Finding>) {
    if let Some(ref stats_val) = ctx.stats {
        let gene_count = stats_val
            .get("total_genes")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);
        if gene_count > 0 {
            findings.push(Finding {
                rule_id: "update.genes_present".to_string(),
                severity: Severity::Info,
                message: format!("Updated gene models: {} genes.", gene_count),
                evidence: format!("total_genes: {}", gene_count),
                citation: None,
            });
        }
    } else {
        findings.push(Finding {
            rule_id: "update.no_gff3".to_string(),
            severity: Severity::Warning,
            message: "No GFF3 found in update directory. PASA update may not have run.".to_string(),
            evidence: "No GFF3 files".to_string(),
            citation: None,
        });
    }
}

fn evaluate_annotate(ctx: &StageContext, findings: &mut Vec<Finding>) {
    // Check for annotation source coverage
    let annotation_sources = [
        "pfam",
        "eggnog",
        "cazyme",
        "merops",
        "interproscan",
        "busco",
        "mmseqs",
    ];
    let mut found = Vec::new();
    let mut missing = Vec::new();

    for src in &annotation_sources {
        if ctx
            .artifacts
            .iter()
            .any(|a| a.name.to_lowercase().contains(src))
            || ctx.notes.iter().any(|n| n.to_lowercase().contains(src))
        {
            found.push(*src);
        } else {
            missing.push(*src);
        }
    }

    if !found.is_empty() {
        findings.push(Finding {
            rule_id: "annotate.sources_found".to_string(),
            severity: Severity::Info,
            message: format!("Annotation sources detected: {}", found.join(", ")),
            evidence: format!("Found: {}", found.join(", ")),
            citation: None,
        });
    }

    if !missing.is_empty() {
        findings.push(Finding {
            rule_id: "annotate.sources_missing".to_string(),
            severity: Severity::Warning,
            message: format!("Annotation sources not detected: {}. Consider running these for more comprehensive annotation.", missing.join(", ")),
            evidence: format!("Missing: {}", missing.join(", ")),
            citation: Some("knowledge:annotate.functional_coverage".to_string()),
        });
    }

    // Gene count from stats
    if let Some(ref stats_val) = ctx.stats {
        let gene_count = stats_val
            .get("total_genes")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);
        if gene_count > 0 {
            findings.push(Finding {
                rule_id: "annotate.gene_count".to_string(),
                severity: Severity::Info,
                message: format!("Annotated gene models: {} genes.", gene_count),
                evidence: format!("total_genes: {}", gene_count),
                citation: None,
            });
        }
    }
}

fn evaluate_submit(ctx: &StageContext, findings: &mut Vec<Finding>) {
    // Check for errorsummary.val
    let val_artifact = ctx.artifacts.iter().find(|a| a.name == "errorsummary.val");

    match val_artifact {
        Some(art) => {
            if let Some(ref preview) = art.preview {
                let error_count = preview.lines().filter(|l| l.contains("ERROR")).count();
                let warning_count = preview.lines().filter(|l| l.contains("WARNING")).count();

                if error_count > 0 {
                    findings.push(Finding {
                        rule_id: "submit.validation_errors".to_string(),
                        severity: Severity::Critical,
                        message: format!("NCBI validation found {} error(s). These must be resolved before submission.", error_count),
                        evidence: format!("{} ERRORs, {} WARNINGs in errorsummary.val", error_count, warning_count),
                        citation: Some("knowledge:submit.ncbi_validation".to_string()),
                    });
                } else if warning_count > 0 {
                    findings.push(Finding {
                        rule_id: "submit.validation_warnings".to_string(),
                        severity: Severity::Warning,
                        message: format!("NCBI validation found {} warning(s) but no errors. Review warnings before submitting.", warning_count),
                        evidence: format!("{} WARNINGs in errorsummary.val", warning_count),
                        citation: Some("knowledge:submit.ncbi_validation".to_string()),
                    });
                } else {
                    findings.push(Finding {
                        rule_id: "submit.validation_clean".to_string(),
                        severity: Severity::Info,
                        message: "NCBI validation passed with no errors or warnings.".to_string(),
                        evidence: "errorsummary.val: clean".to_string(),
                        citation: None,
                    });
                }
            }
        }
        None => {
            findings.push(Finding {
                rule_id: "submit.no_validation".to_string(),
                severity: Severity::Warning,
                message: "No errorsummary.val found. Run `myconote-cli submit --validate-only` to check for NCBI compliance.".to_string(),
                evidence: "errorsummary.val not found".to_string(),
                citation: None,
            });
        }
    }

    // Check for .sqn file
    let has_sqn = ctx.artifacts.iter().any(|a| a.name.ends_with(".sqn"));
    if has_sqn {
        findings.push(Finding {
            rule_id: "submit.sqn_present".to_string(),
            severity: Severity::Info,
            message: "Sequin file (.sqn) found — ready for NCBI submission.".to_string(),
            evidence: ".sqn file present".to_string(),
            citation: None,
        });
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TOML rule loader (extends built-in rules)
// ─────────────────────────────────────────────────────────────────────────────

fn evaluate_toml_rules(ctx: &StageContext) -> Option<Vec<Finding>> {
    let rules_dir = find_assets_rules_dir()?;
    let stage_file = rules_dir.join(format!("{}.toml", ctx.stage));

    if !stage_file.exists() {
        return None;
    }

    let text = std::fs::read_to_string(&stage_file).ok()?;
    let rule_set: RuleSet = toml::from_str(&text).ok()?;

    let mut findings = Vec::new();

    for rule in &rule_set.rules {
        if rule.stage != ctx.stage {
            continue;
        }

        // Extract the field value from stats
        let value = extract_field_value(ctx, &rule.field);
        if let Some(val) = value {
            let triggered = match rule.condition.as_str() {
                "gt" => val > rule.threshold,
                "lt" => val < rule.threshold,
                "gte" => val >= rule.threshold,
                "lte" => val <= rule.threshold,
                "eq" => (val - rule.threshold).abs() < f64::EPSILON,
                _ => false,
            };

            if triggered {
                let sev = match rule.severity.as_str() {
                    "critical" => Severity::Critical,
                    "warning" => Severity::Warning,
                    _ => Severity::Info,
                };
                findings.push(Finding {
                    rule_id: rule.id.clone(),
                    severity: sev,
                    message: rule.message.clone(),
                    evidence: format!("{} = {}", rule.field, val),
                    citation: rule.citation.clone(),
                });
            }
        }
    }

    Some(findings)
}

fn extract_field_value(ctx: &StageContext, field: &str) -> Option<f64> {
    ctx.stats.as_ref().and_then(|v| {
        v.get(field)
            .and_then(|f| f.as_f64().or_else(|| f.as_u64().map(|u| u as f64)))
    })
}

fn find_assets_rules_dir() -> Option<std::path::PathBuf> {
    // Try relative to the binary, then CARGO_MANIFEST_DIR, then cwd
    let candidates = [
        std::env::current_exe()
            .ok()
            .and_then(|p| p.parent().map(|d| d.join("assets/chat/rules"))),
        std::env::var("CARGO_MANIFEST_DIR")
            .ok()
            .map(|d| std::path::PathBuf::from(d).join("assets/chat/rules")),
        Some(std::path::PathBuf::from("assets/chat/rules")),
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
    use crate::chat::context::ArtifactSummary;

    fn make_predict_context(total_genes: u64, total_transcripts: u64) -> StageContext {
        let stats = serde_json::json!({
            "total_genes": total_genes,
            "total_transcripts": total_transcripts,
            "total_cds": total_genes,
            "total_exons": total_genes * 3,
            "gene_lengths": [],
            "chromosome_stats": {},
            "isoforms_per_gene": {},
            "primary_only_genes": null,
            "total_features": total_genes + total_transcripts,
        });

        StageContext {
            stage: "predict".to_string(),
            dir: "/tmp/test".to_string(),
            artifacts: vec![ArtifactSummary {
                name: "consensus.gff3".to_string(),
                size_bytes: 1000,
                line_count: Some(100),
                preview: None,
            }],
            stats: Some(stats),
            notes: vec![],
        }
    }

    #[test]
    fn predict_gene_count_high_fires() {
        let ctx = make_predict_context(25_000, 25_000);
        let findings = evaluate(&ctx);
        assert!(findings
            .iter()
            .any(|f| f.rule_id == "predict.gene_count_high"));
        let high = findings
            .iter()
            .find(|f| f.rule_id == "predict.gene_count_high")
            .unwrap();
        assert_eq!(high.severity, Severity::Warning);
    }

    #[test]
    fn predict_gene_count_low_fires() {
        let ctx = make_predict_context(2_000, 2_000);
        let findings = evaluate(&ctx);
        assert!(findings
            .iter()
            .any(|f| f.rule_id == "predict.gene_count_low"));
    }

    #[test]
    fn predict_normal_no_warnings() {
        let ctx = make_predict_context(8_000, 8_000);
        let findings = evaluate(&ctx);
        assert!(!findings
            .iter()
            .any(|f| f.severity == Severity::Warning || f.severity == Severity::Critical));
        assert!(findings
            .iter()
            .any(|f| f.rule_id == "predict.gene_count_normal"));
    }

    #[test]
    fn predict_no_genes_critical() {
        let ctx = make_predict_context(0, 0);
        let findings = evaluate(&ctx);
        assert!(findings
            .iter()
            .any(|f| f.rule_id == "predict.no_genes" && f.severity == Severity::Critical));
    }

    #[test]
    fn mask_low_masking_fires() {
        let ctx = StageContext {
            stage: "mask".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![],
            stats: None,
            notes: vec!["Soft-masked content: 0.5%".to_string()],
        };
        let findings = evaluate(&ctx);
        assert!(findings
            .iter()
            .any(|f| f.rule_id == "mask.very_low_masking" && f.severity == Severity::Critical));
    }

    #[test]
    fn submit_validation_errors_fires() {
        let ctx = StageContext {
            stage: "submit".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![ArtifactSummary {
                name: "errorsummary.val".to_string(),
                size_bytes: 200,
                line_count: Some(5),
                preview: Some("ERROR: SEQ_FEAT.NoStop\nERROR: SEQ_FEAT.PartialProblem\nWARNING: SEQ_FEAT.NotSpliceConsensus".to_string()),
            }],
            stats: None,
            notes: vec![],
        };
        let findings = evaluate(&ctx);
        let err = findings
            .iter()
            .find(|f| f.rule_id == "submit.validation_errors")
            .unwrap();
        assert_eq!(err.severity, Severity::Critical);
        assert!(err.message.contains("2 error(s)"));
    }

    #[test]
    fn submit_no_validation_file() {
        let ctx = StageContext {
            stage: "submit".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![],
            stats: None,
            notes: vec![],
        };
        let findings = evaluate(&ctx);
        assert!(findings.iter().any(|f| f.rule_id == "submit.no_validation"));
    }

    #[test]
    fn annotate_sources_detection() {
        let ctx = StageContext {
            stage: "annotate".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![
                ArtifactSummary {
                    name: "pfam_results.tsv".to_string(),
                    size_bytes: 1000,
                    line_count: Some(50),
                    preview: None,
                },
                ArtifactSummary {
                    name: "busco_summary.txt".to_string(),
                    size_bytes: 500,
                    line_count: Some(10),
                    preview: None,
                },
            ],
            stats: None,
            notes: vec![],
        };
        let findings = evaluate(&ctx);
        assert!(findings
            .iter()
            .any(|f| f.rule_id == "annotate.sources_found"));
        assert!(findings
            .iter()
            .any(|f| f.rule_id == "annotate.sources_missing"));
    }

    #[test]
    fn severity_display() {
        assert_eq!(format!("{}", Severity::Info), "INFO");
        assert_eq!(format!("{}", Severity::Warning), "WARNING");
        assert_eq!(format!("{}", Severity::Critical), "CRITICAL");
    }
}
