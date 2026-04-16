use super::context::StageContext;
use super::profile::UserProfile;
use super::retrieval::Citation;
use super::rules::Finding;

/// Render the complete system prompt for a stage interpretation.
pub fn render(
    stage: &str,
    ctx: &StageContext,
    findings: &[Finding],
    retrieved: &[Citation],
    profile: &UserProfile,
) -> String {
    let mut prompt = String::with_capacity(8192);

    // System identity
    prompt.push_str(&format!(
        "You are a bioinformatics assistant interpreting the {} stage of a \
         myconote-cli genome-annotation run. Answer concisely, in prose that a lab \
         scientist (not a bioinformatician) can follow.\n\n",
        stage
    ));

    // Citation rules
    prompt.push_str("RULES:\n");
    prompt.push_str("- Every factual claim must cite a source: [knowledge:<id>], [paper:<doi>], [data:<file>:<line>], or [rule:<id>].\n");
    prompt.push_str("- If a claim cannot be cited, omit it.\n");
    prompt.push_str("- Do not invent numbers. Use only values from the context below.\n");
    prompt.push_str("- If a critical file is missing, recommend rerunning the stage.\n");
    prompt.push_str("- End your response with 1-3 recommended next commands from the tool catalog.\n");
    prompt.push_str("- Format commands as: `myconote-cli <command> [args]` with a one-sentence rationale.\n");
    prompt.push_str("- Your interpretation is suggestive, not definitive. Remind the user to verify findings in the context of their specific organism, assembly, and research goals.\n\n");

    // Stage-specific reference ranges
    prompt.push_str("REFERENCE RANGES:\n");
    prompt.push_str(&stage_reference_ranges(stage));
    prompt.push('\n');

    // Deterministic findings
    if !findings.is_empty() {
        prompt.push_str("DETERMINISTIC FINDINGS (from rule engine — incorporate these):\n");
        for f in findings {
            prompt.push_str(&format!(
                "- [{}] {}: {} (evidence: {})",
                f.severity, f.rule_id, f.message, f.evidence
            ));
            if let Some(ref cite) = f.citation {
                prompt.push_str(&format!(" [{}]", cite));
            }
            prompt.push('\n');
        }
        prompt.push('\n');
    }

    // Retrieved knowledge
    if !retrieved.is_empty() {
        prompt.push_str("RETRIEVED KNOWLEDGE:\n");
        for c in retrieved {
            prompt.push_str(&format!("[{}] (score: {:.2})\n{}\n\n", c.id, c.score, c.snippet));
        }
    }

    // User profile
    let profile_summary = profile.summary();
    if profile_summary != "no profile set" {
        prompt.push_str(&format!("USER PROFILE: {}\n\n", profile_summary));
    }

    // Context JSON
    prompt.push_str("CONTEXT (JSON):\n```json\n");
    if let Ok(json) = serde_json::to_string_pretty(ctx) {
        prompt.push_str(&json);
    }
    prompt.push_str("\n```\n");

    prompt
}

/// Render the user message.
pub fn user_message(stage: &str) -> String {
    format!(
        "Interpret the {} stage results. Highlight anything unusual. \
         Recommend 1-3 concrete myconote-cli commands I should run next, \
         with a one-sentence rationale for each.",
        stage
    )
}

fn stage_reference_ranges(stage: &str) -> String {
    match stage {
        "sort" => "- Total contigs: <100 is well-assembled; >1000 suggests fragmentation\n\
                   - Shortest contig: <500 bp may be noise, consider --min-length filter\n\
                   - Assembly naming: NCBI requires clean IDs without special characters\n".to_string(),

        "mask" => "- Fungi: 3-10% masked content is typical\n\
                   - Plants: 50-85% is typical (extensive transposons)\n\
                   - Animals: 30-50% is typical\n\
                   - <1% masked: likely masking failure, will cause over-prediction\n\
                   - >85%: unusual except for large plant genomes\n".to_string(),

        "train" => "- Augustus training needs a close starting species model\n\
                    - Trinity assembly should produce >10,000 transcripts\n\
                    - >20M paired-end reads recommended for training\n".to_string(),

        "predict" => "- Fungal genomes: 5,000-15,000 genes typical\n\
                      - Plants: 20,000-50,000 genes typical\n\
                      - Animals: 15,000-25,000 genes typical\n\
                      - >15,000 genes in fungi usually = over-prediction from unmasked repeats\n\
                      - Augustus vs SNAP disagreement <20% is normal\n".to_string(),

        "update" => "- Gene count may decrease 1-5% (overlapping models merged by PASA)\n\
                     - Average gene length should increase (UTR addition)\n\
                     - >10% gene count increase after update: investigate\n\
                     - >20% decrease: PASA may have over-merged\n".to_string(),

        "annotate" => "- Swiss-Prot: 40-70% of genes should have hits\n\
                       - Pfam: 50-80% of genes should have domains\n\
                       - BUSCO completeness: >90% expected\n\
                       - EggNOG: 60-85% classified\n\
                       - CAZyme: 1-3% for fungi, 1-5% for plants\n".to_string(),

        "submit" => "- Zero errors in errorsummary.val required for submission\n\
                     - Common errors: NoStop (wrong genetic code), InternalStop, PartialProblem\n\
                     - Warnings are acceptable but should be reviewed\n\
                     - .sqn file = ready for NCBI upload\n".to_string(),

        _ => String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chat::rules::Severity;

    #[test]
    fn render_includes_stage() {
        let ctx = StageContext {
            stage: "predict".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![],
            stats: None,
            notes: vec![],
        };
        let prompt = render("predict", &ctx, &[], &[], &UserProfile::default());
        assert!(prompt.contains("predict stage"));
        assert!(prompt.contains("RULES:"));
        assert!(prompt.contains("REFERENCE RANGES:"));
        assert!(prompt.contains("CONTEXT (JSON):"));
    }

    #[test]
    fn render_includes_findings() {
        let ctx = StageContext {
            stage: "predict".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![],
            stats: None,
            notes: vec![],
        };
        let findings = vec![Finding {
            rule_id: "predict.gene_count_high".to_string(),
            severity: Severity::Warning,
            message: "Gene count is high.".to_string(),
            evidence: "total_genes: 25000".to_string(),
            citation: Some("knowledge:predict.gene_count_ranges".to_string()),
        }];
        let prompt = render("predict", &ctx, &findings, &[], &UserProfile::default());
        assert!(prompt.contains("DETERMINISTIC FINDINGS"));
        assert!(prompt.contains("gene_count_high"));
        assert!(prompt.contains("knowledge:predict.gene_count_ranges"));
    }

    #[test]
    fn render_includes_profile() {
        let ctx = StageContext {
            stage: "predict".to_string(),
            dir: "/tmp".to_string(),
            artifacts: vec![],
            stats: None,
            notes: vec![],
        };
        let profile = UserProfile {
            kingdom: Some("fungi".to_string()),
            organism: Some("Aspergillus niger".to_string()),
            typical_genome_size_mb: None,
            preferred_verbosity: None,
        };
        let prompt = render("predict", &ctx, &[], &[], &profile);
        assert!(prompt.contains("USER PROFILE:"));
        assert!(prompt.contains("fungi"));
    }

    #[test]
    fn user_message_mentions_stage() {
        let msg = user_message("annotate");
        assert!(msg.contains("annotate"));
        assert!(msg.contains("commands"));
    }

    #[test]
    fn reference_ranges_exist_for_all_stages() {
        for stage in ["sort", "mask", "train", "predict", "update", "annotate", "submit"] {
            let ranges = stage_reference_ranges(stage);
            assert!(!ranges.is_empty(), "no reference ranges for stage '{}'", stage);
        }
    }
}
