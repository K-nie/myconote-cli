use super::retrieval::Citation;

/// Validate LLM output: strip unresolvable citations.
///
/// Every `[knowledge:...]`, `[paper:...]`, `[data:...]`, or `[rule:...]` tag
/// in the response must match an entry in the retrieval set or findings list.
/// Unresolvable tags are removed from the text.
pub fn validate_citations(response: &str, valid_ids: &[String]) -> ValidationResult {
    let re = regex::Regex::new(r"\[(knowledge|paper|data|rule):([^\]]+)\]").unwrap();

    let mut cleaned = response.to_string();
    let mut removed_count = 0;
    let mut kept_count = 0;

    // Collect all matches first, then process in reverse to preserve positions
    let matches: Vec<_> = re.find_iter(response).collect();

    for m in matches.iter().rev() {
        let tag = m.as_str();
        // Extract the full ID (e.g., "knowledge:predict.gene_count_ranges")
        let id = &tag[1..tag.len() - 1]; // strip [ and ]

        if valid_ids.iter().any(|vid| vid == id || tag.contains(vid.as_str())) {
            kept_count += 1;
        } else {
            // Remove the unresolvable citation
            cleaned.replace_range(m.range(), "");
            removed_count += 1;
        }
    }

    ValidationResult {
        text: cleaned.trim().to_string(),
        citations_kept: kept_count,
        citations_removed: removed_count,
        is_valid: removed_count == 0,
    }
}

/// Build the list of valid citation IDs from findings and retrieved knowledge.
pub fn build_valid_ids(
    findings: &[super::rules::Finding],
    retrieved: &[Citation],
) -> Vec<String> {
    let mut ids = Vec::new();

    for f in findings {
        ids.push(format!("rule:{}", f.rule_id));
        if let Some(ref cite) = f.citation {
            ids.push(cite.clone());
        }
    }

    for c in retrieved {
        ids.push(c.id.clone());
    }

    ids
}

pub struct ValidationResult {
    pub text: String,
    pub citations_kept: usize,
    pub citations_removed: usize,
    pub is_valid: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_citations_kept() {
        let response = "Gene count is high [knowledge:predict.gene_count_ranges]. This is typical.";
        let valid = vec!["knowledge:predict.gene_count_ranges".to_string()];
        let result = validate_citations(response, &valid);
        assert!(result.is_valid);
        assert_eq!(result.citations_kept, 1);
        assert_eq!(result.citations_removed, 0);
        assert!(result.text.contains("[knowledge:predict.gene_count_ranges]"));
    }

    #[test]
    fn invalid_citations_removed() {
        let response = "Gene count is 25000 [knowledge:predict.gene_count_ranges]. \
                        Also see [paper:10.1234/fake] for details.";
        let valid = vec!["knowledge:predict.gene_count_ranges".to_string()];
        let result = validate_citations(response, &valid);
        assert!(!result.is_valid);
        assert_eq!(result.citations_kept, 1);
        assert_eq!(result.citations_removed, 1);
        assert!(result.text.contains("[knowledge:predict.gene_count_ranges]"));
        assert!(!result.text.contains("[paper:10.1234/fake]"));
    }

    #[test]
    fn no_citations_is_valid() {
        let response = "Everything looks normal. Run stats next.";
        let valid = vec![];
        let result = validate_citations(response, &valid);
        assert!(result.is_valid);
        assert_eq!(result.citations_kept, 0);
        assert_eq!(result.citations_removed, 0);
    }

    #[test]
    fn rule_citations_valid() {
        let response = "Low masking detected [rule:mask.very_low_masking].";
        let valid = vec!["rule:mask.very_low_masking".to_string()];
        let result = validate_citations(response, &valid);
        assert!(result.is_valid);
        assert_eq!(result.citations_kept, 1);
    }

    #[test]
    fn build_valid_ids_from_findings() {
        use crate::chat::rules::{Finding, Severity};
        use crate::chat::retrieval::{Citation, SourceType};

        let findings = vec![Finding {
            rule_id: "predict.high".to_string(),
            severity: Severity::Warning,
            message: "high".to_string(),
            evidence: "25000".to_string(),
            citation: Some("knowledge:predict.gene_count_ranges".to_string()),
        }];

        let retrieved = vec![Citation {
            source_type: SourceType::Knowledge,
            id: "knowledge:predict.overprediction_causes".to_string(),
            snippet: "text".to_string(),
            score: 1.5,
        }];

        let ids = build_valid_ids(&findings, &retrieved);
        assert!(ids.contains(&"rule:predict.high".to_string()));
        assert!(ids.contains(&"knowledge:predict.gene_count_ranges".to_string()));
        assert!(ids.contains(&"knowledge:predict.overprediction_causes".to_string()));
    }
}
