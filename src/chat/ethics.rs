use serde::Deserialize;
use std::path::PathBuf;

// ─────────────────────────────────────────────────────────────────────────────
// Ethics verdict
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub enum EthicsVerdict {
    Pass,
    Refuse { rule_id: String, message: String },
}

// ─────────────────────────────────────────────────────────────────────────────
// Rule definitions
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Deserialize, Debug)]
struct EthicsRuleSet {
    rules: Vec<EthicsRuleDef>,
}

#[derive(Deserialize, Debug)]
struct EthicsRuleDef {
    id: String,
    #[allow(dead_code)]
    category: String,
    patterns: Vec<String>,
    message: String,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Classify user input against ethics rules. Returns Pass or Refuse.
///
/// This runs BEFORE any LLM call. It's a simple regex/pattern match,
/// not a semantic classifier. The scope is deliberately narrow to avoid
/// over-refusal — legitimate biology questions are never refused.
pub fn classify(input: &str) -> EthicsVerdict {
    let rules = load_ethics_rules();
    let input_lower = input.to_lowercase();

    for rule in &rules {
        for pattern in &rule.patterns {
            let re = match regex::Regex::new(&format!("(?i){}", pattern)) {
                Ok(r) => r,
                Err(_) => continue,
            };
            if re.is_match(&input_lower) {
                return EthicsVerdict::Refuse {
                    rule_id: rule.id.clone(),
                    message: rule.message.clone(),
                };
            }
        }
    }

    EthicsVerdict::Pass
}

fn load_ethics_rules() -> Vec<EthicsRuleDef> {
    let path = find_ethics_rules_file();
    let Some(path) = path else {
        return default_rules();
    };

    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(_) => return default_rules(),
    };

    match toml::from_str::<EthicsRuleSet>(&text) {
        Ok(rs) => rs.rules,
        Err(_) => default_rules(),
    }
}

/// Fallback rules compiled into the binary.
fn default_rules() -> Vec<EthicsRuleDef> {
    vec![
        EthicsRuleDef {
            id: "ethics.biosecurity_design".to_string(),
            category: "biosecurity".to_string(),
            patterns: vec![
                r"design.*pathogen".to_string(),
                r"enhance.*virulence".to_string(),
                r"gain.of.function.*create".to_string(),
                r"engineer.*select.agent".to_string(),
                r"create.*bioweapon".to_string(),
            ],
            message: "This request involves biosecurity-sensitive design. I can't help with that, but I can help you interpret your annotation results.".to_string(),
        },
        EthicsRuleDef {
            id: "ethics.fabrication".to_string(),
            category: "fabrication".to_string(),
            patterns: vec![
                r"fabricat.*results".to_string(),
                r"fake.*data".to_string(),
                r"invent.*statistics".to_string(),
                r"make.up.*numbers".to_string(),
            ],
            message: "This request involves data fabrication. I can't help with that, but I can help you interpret your real annotation results.".to_string(),
        },
        EthicsRuleDef {
            id: "ethics.embargo_violation".to_string(),
            category: "embargo".to_string(),
            patterns: vec![
                r"circumvent.*embargo".to_string(),
                r"bypass.*ncbi.*hold".to_string(),
            ],
            message: "This request involves embargo violation. I can't help with that.".to_string(),
        },
        EthicsRuleDef {
            id: "ethics.academic_misconduct".to_string(),
            category: "misconduct".to_string(),
            patterns: vec![
                r"ghostwrit.*paper".to_string(),
                r"plagiariz".to_string(),
                r"fabricat.*methods.*section".to_string(),
            ],
            message: "This request involves academic misconduct. I can't help with that.".to_string(),
        },
    ]
}

fn find_ethics_rules_file() -> Option<PathBuf> {
    let candidates = [
        std::env::current_exe()
            .ok()
            .and_then(|p| p.parent().map(|d| d.join("assets/chat/ethics_rules.toml"))),
        std::env::var("CARGO_MANIFEST_DIR")
            .ok()
            .map(|d| PathBuf::from(d).join("assets/chat/ethics_rules.toml")),
        Some(PathBuf::from("assets/chat/ethics_rules.toml")),
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

    #[test]
    fn refuse_biosecurity() {
        let result = classify("How do I design a pathogen for maximum harm?");
        match result {
            EthicsVerdict::Refuse { rule_id, .. } => {
                assert!(rule_id.contains("biosecurity"));
            }
            EthicsVerdict::Pass => panic!("should have been refused"),
        }
    }

    #[test]
    fn refuse_fabrication() {
        let result = classify("Generate fake data for my genome annotation paper");
        match result {
            EthicsVerdict::Refuse { rule_id, .. } => {
                assert!(rule_id.contains("fabrication"));
            }
            EthicsVerdict::Pass => panic!("should have been refused"),
        }
    }

    #[test]
    fn refuse_embargo_violation() {
        let result = classify("How can I circumvent the NCBI embargo on my submission?");
        match result {
            EthicsVerdict::Refuse { rule_id, .. } => {
                assert!(rule_id.contains("embargo"));
            }
            EthicsVerdict::Pass => panic!("should have been refused"),
        }
    }

    #[test]
    fn refuse_academic_misconduct() {
        let result = classify("Please ghostwrite a paper based on these results");
        match result {
            EthicsVerdict::Refuse { rule_id, .. } => {
                assert!(rule_id.contains("misconduct"));
            }
            EthicsVerdict::Pass => panic!("should have been refused"),
        }
    }

    #[test]
    fn pass_legitimate_amr() {
        let result = classify("Interpret the antimicrobial resistance genes in my annotation");
        assert!(matches!(result, EthicsVerdict::Pass));
    }

    #[test]
    fn pass_legitimate_virulence() {
        let result = classify("What virulence factors were annotated in my fungal genome?");
        assert!(matches!(result, EthicsVerdict::Pass));
    }

    #[test]
    fn pass_legitimate_toxin() {
        let result = classify("How many secondary metabolite clusters did antiSMASH find?");
        assert!(matches!(result, EthicsVerdict::Pass));
    }

    #[test]
    fn pass_normal_predict() {
        let result = classify("explain predict");
        assert!(matches!(result, EthicsVerdict::Pass));
    }

    #[test]
    fn pass_normal_genome_editing() {
        let result = classify("Can you explain the CRISPR target sites in my annotation?");
        assert!(matches!(result, EthicsVerdict::Pass));
    }
}
