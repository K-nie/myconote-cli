//! Output validator for LLM-generated interpretations.
//!
//! Catches a real failure mode observed in production: the LLM confidently
//! suggests `myconote-cli <subcommand>` invocations that do not exist
//! (e.g. `myconote-cli repeatmasker`, `myconote-cli gff3stats`). Lines that
//! reference any unknown subcommand are stripped from the response before
//! it reaches the user.
//!
//! The previous citation-validator (which tried to strip fabricated
//! `[knowledge:…]` / `[paper:…]` tags) has been removed — it rarely fired
//! because the model almost never emitted citation tags, and it gave a
//! false sense of rigour.

/// Authoritative list of myconote-cli subcommands as exposed by the CLI
/// router in `src/main.rs`. Any `myconote-cli <word>` reference in LLM
/// output whose word is not in this list is considered hallucinated and
/// the whole line is stripped.
pub const VALID_SUBCOMMANDS: &[&str] = &[
    // Pipeline
    "sort",
    "mask",
    "train",
    "predict",
    "update",
    "annotate",
    "submit",
    "batch",
    "remote",
    // Analysis
    "stats",
    "quant",
    "fetch-rna",
    "de-template",
    "compare",
    "convert",
    "clean",
    "fix", // Utility
    "install",
    "check",
    "setup",
    "species",
    "learn",
    "explain", // Meta
    "help",
];

/// Validate LLM output: strip any line that invokes a non-existent
/// `myconote-cli <subcommand>`. Returns the cleaned text plus counts.
pub fn validate_commands(response: &str) -> ValidationResult {
    let re = regex::Regex::new(r"myconote-cli\s+([A-Za-z][A-Za-z0-9_-]*)").unwrap();

    let mut kept = 0;
    let mut removed = 0;
    let mut out_lines: Vec<String> = Vec::new();

    for line in response.lines() {
        let mentions: Vec<&str> = re
            .captures_iter(line)
            .map(|c| c.get(1).unwrap().as_str())
            .collect();

        if mentions.is_empty() {
            out_lines.push(line.to_string());
            continue;
        }

        let any_bad = mentions
            .iter()
            .any(|m| !VALID_SUBCOMMANDS.contains(&m.to_ascii_lowercase().as_str()));

        if any_bad {
            removed += 1;
        } else {
            kept += mentions.len();
            out_lines.push(line.to_string());
        }
    }

    ValidationResult {
        text: out_lines.join("\n").trim().to_string(),
        commands_kept: kept,
        commands_removed: removed,
        is_valid: removed == 0,
    }
}

pub struct ValidationResult {
    pub text: String,
    pub commands_kept: usize,
    pub commands_removed: usize,
    pub is_valid: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_subcommand_line_is_kept() {
        let r = validate_commands("Next: run `myconote-cli stats predict/consensus.gff3`.");
        assert!(r.is_valid);
        assert_eq!(r.commands_kept, 1);
        assert_eq!(r.commands_removed, 0);
        assert!(r.text.contains("myconote-cli stats"));
    }

    #[test]
    fn hallucinated_subcommand_is_stripped() {
        // Exactly the failure mode observed in production: the LLM
        // invented `repeatmasker` and `gff3stats`.
        let response = "\
You should run `myconote-cli repeatmasker -n 5 ...` to mask repeats.
Also run `myconote-cli gff3stats -i consensus.gff3` for summaries.
For real annotation, `myconote-cli annotate predict/consensus.gff3 --fasta genome.fa`.
Finally, `myconote-cli stats predict/consensus.gff3 --taxon fungi`.";
        let r = validate_commands(response);
        assert!(!r.is_valid);
        assert_eq!(r.commands_removed, 2);
        assert!(!r.text.contains("repeatmasker"));
        assert!(!r.text.contains("gff3stats"));
        assert!(r.text.contains("myconote-cli annotate"));
        assert!(r.text.contains("myconote-cli stats"));
    }

    #[test]
    fn prose_without_commands_passes_through() {
        let r = validate_commands("Gene count is within the fungal norm.");
        assert!(r.is_valid);
        assert_eq!(r.commands_kept, 0);
        assert_eq!(r.commands_removed, 0);
        assert_eq!(r.text, "Gene count is within the fungal norm.");
    }

    #[test]
    fn multiple_commands_on_one_line_all_valid() {
        let r = validate_commands(
            "First `myconote-cli predict ...`, then `myconote-cli annotate ...`.",
        );
        assert!(r.is_valid);
        assert_eq!(r.commands_kept, 2);
    }

    #[test]
    fn case_insensitive_command_match() {
        let r = validate_commands("Try `myconote-cli Stats file.gff3`.");
        assert!(r.is_valid);
        assert_eq!(r.commands_kept, 1);
    }

    #[test]
    fn removed_subcommands_are_stripped() {
        // phylogeny and place were real subcommands in v0.2.0 but were
        // removed when the tool's scope narrowed to fungal annotation.
        // If someone re-adds either to VALID_SUBCOMMANDS without wiring up
        // the module, this test fails and flags the inconsistency.
        for removed in ["phylogeny", "place"] {
            assert!(
                !VALID_SUBCOMMANDS.contains(&removed),
                "'{removed}' should no longer be in VALID_SUBCOMMANDS",
            );
            let line = format!("Run `myconote-cli {removed} input.fa` to continue.");
            let r = validate_commands(&line);
            assert!(
                !r.is_valid,
                "line invoking removed subcommand '{removed}' should be invalid",
            );
            assert_eq!(r.commands_removed, 1);
            assert!(!r.text.contains(removed));
        }
    }
}
