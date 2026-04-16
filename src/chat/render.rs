use crate::learn::{C_BOLD, C_CYAN, C_DIM, C_GREEN, C_RED, C_RESET, C_YELLOW};
use super::commands::CommandRecommendation;
use super::rules::{Finding, Severity};

/// Verbosity level for output rendering.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Verbosity {
    Short,
    Verbose,
    Trace,
}

/// Render findings to stdout.
pub fn render_findings(findings: &[Finding], verbosity: Verbosity) {
    if findings.is_empty() {
        println!("  {}No findings.{}", C_DIM, C_RESET);
        return;
    }

    println!("\n  {}{}Findings{}", C_BOLD, C_CYAN, C_RESET);
    println!("  {}{}{}", C_DIM, "─".repeat(60), C_RESET);

    for f in findings {
        let (icon, color) = match f.severity {
            Severity::Critical => ("!", C_RED),
            Severity::Warning => ("~", C_YELLOW),
            Severity::Info => ("i", C_GREEN),
        };

        println!("  {}{}{}{} {}", color, C_BOLD, icon, C_RESET, f.message);

        if verbosity == Verbosity::Verbose || verbosity == Verbosity::Trace {
            println!("    {}Rule: {}{}", C_DIM, f.rule_id, C_RESET);
            println!("    {}Evidence: {}{}", C_DIM, f.evidence, C_RESET);
            if let Some(ref cite) = f.citation {
                println!("    {}Citation: {}{}", C_DIM, cite, C_RESET);
            }
        }
    }
    println!();
}

/// Render command recommendations to stdout.
pub fn render_commands(recs: &[CommandRecommendation]) {
    if recs.is_empty() {
        return;
    }

    println!("  {}{}Recommended commands{}", C_BOLD, C_CYAN, C_RESET);
    println!("  {}{}{}", C_DIM, "─".repeat(60), C_RESET);

    for (i, rec) in recs.iter().enumerate() {
        println!("  {}{}. {}{}{}", C_BOLD, i + 1, C_GREEN, rec.command, C_RESET);
        println!("     {}{}{}", C_DIM, rec.rationale, C_RESET);
    }
    println!();
}

/// Render the LLM response to stdout.
pub fn render_llm_response(response: &str) {
    println!("\n  {}{}Interpretation{}", C_BOLD, C_CYAN, C_RESET);
    println!("  {}{}{}", C_DIM, "─".repeat(60), C_RESET);

    for line in response.lines() {
        println!("  {}", line);
    }
    println!();
}

/// Render trace info (prompt, retrieval scores) for --trace mode.
pub fn render_trace(
    prompt: &str,
    retrieved_count: usize,
    citations_kept: usize,
    citations_removed: usize,
) {
    println!("\n  {}{}Trace{}", C_BOLD, C_YELLOW, C_RESET);
    println!("  {}{}{}", C_DIM, "─".repeat(60), C_RESET);
    println!("  Prompt length: {} chars", prompt.len());
    println!("  Retrieved snippets: {}", retrieved_count);
    println!("  Citations kept: {}", citations_kept);
    println!("  Citations removed: {}", citations_removed);
    println!("\n  {}Full prompt:{}", C_DIM, C_RESET);
    for line in prompt.lines() {
        println!("  {}{}{}", C_DIM, line, C_RESET);
    }
    println!();
}

/// Render the "LLM unavailable" fallback message.
pub fn render_llm_unavailable() {
    println!(
        "\n  {}LLM unavailable{} — showing rule-based findings only.",
        C_YELLOW, C_RESET
    );
    println!(
        "  Start Ollama for a full interpretation: {}ollama serve{}",
        C_BOLD, C_RESET
    );
}

/// Render a "rules-only mode" header.
pub fn render_rules_only_header() {
    println!(
        "\n  {}Rules-only mode{} — no LLM call.",
        C_YELLOW, C_RESET
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verbosity_equality() {
        assert_eq!(Verbosity::Short, Verbosity::Short);
        assert_ne!(Verbosity::Short, Verbosity::Verbose);
    }
}
