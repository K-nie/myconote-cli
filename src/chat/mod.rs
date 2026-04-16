/// LLM-powered result interpreter — `myconote-cli explain`
///
/// A local-first chatbot that interprets genome annotation pipeline outputs.
/// Works in two modes:
///   - Stage mode: `explain <stage> [--dir <path>]` — interprets a pipeline stage
///   - Paste mode: `explain --paste` — analyzes pasted output from stdin
///
/// Can run with or without an LLM:
///   - Full mode (default): deterministic rules + Ollama LLM interpretation
///   - Rules-only mode (`--no-llm`): no Ollama needed, prints findings + commands
pub mod backend;
pub mod config;
pub mod context;
pub mod rules;

use crate::learn::{C_BOLD, C_CYAN, C_DIM, C_GREEN, C_RESET, C_YELLOW};
use crate::utils::error::{MycoNoteError, Result};
use config::{ChatConfig, CliOverrides};

/// Entry point for `myconote-cli explain <args>`.
pub fn run_explain(args: &[String]) -> Result<()> {
    // ── Help ──
    if args.is_empty() || args.iter().any(|a| a == "--help" || a == "-h") {
        print_explain_help();
        return Ok(());
    }

    // ── Parse CLI flags ──
    let (stage_or_paste, dir, cli) = parse_args(args)?;

    // ── Load config (TOML → env → CLI) ──
    let cfg = ChatConfig::load(&cli)?;

    if !cfg.enabled {
        println!(
            "  {}explain is disabled.{} Set [chat] enabled = true in ~/.myconote/config.toml",
            C_DIM, C_RESET
        );
        return Ok(());
    }

    // ── Determine mode ──
    let is_paste = stage_or_paste == "__paste__";

    if is_paste {
        println!(
            "  {}{}Paste mode{} — reading from stdin...",
            C_BOLD, C_CYAN, C_RESET
        );
        println!(
            "  {}(paste mode will be available in a future commit){}", C_DIM, C_RESET
        );
        return Ok(());
    }

    // ── Stage mode ──
    let stage_name = &stage_or_paste;
    let valid_stages = ["sort", "mask", "train", "predict", "update", "annotate", "submit"];
    if !valid_stages.contains(&stage_name.as_str()) {
        return Err(MycoNoteError::ChatContext(format!(
            "unknown stage '{}'. Expected one of: {}",
            stage_name,
            valid_stages.join(", ")
        )));
    }

    // Resolve directory
    let work_dir = resolve_dir(stage_name, dir.as_deref())?;

    println!(
        "\n  {}{}myconote explain{} — stage: {}{}{}, dir: {}",
        C_BOLD, C_CYAN, C_RESET, C_GREEN, stage_name, C_RESET, work_dir.display()
    );

    if cfg.no_llm {
        println!(
            "  {}Rules-only mode{} — no LLM call will be made.\n",
            C_YELLOW, C_RESET
        );
        println!("  {}(rule engine will be available in a future commit){}", C_DIM, C_RESET);
        return Ok(());
    }

    // ── Try LLM ──
    println!(
        "  Model: {}{}{} @ {}{}{}",
        C_BOLD, cfg.model, C_RESET, C_DIM, cfg.endpoint, C_RESET
    );

    let backend = backend::ollama::OllamaBackend::new(&cfg.endpoint, &cfg.model, cfg.timeout_s)?;

    if !backend.is_available() {
        println!(
            "\n  {}LLM unavailable{} — showing rule-based findings only.",
            C_YELLOW, C_RESET
        );
        println!(
            "  Start Ollama for a full interpretation: {}ollama serve{}",
            C_BOLD, C_RESET
        );
        println!("  {}(rule engine will be available in a future commit){}", C_DIM, C_RESET);
        return Ok(());
    }

    if cfg.dry_run {
        println!(
            "\n  {}--dry-run{}: would assemble prompt for stage '{}' from {}",
            C_YELLOW, C_RESET, stage_name, work_dir.display()
        );
        println!("  {}(prompt assembly will be available in a future commit){}", C_DIM, C_RESET);
        return Ok(());
    }

    // Placeholder: full pipeline (context → rules → retrieval → prompt → LLM → validate → render → bundle)
    println!(
        "\n  {}(full pipeline will be built in subsequent commits){}", C_DIM, C_RESET
    );

    Ok(())
}

/// Parse args into (stage_or_paste, optional_dir, cli_overrides).
fn parse_args(args: &[String]) -> Result<(String, Option<String>, CliOverrides)> {
    let mut stage = None;
    let mut dir = None;
    let mut cli = CliOverrides::default();
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--paste" => {
                stage = Some("__paste__".to_string());
                i += 1;
            }
            "--dir" => {
                if i + 1 >= args.len() {
                    return Err(MycoNoteError::ChatConfig("--dir requires a path argument".to_string()));
                }
                dir = Some(args[i + 1].clone());
                i += 2;
            }
            "--model" => {
                if i + 1 >= args.len() {
                    return Err(MycoNoteError::ChatConfig("--model requires a value".to_string()));
                }
                cli.model = Some(args[i + 1].clone());
                i += 2;
            }
            "--endpoint" => {
                if i + 1 >= args.len() {
                    return Err(MycoNoteError::ChatConfig("--endpoint requires a URL".to_string()));
                }
                cli.endpoint = Some(args[i + 1].clone());
                i += 2;
            }
            "--no-llm" => {
                cli.no_llm = true;
                i += 1;
            }
            "--verbose" => {
                cli.verbose = true;
                i += 1;
            }
            "--trace" => {
                cli.trace = true;
                i += 1;
            }
            "--dry-run" => {
                cli.dry_run = true;
                i += 1;
            }
            other if !other.starts_with('-') && stage.is_none() => {
                stage = Some(other.to_lowercase());
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }

    let stage = stage.ok_or_else(|| {
        MycoNoteError::ChatContext(
            "missing stage. Usage: myconote-cli explain <stage> [options]\n  \
             Stages: sort | mask | train | predict | update | annotate | submit\n  \
             Or use: myconote-cli explain --paste".to_string()
        )
    })?;

    Ok((stage, dir, cli))
}

/// Resolve the working directory for a stage.
///
/// Priority: explicit `--dir` > auto-discover `<NN>_<stage>/` in cwd > cwd.
fn resolve_dir(stage: &str, explicit: Option<&str>) -> Result<std::path::PathBuf> {
    use std::path::PathBuf;

    if let Some(dir) = explicit {
        let p = PathBuf::from(dir);
        if !p.exists() {
            return Err(MycoNoteError::ChatContext(format!(
                "directory not found: {}",
                p.display()
            )));
        }
        return Ok(p);
    }

    // Auto-discover: look for <NN>_<stage>/ in cwd
    let cwd = std::env::current_dir()
        .map_err(|e| MycoNoteError::ChatContext(format!("cannot read cwd: {}", e)))?;

    let pattern = regex::Regex::new(&format!(r"^\d+_{}", regex::escape(stage)))
        .map_err(|e| MycoNoteError::ChatContext(format!("regex error: {}", e)))?;

    if let Ok(entries) = std::fs::read_dir(&cwd) {
        for entry in entries.flatten() {
            if entry.file_type().map(|t| t.is_dir()).unwrap_or(false) {
                if let Some(name) = entry.file_name().to_str() {
                    if pattern.is_match(name) {
                        return Ok(entry.path());
                    }
                }
            }
        }
    }

    // Fall back to cwd
    Ok(cwd)
}

fn print_explain_help() {
    println!();
    println!(
        "  {}{}myconote explain{} — LLM-powered result interpreter",
        C_BOLD, C_CYAN, C_RESET
    );
    println!(
        "  {}Interprets pipeline outputs, recommends next commands, cites Q1 papers.{}",
        C_DIM, C_RESET
    );
    println!();
    println!("  {}Usage:{}",  C_BOLD, C_RESET);
    println!("    myconote-cli explain <stage> [options]");
    println!("    myconote-cli explain --paste [options]");
    println!();
    println!("  {}Stages:{}",  C_BOLD, C_RESET);
    println!("    sort       Sort + rename genome contigs");
    println!("    mask       Repeat masking results");
    println!("    train      Gene predictor training reports");
    println!("    predict    Gene prediction summary + GFF3");
    println!("    update     PASA UTR update results");
    println!("    annotate   Functional annotation coverage");
    println!("    submit     NCBI validation + submission files");
    println!();
    println!("  {}Options:{}",  C_BOLD, C_RESET);
    println!("    --dir <path>       Working directory (default: auto-discover)");
    println!("    --model <name>     Ollama model (default: llama3.1)");
    println!("    --endpoint <url>   Ollama endpoint (default: http://localhost:11434)");
    println!("    --no-llm           Rules-only mode: no Ollama needed");
    println!("    --verbose          Show all findings with evidence");
    println!("    --trace            Show prompt, retrieval scores, validator decisions");
    println!("    --dry-run          Print assembled prompt without calling LLM");
    println!("    --paste            Read from stdin instead of scanning a directory");
    println!();
    println!("  {}Examples:{}",  C_BOLD, C_RESET);
    println!("    myconote-cli explain predict");
    println!("    myconote-cli explain predict --no-llm");
    println!("    myconote-cli explain annotate --dir 05_annotate/ --verbose");
    println!("    echo \"ERROR: SEQ_FEAT.NoStop\" | myconote-cli explain --paste");
    println!();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_stage_only() {
        let args: Vec<String> = vec!["predict".to_string()];
        let (stage, dir, cli) = parse_args(&args).unwrap();
        assert_eq!(stage, "predict");
        assert!(dir.is_none());
        assert!(!cli.no_llm);
    }

    #[test]
    fn parse_stage_with_flags() {
        let args: Vec<String> = vec![
            "annotate".to_string(),
            "--dir".to_string(),
            "/tmp/test".to_string(),
            "--no-llm".to_string(),
            "--verbose".to_string(),
        ];
        let (stage, dir, cli) = parse_args(&args).unwrap();
        assert_eq!(stage, "annotate");
        assert_eq!(dir.unwrap(), "/tmp/test");
        assert!(cli.no_llm);
        assert!(cli.verbose);
    }

    #[test]
    fn parse_paste_mode() {
        let args: Vec<String> = vec!["--paste".to_string(), "--no-llm".to_string()];
        let (stage, _dir, cli) = parse_args(&args).unwrap();
        assert_eq!(stage, "__paste__");
        assert!(cli.no_llm);
    }

    #[test]
    fn parse_missing_stage_errors() {
        let args: Vec<String> = vec!["--verbose".to_string()];
        let result = parse_args(&args);
        assert!(result.is_err());
    }

    #[test]
    fn unknown_stage_errors() {
        let result = run_explain(&[String::from("bogus")]);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("unknown stage 'bogus'"));
        assert!(err.contains("sort"));
    }
}
