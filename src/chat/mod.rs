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
pub mod bundle;
pub mod commands;
pub mod config;
pub mod context;
pub mod ethics;
pub mod history;
pub mod paste;
pub mod profile;
pub mod prompts;
pub mod render;
pub mod retrieval;
pub mod rules;
pub mod validator;

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

    let verbosity = if cfg.trace {
        render::Verbosity::Trace
    } else if cfg.verbose {
        render::Verbosity::Verbose
    } else {
        render::Verbosity::Short
    };

    // ── Ethics check (runs on any user-provided text) ──
    let all_input = format!("{} {}", stage_or_paste, dir.as_deref().unwrap_or(""));
    match ethics::classify(&all_input) {
        ethics::EthicsVerdict::Refuse {
            rule_id: _,
            message,
        } => {
            println!("\n  {}{}{}", C_YELLOW, message, C_RESET);
            return Ok(());
        }
        ethics::EthicsVerdict::Pass => {}
    }

    // ── Paste mode ──
    if stage_or_paste == "__paste__" {
        return run_paste_mode(&cfg, verbosity);
    }

    // ── Stage mode ──
    let stage_name = &stage_or_paste;
    let stage = context::Stage::from_str(stage_name)?;
    let work_dir = resolve_dir(stage_name, dir.as_deref())?;

    println!(
        "\n  {}{}myconote explain{} — stage: {}{}{}, dir: {}",
        C_BOLD,
        C_CYAN,
        C_RESET,
        C_GREEN,
        stage_name,
        C_RESET,
        work_dir.display()
    );

    // Hint for first-time users: suggest the relevant learn lesson
    if is_first_run() {
        if let Some(lesson) = learn_lesson_for_stage(stage_name) {
            println!(
                "  {}New to this stage? Run `myconote-cli learn {}` for a walkthrough.{}",
                C_DIM, lesson, C_RESET
            );
        }
    }

    // ── Build context ──
    let ctx = context::build_context(stage, &work_dir)?;

    // ── Run rule engine ──
    let findings = rules::evaluate(&ctx);

    // ── Retrieve knowledge + papers ──
    let query = build_retrieval_query(&ctx, &findings);
    let mut retrieved = retrieval::retrieve_knowledge(stage_name, &query, 5);

    // Retrieve from paper corpus (if available)
    let corpus_dir = cfg.corpus_dir.as_deref();
    let paper_citations = retrieval::retrieve_papers(&query, 3, corpus_dir);
    if paper_citations.is_empty() && retrieval::corpus_status().is_none() {
        // Only mention on first run (when trace is enabled)
        if verbosity == render::Verbosity::Trace {
            println!(
                "  {}Paper citations unavailable — run `myconote-cli setup --chat-corpus` to enable.{}",
                C_DIM, C_RESET
            );
        }
    }
    retrieved.extend(paper_citations);

    // ── Generate command recommendations ──
    let recs = commands::recommend(&ctx, &findings);

    // ── Load user profile ──
    let user_profile = profile::UserProfile::load();

    // ── Rules-only mode or LLM unavailable ──
    if cfg.no_llm {
        render::render_rules_only_header();
        render::render_findings(&findings, verbosity);
        render::render_commands(&recs);
        render::render_disclaimer();
        write_bundle(
            stage_name, &cfg, &ctx, &findings, &retrieved, &recs, None, None, 0, 0,
        )?;
        record_history(stage_name, &ctx, &findings, false, &cfg)?;
        return Ok(());
    }

    // ── Assemble prompt ──
    let system_prompt = prompts::render(stage_name, &ctx, &findings, &retrieved, &user_profile);
    let user_msg = prompts::user_message(stage_name);

    // ── Dry-run ──
    if cfg.dry_run {
        println!(
            "\n  {}--dry-run{}: showing assembled prompt (no LLM call)",
            C_YELLOW, C_RESET
        );
        render::render_trace(&system_prompt, retrieved.len(), 0, 0);
        return Ok(());
    }

    // ── Try LLM ──
    println!(
        "  Model: {}{}{} @ {}{}{}",
        C_BOLD, cfg.model, C_RESET, C_DIM, cfg.endpoint, C_RESET
    );

    let llm_backend =
        backend::ollama::OllamaBackend::new(&cfg.endpoint, &cfg.model, cfg.timeout_s)?;

    if !llm_backend.is_available() {
        render::render_llm_unavailable();
        render::render_findings(&findings, verbosity);
        render::render_commands(&recs);
        render::render_disclaimer();
        write_bundle(
            stage_name,
            &cfg,
            &ctx,
            &findings,
            &retrieved,
            &recs,
            Some(&system_prompt),
            None,
            0,
            0,
        )?;
        record_history(stage_name, &ctx, &findings, false, &cfg)?;
        return Ok(());
    }

    // ── Call LLM ──
    use backend::ChatBackend;
    let messages = vec![
        backend::Message::system(&system_prompt),
        backend::Message::user(&user_msg),
    ];

    let response = match llm_backend.chat(&messages) {
        Ok(msg) => msg.content,
        Err(e) => {
            println!("\n  {}LLM error: {}{}", C_YELLOW, e, C_RESET);
            render::render_findings(&findings, verbosity);
            render::render_commands(&recs);
            render::render_disclaimer();
            write_bundle(
                stage_name,
                &cfg,
                &ctx,
                &findings,
                &retrieved,
                &recs,
                Some(&system_prompt),
                None,
                0,
                0,
            )?;
            record_history(stage_name, &ctx, &findings, false, &cfg)?;
            return Ok(());
        }
    };

    // ── Validate output (strip lines with hallucinated subcommands) ──
    let validation = validator::validate_commands(&response);

    // ── Render ──
    render::render_findings(&findings, verbosity);
    render::render_llm_response(&validation.text);
    render::render_commands(&recs);

    if verbosity == render::Verbosity::Trace {
        render::render_trace(
            &system_prompt,
            retrieved.len(),
            validation.commands_kept,
            validation.commands_removed,
        );
    }

    render::render_disclaimer();

    // ── Write bundle ──
    write_bundle(
        stage_name,
        &cfg,
        &ctx,
        &findings,
        &retrieved,
        &recs,
        Some(&system_prompt),
        Some(&validation.text),
        validation.commands_kept,
        validation.commands_removed,
    )?;

    // ── Record history ──
    record_history(stage_name, &ctx, &findings, true, &cfg)?;

    Ok(())
}

/// Paste mode: read stdin, detect format, run rules, optionally call LLM.
fn run_paste_mode(_cfg: &ChatConfig, verbosity: render::Verbosity) -> Result<()> {
    println!(
        "\n  {}{}Paste mode{} — reading from stdin...",
        C_BOLD, C_CYAN, C_RESET
    );

    let (input, format) = paste::read_stdin();
    if input.trim().is_empty() {
        println!("  {}No input received.{}", C_DIM, C_RESET);
        return Ok(());
    }

    // Ethics check on pasted content
    match ethics::classify(&input) {
        ethics::EthicsVerdict::Refuse {
            rule_id: _,
            message,
        } => {
            println!("\n  {}{}{}", C_YELLOW, message, C_RESET);
            return Ok(());
        }
        ethics::EthicsVerdict::Pass => {}
    }

    println!("  Detected format: {}{}{}", C_GREEN, format, C_RESET);
    println!(
        "  Input: {} lines, {} bytes",
        input.lines().count(),
        input.len()
    );

    // For paste mode, determine the most likely stage from the format
    let inferred_stage = match format {
        paste::PasteFormat::Gff3 => "predict",
        paste::PasteFormat::FastaHeader => "sort",
        paste::PasteFormat::ValidationError => "submit",
        paste::PasteFormat::LogOutput => "predict",
        paste::PasteFormat::Tsv => "annotate",
        paste::PasteFormat::Unknown => "predict",
    };

    // Build a minimal context from the pasted text
    let ctx = context::StageContext {
        stage: inferred_stage.to_string(),
        dir: "(stdin)".to_string(),
        artifacts: vec![context::ArtifactSummary {
            name: "(pasted input)".to_string(),
            size_bytes: input.len() as u64,
            line_count: Some(input.lines().count()),
            preview: Some(input.chars().take(50_000).collect()),
        }],
        stats: None,
        notes: vec![format!("Paste mode: detected as {}", format)],
    };

    let findings = rules::evaluate(&ctx);
    let recs = commands::recommend(&ctx, &findings);

    render::render_findings(&findings, verbosity);
    render::render_commands(&recs);
    render::render_disclaimer();

    Ok(())
}

fn build_retrieval_query(ctx: &context::StageContext, findings: &[rules::Finding]) -> String {
    let mut parts = vec![ctx.stage.clone()];
    for f in findings {
        parts.push(f.message.clone());
    }
    for note in &ctx.notes {
        parts.push(note.clone());
    }
    parts.join(" ")
}

fn write_bundle(
    stage_name: &str,
    cfg: &ChatConfig,
    ctx: &context::StageContext,
    findings: &[rules::Finding],
    retrieved: &[retrieval::Citation],
    recs: &[commands::CommandRecommendation],
    prompt: Option<&str>,
    response: Option<&str>,
    commands_kept: usize,
    commands_removed: usize,
) -> Result<()> {
    let b = bundle::ExplainBundle {
        version: env!("CARGO_PKG_VERSION").to_string(),
        timestamp: bundle::timestamp(),
        stage: stage_name.to_string(),
        model: if cfg.no_llm {
            None
        } else {
            Some(cfg.model.clone())
        },
        endpoint: if cfg.no_llm {
            None
        } else {
            Some(cfg.endpoint.clone())
        },
        prompt: prompt.map(|s| s.to_string()),
        context: ctx.clone(),
        findings: findings.to_vec(),
        retrieved: retrieved.to_vec(),
        commands: recs.to_vec(),
        response: response.map(|s| s.to_string()),
        commands_kept,
        commands_removed,
    };
    match b.write() {
        Ok(dir) => {
            println!("  {}Bundle: {}{}", C_DIM, dir.display(), C_RESET);
        }
        Err(e) => {
            eprintln!("  {}Bundle write failed: {}{}", C_DIM, e, C_RESET);
        }
    }
    Ok(())
}

fn record_history(
    stage_name: &str,
    ctx: &context::StageContext,
    findings: &[rules::Finding],
    llm_used: bool,
    cfg: &ChatConfig,
) -> Result<()> {
    let entry = history::HistoryEntry {
        timestamp: bundle::timestamp(),
        stage: stage_name.to_string(),
        dir: ctx.dir.clone(),
        findings_count: findings.len(),
        llm_used,
        model: if llm_used {
            Some(cfg.model.clone())
        } else {
            None
        },
    };
    // Non-fatal: don't fail the explain call if history can't be written
    let _ = history::append(&entry);
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
                    return Err(MycoNoteError::ChatConfig(
                        "--dir requires a path argument".to_string(),
                    ));
                }
                dir = Some(args[i + 1].clone());
                i += 2;
            }
            "--model" => {
                if i + 1 >= args.len() {
                    return Err(MycoNoteError::ChatConfig(
                        "--model requires a value".to_string(),
                    ));
                }
                cli.model = Some(args[i + 1].clone());
                i += 2;
            }
            "--endpoint" => {
                if i + 1 >= args.len() {
                    return Err(MycoNoteError::ChatConfig(
                        "--endpoint requires a URL".to_string(),
                    ));
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
             Or use: myconote-cli explain --paste"
                .to_string(),
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
    println!("  {}Usage:{}", C_BOLD, C_RESET);
    println!("    myconote-cli explain <stage> [options]");
    println!("    myconote-cli explain --paste [options]");
    println!();
    println!("  {}Stages:{}", C_BOLD, C_RESET);
    println!("    sort       Sort + rename genome contigs");
    println!("    mask       Repeat masking results");
    println!("    train      Gene predictor training reports");
    println!("    predict    Gene prediction summary + GFF3");
    println!("    update     PASA UTR update results");
    println!("    annotate   Functional annotation coverage");
    println!("    submit     NCBI validation + submission files");
    println!();
    println!("  {}Options:{}", C_BOLD, C_RESET);
    println!("    --dir <path>       Working directory (default: auto-discover)");
    println!("    --model <name>     Ollama model (default: llama3.1)");
    println!("    --endpoint <url>   Ollama endpoint (default: http://localhost:11434)");
    println!("    --no-llm           Rules-only mode: no Ollama needed");
    println!("    --verbose          Show all findings with evidence");
    println!("    --trace            Show prompt, retrieval scores, validator decisions");
    println!("    --dry-run          Print assembled prompt without calling LLM");
    println!("    --paste            Read from stdin instead of scanning a directory");
    println!();
    println!("  {}Examples:{}", C_BOLD, C_RESET);
    println!("    myconote-cli explain predict");
    println!("    myconote-cli explain predict --no-llm");
    println!("    myconote-cli explain annotate --dir 05_annotate/ --verbose");
    println!("    echo \"ERROR: SEQ_FEAT.NoStop\" | myconote-cli explain --paste");
    println!();
    println!(
        "  {}Note:{} All findings and interpretations are suggestive, not definitive.",
        C_YELLOW, C_RESET
    );
    println!("  Always verify results in the context of your organism, assembly, and");
    println!("  research goals before drawing scientific conclusions.");
    println!();
}

/// Check if this is the user's first time running explain (no history directory).
fn is_first_run() -> bool {
    let home = match std::env::var("HOME") {
        Ok(h) => h,
        Err(_) => return true,
    };
    let history_dir = std::path::PathBuf::from(home)
        .join(".myconote")
        .join("chat_history");
    !history_dir.exists()
}

/// Map a pipeline stage to the corresponding `learn` lesson number.
fn learn_lesson_for_stage(stage: &str) -> Option<&'static str> {
    match stage {
        "sort" | "mask" => Some("3"),
        "train" => Some("3"),
        "predict" => Some("4"),
        "update" => Some("4"),
        "annotate" => Some("5"),
        "submit" => Some("7"),
        _ => None,
    }
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
