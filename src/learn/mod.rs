/// Interactive tutorial system — "myconote learn"
///
/// A self-paced, swirl-style interactive learning environment that teaches
/// users how to use myconote-cli by walking them through real commands
/// with explanations, quizzes, and hands-on exercises.
///
/// Inspired by R's swirl package: the user types answers directly into
/// the terminal and gets instant feedback, with the ability to skip
/// questions, get hints, and resume where they left off.
///
/// Usage:
///   myconote-cli learn                   # list available lessons
///   myconote-cli learn 1                 # start lesson 1
///   myconote-cli learn --reset           # reset all progress
///   myconote-cli learn --resume          # resume last lesson
pub mod engine;
pub mod lessons;
pub mod progress;

use crate::utils::error::Result;

// ─────────────────────────────────────────────────────────────────────────────
// ANSI colours for the tutorial UI
// ─────────────────────────────────────────────────────────────────────────────

pub const C_RESET: &str = "\x1b[0m";
pub const C_BOLD: &str = "\x1b[1m";
pub const C_DIM: &str = "\x1b[2m";
pub const C_ITALIC: &str = "\x1b[3m";
pub const C_GREEN: &str = "\x1b[32m";
pub const C_YELLOW: &str = "\x1b[33m";
pub const C_CYAN: &str = "\x1b[36m";
pub const C_RED: &str = "\x1b[31m";
pub const C_MAGENTA: &str = "\x1b[35m";
pub const C_BLUE: &str = "\x1b[34m";
pub const C_BG_DARK: &str = "\x1b[48;5;236m";

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_learn(args: &[String]) -> Result<()> {
    if args.is_empty() || args.iter().any(|a| a == "--help" || a == "-h") {
        print_learn_help();
        return Ok(());
    }

    if args.iter().any(|a| a == "--reset") {
        progress::reset_progress()?;
        println!("  {}✓{} Tutorial progress reset.", C_GREEN, C_RESET);
        return Ok(());
    }

    if args.iter().any(|a| a == "--list" || a == "list") {
        list_lessons();
        return Ok(());
    }

    if args.iter().any(|a| a == "--resume") {
        let next = progress::get_resume_point()?;
        println!("  Resuming from lesson {}...\n", next);
        return engine::run_lesson(next);
    }

    // Parse lesson number
    if let Some(num_str) = args.first() {
        if let Ok(n) = num_str.parse::<usize>() {
            return engine::run_lesson(n);
        }
        // Try lesson name match
        if let Some(idx) = lessons::find_lesson_by_name(num_str) {
            return engine::run_lesson(idx);
        }
    }

    list_lessons();
    Ok(())
}

fn print_learn_help() {
    println!();
    println!(
        "  {}{}myconote learn{} — Interactive Tutorial System",
        C_BOLD, C_CYAN, C_RESET
    );
    println!(
        "  {}Inspired by R's swirl: learn by doing, right in your terminal.{}",
        C_DIM, C_RESET
    );
    println!();
    println!("  {}Usage:{}", C_BOLD, C_RESET);
    println!("    myconote-cli learn              List all available lessons");
    println!("    myconote-cli learn 1            Start lesson 1");
    println!("    myconote-cli learn basics       Start a lesson by name");
    println!("    myconote-cli learn --resume     Resume where you left off");
    println!("    myconote-cli learn --reset      Reset all progress");
    println!();
    println!("  {}During a lesson:{}", C_BOLD, C_RESET);
    println!("    Type your answer and press Enter");
    println!(
        "    {}skip{}      Skip the current question",
        C_YELLOW, C_RESET
    );
    println!("    {}hint{}      Get a hint", C_YELLOW, C_RESET);
    println!(
        "    {}info{}      Show background information",
        C_YELLOW, C_RESET
    );
    println!(
        "    {}quit{}      Save progress and exit",
        C_YELLOW, C_RESET
    );
    println!("    {}bye{}       Same as quit", C_YELLOW, C_RESET);
    println!();

    list_lessons();
}

fn list_lessons() {
    let all = lessons::all_lessons();
    let completed = progress::load_progress().unwrap_or_default();

    println!();
    println!("  {}{}Available Lessons{}", C_BOLD, C_CYAN, C_RESET);
    println!("  {}", "─".repeat(60));

    for (i, lesson) in all.iter().enumerate() {
        let num = i + 1;
        let status = if completed.contains(&num) {
            format!("{}✓{}", C_GREEN, C_RESET)
        } else {
            format!("{}○{}", C_DIM, C_RESET)
        };
        let duration = format!("{}~{} min{}", C_DIM, lesson.est_minutes, C_RESET);
        println!(
            "  {} {}{}. {}{:<40} {}",
            status, C_BOLD, num, C_RESET, lesson.title, duration
        );
        println!("     {}{}{}", C_DIM, lesson.description, C_RESET);
    }

    let done = completed.len();
    let total = all.len();
    println!();
    if done == 0 {
        println!("  Start with: {}myconote-cli learn 1{}", C_CYAN, C_RESET);
    } else if done < total {
        let next = (1..=total).find(|n| !completed.contains(n)).unwrap_or(1);
        println!("  Progress: {}/{} lessons completed", done, total);
        println!(
            "  Continue: {}myconote-cli learn {}{}",
            C_CYAN, next, C_RESET
        );
    } else {
        println!(
            "  {}Congratulations! All {} lessons completed!{}",
            C_GREEN, total, C_RESET
        );
    }
    println!();
}
