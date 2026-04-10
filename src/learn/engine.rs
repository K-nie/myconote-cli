use super::lessons::{all_lessons, Lesson, LessonItem, QuestionKind};
use super::progress as prog;
/// Tutorial engine — the interactive lesson runner
///
/// Handles:
///   - Stepping through lesson items (text, questions, exercises)
///   - Reading user input from stdin
///   - Validating answers with fuzzy matching
///   - Providing hints and info panels
///   - Tracking progress and allowing skip/quit
///   - Running actual myconote-cli commands as demonstrations
use super::*;
use crate::utils::error::{MycoNoteError, Result};
use std::io::{self, BufRead, Write};

// ─────────────────────────────────────────────────────────────────────────────
// Lesson runner
// ─────────────────────────────────────────────────────────────────────────────

pub fn run_lesson(lesson_num: usize) -> Result<()> {
    let all = all_lessons();
    if lesson_num == 0 || lesson_num > all.len() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "Lesson {} does not exist. There are {} lessons available.",
            lesson_num,
            all.len()
        )));
    }

    let lesson = &all[lesson_num - 1];
    print_lesson_header(lesson, lesson_num, all.len());

    let mut correct_count = 0usize;
    let mut total_questions = 0usize;
    let mut skipped = 0usize;
    let stdin = io::stdin();

    for (i, item) in lesson.items.iter().enumerate() {
        match item {
            LessonItem::Text(text) => {
                print_text_block(text);
            }

            LessonItem::Info(title, body) => {
                print_info_panel(title, body);
            }

            LessonItem::CodeExample(description, code) => {
                print_code_example(description, code);
            }

            LessonItem::Demo(description, command) => {
                print_demo(description, command);
            }

            LessonItem::Question {
                prompt,
                kind,
                hint,
                explanation,
            } => {
                total_questions += 1;
                let result = ask_question(
                    prompt,
                    kind,
                    hint.as_deref(),
                    explanation.as_deref(),
                    &stdin,
                    i + 1,
                    lesson.items.len(),
                )?;

                match result {
                    AnswerResult::Correct => {
                        correct_count += 1;
                        print_correct(explanation.as_deref());
                    }
                    AnswerResult::Skipped => {
                        skipped += 1;
                        print_skipped(explanation.as_deref());
                    }
                    AnswerResult::Quit => {
                        print_quit_message(lesson_num, correct_count, total_questions);
                        // Save progress to resume later
                        prog::save_resume_point(lesson_num)?;
                        return Ok(());
                    }
                }
            }

            LessonItem::TryIt {
                instruction,
                command_template,
                hint,
            } => {
                print_try_it(instruction, command_template, hint.as_deref(), &stdin)?;
            }

            LessonItem::Checkpoint(msg) => {
                print_checkpoint(msg, i + 1, lesson.items.len());
            }
        }
    }

    // Lesson complete
    print_lesson_complete(
        lesson_num,
        correct_count,
        total_questions,
        skipped,
        all.len(),
    );
    prog::mark_lesson_complete(lesson_num)?;

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Answer handling
// ─────────────────────────────────────────────────────────────────────────────

enum AnswerResult {
    Correct,
    Skipped,
    Quit,
}

fn ask_question(
    prompt: &str,
    kind: &QuestionKind,
    hint: Option<&str>,
    _explanation: Option<&str>,
    stdin: &io::Stdin,
    step: usize,
    total: usize,
) -> Result<AnswerResult> {
    println!();
    println!(
        "  {}{}Question [{}/{}]{}",
        C_BOLD, C_YELLOW, step, total, C_RESET
    );
    println!("  {}", prompt);

    match kind {
        QuestionKind::FreeText {
            answer,
            accept_regex,
        } => ask_free_text(answer, accept_regex.as_deref(), hint, stdin),
        QuestionKind::MultipleChoice {
            choices,
            correct_index,
        } => ask_multiple_choice(choices, *correct_index, hint, stdin),
        QuestionKind::TrueFalse { answer } => ask_true_false(*answer, hint, stdin),
        QuestionKind::FillBlank { template, answer } => {
            println!("  {}Complete the command:{}", C_DIM, C_RESET);
            println!("  {}{}{}", C_CYAN, template, C_RESET);
            ask_free_text(answer, None, hint, stdin)
        }
        QuestionKind::OrderSteps {
            steps,
            correct_order,
        } => ask_order_steps(steps, correct_order, hint, stdin),
    }
}

fn ask_free_text(
    expected: &str,
    accept_regex: Option<&str>,
    hint: Option<&str>,
    stdin: &io::Stdin,
) -> Result<AnswerResult> {
    let mut attempts = 0;
    loop {
        print!("\n  {}> {}", C_GREEN, C_RESET);
        io::stdout().flush().ok();

        let mut input = String::new();
        stdin
            .lock()
            .read_line(&mut input)
            .map_err(MycoNoteError::Io)?;
        let trimmed = input.trim();

        match trimmed.to_lowercase().as_str() {
            "skip" | "s" => return Ok(AnswerResult::Skipped),
            "quit" | "bye" | "q" | "exit" => return Ok(AnswerResult::Quit),
            "hint" | "h" => {
                if let Some(h) = hint {
                    println!("  {}Hint:{} {}", C_YELLOW, C_RESET, h);
                } else {
                    println!("  {}No hint available for this question.{}", C_DIM, C_RESET);
                }
                continue;
            }
            "info" | "i" => {
                println!(
                    "  {}The expected answer relates to: {}{}",
                    C_DIM, expected, C_RESET
                );
                continue;
            }
            _ => {}
        }

        // Check answer
        let correct = if let Some(re_str) = accept_regex {
            if let Ok(re) = regex::Regex::new(re_str) {
                re.is_match(trimmed)
            } else {
                fuzzy_match(trimmed, expected)
            }
        } else {
            fuzzy_match(trimmed, expected)
        };

        if correct {
            return Ok(AnswerResult::Correct);
        }

        attempts += 1;
        if attempts >= 3 {
            println!(
                "  {}The answer is:{} {}{}{}",
                C_RED, C_RESET, C_BOLD, expected, C_RESET
            );
            println!("  {}Don't worry — this is how we learn!{}", C_DIM, C_RESET);
            return Ok(AnswerResult::Correct); // count as learned after reveal
        }

        println!(
            "  {}Not quite.{} Try again, or type {}hint{} / {}skip{}.",
            C_RED, C_RESET, C_YELLOW, C_RESET, C_YELLOW, C_RESET
        );
    }
}

fn ask_multiple_choice(
    choices: &[String],
    correct_idx: usize,
    hint: Option<&str>,
    stdin: &io::Stdin,
) -> Result<AnswerResult> {
    // Display choices
    for (i, choice) in choices.iter().enumerate() {
        let letter = (b'a' + i as u8) as char;
        println!("    {}{}{}) {}{}", C_BOLD, letter, C_RESET, choice, C_RESET);
    }

    let mut attempts = 0;
    loop {
        print!("\n  {}> {}", C_GREEN, C_RESET);
        io::stdout().flush().ok();

        let mut input = String::new();
        stdin
            .lock()
            .read_line(&mut input)
            .map_err(MycoNoteError::Io)?;
        let trimmed = input.trim().to_lowercase();

        match trimmed.as_str() {
            "skip" | "s" => return Ok(AnswerResult::Skipped),
            "quit" | "bye" | "q" | "exit" => return Ok(AnswerResult::Quit),
            "hint" | "h" => {
                if let Some(h) = hint {
                    println!("  {}Hint:{} {}", C_YELLOW, C_RESET, h);
                } else {
                    println!("  {}No hint available.{}", C_DIM, C_RESET);
                }
                continue;
            }
            _ => {}
        }

        // Accept letter (a, b, c...) or number (1, 2, 3...)
        let selected = if trimmed.len() == 1 {
            let ch = trimmed.chars().next().unwrap_or(' ');
            if ch.is_ascii_lowercase() {
                Some((ch as u8 - b'a') as usize)
            } else if ch.is_ascii_digit() {
                ch.to_digit(10).map(|d| d as usize - 1)
            } else {
                None
            }
        } else if let Ok(n) = trimmed.parse::<usize>() {
            Some(n - 1)
        } else {
            // Try fuzzy match against choice text
            choices.iter().position(|c| fuzzy_match(&trimmed, c))
        };

        if let Some(idx) = selected {
            if idx == correct_idx {
                return Ok(AnswerResult::Correct);
            }
        }

        attempts += 1;
        if attempts >= 3 {
            let letter = (b'a' + correct_idx as u8) as char;
            println!(
                "  {}The answer is:{} {}{}) {}{}",
                C_RED, C_RESET, C_BOLD, letter, choices[correct_idx], C_RESET
            );
            return Ok(AnswerResult::Correct);
        }

        println!(
            "  {}Not quite.{} Enter a letter (a/b/c/...) or type {}hint{}.",
            C_RED, C_RESET, C_YELLOW, C_RESET
        );
    }
}

fn ask_true_false(expected: bool, hint: Option<&str>, stdin: &io::Stdin) -> Result<AnswerResult> {
    println!("    {}a){} True", C_BOLD, C_RESET);
    println!("    {}b){} False", C_BOLD, C_RESET);

    loop {
        print!("\n  {}> {}", C_GREEN, C_RESET);
        io::stdout().flush().ok();

        let mut input = String::new();
        stdin
            .lock()
            .read_line(&mut input)
            .map_err(MycoNoteError::Io)?;
        let trimmed = input.trim().to_lowercase();

        match trimmed.as_str() {
            "skip" | "s" => return Ok(AnswerResult::Skipped),
            "quit" | "bye" | "q" | "exit" => return Ok(AnswerResult::Quit),
            "hint" | "h" => {
                if let Some(h) = hint {
                    println!("  {}Hint:{} {}", C_YELLOW, C_RESET, h);
                }
                continue;
            }
            _ => {}
        }

        let answer = match trimmed.as_str() {
            "true" | "t" | "a" | "1" | "yes" | "y" => Some(true),
            "false" | "f" | "b" | "2" | "no" | "n" => Some(false),
            _ => None,
        };

        if let Some(ans) = answer {
            if ans == expected {
                return Ok(AnswerResult::Correct);
            }
            let label = if expected { "True" } else { "False" };
            println!(
                "  {}The answer is:{} {}{}{}",
                C_RED, C_RESET, C_BOLD, label, C_RESET
            );
            return Ok(AnswerResult::Correct); // learned it
        }

        println!(
            "  Please answer {}true{} or {}false{}.",
            C_BOLD, C_RESET, C_BOLD, C_RESET
        );
    }
}

fn ask_order_steps(
    steps: &[String],
    correct_order: &[usize],
    hint: Option<&str>,
    stdin: &io::Stdin,
) -> Result<AnswerResult> {
    println!(
        "  {}Put these steps in the correct pipeline order:{}",
        C_DIM, C_RESET
    );
    for (i, step) in steps.iter().enumerate() {
        println!("    {}{}{}) {}", C_BOLD, i + 1, C_RESET, step);
    }
    println!(
        "  {}Enter the numbers in order (e.g. 3,1,2,4):{}",
        C_DIM, C_RESET
    );

    let mut attempts = 0;
    loop {
        print!("\n  {}> {}", C_GREEN, C_RESET);
        io::stdout().flush().ok();

        let mut input = String::new();
        stdin
            .lock()
            .read_line(&mut input)
            .map_err(MycoNoteError::Io)?;
        let trimmed = input.trim().to_lowercase();

        match trimmed.as_str() {
            "skip" | "s" => return Ok(AnswerResult::Skipped),
            "quit" | "bye" | "q" | "exit" => return Ok(AnswerResult::Quit),
            "hint" | "h" => {
                if let Some(h) = hint {
                    println!("  {}Hint:{} {}", C_YELLOW, C_RESET, h);
                }
                continue;
            }
            _ => {}
        }

        let user_order: Vec<usize> = trimmed
            .split(|c: char| c == ',' || c == ' ' || c == '-' || c == '>')
            .filter_map(|s| s.trim().parse::<usize>().ok())
            .collect();

        if user_order == correct_order.iter().map(|x| x + 1).collect::<Vec<_>>()
            || user_order == correct_order.to_vec()
        {
            return Ok(AnswerResult::Correct);
        }

        attempts += 1;
        if attempts >= 3 {
            let answer: Vec<String> = correct_order
                .iter()
                .map(|&i| format!("{}", i + 1))
                .collect();
            println!(
                "  {}The correct order is:{} {}{}{}",
                C_RED,
                C_RESET,
                C_BOLD,
                answer.join(", "),
                C_RESET
            );
            return Ok(AnswerResult::Correct);
        }

        println!("  {}Not quite.{} Try again.", C_RED, C_RESET);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Fuzzy matching
// ─────────────────────────────────────────────────────────────────────────────

fn fuzzy_match(input: &str, expected: &str) -> bool {
    let i = input.trim().to_lowercase();
    let e = expected.trim().to_lowercase();

    // Exact match
    if i == e {
        return true;
    }

    // Contains match (for long answers)
    if e.len() > 10 && i.contains(&e) {
        return true;
    }
    if i.len() > 5 && e.contains(&i) {
        return true;
    }

    // Strip common prefixes/suffixes
    let i_clean = i.replace("myconote-cli ", "").replace("myconote ", "");
    let e_clean = e.replace("myconote-cli ", "").replace("myconote ", "");
    if i_clean == e_clean {
        return true;
    }

    // Allow minor differences (1-2 chars for short strings)
    if e.len() <= 15 {
        let dist = levenshtein(&i, &e);
        if dist <= 1 {
            return true;
        }
    }

    false
}

fn levenshtein(a: &str, b: &str) -> usize {
    let a_chars: Vec<char> = a.chars().collect();
    let b_chars: Vec<char> = b.chars().collect();
    let m = a_chars.len();
    let n = b_chars.len();

    let mut prev = (0..=n).collect::<Vec<_>>();
    let mut curr = vec![0; n + 1];

    for i in 1..=m {
        curr[0] = i;
        for j in 1..=n {
            let cost = if a_chars[i - 1] == b_chars[j - 1] {
                0
            } else {
                1
            };
            curr[j] = (prev[j] + 1).min(curr[j - 1] + 1).min(prev[j - 1] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[n]
}

// ─────────────────────────────────────────────────────────────────────────────
// Pretty-printing helpers
// ─────────────────────────────────────────────────────────────────────────────

fn print_lesson_header(lesson: &Lesson, num: usize, total: usize) {
    println!();
    println!("  {}{}", C_DIM, "═".repeat(60));
    println!(
        "  {}{}Lesson {}/{}: {}{}",
        C_BOLD, C_CYAN, num, total, lesson.title, C_RESET
    );
    println!("  {}{}", C_DIM, "═".repeat(60));
    println!("  {}{}{}", C_DIM, lesson.description, C_RESET);
    println!();
    println!("  {}Type your answers. Special commands:{}", C_DIM, C_RESET);
    println!(
        "    {}skip{}  {}hint{}  {}info{}  {}quit{}",
        C_YELLOW, C_RESET, C_YELLOW, C_RESET, C_YELLOW, C_RESET, C_YELLOW, C_RESET
    );
    println!("  {}{}", C_DIM, "─".repeat(60));
    println!("{}", C_RESET);
}

fn print_text_block(text: &str) {
    println!();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            println!();
        } else {
            println!("  {}", trimmed);
        }
    }
}

fn print_info_panel(title: &str, body: &str) {
    println!();
    println!(
        "  {}{}  {} {}{}",
        C_BG_DARK, C_CYAN, title, C_RESET, C_RESET
    );
    for line in body.lines() {
        println!("  {}{}  {}{}", C_BG_DARK, C_DIM, line.trim(), C_RESET);
    }
    println!();
}

fn print_code_example(description: &str, code: &str) {
    println!();
    println!("  {}", description);
    println!();
    for line in code.lines() {
        println!("    {}$ {}{}", C_GREEN, line.trim(), C_RESET);
    }
    println!();
}

fn print_demo(description: &str, command: &str) {
    println!();
    println!("  {}{}{}", C_MAGENTA, description, C_RESET);
    println!("    {}$ {}{}", C_GREEN, command, C_RESET);
    println!();
    println!(
        "  {}(In a real session you would run this command now.){}",
        C_DIM, C_RESET
    );
    println!("  {}Press Enter to continue...{}", C_DIM, C_RESET);

    let _ = io::stdin().lock().lines().next();
}

fn print_try_it(
    instruction: &str,
    command_template: &str,
    hint: Option<&str>,
    stdin: &io::Stdin,
) -> Result<()> {
    println!();
    println!("  {}{}Try it yourself!{}", C_BOLD, C_MAGENTA, C_RESET);
    println!("  {}", instruction);
    println!();
    println!("  {}Template:{} {}", C_DIM, C_RESET, command_template);
    if let Some(h) = hint {
        println!("  {}Hint:{} {}", C_YELLOW, C_RESET, h);
    }

    print!("\n  {}> {}", C_GREEN, C_RESET);
    io::stdout().flush().ok();

    let mut input = String::new();
    stdin
        .lock()
        .read_line(&mut input)
        .map_err(MycoNoteError::Io)?;
    let trimmed = input.trim();

    if !trimmed.is_empty() && !matches!(trimmed, "skip" | "quit" | "bye") {
        println!("  {}Great! You typed:{} {}", C_GREEN, C_RESET, trimmed);
    }

    Ok(())
}

fn print_checkpoint(msg: &str, step: usize, total: usize) {
    let pct = (step as f64 / total as f64 * 100.0) as usize;
    let bar_len = 30;
    let filled = bar_len * pct / 100;
    let empty = bar_len - filled;
    let bar = format!("{}{}", "█".repeat(filled), "░".repeat(empty));

    println!();
    println!("  {}{}Checkpoint{} — {}", C_BOLD, C_BLUE, C_RESET, msg);
    println!("  [{}{}{}] {}%", C_CYAN, bar, C_RESET, pct);
    println!();
}

fn print_correct(explanation: Option<&str>) {
    println!("  {}✓ Correct!{}", C_GREEN, C_RESET);
    if let Some(exp) = explanation {
        println!("  {}{}{}", C_DIM, exp, C_RESET);
    }
}

fn print_skipped(explanation: Option<&str>) {
    println!("  {}→ Skipped.{}", C_YELLOW, C_RESET);
    if let Some(exp) = explanation {
        println!("  {}Answer: {}{}", C_DIM, exp, C_RESET);
    }
}

fn print_quit_message(lesson_num: usize, correct: usize, total: usize) {
    println!();
    println!("  {}Saving progress...{}", C_DIM, C_RESET);
    println!(
        "  You answered {}/{} questions before leaving.",
        correct, total
    );
    println!(
        "  Resume later: {}myconote-cli learn {}{}",
        C_CYAN, lesson_num, C_RESET
    );
    println!();
}

fn print_lesson_complete(
    lesson_num: usize,
    correct: usize,
    total: usize,
    skipped: usize,
    total_lessons: usize,
) {
    println!();
    println!("  {}{}", C_DIM, "═".repeat(60));
    println!(
        "  {}{}Lesson {} Complete!{}",
        C_BOLD, C_GREEN, lesson_num, C_RESET
    );
    println!("  {}{}", C_DIM, "═".repeat(60));
    println!();

    if total > 0 {
        let pct = correct as f64 / total as f64 * 100.0;
        println!("  Score: {}/{} ({:.0}%)", correct, total, pct);
        if skipped > 0 {
            println!("  Skipped: {}", skipped);
        }

        if pct >= 90.0 {
            println!("  {}Excellent work!{}", C_GREEN, C_RESET);
        } else if pct >= 70.0 {
            println!(
                "  {}Good job! Review the tricky parts.{}",
                C_YELLOW, C_RESET
            );
        } else {
            println!(
                "  {}Consider re-running this lesson to reinforce the concepts.{}",
                C_YELLOW, C_RESET
            );
        }
    }

    if lesson_num < total_lessons {
        println!();
        println!(
            "  Next lesson: {}myconote-cli learn {}{}",
            C_CYAN,
            lesson_num + 1,
            C_RESET
        );
    } else {
        println!();
        println!(
            "  {}You've completed all lessons! You're ready to annotate genomes.{}",
            C_GREEN, C_RESET
        );
    }
    println!();
}
