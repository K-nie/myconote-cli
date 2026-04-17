/// TTY-aware progress dashboard for batch genome annotation.
///
/// - Interactive terminal: live multi-progress bars via `indicatif`
/// - HPC batch job (no TTY): timestamped log lines, one per event
///
/// Auto-detects the environment — no configuration needed.
use std::collections::HashMap;
use std::io::IsTerminal;
use std::time::Instant;

use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

use super::state::BatchState;

// ANSI helpers
const C_GREEN: &str = "\x1b[32m";
const C_RED: &str = "\x1b[31m";
const C_CYAN: &str = "\x1b[36m";
const C_DIM: &str = "\x1b[2m";
const C_BOLD: &str = "\x1b[1m";
const C_YELLOW: &str = "\x1b[33m";
const C_RESET: &str = "\x1b[0m";

const STAGES: &[&str] = &[
    "sort", "mask", "train", "predict", "update", "annotate", "submit",
];

/// Dashboard abstraction — interactive or log-based.
pub struct Dashboard {
    interactive: bool,
    multi: Option<MultiProgress>,
    bars: HashMap<String, ProgressBar>,
    overall: Option<ProgressBar>,
    start: Instant,
    total_genomes: usize,
}

impl Dashboard {
    /// Create a new dashboard. Auto-detects TTY.
    pub fn new(total_genomes: usize) -> Self {
        let interactive = std::io::stdout().is_terminal();

        let (multi, overall) = if interactive {
            let mp = MultiProgress::new();

            // Header bar
            let header = mp.add(ProgressBar::new(total_genomes as u64));
            header.set_style(
                ProgressStyle::with_template(
                    "  {msg}\n  [{bar:40.cyan/blue}] {pos}/{len} genomes ({eta})",
                )
                .unwrap_or_else(|_| ProgressStyle::default_bar())
                .progress_chars("\u{2588}\u{2593}\u{2591}"),
            );
            header.set_message(format!(
                "{}{}myconote batch{} \u{2014} 0 / {} genomes complete",
                C_BOLD, C_CYAN, C_RESET, total_genomes
            ));
            (Some(mp), Some(header))
        } else {
            (None, None)
        };

        Self {
            interactive,
            multi,
            bars: HashMap::new(),
            overall,
            start: Instant::now(),
            total_genomes,
        }
    }

    /// Register a genome in the dashboard.
    pub fn add_genome(&mut self, name: &str) {
        if self.interactive {
            if let Some(ref mp) = self.multi {
                let pb = mp.add(ProgressBar::new_spinner());
                pb.set_style(
                    ProgressStyle::with_template("  {spinner:.cyan} {msg}")
                        .unwrap_or_else(|_| ProgressStyle::default_spinner())
                        .tick_strings(&["\u{00b7}", "\u{00b7}", "\u{00b7}", "\u{00b7}"]),
                );
                pb.set_message(format!("{:<30} queued", name));
                self.bars.insert(name.to_string(), pb);
            }
        }
    }

    /// Report that a genome is starting a stage.
    pub fn stage_start(&self, genome: &str, stage: &str) {
        if self.interactive {
            if let Some(pb) = self.bars.get(genome) {
                pb.set_style(
                    ProgressStyle::with_template("  {spinner:.cyan} {msg}")
                        .unwrap_or_else(|_| ProgressStyle::default_spinner())
                        .tick_strings(&[
                            "\u{2800}\u{2801}",
                            "\u{2800}\u{2809}",
                            "\u{2800}\u{2839}",
                            "\u{2800}\u{2838}",
                            "\u{2800}\u{283c}",
                            "\u{2800}\u{2834}",
                            "\u{2800}\u{2826}",
                            "\u{2800}\u{2827}",
                            "\u{2800}\u{2807}",
                            "\u{2800}\u{280f}",
                        ]),
                );
                pb.enable_steady_tick(std::time::Duration::from_millis(80));
                let progress = render_stage_progress(genome, stage, &[]);
                pb.set_message(progress);
            }
        } else {
            log_line(genome, stage, "started");
        }
    }

    /// Report that a stage completed for a genome.
    pub fn stage_done(&self, genome: &str, stage: &str, completed: &[&str]) {
        if self.interactive {
            if let Some(pb) = self.bars.get(genome) {
                let progress = render_stage_progress(genome, stage, completed);
                pb.set_message(progress);
            }
        } else {
            log_line(genome, stage, "done");
        }
    }

    /// Report that a genome completed all stages.
    pub fn genome_done(&self, genome: &str, elapsed_s: f64) {
        if self.interactive {
            if let Some(pb) = self.bars.get(genome) {
                pb.set_style(
                    ProgressStyle::with_template("  \u{2713} {msg}")
                        .unwrap_or_else(|_| ProgressStyle::default_spinner()),
                );
                let time_str = format_duration(elapsed_s);
                let stages_str = STAGES
                    .iter()
                    .map(|s| format!("{}{}{}", C_DIM, s, C_RESET))
                    .collect::<Vec<_>>()
                    .join(" \u{2192} ");
                pb.finish_with_message(format!(
                    "{:<30} {}  {}{}{}",
                    genome, stages_str, C_DIM, time_str, C_RESET
                ));
            }
            if let Some(ref overall) = self.overall {
                overall.inc(1);
                let (done, _, _) = self.current_counts();
                overall.set_message(format!(
                    "{}{}myconote batch{} \u{2014} {} / {} genomes complete",
                    C_BOLD, C_CYAN, C_RESET, done, self.total_genomes
                ));
            }
        } else {
            log_line(genome, "ALL", &format!("done ({:.0}s)", elapsed_s));
        }
    }

    /// Report that a genome failed at a specific stage.
    pub fn genome_failed(&self, genome: &str, stage: &str, error: &str) {
        if self.interactive {
            if let Some(pb) = self.bars.get(genome) {
                pb.set_style(
                    ProgressStyle::with_template("  \u{2717} {msg}")
                        .unwrap_or_else(|_| ProgressStyle::default_spinner()),
                );
                pb.finish_with_message(format!(
                    "{}{:<30}{} {}{} failed{}: {}",
                    C_RED, genome, C_RESET, C_RED, stage, C_RESET, error
                ));
            }
            if let Some(ref overall) = self.overall {
                overall.inc(1);
            }
        } else {
            log_line(genome, stage, &format!("FAILED \u{2014} {}", error));
        }
    }

    /// Print the final summary table after all genomes are processed.
    pub fn print_summary(&self, state: &BatchState) {
        let (done, failed, _pending) = state.summary();
        let elapsed = self.start.elapsed().as_secs_f64();

        println!();
        if failed == 0 {
            println!(
                "  {}\u{2500}\u{2500} Batch complete: {}{}{} succeeded {}({}){}",
                C_DIM,
                C_GREEN,
                done,
                C_RESET,
                C_DIM,
                format_duration(elapsed),
                C_RESET
            );
        } else {
            println!(
                "  {}\u{2500}\u{2500} Batch complete: {}{}{} succeeded, {}{} failed{} {}({}){}",
                C_DIM,
                C_GREEN,
                done,
                C_RESET,
                C_RED,
                failed,
                C_RESET,
                C_DIM,
                format_duration(elapsed),
                C_RESET
            );
        }
        println!();

        // Table header
        println!(
            "  {}{:<30} {:>8} {:>9} {:>8}  {}{}",
            C_BOLD, "Genome", "Genes", "Masked%", "Time", "Status", C_RESET
        );
        println!("  {}{}{}", C_DIM, "\u{2500}".repeat(75), C_RESET);

        let mut names: Vec<&String> = state.genomes.keys().collect();
        names.sort();

        for name in &names {
            let gs = &state.genomes[*name];
            let time_str = format_duration(gs.elapsed_s);

            if gs.failed {
                let err = gs.error.as_deref().unwrap_or("unknown error");
                println!(
                    "  {:<30} {:>8} {:>9} {:>8}  {}\u{2717} {}{}",
                    name, "\u{2014}", "\u{2014}", time_str, C_RED, err, C_RESET
                );
            } else {
                let genes = gs
                    .outputs
                    .get("gene_count")
                    .map(|s| format_number(s))
                    .unwrap_or_else(|| "\u{2014}".to_string());
                let masked = gs
                    .outputs
                    .get("masked_pct")
                    .map(|s| format!("{}%", s))
                    .unwrap_or_else(|| "\u{2014}".to_string());
                println!(
                    "  {:<30} {:>8} {:>9} {:>8}  {}\u{2713}{}",
                    name, genes, masked, time_str, C_GREEN, C_RESET
                );
            }
        }

        println!("  {}{}{}", C_DIM, "\u{2500}".repeat(75), C_RESET);
        println!();
        println!("  Logs: {}{}/{}", C_CYAN, state.batch_dir, C_RESET);

        if done >= 2 {
            println!(
                "  Next: {}myconote-cli compare {}/**/annotate_out/*.gff3{}",
                C_YELLOW, state.batch_dir, C_RESET
            );
        }
        println!();
    }

    fn current_counts(&self) -> (usize, usize, usize) {
        // Approximate from progress bar state
        if let Some(ref overall) = self.overall {
            let done = overall.position() as usize;
            (done, 0, self.total_genomes - done)
        } else {
            (0, 0, self.total_genomes)
        }
    }
}

/// Render the inline stage progress for interactive mode.
fn render_stage_progress(genome: &str, current: &str, _completed: &[&str]) -> String {
    let mut parts = Vec::new();

    for &stage in STAGES {
        if stage == current {
            parts.push(format!("{}{}\u{25b6} {}{}", C_BOLD, C_CYAN, stage, C_RESET));
        } else {
            // Both completed (before current) and pending (after current) are dimmed.
            parts.push(format!("{}{}{}", C_DIM, stage, C_RESET));
        }
    }

    format!("{:<30} {}", genome, parts.join(" \u{2192} "))
}

/// Print a structured log line for HPC/non-interactive mode.
fn log_line(genome: &str, stage: &str, status: &str) {
    // Use same chrono_now from state.rs pattern
    let ts = timestamp_now();
    println!("[{}] {:<30} | {:<10} | {}", ts, genome, stage, status);
}

/// Simple timestamp for log lines.
fn timestamp_now() -> String {
    use std::time::SystemTime;
    let dur = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = dur.as_secs();
    let time_of_day = secs % 86400;
    let hours = time_of_day / 3600;
    let mins = (time_of_day % 3600) / 60;
    let s = time_of_day % 60;

    let days = secs / 86400;
    let mut y = 1970u64;
    let mut remaining = days;
    loop {
        let diy = if (y % 4 == 0 && y % 100 != 0) || y % 400 == 0 {
            366
        } else {
            365
        };
        if remaining < diy {
            break;
        }
        remaining -= diy;
        y += 1;
    }
    let leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
    let month_days: &[u64] = if leap {
        &[31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        &[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };
    let mut m = 1u64;
    for &md in month_days {
        if remaining < md {
            break;
        }
        remaining -= md;
        m += 1;
    }
    let d = remaining + 1;
    format!(
        "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
        y, m, d, hours, mins, s
    )
}

/// Format seconds into human-readable duration.
pub fn format_duration(secs: f64) -> String {
    if secs < 60.0 {
        format!("{:.0}s", secs)
    } else if secs < 3600.0 {
        format!("{:.0}m{:.0}s", secs / 60.0, secs % 60.0)
    } else {
        let h = (secs / 3600.0).floor();
        let m = ((secs % 3600.0) / 60.0).floor();
        format!("{:.0}h{:.0}m", h, m)
    }
}

/// Format a numeric string with commas: "12345" -> "12,345"
fn format_number(s: &str) -> String {
    let n: i64 = s.parse().unwrap_or(0);
    if n == 0 {
        return s.to_string();
    }
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_duration_seconds() {
        assert_eq!(format_duration(45.0), "45s");
    }

    #[test]
    fn test_format_duration_minutes() {
        assert_eq!(format_duration(125.0), "2m5s");
    }

    #[test]
    fn test_format_duration_hours() {
        assert_eq!(format_duration(7380.0), "2h3m");
    }

    #[test]
    fn test_format_number() {
        assert_eq!(format_number("12345"), "12,345");
        assert_eq!(format_number("999"), "999");
        assert_eq!(format_number("1000000"), "1,000,000");
    }

    #[test]
    fn test_timestamp_now_format() {
        let ts = timestamp_now();
        assert_eq!(ts.len(), 19);
        assert_eq!(&ts[4..5], "-");
    }
}
