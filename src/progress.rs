/// Shared progress bar utilities using `indicatif`
///
/// Provides consistent, pre-styled progress indicators used across
/// the mask, predict, and annotate pipelines.

use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
use std::time::Duration;

// ─────────────────────────────────────────────────────────────────────────────
// Style templates
// ─────────────────────────────────────────────────────────────────────────────

const SPINNER_TEMPLATE: &str = "  {spinner:.cyan} {msg}";
const BAR_TEMPLATE:     &str = "  {msg}\n  [{bar:40.cyan/blue}] {pos}/{len} ({eta})";
#[allow(dead_code)]
const BYTES_TEMPLATE:   &str = "  {msg}\n  [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({bytes_per_sec}, {eta})";

// ─────────────────────────────────────────────────────────────────────────────
// Spinner — for tasks of unknown duration
// ─────────────────────────────────────────────────────────────────────────────

/// Create a spinner for a long-running step with unknown duration.
/// Call `.finish_with_message("done")` when complete.
pub fn spinner(msg: impl Into<String>) -> ProgressBar {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template(SPINNER_TEMPLATE)
            .unwrap_or_else(|_| ProgressStyle::default_spinner())
            .tick_strings(&["⠋","⠙","⠹","⠸","⠼","⠴","⠦","⠧","⠇","⠏"])
    );
    pb.set_message(msg.into());
    pb.enable_steady_tick(Duration::from_millis(80));
    pb
}

/// Finish a spinner with a green checkmark.
pub fn finish_spinner(pb: &ProgressBar, msg: impl Into<String>) {
    pb.set_style(
        ProgressStyle::with_template("  ✓ {msg}")
            .unwrap_or_else(|_| ProgressStyle::default_spinner())
    );
    pb.finish_with_message(msg.into());
}

/// Finish a spinner indicating a warning (non-fatal skip).
pub fn warn_spinner(pb: &ProgressBar, msg: impl Into<String>) {
    pb.set_style(
        ProgressStyle::with_template("  ⚠ {msg}")
            .unwrap_or_else(|_| ProgressStyle::default_spinner())
    );
    pb.finish_with_message(msg.into());
}

// ─────────────────────────────────────────────────────────────────────────────
// Count bar — for iterating over a known number of items
// ─────────────────────────────────────────────────────────────────────────────

/// Create a counting progress bar (items processed / total).
pub fn count_bar(total: u64, msg: impl Into<String>) -> ProgressBar {
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::with_template(BAR_TEMPLATE)
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars("█▓░")
    );
    pb.set_message(msg.into());
    pb
}

// ─────────────────────────────────────────────────────────────────────────────
// Multi-progress — show several bars at once
// ─────────────────────────────────────────────────────────────────────────────

pub fn multi() -> MultiProgress {
    MultiProgress::new()
}

/// Add a spinner to a MultiProgress group.
pub fn add_spinner(multi: &MultiProgress, msg: impl Into<String>) -> ProgressBar {
    let pb = multi.add(ProgressBar::new_spinner());
    pb.set_style(
        ProgressStyle::with_template(SPINNER_TEMPLATE)
            .unwrap_or_else(|_| ProgressStyle::default_spinner())
            .tick_strings(&["⠋","⠙","⠹","⠸","⠼","⠴","⠦","⠧","⠇","⠏"])
    );
    pb.set_message(msg.into());
    pb.enable_steady_tick(Duration::from_millis(80));
    pb
}

// ─────────────────────────────────────────────────────────────────────────────
// Step label — for numbered pipeline steps without a progress count
// ─────────────────────────────────────────────────────────────────────────────

/// Print a numbered pipeline step header.
pub fn step(n: usize, total: usize, msg: &str) {
    println!("  [{}/{}] {}", n, total, msg);
}
