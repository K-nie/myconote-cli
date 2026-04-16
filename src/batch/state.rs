/// Persistent state for batch runs — enables resume after interruption.
///
/// Writes `status.json` after every stage completion so that `batch --resume`
/// can skip completed work and pick up where it left off.
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::utils::error::{MycoNoteError, Result};

/// Overall batch run state.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct BatchState {
    /// When the batch run started (ISO 8601).
    pub started_at: String,
    /// Batch output directory.
    pub batch_dir: String,
    /// Per-genome state keyed by genome filename.
    pub genomes: HashMap<String, GenomeState>,
    /// Global settings used for this run.
    pub settings: BatchSettings,
}

/// Settings snapshot — stored so resume uses the same config.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct BatchSettings {
    pub kingdom: String,
    pub threads: usize,
    pub max_parallel: usize,
    pub stages: Vec<String>,
}

/// Per-genome progress.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct GenomeState {
    /// Stages that completed successfully, in order.
    pub completed_stages: Vec<String>,
    /// The stage currently running (if interrupted mid-stage).
    pub current_stage: Option<String>,
    /// Whether this genome failed.
    pub failed: bool,
    /// Error message if failed.
    pub error: Option<String>,
    /// Elapsed time in seconds.
    pub elapsed_s: f64,
    /// Key output paths produced so far.
    pub outputs: HashMap<String, String>,
}

impl GenomeState {
    pub fn new() -> Self {
        Self {
            completed_stages: Vec::new(),
            current_stage: None,
            failed: false,
            error: None,
            elapsed_s: 0.0,
            outputs: HashMap::new(),
        }
    }

    /// Check whether a given stage is already done.
    pub fn is_stage_done(&self, stage: &str) -> bool {
        self.completed_stages.iter().any(|s| s == stage)
    }

    /// Mark a stage as completed and record its primary output path.
    pub fn complete_stage(&mut self, stage: &str, output_path: Option<&str>) {
        self.current_stage = None;
        if !self.is_stage_done(stage) {
            self.completed_stages.push(stage.to_string());
        }
        if let Some(p) = output_path {
            self.outputs.insert(stage.to_string(), p.to_string());
        }
    }

    /// Mark a stage as started.
    pub fn start_stage(&mut self, stage: &str) {
        self.current_stage = Some(stage.to_string());
    }

    /// Mark this genome as failed.
    pub fn mark_failed(&mut self, stage: &str, error: &str) {
        self.current_stage = None;
        self.failed = true;
        self.error = Some(format!("{}: {}", stage, error));
    }
}

impl BatchState {
    /// Create a new state for a fresh batch run.
    pub fn new(batch_dir: &Path, settings: BatchSettings) -> Self {
        Self {
            started_at: chrono_now(),
            batch_dir: batch_dir.display().to_string(),
            genomes: HashMap::new(),
            settings,
        }
    }

    /// Load state from a `status.json` file.
    pub fn load(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| MycoNoteError::BatchError(format!("cannot read state file: {}", e)))?;
        let state: Self = serde_json::from_str(&text)
            .map_err(|e| MycoNoteError::BatchError(format!("invalid state JSON: {}", e)))?;
        Ok(state)
    }

    /// Persist state to `<batch_dir>/status.json`.
    pub fn save(&self) -> Result<()> {
        let path = PathBuf::from(&self.batch_dir).join("status.json");
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| MycoNoteError::BatchError(format!("cannot serialize state: {}", e)))?;
        std::fs::write(&path, json).map_err(|e| {
            MycoNoteError::BatchError(format!("cannot write {}: {}", path.display(), e))
        })?;
        Ok(())
    }

    /// Get or create state for a genome.
    pub fn genome_mut(&mut self, name: &str) -> &mut GenomeState {
        self.genomes
            .entry(name.to_string())
            .or_insert_with(GenomeState::new)
    }

    /// Count genomes by status.
    pub fn summary(&self) -> (usize, usize, usize) {
        let done = self
            .genomes
            .values()
            .filter(|g| !g.failed && !g.completed_stages.is_empty())
            .count();
        let failed = self.genomes.values().filter(|g| g.failed).count();
        let pending = self.genomes.len() - done - failed;
        (done, failed, pending)
    }
}

/// Simple ISO 8601 timestamp without pulling in chrono crate.
fn chrono_now() -> String {
    use std::time::SystemTime;
    let dur = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = dur.as_secs();
    // Format as YYYY-MM-DDTHH:MM:SSZ (approximate — good enough for batch IDs)
    let days = secs / 86400;
    let time_of_day = secs % 86400;
    let hours = time_of_day / 3600;
    let mins = (time_of_day % 3600) / 60;
    let s = time_of_day % 60;

    // Days since epoch to Y/M/D (simplified leap year calculation)
    let mut y = 1970u64;
    let mut remaining = days;
    loop {
        let days_in_year = if is_leap(y) { 366 } else { 365 };
        if remaining < days_in_year {
            break;
        }
        remaining -= days_in_year;
        y += 1;
    }
    let month_days: &[u64] = if is_leap(y) {
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
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        y, m, d, hours, mins, s
    )
}

fn is_leap(y: u64) -> bool {
    (y % 4 == 0 && y % 100 != 0) || y % 400 == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genome_state_lifecycle() {
        let mut gs = GenomeState::new();
        assert!(!gs.is_stage_done("sort"));

        gs.start_stage("sort");
        assert_eq!(gs.current_stage.as_deref(), Some("sort"));

        gs.complete_stage("sort", Some("/out/sorted.fa"));
        assert!(gs.is_stage_done("sort"));
        assert!(gs.current_stage.is_none());
        assert_eq!(gs.outputs.get("sort").unwrap(), "/out/sorted.fa");
    }

    #[test]
    fn genome_state_failure() {
        let mut gs = GenomeState::new();
        gs.start_stage("predict");
        gs.mark_failed("predict", "Augustus crashed");
        assert!(gs.failed);
        assert!(gs.error.as_ref().unwrap().contains("Augustus crashed"));
        assert!(gs.current_stage.is_none());
    }

    #[test]
    fn batch_state_summary() {
        let settings = BatchSettings {
            kingdom: "fungi".to_string(),
            threads: 4,
            max_parallel: 2,
            stages: vec!["sort".to_string(), "mask".to_string()],
        };
        let mut state = BatchState::new(Path::new("/tmp/batch"), settings);

        state.genome_mut("a.fa").complete_stage("sort", None);
        state.genome_mut("a.fa").complete_stage("mask", None);
        state.genome_mut("b.fa").mark_failed("sort", "bad fasta");
        state.genome_mut("c.fa"); // pending — no stages done

        let (done, failed, pending) = state.summary();
        assert_eq!(done, 1);
        assert_eq!(failed, 1);
        assert_eq!(pending, 1);
    }

    #[test]
    fn state_roundtrip() {
        let settings = BatchSettings {
            kingdom: "fungi".to_string(),
            threads: 4,
            max_parallel: 2,
            stages: vec!["sort".to_string()],
        };
        let mut state = BatchState::new(Path::new("/tmp/test_batch"), settings);
        state
            .genome_mut("genome.fa")
            .complete_stage("sort", Some("/out/sorted.fa"));

        let json = serde_json::to_string_pretty(&state).unwrap();
        let loaded: BatchState = serde_json::from_str(&json).unwrap();
        assert_eq!(loaded.genomes.len(), 1);
        assert!(loaded.genomes["genome.fa"].is_stage_done("sort"));
    }

    #[test]
    fn chrono_now_format() {
        let ts = chrono_now();
        // Should look like 20XX-MM-DDTHH:MM:SSZ
        assert!(ts.ends_with('Z'));
        assert_eq!(ts.len(), 20);
        assert_eq!(&ts[4..5], "-");
        assert_eq!(&ts[7..8], "-");
        assert_eq!(&ts[10..11], "T");
    }
}
