use crate::utils::error::{MycoNoteError, Result};
use serde::Serialize;
use std::path::PathBuf;

use super::commands::CommandRecommendation;
use super::context::StageContext;
use super::retrieval::Citation;
use super::rules::Finding;

/// A reproducibility bundle capturing the full explain call.
#[derive(Serialize)]
pub struct ExplainBundle {
    pub version: String,
    pub timestamp: String,
    pub stage: String,
    pub model: Option<String>,
    pub endpoint: Option<String>,
    pub prompt: Option<String>,
    pub context: StageContext,
    pub findings: Vec<Finding>,
    pub retrieved: Vec<Citation>,
    pub commands: Vec<CommandRecommendation>,
    pub response: Option<String>,
    pub commands_kept: usize,
    pub commands_removed: usize,
}

impl ExplainBundle {
    /// Write the bundle to `./explain_<stage>_<timestamp>/`.
    pub fn write(&self) -> Result<PathBuf> {
        let dir_name = format!(
            "explain_{}_{}",
            self.stage,
            self.timestamp
                .replace([':', ' ', 'T'], "_")
                .replace('+', "")
        );
        let dir = PathBuf::from(&dir_name);

        std::fs::create_dir_all(&dir).map_err(|e| {
            MycoNoteError::ChatContext(format!("cannot create bundle dir {}: {}", dir.display(), e))
        })?;

        // Write the full bundle as JSON
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| MycoNoteError::ChatContext(format!("cannot serialize bundle: {}", e)))?;
        std::fs::write(dir.join("bundle.json"), &json)
            .map_err(|e| MycoNoteError::ChatContext(format!("cannot write bundle.json: {}", e)))?;

        // Write the prompt separately for easy inspection
        if let Some(ref prompt) = self.prompt {
            let _ = std::fs::write(dir.join("prompt.txt"), prompt);
        }

        // Write the response separately
        if let Some(ref response) = self.response {
            let _ = std::fs::write(dir.join("response.txt"), response);
        }

        // Write findings as a human-readable text file
        let mut findings_text = String::new();
        for f in &self.findings {
            findings_text.push_str(&format!(
                "[{}] {}: {}\n  Evidence: {}\n",
                f.severity, f.rule_id, f.message, f.evidence
            ));
            if let Some(ref cite) = f.citation {
                findings_text.push_str(&format!("  Citation: {}\n", cite));
            }
            findings_text.push('\n');
        }
        let _ = std::fs::write(dir.join("findings.txt"), &findings_text);

        Ok(dir)
    }
}

/// Generate an ISO 8601 timestamp.
pub fn timestamp() -> String {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    // Simple UTC timestamp without chrono dependency
    let secs_per_day = 86400u64;
    let days = now / secs_per_day;
    let time_of_day = now % secs_per_day;

    let hours = time_of_day / 3600;
    let minutes = (time_of_day % 3600) / 60;
    let seconds = time_of_day % 60;

    // Approximate date calculation (good enough for filenames)
    let mut year = 1970u64;
    let mut remaining_days = days;

    loop {
        let days_in_year = if year % 4 == 0 && (year % 100 != 0 || year % 400 == 0) {
            366
        } else {
            365
        };
        if remaining_days < days_in_year {
            break;
        }
        remaining_days -= days_in_year;
        year += 1;
    }

    let month_days = [
        31,
        28 + if year % 4 == 0 && (year % 100 != 0 || year % 400 == 0) {
            1
        } else {
            0
        },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    let mut month = 1u64;
    for md in &month_days {
        if remaining_days < *md as u64 {
            break;
        }
        remaining_days -= *md as u64;
        month += 1;
    }
    let day = remaining_days + 1;

    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        year, month, day, hours, minutes, seconds
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn timestamp_format() {
        let ts = timestamp();
        assert!(ts.contains('T'));
        assert!(ts.ends_with('Z'));
        assert!(ts.len() == 20);
    }

    #[test]
    fn bundle_write_creates_dir() {
        let dir = tempfile::tempdir().unwrap();
        let old_cwd = std::env::current_dir().unwrap();
        std::env::set_current_dir(dir.path()).unwrap();

        let bundle = ExplainBundle {
            version: "0.1.0".to_string(),
            timestamp: "2024-01-01T00_00_00Z".to_string(),
            stage: "predict".to_string(),
            model: Some("llama3.1".to_string()),
            endpoint: Some("http://localhost:11434".to_string()),
            prompt: Some("test prompt".to_string()),
            context: StageContext {
                stage: "predict".to_string(),
                dir: "/tmp".to_string(),
                artifacts: vec![],
                stats: None,
                notes: vec![],
            },
            findings: vec![],
            retrieved: vec![],
            commands: vec![],
            response: Some("test response".to_string()),
            commands_kept: 0,
            commands_removed: 0,
        };

        let result = bundle.write();
        std::env::set_current_dir(old_cwd).unwrap();

        assert!(result.is_ok());
        let path = result.unwrap();
        assert!(dir.path().join(&path).join("bundle.json").exists());
        assert!(dir.path().join(&path).join("prompt.txt").exists());
        assert!(dir.path().join(&path).join("response.txt").exists());
    }
}
