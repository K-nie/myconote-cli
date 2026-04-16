use crate::utils::error::{MycoNoteError, Result};
use serde::Serialize;
use std::path::PathBuf;

#[derive(Serialize)]
pub struct HistoryEntry {
    pub timestamp: String,
    pub stage: String,
    pub dir: String,
    pub findings_count: usize,
    pub llm_used: bool,
    pub model: Option<String>,
}

/// Append an entry to the history JSONL file.
pub fn append(entry: &HistoryEntry) -> Result<()> {
    let dir = history_dir()?;
    let date = &entry.timestamp[..10]; // YYYY-MM-DD
    let file = dir.join(format!("{}.jsonl", date));

    let json = serde_json::to_string(entry)
        .map_err(|e| MycoNoteError::ChatConfig(format!("cannot serialize history: {}", e)))?;

    use std::io::Write;
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&file)
        .map_err(|e| MycoNoteError::ChatConfig(format!("cannot open {}: {}", file.display(), e)))?;

    writeln!(f, "{}", json)
        .map_err(|e| MycoNoteError::ChatConfig(format!("cannot write to {}: {}", file.display(), e)))?;

    Ok(())
}

fn history_dir() -> Result<PathBuf> {
    let base = super::config::myconote_dir()?;
    let dir = base.join("chat_history");
    if !dir.exists() {
        std::fs::create_dir_all(&dir)
            .map_err(|e| MycoNoteError::ChatConfig(format!("cannot create {}: {}", dir.display(), e)))?;
    }
    Ok(dir)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn history_entry_serializes() {
        let entry = HistoryEntry {
            timestamp: "2024-01-15T10:30:00Z".to_string(),
            stage: "predict".to_string(),
            dir: "/tmp/test".to_string(),
            findings_count: 3,
            llm_used: true,
            model: Some("llama3.1".to_string()),
        };
        let json = serde_json::to_string(&entry).unwrap();
        assert!(json.contains("predict"));
        assert!(json.contains("llama3.1"));
    }
}
