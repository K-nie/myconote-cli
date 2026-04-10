/// Tutorial progress tracker
///
/// Saves lesson completion status to ~/.myconote/learn_progress.json
/// so users can resume where they left off across sessions.
use crate::utils::error::{MycoNoteError, Result};
use std::path::PathBuf;

fn progress_dir() -> PathBuf {
    std::env::var("HOME")
        .or_else(|_| std::env::var("USERPROFILE"))
        .map(PathBuf::from)
        .unwrap_or_else(|_| std::env::temp_dir())
        .join(".myconote")
}

fn progress_file() -> PathBuf {
    progress_dir().join("learn_progress.json")
}

#[derive(serde::Serialize, serde::Deserialize, Default)]
struct ProgressData {
    completed_lessons: Vec<usize>,
    resume_lesson: Option<usize>,
}

pub fn load_progress() -> Result<Vec<usize>> {
    let path = progress_file();
    if !path.exists() {
        return Ok(Vec::new());
    }
    let content = std::fs::read_to_string(&path).map_err(MycoNoteError::Io)?;
    let data: ProgressData = serde_json::from_str(&content).unwrap_or_default();
    Ok(data.completed_lessons)
}

fn load_full() -> ProgressData {
    let path = progress_file();
    if !path.exists() {
        return ProgressData::default();
    }
    std::fs::read_to_string(&path)
        .ok()
        .and_then(|c| serde_json::from_str(&c).ok())
        .unwrap_or_default()
}

fn save_full(data: &ProgressData) -> Result<()> {
    let dir = progress_dir();
    std::fs::create_dir_all(&dir).map_err(MycoNoteError::Io)?;
    let json = serde_json::to_string_pretty(data)
        .map_err(|e| MycoNoteError::InvalidFormat(format!("JSON: {}", e)))?;
    std::fs::write(progress_file(), json).map_err(MycoNoteError::Io)?;
    Ok(())
}

pub fn mark_lesson_complete(lesson_num: usize) -> Result<()> {
    let mut data = load_full();
    if !data.completed_lessons.contains(&lesson_num) {
        data.completed_lessons.push(lesson_num);
        data.completed_lessons.sort();
    }
    // Auto-advance resume point
    data.resume_lesson = Some(lesson_num + 1);
    save_full(&data)
}

pub fn save_resume_point(lesson_num: usize) -> Result<()> {
    let mut data = load_full();
    data.resume_lesson = Some(lesson_num);
    save_full(&data)
}

pub fn get_resume_point() -> Result<usize> {
    let data = load_full();
    Ok(data.resume_lesson.unwrap_or(1))
}

pub fn reset_progress() -> Result<()> {
    let path = progress_file();
    if path.exists() {
        std::fs::remove_file(&path).map_err(MycoNoteError::Io)?;
    }
    Ok(())
}
