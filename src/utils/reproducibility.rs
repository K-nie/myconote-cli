/// Reproducibility report generator
///
/// Generates a comprehensive workflow summary documenting:
///   - myconote-cli version and build info
///   - All external tool versions used in the pipeline
///   - Database versions and download dates
///   - Input files and their checksums
///   - Command-line parameters used
///   - Runtime timestamps and durations
///   - System information (OS, CPU, memory)
///
/// Output: JSON and human-readable text reports.
use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

// ─────────────────────────────────────────────────────────────────────────────
// Report structure
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, serde::Serialize)]
pub struct ReproducibilityReport {
    /// myconote-cli version
    pub myconote_version: String,
    /// Timestamp when the pipeline was started (ISO 8601)
    pub timestamp: String,
    /// Pipeline steps executed
    pub pipeline_steps: Vec<String>,
    /// External tool versions
    pub tool_versions: HashMap<String, String>,
    /// Database versions and paths
    pub databases: Vec<DatabaseInfo>,
    /// Input files with checksums
    pub input_files: Vec<FileInfo>,
    /// Output files
    pub output_files: Vec<String>,
    /// Command-line parameters
    pub parameters: HashMap<String, String>,
    /// System information
    pub system_info: SystemInfo,
    /// Runtime duration in seconds
    pub runtime_seconds: Option<f64>,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct DatabaseInfo {
    pub name: String,
    pub path: String,
    pub version: Option<String>,
    pub download_date: Option<String>,
    pub size_bytes: Option<u64>,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct FileInfo {
    pub path: String,
    pub size_bytes: u64,
    pub md5: Option<String>,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct SystemInfo {
    pub os: String,
    pub arch: String,
    pub hostname: String,
    pub cpus: usize,
}

impl Default for ReproducibilityReport {
    fn default() -> Self {
        Self {
            myconote_version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp: current_iso8601(),
            pipeline_steps: Vec::new(),
            tool_versions: HashMap::new(),
            databases: Vec::new(),
            input_files: Vec::new(),
            output_files: Vec::new(),
            parameters: HashMap::new(),
            system_info: detect_system_info(),
            runtime_seconds: None,
        }
    }
}

impl ReproducibilityReport {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a pipeline step
    pub fn add_step(&mut self, step: &str) {
        self.pipeline_steps.push(step.to_string());
    }

    /// Record an external tool's version
    pub fn record_tool_version(&mut self, tool_name: &str) {
        if let Some(version) = get_tool_version(tool_name) {
            self.tool_versions.insert(tool_name.to_string(), version);
        } else {
            self.tool_versions
                .insert(tool_name.to_string(), "not found".to_string());
        }
    }

    /// Record multiple tools at once
    pub fn record_tools(&mut self, tools: &[&str]) {
        for tool in tools {
            self.record_tool_version(tool);
        }
    }

    /// Add a database entry
    pub fn add_database(&mut self, name: &str, path: &Path) {
        let size_bytes = std::fs::metadata(path).map(|m| m.len()).ok();
        let download_date = std::fs::metadata(path)
            .and_then(|m| m.modified())
            .ok()
            .map(|t| {
                let secs = t.duration_since(UNIX_EPOCH).unwrap_or_default().as_secs();
                timestamp_to_iso8601(secs)
            });

        self.databases.push(DatabaseInfo {
            name: name.to_string(),
            path: path.display().to_string(),
            version: None,
            download_date,
            size_bytes,
        });
    }

    /// Add an input file with optional MD5 checksum
    pub fn add_input_file(&mut self, path: &Path, compute_md5: bool) {
        let size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        let md5 = if compute_md5 {
            compute_file_md5(path)
        } else {
            None
        };

        self.input_files.push(FileInfo {
            path: path.display().to_string(),
            size_bytes: size,
            md5,
        });
    }

    /// Add a parameter
    pub fn add_param(&mut self, key: &str, value: &str) {
        self.parameters.insert(key.to_string(), value.to_string());
    }

    /// Set runtime duration
    pub fn set_runtime(&mut self, seconds: f64) {
        self.runtime_seconds = Some(seconds);
    }

    /// Write report as JSON
    pub fn write_json(&self, path: &Path) -> Result<()> {
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| MycoNoteError::InvalidFormat(format!("JSON serialization: {}", e)))?;
        std::fs::write(path, json).map_err(MycoNoteError::Io)?;
        Ok(())
    }

    /// Write human-readable report
    pub fn write_text(&self, path: &Path) -> Result<()> {
        let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;

        writeln!(
            f,
            "================================================================"
        )
        .map_err(MycoNoteError::Io)?;
        writeln!(f, "  myconote-cli Reproducibility Report").map_err(MycoNoteError::Io)?;
        writeln!(
            f,
            "================================================================"
        )
        .map_err(MycoNoteError::Io)?;
        writeln!(f).map_err(MycoNoteError::Io)?;
        writeln!(f, "Version   : {}", self.myconote_version).map_err(MycoNoteError::Io)?;
        writeln!(f, "Timestamp : {}", self.timestamp).map_err(MycoNoteError::Io)?;
        writeln!(
            f,
            "System    : {} / {}",
            self.system_info.os, self.system_info.arch
        )
        .map_err(MycoNoteError::Io)?;
        writeln!(f, "Hostname  : {}", self.system_info.hostname).map_err(MycoNoteError::Io)?;
        writeln!(f, "CPUs      : {}", self.system_info.cpus).map_err(MycoNoteError::Io)?;

        if let Some(secs) = self.runtime_seconds {
            let mins = secs / 60.0;
            if mins > 60.0 {
                writeln!(f, "Runtime   : {:.1} hours", mins / 60.0).map_err(MycoNoteError::Io)?;
            } else {
                writeln!(f, "Runtime   : {:.1} minutes", mins).map_err(MycoNoteError::Io)?;
            }
        }

        // Pipeline steps
        if !self.pipeline_steps.is_empty() {
            writeln!(f).map_err(MycoNoteError::Io)?;
            writeln!(f, "Pipeline Steps:").map_err(MycoNoteError::Io)?;
            for (i, step) in self.pipeline_steps.iter().enumerate() {
                writeln!(f, "  {}. {}", i + 1, step).map_err(MycoNoteError::Io)?;
            }
        }

        // Parameters
        if !self.parameters.is_empty() {
            writeln!(f).map_err(MycoNoteError::Io)?;
            writeln!(f, "Parameters:").map_err(MycoNoteError::Io)?;
            let mut params: Vec<_> = self.parameters.iter().collect();
            params.sort_by_key(|(k, _)| (*k).clone());
            for (k, v) in params {
                writeln!(f, "  {:<24} {}", k, v).map_err(MycoNoteError::Io)?;
            }
        }

        // Tool versions
        if !self.tool_versions.is_empty() {
            writeln!(f).map_err(MycoNoteError::Io)?;
            writeln!(f, "External Tool Versions:").map_err(MycoNoteError::Io)?;
            let mut tools: Vec<_> = self.tool_versions.iter().collect();
            tools.sort_by_key(|(k, _)| (*k).clone());
            for (tool, ver) in tools {
                writeln!(f, "  {:<24} {}", tool, ver).map_err(MycoNoteError::Io)?;
            }
        }

        // Databases
        if !self.databases.is_empty() {
            writeln!(f).map_err(MycoNoteError::Io)?;
            writeln!(f, "Databases:").map_err(MycoNoteError::Io)?;
            for db in &self.databases {
                writeln!(f, "  {:<20} {}", db.name, db.path).map_err(MycoNoteError::Io)?;
                if let Some(ref date) = db.download_date {
                    writeln!(f, "    Downloaded: {}", date).map_err(MycoNoteError::Io)?;
                }
                if let Some(size) = db.size_bytes {
                    writeln!(f, "    Size: {} MB", size / 1_000_000).map_err(MycoNoteError::Io)?;
                }
            }
        }

        // Input files
        if !self.input_files.is_empty() {
            writeln!(f).map_err(MycoNoteError::Io)?;
            writeln!(f, "Input Files:").map_err(MycoNoteError::Io)?;
            for fi in &self.input_files {
                writeln!(f, "  {} ({} bytes)", fi.path, fi.size_bytes)
                    .map_err(MycoNoteError::Io)?;
                if let Some(ref md5) = fi.md5 {
                    writeln!(f, "    MD5: {}", md5).map_err(MycoNoteError::Io)?;
                }
            }
        }

        writeln!(f).map_err(MycoNoteError::Io)?;
        writeln!(
            f,
            "================================================================"
        )
        .map_err(MycoNoteError::Io)?;

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn current_iso8601() -> String {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    timestamp_to_iso8601(now)
}

fn timestamp_to_iso8601(secs: u64) -> String {
    // Simple UTC timestamp without chrono dependency
    let days = secs / 86400;
    let time_secs = secs % 86400;
    let hours = time_secs / 3600;
    let minutes = (time_secs % 3600) / 60;
    let seconds = time_secs % 60;

    // Approximate year/month/day (good enough for logging)
    let mut y = 1970u64;
    let mut remaining_days = days;
    loop {
        let days_in_year = if is_leap(y) { 366 } else { 365 };
        if remaining_days < days_in_year {
            break;
        }
        remaining_days -= days_in_year;
        y += 1;
    }
    let months_days = if is_leap(y) {
        [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };
    let mut m = 1u64;
    for md in &months_days {
        if remaining_days < *md as u64 {
            break;
        }
        remaining_days -= *md as u64;
        m += 1;
    }
    let d = remaining_days + 1;

    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        y, m, d, hours, minutes, seconds
    )
}

fn is_leap(y: u64) -> bool {
    (y % 4 == 0 && y % 100 != 0) || y % 400 == 0
}

fn detect_system_info() -> SystemInfo {
    let os = std::env::consts::OS.to_string();
    let arch = std::env::consts::ARCH.to_string();

    let hostname = Command::new("hostname")
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|_| "unknown".to_string());

    let cpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);

    SystemInfo {
        os,
        arch,
        hostname,
        cpus,
    }
}

/// Get the version string of an external tool.
fn get_tool_version(tool: &str) -> Option<String> {
    let version_args: &[&str] = match tool {
        "augustus" | "RepeatMasker" | "RepeatModeler" | "samtools" | "minimap2" | "Trinity"
        | "antismash" | "busco" | "tRNAscan-SE" | "emapper.py" | "signalp6" | "biolib" => {
            &["--version"]
        }
        "mmseqs" | "diamond" => &["version"],
        "blastp" => &["-version"],
        "hmmscan" => &["-h"],
        "snap" | "glimmerhmm" => &["--help"],
        "iqtree" | "iqtree2" => &["--version"],
        _ => &["--version"],
    };

    let output = Command::new(tool).args(version_args).output().ok()?;

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    // Take first non-empty line from stdout or stderr
    let version = stdout
        .lines()
        .chain(stderr.lines())
        .find(|l| !l.trim().is_empty())
        .unwrap_or("unknown")
        .trim()
        .to_string();

    // Truncate long version strings
    if version.len() > 120 {
        Some(format!("{}...", &version[..117]))
    } else {
        Some(version)
    }
}

/// Compute MD5 checksum of a file using the system's md5/md5sum command.
fn compute_file_md5(path: &Path) -> Option<String> {
    let path_str = path.to_str()?;

    // Try md5sum (Linux) first, then md5 (macOS)
    let output = Command::new("md5sum")
        .arg(path_str)
        .output()
        .or_else(|_| Command::new("md5").arg("-q").arg(path_str).output())
        .ok()?;

    if !output.status.success() {
        return None;
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    // md5sum format: "hash  filename", md5 -q format: "hash"
    stdout.split_whitespace().next().map(|s| s.to_string())
}
