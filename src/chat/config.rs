use crate::utils::error::{MycoNoteError, Result};
use serde::Deserialize;
use std::path::PathBuf;

const DEFAULT_ENDPOINT: &str = "http://localhost:11434";
const DEFAULT_MODEL: &str = "llama3.1";
const DEFAULT_TIMEOUT_S: u64 = 120;

/// Configuration for the `explain` subcommand.
#[derive(Debug, Clone)]
pub struct ChatConfig {
    pub enabled: bool,
    pub endpoint: String,
    pub model: String,
    pub api_key: Option<String>,
    pub provider: String,
    pub timeout_s: u64,
    pub corpus_dir: Option<PathBuf>,
    pub no_llm: bool,
    pub verbose: bool,
    pub trace: bool,
    pub dry_run: bool,
}

impl Default for ChatConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            endpoint: DEFAULT_ENDPOINT.to_string(),
            model: DEFAULT_MODEL.to_string(),
            api_key: None,
            provider: "ollama".to_string(),
            timeout_s: DEFAULT_TIMEOUT_S,
            corpus_dir: None,
            no_llm: false,
            verbose: false,
            trace: false,
            dry_run: false,
        }
    }
}

/// TOML representation of the `[chat]` section in `~/.myconote/config.toml`.
#[derive(Deserialize, Default)]
struct ChatToml {
    enabled: Option<bool>,
    endpoint: Option<String>,
    model: Option<String>,
    api_key: Option<String>,
    provider: Option<String>,
    timeout_s: Option<u64>,
    corpus_dir: Option<String>,
}

/// Top-level config file shape (only `[chat]` used here).
#[derive(Deserialize, Default)]
struct ConfigFile {
    chat: Option<ChatToml>,
}

impl ChatConfig {
    /// Load configuration with precedence: CLI flags > env vars > TOML file > defaults.
    pub fn load(cli: &CliOverrides) -> Result<Self> {
        let mut cfg = Self::default();

        // Layer 1: TOML file
        if let Some(toml_cfg) = load_toml_config() {
            if let Some(chat) = toml_cfg.chat {
                if let Some(v) = chat.enabled { cfg.enabled = v; }
                if let Some(v) = chat.endpoint { cfg.endpoint = v; }
                if let Some(v) = chat.model { cfg.model = v; }
                if let Some(v) = chat.api_key { cfg.api_key = Some(v); }
                if let Some(v) = chat.provider { cfg.provider = v; }
                if let Some(v) = chat.timeout_s { cfg.timeout_s = v; }
                if let Some(v) = chat.corpus_dir { cfg.corpus_dir = Some(PathBuf::from(v)); }
            }
        }

        // Layer 2: Environment variables
        if let Ok(v) = std::env::var("MYCONOTE_CHAT_ENDPOINT") {
            cfg.endpoint = v;
        }
        if let Ok(v) = std::env::var("MYCONOTE_CHAT_MODEL") {
            cfg.model = v;
        }
        if let Ok(v) = std::env::var("MYCONOTE_CHAT_API_KEY") {
            cfg.api_key = Some(v);
        }

        // Layer 3: CLI flags
        if let Some(ref v) = cli.endpoint { cfg.endpoint = v.clone(); }
        if let Some(ref v) = cli.model { cfg.model = v.clone(); }
        cfg.no_llm = cli.no_llm;
        cfg.verbose = cli.verbose;
        cfg.trace = cli.trace;
        cfg.dry_run = cli.dry_run;

        Ok(cfg)
    }
}

/// CLI flag overrides parsed in `mod.rs`.
#[derive(Debug, Default)]
pub struct CliOverrides {
    pub endpoint: Option<String>,
    pub model: Option<String>,
    pub no_llm: bool,
    pub verbose: bool,
    pub trace: bool,
    pub dry_run: bool,
}

/// Try to load `~/.myconote/config.toml`. Returns `None` if the file
/// doesn't exist or can't be parsed (non-fatal).
fn load_toml_config() -> Option<ConfigFile> {
    let home = std::env::var("HOME").ok()?;
    let path = PathBuf::from(home).join(".myconote").join("config.toml");
    let text = std::fs::read_to_string(&path).ok()?;
    toml::from_str(&text).ok()
}

/// Return the `~/.myconote/` directory, creating it if necessary.
pub fn myconote_dir() -> Result<PathBuf> {
    let home = std::env::var("HOME")
        .map_err(|_| MycoNoteError::ChatConfig("HOME environment variable not set".to_string()))?;
    let dir = PathBuf::from(home).join(".myconote");
    if !dir.exists() {
        std::fs::create_dir_all(&dir)
            .map_err(|e| MycoNoteError::ChatConfig(format!("failed to create {}: {}", dir.display(), e)))?;
    }
    Ok(dir)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_config_values() {
        let cfg = ChatConfig::default();
        assert!(cfg.enabled);
        assert_eq!(cfg.endpoint, "http://localhost:11434");
        assert_eq!(cfg.model, "llama3.1");
        assert_eq!(cfg.provider, "ollama");
        assert_eq!(cfg.timeout_s, 120);
        assert!(!cfg.no_llm);
        assert!(!cfg.verbose);
        assert!(!cfg.trace);
        assert!(!cfg.dry_run);
    }

    #[test]
    fn toml_parse() {
        let toml_str = r#"
[chat]
enabled = false
endpoint = "http://example.com:11434"
model = "mistral"
timeout_s = 60
"#;
        let config: ConfigFile = toml::from_str(toml_str).unwrap();
        let chat = config.chat.unwrap();
        assert_eq!(chat.enabled, Some(false));
        assert_eq!(chat.endpoint.unwrap(), "http://example.com:11434");
        assert_eq!(chat.model.unwrap(), "mistral");
        assert_eq!(chat.timeout_s, Some(60));
    }

    #[test]
    fn cli_overrides_take_precedence() {
        let cli = CliOverrides {
            endpoint: Some("http://custom:1234".to_string()),
            model: Some("phi3".to_string()),
            no_llm: true,
            verbose: true,
            trace: false,
            dry_run: true,
        };
        let cfg = ChatConfig::load(&cli).unwrap();
        assert_eq!(cfg.endpoint, "http://custom:1234");
        assert_eq!(cfg.model, "phi3");
        assert!(cfg.no_llm);
        assert!(cfg.verbose);
        assert!(cfg.dry_run);
    }
}
