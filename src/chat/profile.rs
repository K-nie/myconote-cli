use crate::utils::error::{MycoNoteError, Result};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct UserProfile {
    pub kingdom: Option<String>,
    pub organism: Option<String>,
    pub typical_genome_size_mb: Option<f64>,
    pub preferred_verbosity: Option<String>,
}

impl UserProfile {
    pub fn load() -> Self {
        let path = profile_path();
        match path {
            Some(p) if p.exists() => std::fs::read_to_string(&p)
                .ok()
                .and_then(|t| toml::from_str(&t).ok())
                .unwrap_or_default(),
            _ => Self::default(),
        }
    }

    pub fn save(&self) -> Result<()> {
        let dir = super::config::myconote_dir()?;
        let path = dir.join("profile.toml");
        let text = toml::to_string_pretty(self).map_err(|e| {
            MycoNoteError::ChatConfig(format!("failed to serialize profile: {}", e))
        })?;
        std::fs::write(&path, text).map_err(|e| {
            MycoNoteError::ChatConfig(format!("failed to write {}: {}", path.display(), e))
        })?;
        Ok(())
    }

    pub fn summary(&self) -> String {
        let mut parts = Vec::new();
        if let Some(ref k) = self.kingdom {
            parts.push(format!("kingdom: {}", k));
        }
        if let Some(ref o) = self.organism {
            parts.push(format!("organism: {}", o));
        }
        if let Some(s) = self.typical_genome_size_mb {
            parts.push(format!("genome: {:.0} Mb", s));
        }
        if parts.is_empty() {
            "no profile set".to_string()
        } else {
            parts.join(", ")
        }
    }
}

fn profile_path() -> Option<PathBuf> {
    let home = std::env::var("HOME").ok()?;
    Some(PathBuf::from(home).join(".myconote").join("profile.toml"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_profile() {
        let p = UserProfile::default();
        assert!(p.kingdom.is_none());
        assert_eq!(p.summary(), "no profile set");
    }

    #[test]
    fn profile_summary() {
        let p = UserProfile {
            kingdom: Some("fungi".to_string()),
            organism: Some("Aspergillus niger".to_string()),
            typical_genome_size_mb: Some(35.0),
            preferred_verbosity: None,
        };
        let s = p.summary();
        assert!(s.contains("fungi"));
        assert!(s.contains("Aspergillus"));
        assert!(s.contains("35 Mb"));
    }
}
