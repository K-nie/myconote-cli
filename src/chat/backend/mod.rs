use crate::utils::error::Result;
use serde::{Deserialize, Serialize};

/// A single message in a chat conversation.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Message {
    pub role: String,
    pub content: String,
}

impl Message {
    pub fn system(content: impl Into<String>) -> Self {
        Self {
            role: "system".to_string(),
            content: content.into(),
        }
    }

    pub fn user(content: impl Into<String>) -> Self {
        Self {
            role: "user".to_string(),
            content: content.into(),
        }
    }

    pub fn assistant(content: impl Into<String>) -> Self {
        Self {
            role: "assistant".to_string(),
            content: content.into(),
        }
    }
}

/// Backend trait for LLM providers.
///
/// v1 ships only `OllamaBackend`; the trait exists so future providers
/// (Anthropic, OpenAI) can be added without touching the orchestration code.
pub trait ChatBackend {
    fn chat(&self, messages: &[Message]) -> Result<Message>;
}

pub mod ollama;
