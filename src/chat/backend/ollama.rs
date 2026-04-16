use super::{ChatBackend, Message};
use crate::utils::error::{MycoNoteError, Result};
use serde::{Deserialize, Serialize};
use std::time::Duration;

/// Ollama chat backend — POSTs to the local `/api/chat` endpoint.
pub struct OllamaBackend {
    client: reqwest::blocking::Client,
    endpoint: String,
    model: String,
}

#[derive(Serialize)]
struct OllamaRequest<'a> {
    model: &'a str,
    messages: &'a [Message],
    stream: bool,
}

#[derive(Deserialize)]
struct OllamaResponse {
    message: OllamaMessage,
}

#[derive(Deserialize)]
struct OllamaMessage {
    role: String,
    content: String,
}

impl OllamaBackend {
    pub fn new(endpoint: &str, model: &str, timeout_s: u64) -> Result<Self> {
        let client = reqwest::blocking::Client::builder()
            .user_agent("myconote-cli/0.1 (genome annotation tool)")
            .timeout(Duration::from_secs(timeout_s))
            .build()
            .map_err(|e| {
                MycoNoteError::ChatBackend(format!("failed to build HTTP client: {}", e))
            })?;

        Ok(Self {
            client,
            endpoint: endpoint.trim_end_matches('/').to_string(),
            model: model.to_string(),
        })
    }

    /// Check whether the Ollama server is reachable.
    pub fn is_available(&self) -> bool {
        let url = format!("{}/api/tags", self.endpoint);
        self.client.get(&url).send().is_ok()
    }
}

impl ChatBackend for OllamaBackend {
    fn chat(&self, messages: &[Message]) -> Result<Message> {
        let url = format!("{}/api/chat", self.endpoint);

        let request = OllamaRequest {
            model: &self.model,
            messages,
            stream: false,
        };

        let response = self.client.post(&url).json(&request).send().map_err(|e| {
            if e.is_connect() {
                MycoNoteError::ChatBackend(format!(
                    "connection refused to {} — is `ollama serve` running?",
                    self.endpoint
                ))
            } else if e.is_timeout() {
                MycoNoteError::ChatBackend(format!(
                    "request timed out after {}s — try a smaller model or increase --timeout",
                    e
                ))
            } else {
                MycoNoteError::ChatBackend(format!("HTTP error: {}", e))
            }
        })?;

        let status = response.status();
        if !status.is_success() {
            let body = response.text().unwrap_or_default();
            return Err(MycoNoteError::ChatBackend(format!(
                "Ollama returned HTTP {}: {}",
                status, body
            )));
        }

        let resp: OllamaResponse = response.json().map_err(|e| {
            MycoNoteError::ChatBackend(format!("failed to parse Ollama response: {}", e))
        })?;

        Ok(Message {
            role: resp.message.role,
            content: resp.message.content,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn request_serialization() {
        let messages = vec![Message::user("hello")];
        let req = OllamaRequest {
            model: "llama3.1",
            messages: &messages,
            stream: false,
        };
        let json = serde_json::to_value(&req).unwrap();
        assert_eq!(json["model"], "llama3.1");
        assert_eq!(json["stream"], false);
        assert_eq!(json["messages"][0]["role"], "user");
        assert_eq!(json["messages"][0]["content"], "hello");
    }

    #[test]
    fn message_constructors() {
        let s = Message::system("sys");
        assert_eq!(s.role, "system");
        assert_eq!(s.content, "sys");

        let u = Message::user("usr");
        assert_eq!(u.role, "user");

        let a = Message::assistant("asst");
        assert_eq!(a.role, "assistant");
    }
}
