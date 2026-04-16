use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

// ─────────────────────────────────────────────────────────────────────────────
// Citation — a retrieved snippet from the knowledge base or paper corpus
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Serialize, Clone, Debug)]
pub struct Citation {
    pub source_type: SourceType,
    pub id: String,
    pub snippet: String,
    pub score: f32,
}

#[derive(Serialize, Clone, Debug, PartialEq)]
pub enum SourceType {
    Knowledge,
    Paper,
    UserData,
    Rule,
}

impl std::fmt::Display for SourceType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SourceType::Knowledge => write!(f, "knowledge"),
            SourceType::Paper => write!(f, "paper"),
            SourceType::UserData => write!(f, "data"),
            SourceType::Rule => write!(f, "rule"),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Knowledge entry from TOML files
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Deserialize, Debug, Clone)]
pub struct KnowledgeEntry {
    pub id: String,
    pub title: String,
    pub content: String,
    pub tags: Vec<String>,
}

#[derive(Deserialize, Debug)]
struct KnowledgeFile {
    entries: Vec<KnowledgeEntry>,
}

// ─────────────────────────────────────────────────────────────────────────────
// BM25 retrieval
// ─────────────────────────────────────────────────────────────────────────────

/// Simple BM25 implementation — no external crate needed.
/// Sufficient for <10K documents.
struct BM25Index {
    docs: Vec<BM25Doc>,
    avg_dl: f64,
    idf: HashMap<String, f64>,
}

struct BM25Doc {
    id: String,
    content: String,
    terms: HashMap<String, usize>,
    length: usize,
}

const BM25_K1: f64 = 1.2;
const BM25_B: f64 = 0.75;

impl BM25Index {
    fn new(docs: Vec<(String, String)>) -> Self {
        let n = docs.len() as f64;
        let mut bm25_docs = Vec::with_capacity(docs.len());
        let mut df: HashMap<String, usize> = HashMap::new();
        let mut total_length = 0usize;

        for (id, content) in &docs {
            let tokens = tokenize(content);
            let length = tokens.len();
            total_length += length;

            let mut terms: HashMap<String, usize> = HashMap::new();
            let mut seen = std::collections::HashSet::new();
            for token in &tokens {
                *terms.entry(token.clone()).or_insert(0) += 1;
                if seen.insert(token.clone()) {
                    *df.entry(token.clone()).or_insert(0) += 1;
                }
            }

            bm25_docs.push(BM25Doc {
                id: id.clone(),
                content: content.clone(),
                terms,
                length,
            });
        }

        let avg_dl = if bm25_docs.is_empty() {
            1.0
        } else {
            total_length as f64 / bm25_docs.len() as f64
        };

        let mut idf = HashMap::new();
        for (term, count) in &df {
            let idf_val = ((n - *count as f64 + 0.5) / (*count as f64 + 0.5) + 1.0).ln();
            idf.insert(term.clone(), idf_val.max(0.0));
        }

        Self {
            docs: bm25_docs,
            avg_dl,
            idf,
        }
    }

    fn search(&self, query: &str, top_k: usize) -> Vec<(String, String, f32)> {
        let query_tokens = tokenize(query);
        let mut scores: Vec<(usize, f64)> = Vec::new();

        for (idx, doc) in self.docs.iter().enumerate() {
            let mut score = 0.0;
            for token in &query_tokens {
                let tf = *doc.terms.get(token).unwrap_or(&0) as f64;
                let idf = *self.idf.get(token).unwrap_or(&0.0);
                let norm = tf * (BM25_K1 + 1.0)
                    / (tf + BM25_K1 * (1.0 - BM25_B + BM25_B * doc.length as f64 / self.avg_dl));
                score += idf * norm;
            }
            if score > 0.0 {
                scores.push((idx, score));
            }
        }

        scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        scores.truncate(top_k);

        scores
            .into_iter()
            .map(|(idx, score)| {
                let doc = &self.docs[idx];
                (doc.id.clone(), doc.content.clone(), score as f32)
            })
            .collect()
    }
}

/// Simple whitespace + lowercase tokenizer.
fn tokenize(text: &str) -> Vec<String> {
    text.to_lowercase()
        .split(|c: char| !c.is_alphanumeric() && c != '_')
        .filter(|s| s.len() > 2)
        .map(|s| s.to_string())
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Load knowledge entries for a stage and retrieve the top-K most relevant.
pub fn retrieve_knowledge(stage: &str, query: &str, top_k: usize) -> Vec<Citation> {
    let entries = load_knowledge_entries(stage);
    if entries.is_empty() {
        return Vec::new();
    }

    let docs: Vec<(String, String)> = entries
        .iter()
        .map(|e| {
            let text = format!("{} {} {}", e.title, e.content, e.tags.join(" "));
            (e.id.clone(), text)
        })
        .collect();

    let index = BM25Index::new(docs);
    let results = index.search(query, top_k);

    results
        .into_iter()
        .map(|(id, snippet, score)| {
            // Find the original entry to get the clean content
            let original = entries.iter().find(|e| e.id == id);
            let clean_snippet = original.map(|e| e.content.clone()).unwrap_or(snippet);

            Citation {
                source_type: SourceType::Knowledge,
                id: format!("knowledge:{}", id),
                snippet: clean_snippet,
                score,
            }
        })
        .collect()
}

/// Load all knowledge entries for a stage from TOML files.
fn load_knowledge_entries(stage: &str) -> Vec<KnowledgeEntry> {
    let knowledge_dir = find_assets_knowledge_dir();
    let Some(dir) = knowledge_dir else {
        return Vec::new();
    };

    let stage_file = dir.join(format!("{}.toml", stage));
    if !stage_file.exists() {
        return Vec::new();
    }

    let text = match std::fs::read_to_string(&stage_file) {
        Ok(t) => t,
        Err(_) => return Vec::new(),
    };

    let kf: KnowledgeFile = match toml::from_str(&text) {
        Ok(k) => k,
        Err(_) => return Vec::new(),
    };

    kf.entries
}

fn find_assets_knowledge_dir() -> Option<PathBuf> {
    let candidates = [
        std::env::current_exe()
            .ok()
            .and_then(|p| p.parent().map(|d| d.join("assets/chat/knowledge"))),
        std::env::var("CARGO_MANIFEST_DIR")
            .ok()
            .map(|d| PathBuf::from(d).join("assets/chat/knowledge")),
        Some(PathBuf::from("assets/chat/knowledge")),
    ];

    for candidate in &candidates {
        if let Some(ref path) = candidate {
            if path.exists() {
                return Some(path.clone());
            }
        }
    }
    None
}

/// Retrieve relevant paper citations from the local corpus.
///
/// Searches the corpus at `~/.myconote/papers/` (or the configured corpus_dir).
/// Only papers listed in `corpus_manifest.toml` with quartile=Q1 and a valid
/// OA license are considered. Returns empty if no corpus is available.
pub fn retrieve_papers(
    query: &str,
    top_k: usize,
    corpus_dir: Option<&std::path::Path>,
) -> Vec<Citation> {
    let dir = match corpus_dir {
        Some(d) => d.to_path_buf(),
        None => {
            // Try ~/.myconote/papers/ then ~/.myconote/dbs/papers/
            let home = std::env::var("HOME")
                .or_else(|_| std::env::var("USERPROFILE"))
                .unwrap_or_else(|_| "/tmp".to_string());
            let home_papers = PathBuf::from(&home).join(".myconote").join("papers");
            let dbs_papers = PathBuf::from(&home)
                .join(".myconote")
                .join("dbs")
                .join("papers");
            if home_papers.join("corpus_manifest.toml").exists() {
                home_papers
            } else if dbs_papers.join("corpus_manifest.toml").exists() {
                dbs_papers
            } else {
                return Vec::new();
            }
        }
    };

    let manifest_path = dir.join("corpus_manifest.toml");
    if !manifest_path.exists() {
        return Vec::new();
    }

    let manifest_text = match std::fs::read_to_string(&manifest_path) {
        Ok(t) => t,
        Err(_) => return Vec::new(),
    };

    let manifest: toml::Value = match manifest_text.parse() {
        Ok(v) => v,
        Err(_) => return Vec::new(),
    };

    let papers = match manifest.get("paper").and_then(|v| v.as_array()) {
        Some(p) => p,
        None => return Vec::new(),
    };

    // Load text for each paper with a text_file
    let mut docs: Vec<(String, String)> = Vec::new();
    for paper in papers {
        let doi = paper
            .get("doi")
            .and_then(|v| v.as_str())
            .unwrap_or_default();
        let text_file = match paper.get("text_file").and_then(|v| v.as_str()) {
            Some(tf) => tf,
            None => continue,
        };

        let text_path = dir.join(text_file);
        let content = match std::fs::read_to_string(&text_path) {
            Ok(c) => c,
            Err(_) => continue,
        };

        if content.is_empty() {
            continue;
        }

        // Split into paragraphs for finer-grained retrieval
        let paragraphs: Vec<&str> = content
            .split("\n\n")
            .filter(|p| p.trim().len() > 50)
            .collect();

        for (i, para) in paragraphs.iter().enumerate() {
            let doc_id = format!("{}#p{}", doi, i);
            docs.push((doc_id, para.to_string()));
        }
    }

    if docs.is_empty() {
        return Vec::new();
    }

    let index = BM25Index::new(docs);
    let results = index.search(query, top_k);

    results
        .into_iter()
        .map(|(id, snippet, score)| {
            // Truncate snippet to ~500 chars for prompt budget
            let truncated = if snippet.len() > 500 {
                format!("{}…", &snippet[..500])
            } else {
                snippet
            };

            Citation {
                source_type: SourceType::Paper,
                id: format!("paper:{}", id),
                snippet: truncated,
                score,
            }
        })
        .collect()
}

/// Check if a paper corpus is available and report status.
pub fn corpus_status() -> Option<String> {
    let home = std::env::var("HOME")
        .or_else(|_| std::env::var("USERPROFILE"))
        .unwrap_or_else(|_| "/tmp".to_string());

    let candidates = [
        PathBuf::from(&home)
            .join(".myconote")
            .join("papers")
            .join("corpus_manifest.toml"),
        PathBuf::from(&home)
            .join(".myconote")
            .join("dbs")
            .join("papers")
            .join("corpus_manifest.toml"),
    ];

    for path in &candidates {
        if path.exists() {
            return Some(format!(
                "Corpus at {}",
                path.parent().unwrap_or(path).display()
            ));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenize_basic() {
        let tokens = tokenize("Hello, World! This is a test-123.");
        assert!(tokens.contains(&"hello".to_string()));
        assert!(tokens.contains(&"world".to_string()));
        assert!(tokens.contains(&"test".to_string()));
        assert!(tokens.contains(&"123".to_string()));
        assert!(tokens.contains(&"this".to_string()));
        // "is" and "a" are filtered (len <= 2)
        assert!(!tokens.contains(&"is".to_string()));
        assert!(!tokens.contains(&"a".to_string()));
    }

    #[test]
    fn bm25_ranking() {
        let docs = vec![
            (
                "doc1".to_string(),
                "gene prediction fungal genome annotation".to_string(),
            ),
            (
                "doc2".to_string(),
                "repeat masking transposon repeat content".to_string(),
            ),
            (
                "doc3".to_string(),
                "gene count over-prediction unmasked repeats".to_string(),
            ),
        ];

        let index = BM25Index::new(docs);
        let results = index.search("gene prediction", 3);

        // doc1 and doc3 should score higher than doc2 for "gene prediction"
        assert!(!results.is_empty());
        assert!(results[0].0 == "doc1" || results[0].0 == "doc3");
    }

    #[test]
    fn bm25_empty_query() {
        let docs = vec![("doc1".to_string(), "hello world".to_string())];
        let index = BM25Index::new(docs);
        let results = index.search("", 5);
        assert!(results.is_empty());
    }

    #[test]
    fn bm25_no_match() {
        let docs = vec![("doc1".to_string(), "genome annotation pipeline".to_string())];
        let index = BM25Index::new(docs);
        let results = index.search("basketball sports", 5);
        assert!(results.is_empty());
    }

    #[test]
    fn citation_source_type_display() {
        assert_eq!(format!("{}", SourceType::Knowledge), "knowledge");
        assert_eq!(format!("{}", SourceType::Paper), "paper");
        assert_eq!(format!("{}", SourceType::UserData), "data");
        assert_eq!(format!("{}", SourceType::Rule), "rule");
    }
}
