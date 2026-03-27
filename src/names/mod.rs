/// Gene name resolver
///
/// Maps gene IDs (PromBase, FungiDB, NCBI, etc.) to human-readable gene names.
///
/// Resolution priority:
///   1. Local override file (`--names id_to_name.tsv`)
///   2. Local cache (`~/.myconote/names_cache.tsv`)
///   3. Remote API (NCBI Gene → UniProt → FungiDB/VEuPathDB)
///
/// The cache is a plain TSV: `gene_id\tgene_name\tsource\ttimestamp`
/// Once an ID is resolved it is never re-fetched (unless the cache is cleared).

pub mod ncbi;
pub mod uniprot;
pub mod fungidb;

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

// ─────────────────────────────────────────────────────────────────────────────
// Public types
// ─────────────────────────────────────────────────────────────────────────────

/// Where a resolved name came from — useful for debugging / provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NameSource {
    LocalFile,
    Cache,
    Ncbi,
    UniProt,
    FungiDb,
    Unknown,
}

impl std::fmt::Display for NameSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NameSource::LocalFile => write!(f, "local"),
            NameSource::Cache    => write!(f, "cache"),
            NameSource::Ncbi     => write!(f, "ncbi"),
            NameSource::UniProt  => write!(f, "uniprot"),
            NameSource::FungiDb  => write!(f, "fungidb"),
            NameSource::Unknown  => write!(f, "unknown"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ResolvedName {
    pub gene_id:   String,
    pub gene_name: String,
    pub source:    NameSource,
}

// ─────────────────────────────────────────────────────────────────────────────
// NameResolver
// ─────────────────────────────────────────────────────────────────────────────

/// Main entry point for gene name resolution.
pub struct NameResolver {
    /// ID → resolved name (all sources merged here)
    map: HashMap<String, ResolvedName>,
    cache_path: PathBuf,
    dirty: bool, // cache needs flushing
}

impl NameResolver {
    /// Create a resolver, loading the local cache if it exists.
    pub fn new() -> Self {
        let cache_path = cache_file_path();
        let mut resolver = NameResolver {
            map:        HashMap::new(),
            cache_path: cache_path.clone(),
            dirty:      false,
        };
        if cache_path.exists() {
            resolver.load_tsv(&cache_path, NameSource::Cache);
        }
        resolver
    }

    /// Load an additional local override file (highest priority).
    /// Format: `gene_id\tgene_name` (tab-separated, one per line, `#` comments ok)
    pub fn load_local_file<P: AsRef<Path>>(&mut self, path: P) {
        self.load_tsv(path.as_ref(), NameSource::LocalFile);
    }

    /// Return the resolved name for `id`, or `None` if not found in loaded data.
    pub fn get(&self, id: &str) -> Option<&str> {
        self.map.get(id).map(|r| r.gene_name.as_str())
    }

    /// Return the best display label for `id`:
    /// - the resolved gene name if available
    /// - otherwise the original ID
    pub fn label(&self, id: &str) -> String {
        self.get(id)
            .map(str::to_string)
            .unwrap_or_else(|| id.to_string())
    }

    /// Resolve a batch of IDs that are not yet in the map.
    /// Tries NCBI, then UniProt, then FungiDB in order.
    /// Results are stored in the map and the cache is written.
    pub fn fetch_missing(&mut self, ids: &[String], organism_taxon: Option<u32>) {
        let missing: Vec<String> = ids.iter()
            .filter(|id| !self.map.contains_key(id.as_str()))
            .cloned()
            .collect();

        if missing.is_empty() { return; }

        eprintln!("  🔍 Fetching names for {} unresolved IDs...", missing.len());

        // Try NCBI Gene
        let ncbi_results = ncbi::fetch_gene_names(&missing, organism_taxon);
        for (id, name) in &ncbi_results {
            self.insert(id, name, NameSource::Ncbi);
        }

        // Remaining unresolved → UniProt
        let still_missing: Vec<String> = missing.iter()
            .filter(|id| !self.map.contains_key(id.as_str()))
            .cloned()
            .collect();

        if !still_missing.is_empty() {
            let uniprot_results = uniprot::fetch_gene_names(&still_missing);
            for (id, name) in &uniprot_results {
                self.insert(id, name, NameSource::UniProt);
            }
        }

        // Remaining → FungiDB
        let still_missing2: Vec<String> = missing.iter()
            .filter(|id| !self.map.contains_key(id.as_str()))
            .cloned()
            .collect();

        if !still_missing2.is_empty() {
            let fungidb_results = fungidb::fetch_gene_names(&still_missing2);
            for (id, name) in &fungidb_results {
                self.insert(id, name, NameSource::FungiDb);
            }
        }

        if self.dirty {
            self.flush_cache();
        }

        let resolved = ids.iter().filter(|id| self.map.contains_key(id.as_str())).count();
        eprintln!("  ✓ Resolved {}/{} gene names", resolved, ids.len());
    }

    /// Insert a resolved name (skips if a higher-priority source already exists).
    fn insert(&mut self, id: &str, name: &str, source: NameSource) {
        // Local file > cache > remote — don't overwrite higher-priority entries
        if let Some(existing) = self.map.get(id) {
            if existing.source == NameSource::LocalFile { return; }
            if existing.source == NameSource::Cache && source != NameSource::LocalFile { return; }
        }
        self.map.insert(id.to_string(), ResolvedName {
            gene_id:   id.to_string(),
            gene_name: name.to_string(),
            source,
        });
        self.dirty = true;
    }

    /// Write the current map to the local cache file.
    pub fn flush_cache(&mut self) {
        if let Some(parent) = self.cache_path.parent() {
            let _ = fs::create_dir_all(parent);
        }
        match File::create(&self.cache_path) {
            Ok(mut file) => {
                let _ = writeln!(file, "# myconote gene name cache — do not edit manually");
                for r in self.map.values() {
                    let _ = writeln!(file, "{}\t{}\t{}", r.gene_id, r.gene_name, r.source);
                }
                self.dirty = false;
            }
            Err(e) => eprintln!("  ⚠  Could not write name cache: {}", e),
        }
    }

    /// Export the current map as a TSV string (for embedding in HTML).
    pub fn to_tsv(&self) -> String {
        self.map.values()
            .map(|r| format!("{}\t{}", r.gene_id, r.gene_name))
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// Number of resolved names currently loaded.
    pub fn len(&self) -> usize { self.map.len() }

    pub fn is_empty(&self) -> bool { self.map.is_empty() }

    // ── Private helpers ───────────────────────────────────────────────────

    fn load_tsv(&mut self, path: &Path, source: NameSource) {
        let Ok(file) = File::open(path) else { return };
        for line in BufReader::new(file).lines().flatten() {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            let mut cols = t.splitn(3, '\t');
            if let (Some(id), Some(name)) = (cols.next(), cols.next()) {
                let s = source.clone();
                self.insert(id.trim(), name.trim(), s);
            }
        }
        // Loading from file doesn't make the cache dirty
        if source != NameSource::LocalFile {
            self.dirty = false;
        }
    }
}

impl Drop for NameResolver {
    fn drop(&mut self) {
        if self.dirty { self.flush_cache(); }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cache path
// ─────────────────────────────────────────────────────────────────────────────

fn cache_file_path() -> PathBuf {
    // Respect XDG_CACHE_HOME if set, otherwise ~/.myconote/
    let base = std::env::var("XDG_CACHE_HOME")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            dirs_next()
                .unwrap_or_else(|| PathBuf::from("."))
                .join(".myconote")
        });
    base.join("names_cache.tsv")
}

fn dirs_next() -> Option<PathBuf> {
    std::env::var("HOME").ok().map(PathBuf::from)
}
