//! ENA REST client + FASTQ downloader.
//!
//! Default backend for `fetch-rna`. Hits the filereport endpoint,
//! parses the TSV response, and streams FASTQ files over HTTPS with
//! on-the-fly MD5 verification. No credentials, no `vdb-config`, no
//! sra-toolkit — most public runs are mirrored on the ENA side and
//! resolvable via a single HTTPS GET.
//!
//! TSV schema (verified empirically against the live ENA API on
//! 2026-04-24):
//!
//! ```text
//! run_accession\tfastq_ftp\tfastq_md5\tlibrary_layout
//! SRR453566\tftp.sra.ebi.ac.uk/...R1.fastq.gz;ftp.sra.ebi.ac.uk/...R2.fastq.gz\tmd5a;md5b\tPAIRED
//! SRR7082700\tftp.sra.ebi.ac.uk/...SRR7082700.fastq.gz\t4a6db...\tSINGLE
//! ```
//!
//! Paired rows separate R1 / R2 URLs and MD5s with a `;`. The URL
//! column is returned without a scheme; we prepend `https://`. Study
//! accessions (PRJ*/SRP/ERP/DRP) are transparently expanded by ENA
//! into one row per run, so the caller does not need to pre-resolve.

use crate::utils::error::{MycoNoteError, Result};
use md5::{Digest, Md5};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

pub const ENA_API: &str = "https://www.ebi.ac.uk/ena/portal/api/filereport";
pub const FIELDS: &str = "run_accession,fastq_ftp,fastq_md5,library_layout";

/// One row of the ENA filereport TSV, normalized into Rust.
#[derive(Debug, Clone, PartialEq)]
pub struct FileReport {
    pub run_accession: String,
    /// URLs exactly as ENA returns them (no scheme prefix). One entry
    /// for single-end, two for paired. Downloaders prepend `https://`.
    pub fastq_urls: Vec<String>,
    /// Per-URL MD5 checksums. Same length and order as `fastq_urls`.
    pub fastq_md5s: Vec<String>,
    pub paired: bool,
}

impl FileReport {
    /// `true` when both R1 and R2 URLs are present (normal paired-end
    /// case). Some paired-lib runs have only a single URL because
    /// upstream only deposited the merged read — we treat those as
    /// single-end for downstream purposes.
    pub fn has_r2(&self) -> bool {
        self.paired && self.fastq_urls.len() >= 2
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TSV parser (pure — no network)
// ─────────────────────────────────────────────────────────────────────────────

/// Parse an ENA filereport TSV into a list of `FileReport` rows.
/// Rejects malformed headers; skips rows with empty `fastq_ftp`
/// (ENA occasionally publishes metadata before FASTQs are mirrored).
pub fn parse_filereport_tsv(text: &str) -> Result<Vec<FileReport>> {
    let mut out = Vec::new();
    let mut col: HashMap<String, usize> = HashMap::new();
    let mut saw_header = false;

    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();

        if !saw_header {
            for (idx, h) in cells.iter().enumerate() {
                col.insert(h.to_string(), idx);
            }
            for required in ["run_accession", "fastq_ftp", "fastq_md5", "library_layout"] {
                if !col.contains_key(required) {
                    return Err(MycoNoteError::QuantTool {
                        tool: "ena".to_string(),
                        message: format!(
                            "filereport header missing '{required}'; got columns: {:?}",
                            cells
                        ),
                    });
                }
            }
            saw_header = true;
            continue;
        }

        let get = |name: &str| -> &str {
            col.get(name)
                .and_then(|i| cells.get(*i))
                .copied()
                .unwrap_or("")
        };

        let fastq_ftp = get("fastq_ftp");
        if fastq_ftp.trim().is_empty() {
            // Metadata-only row — no public URLs to download. Not an
            // error: callers report "N runs had no FASTQs" at the end.
            continue;
        }

        let urls: Vec<String> = fastq_ftp.split(';').map(String::from).collect();
        let md5s: Vec<String> = get("fastq_md5").split(';').map(String::from).collect();

        if urls.len() != md5s.len() {
            return Err(MycoNoteError::QuantTool {
                tool: "ena".to_string(),
                message: format!(
                    "run {}: fastq_ftp has {} urls but fastq_md5 has {} entries",
                    get("run_accession"),
                    urls.len(),
                    md5s.len()
                ),
            });
        }

        // Reject empty MD5 entries — a row with a URL but no checksum
        // defeats the verification and is almost always a partial
        // ENA response we shouldn't consume.
        if md5s.iter().any(|m| m.trim().is_empty()) {
            return Err(MycoNoteError::QuantTool {
                tool: "ena".to_string(),
                message: format!(
                    "run {}: one or more fastq_md5 entries is empty",
                    get("run_accession")
                ),
            });
        }

        let paired = get("library_layout").eq_ignore_ascii_case("PAIRED");
        out.push(FileReport {
            run_accession: get("run_accession").to_string(),
            fastq_urls: urls,
            fastq_md5s: md5s,
            paired,
        });
    }

    if !saw_header {
        return Err(MycoNoteError::QuantTool {
            tool: "ena".to_string(),
            message: "filereport response is empty".to_string(),
        });
    }
    Ok(out)
}

// ─────────────────────────────────────────────────────────────────────────────
// Network — filereport query
// ─────────────────────────────────────────────────────────────────────────────

/// Query ENA's filereport endpoint for an accession and parse the
/// response. One accession can resolve to many runs for study /
/// project / sample IDs.
pub fn fetch_filereport(accession: &str) -> Result<Vec<FileReport>> {
    let url = format!("{ENA_API}?accession={accession}&result=read_run&fields={FIELDS}&format=tsv");
    let client = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(60))
        .build()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "ena".to_string(),
            message: format!("client build: {e}"),
        })?;

    let resp = client
        .get(&url)
        .header(
            "User-Agent",
            concat!("myconote-cli/", env!("CARGO_PKG_VERSION")),
        )
        .send()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "ena".to_string(),
            message: format!("GET {url}: {e}"),
        })?;

    let status = resp.status();
    let text = resp.text().map_err(|e| MycoNoteError::QuantTool {
        tool: "ena".to_string(),
        message: format!("read body {url}: {e}"),
    })?;

    if !status.is_success() {
        return Err(MycoNoteError::QuantTool {
            tool: "ena".to_string(),
            message: format!(
                "ENA returned {status} for accession '{accession}'. \
                 If the accession is an SRA-only archive, retry with `--backend sra`."
            ),
        });
    }

    parse_filereport_tsv(&text)
}

// ─────────────────────────────────────────────────────────────────────────────
// Network — download with MD5 verification
// ─────────────────────────────────────────────────────────────────────────────

/// Download `url` (with or without an `https://` prefix) to `dest`,
/// computing MD5 on the fly. Retries on transient failures up to
/// `retries` times. On MD5 mismatch the partial file is deleted and
/// the failure propagates; the caller can then retry.
pub fn download_with_md5(url: &str, dest: &Path, expected_md5: &str, retries: usize) -> Result<()> {
    let full_url = if url.starts_with("http://") || url.starts_with("https://") {
        url.to_string()
    } else {
        format!("https://{url}")
    };

    let tries = retries.max(1);
    let mut last_err: Option<String> = None;
    for attempt in 1..=tries {
        match try_download_once(&full_url, dest, expected_md5) {
            Ok(()) => return Ok(()),
            Err(e) => {
                let msg = format!("{e}");
                if attempt < tries {
                    eprintln!("  ⚠ attempt {attempt}/{tries} for {full_url}: {msg}");
                }
                last_err = Some(msg);
            }
        }
    }
    Err(MycoNoteError::QuantTool {
        tool: "ena".to_string(),
        message: format!(
            "exhausted {tries} attempt(s) on {full_url}: {}",
            last_err.unwrap_or_default()
        ),
    })
}

fn try_download_once(url: &str, dest: &Path, expected_md5: &str) -> Result<()> {
    let client = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(1800))
        .build()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "ena".to_string(),
            message: format!("client build: {e}"),
        })?;

    let mut resp = client
        .get(url)
        .header(
            "User-Agent",
            concat!("myconote-cli/", env!("CARGO_PKG_VERSION")),
        )
        .send()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "ena".to_string(),
            message: format!("GET {url}: {e}"),
        })?;

    if !resp.status().is_success() {
        return Err(MycoNoteError::QuantTool {
            tool: "ena".to_string(),
            message: format!("GET {url} returned {}", resp.status()),
        });
    }

    // Ensure parent dir exists even when the caller passed a nested
    // target path.
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut out = File::create(dest)?;
    let mut hasher = Md5::new();
    let mut buf = [0u8; 65_536];

    loop {
        let n = resp.read(&mut buf).map_err(|e| MycoNoteError::QuantTool {
            tool: "ena".to_string(),
            message: format!("read {url}: {e}"),
        })?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
        out.write_all(&buf[..n])?;
    }
    out.flush()?;

    let got = format!("{:x}", hasher.finalize());
    if !got.eq_ignore_ascii_case(expected_md5) {
        let _ = std::fs::remove_file(dest);
        return Err(MycoNoteError::QuantTool {
            tool: "ena".to_string(),
            message: format!(
                "md5 mismatch for {url}: expected {expected_md5}, got {got}. \
                 Partial download removed."
            ),
        });
    }

    Ok(())
}

/// Compute MD5 of a file. Handy for verifying already-downloaded
/// files against an ENA checksum without re-fetching.
pub fn md5_file(path: &Path) -> Result<String> {
    let mut f = File::open(path)?;
    let mut hasher = Md5::new();
    let mut buf = [0u8; 65_536];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    // Real ENA filereport responses, captured on 2026-04-24. Kept
    // inline so the parser is tested against exact byte-for-byte
    // output rather than a paraphrased fixture.

    const PAIRED_TSV: &str = "\
run_accession\tfastq_ftp\tfastq_md5\tlibrary_layout
SRR453566\tftp.sra.ebi.ac.uk/vol1/fastq/SRR453/SRR453566/SRR453566_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR453/SRR453566/SRR453566_2.fastq.gz\t807c55b7badde864a573745b4543dc35;0728a92e8dae2d8d5699c3f283a2af60\tPAIRED
";

    const SINGLE_TSV: &str = "\
run_accession\tfastq_ftp\tfastq_md5\tlibrary_layout
SRR7082700\tftp.sra.ebi.ac.uk/vol1/fastq/SRR708/000/SRR7082700/SRR7082700.fastq.gz\t4a6db949b3b88f3216a1c91d6e40c927\tSINGLE
";

    const STUDY_MULTI_TSV: &str = "\
run_accession\tfastq_ftp\tfastq_md5\tlibrary_layout
SRR7082700\tftp.sra.ebi.ac.uk/vol1/fastq/SRR708/000/SRR7082700/SRR7082700.fastq.gz\t4a6db949b3b88f3216a1c91d6e40c927\tSINGLE
SRR7082694\tftp.sra.ebi.ac.uk/vol1/fastq/SRR708/004/SRR7082694/SRR7082694.fastq.gz\te71aa37f8cb417ea213b2d449cce048d\tSINGLE
SRR7082695\tftp.sra.ebi.ac.uk/vol1/fastq/SRR708/005/SRR7082695/SRR7082695.fastq.gz\t7ebf89c5a62896af17cd7e3134c44fbe\tSINGLE
";

    #[test]
    fn parse_paired_row() {
        let rows = parse_filereport_tsv(PAIRED_TSV).unwrap();
        assert_eq!(rows.len(), 1);
        let r = &rows[0];
        assert_eq!(r.run_accession, "SRR453566");
        assert!(r.paired);
        assert_eq!(r.fastq_urls.len(), 2);
        assert_eq!(r.fastq_md5s.len(), 2);
        assert!(r.has_r2());
    }

    #[test]
    fn parse_single_end_row() {
        let rows = parse_filereport_tsv(SINGLE_TSV).unwrap();
        assert_eq!(rows.len(), 1);
        let r = &rows[0];
        assert_eq!(r.run_accession, "SRR7082700");
        assert!(!r.paired);
        assert_eq!(r.fastq_urls.len(), 1);
        assert!(!r.has_r2());
    }

    #[test]
    fn parse_study_expands_to_multiple_rows() {
        let rows = parse_filereport_tsv(STUDY_MULTI_TSV).unwrap();
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].run_accession, "SRR7082700");
        assert_eq!(rows[1].run_accession, "SRR7082694");
        assert_eq!(rows[2].run_accession, "SRR7082695");
    }

    #[test]
    fn parse_rejects_missing_header_column() {
        let bad = "run_accession\tother\nSRR1\tx\n";
        let err = parse_filereport_tsv(bad).unwrap_err();
        assert!(format!("{err}").contains("missing 'fastq_ftp'"));
    }

    #[test]
    fn parse_rejects_url_md5_mismatch() {
        let bad = "run_accession\tfastq_ftp\tfastq_md5\tlibrary_layout\n\
                   SRR1\ta.fq.gz;b.fq.gz\tmd5a\tPAIRED\n";
        let err = parse_filereport_tsv(bad).unwrap_err();
        assert!(format!("{err}").contains("2 urls but fastq_md5 has 1"));
    }

    #[test]
    fn parse_rejects_empty_md5() {
        let bad = "run_accession\tfastq_ftp\tfastq_md5\tlibrary_layout\n\
                   SRR1\ta.fq.gz\t\tSINGLE\n";
        let err = parse_filereport_tsv(bad).unwrap_err();
        assert!(format!("{err}").contains("md5"));
    }

    #[test]
    fn parse_skips_rows_with_empty_fastq_ftp() {
        let tsv = "run_accession\tfastq_ftp\tfastq_md5\tlibrary_layout\n\
                   SRR1\t\t\tSINGLE\n\
                   SRR2\tftp.sra.ebi.ac.uk/real.fq.gz\tabc123\tSINGLE\n";
        let rows = parse_filereport_tsv(tsv).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].run_accession, "SRR2");
    }

    #[test]
    fn parse_empty_tsv_errors() {
        let err = parse_filereport_tsv("").unwrap_err();
        assert!(format!("{err}").contains("empty"));
    }

    #[test]
    fn md5_file_computes_known_value() {
        let tmp = TempDir::new().unwrap();
        let p = tmp.path().join("hello.txt");
        std::fs::write(&p, b"hello").unwrap();
        // MD5("hello") — precomputed.
        assert_eq!(md5_file(&p).unwrap(), "5d41402abc4b2a76b9719d911017c592");
    }

    // ── Network-dependent tests (ignored unless explicitly run) ──────────────

    #[test]
    #[ignore]
    fn fetch_filereport_live_single_run() {
        // Small, public S. cerevisiae run. ~100 MB FASTQ total, but we
        // don't download anything here — just the 1-row TSV response.
        let rows = fetch_filereport("SRR453566").expect("ENA reachable");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].run_accession, "SRR453566");
        assert!(rows[0].paired);
    }

    #[test]
    #[ignore]
    fn fetch_filereport_live_unknown_accession_errors() {
        // ENA returns an empty body for unknown accessions; our parser
        // surfaces that as the "empty" error.
        let err = fetch_filereport("XXXXXXXXXX").unwrap_err();
        assert!(format!("{err}").contains("empty") || format!("{err}").contains("returned"));
    }
}
