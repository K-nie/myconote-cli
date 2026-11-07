use super::AnnotateConfig;
/// Pfam domain annotation via HMMER
///
/// Searches predicted proteins against the Pfam-A HMM database using either:
///   - hmmsearch (preferred): searches each HMM profile against all proteins
///     at once. Much faster for large protein sets (5,000+) because HMMER's
///     internal SSV/MSV acceleration works on the full sequence database.
///   - hmmscan (fallback): searches each protein against all HMM profiles.
///     Slower but available when hmmsearch can't be used.
///
/// Performance note (HMMER docs):
///   "hmmsearch is generally faster than hmmscan for searching a large
///    sequence database with a profile database, because the overhead of
///    loading each profile is amortized across all sequences."
///
/// For a typical fungal genome (5,000 proteins × 20,795 Pfam-A profiles):
///   hmmscan:   ~45 min (searches 5K seqs × 20K profiles one-by-one)
///   hmmsearch: ~8 min  (loads 20K profiles × searches all 5K seqs at once)
///
/// Install: conda install -c bioconda hmmer
/// DB setup: myconote annotate --download-dbs   (downloads Pfam-A.hmm)
use crate::utils::error::{MycoNoteError, Result};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─────────────────────────────────────────────────────────────────────────────
// Hit record
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct PfamHit {
    pub query_id: String,
    pub domain_id: String,   // e.g. "PF00001.21"
    pub domain_name: String, // e.g. "7tm_1"
    pub evalue: f64,
    pub bitscore: f64,
    pub ali_from: u64,
    pub ali_to: u64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

pub fn run(
    query_fa: &Path,
    pfam_hmm: &Path,
    output_tsv: &Path,
    config: &AnnotateConfig,
) -> Result<Vec<PfamHit>> {
    // Press the database if .h3i index is missing
    let h3i = pfam_hmm.with_extension("hmm.h3i");
    if !h3i.exists() {
        press_pfam(pfam_hmm)?;
    }

    // Try hmmsearch first (much faster), fall back to hmmscan
    let hmmsearch_available = which::which("hmmsearch").is_ok();
    let hmmscan_available = which::which("hmmscan").is_ok();

    if !hmmsearch_available && !hmmscan_available {
        return Err(MycoNoteError::UnsupportedFormat(
            "Neither hmmsearch nor hmmscan found in PATH.\n\
             Install with: conda install -c bioconda hmmer"
                .to_string(),
        ));
    }

    let domtblout = output_tsv.with_extension("domtblout");

    let hits = if hmmsearch_available {
        run_hmmsearch(query_fa, pfam_hmm, &domtblout, config)?
    } else {
        run_hmmscan(query_fa, pfam_hmm, &domtblout, config)?
    };

    // Write a clean TSV summary
    write_pfam_tsv(output_tsv, &hits)?;

    Ok(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// hmmsearch (preferred — faster for large protein sets)
// ─────────────────────────────────────────────────────────────────────────────

fn run_hmmsearch(
    query_fa: &Path,
    pfam_hmm: &Path,
    domtblout: &Path,
    config: &AnnotateConfig,
) -> Result<Vec<PfamHit>> {
    let hmmsearch = which::which("hmmsearch")
        .map_err(|_| MycoNoteError::UnsupportedFormat("hmmsearch not found".to_string()))?;

    let (n_shards, cpu_per_shard) = shard_plan(config.threads);
    let evalue_s = config.evalue.to_string();
    let query_s = query_fa.to_string_lossy().into_owned();

    if n_shards <= 1 {
        println!("  Running hmmsearch (Pfam-A) — faster than hmmscan for large proteomes");
        run_one_hmmsearch(
            &hmmsearch,
            pfam_hmm,
            &query_s,
            domtblout,
            cpu_per_shard,
            &evalue_s,
        )?;
    } else {
        // `hmmsearch --cpu N` stops scaling after a few threads (one master
        // thread feeds the workers), so split the PROFILE database into shards
        // and run one hmmsearch per shard in parallel. Splitting profiles, not
        // proteins, keeps every E-value identical to a single run: per-sequence
        // E-values use the number of target proteins (unchanged in every shard)
        // and per-domain E-values depend only on that profile's own hits.
        let fallback_dir = domtblout.parent().unwrap_or_else(|| Path::new("."));
        let shards = match shard_hmm_db(pfam_hmm, n_shards, fallback_dir) {
            Ok(shards) => shards,
            Err(e) => {
                eprintln!("  ⚠  Could not shard Pfam database ({e}); running one hmmsearch.");
                run_one_hmmsearch(
                    &hmmsearch,
                    pfam_hmm,
                    &query_s,
                    domtblout,
                    config.threads,
                    &evalue_s,
                )?;
                return finish_hmmsearch(domtblout);
            }
        };
        println!(
            "  Running hmmsearch (Pfam-A) in {} parallel shards × {} CPU",
            shards.len(),
            cpu_per_shard
        );
        let shard_outs: Vec<PathBuf> = (0..shards.len())
            .map(|i| domtblout.with_extension(format!("shard{i}.domtblout")))
            .collect();
        let results: Vec<Result<()>> = std::thread::scope(|scope| {
            let handles: Vec<_> = shards
                .iter()
                .zip(&shard_outs)
                .map(|(shard, out)| {
                    let (hmmsearch, query_s, evalue_s) = (&hmmsearch, &query_s, &evalue_s);
                    scope.spawn(move || {
                        run_one_hmmsearch(hmmsearch, shard, query_s, out, cpu_per_shard, evalue_s)
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| {
                    h.join().unwrap_or_else(|_| {
                        Err(MycoNoteError::InvalidFormat(
                            "hmmsearch shard thread panicked".to_string(),
                        ))
                    })
                })
                .collect()
        });
        for r in results {
            r?;
        }
        concat_files(&shard_outs, domtblout)?;
        for out in &shard_outs {
            let _ = std::fs::remove_file(out);
        }
    }

    finish_hmmsearch(domtblout)
}

fn finish_hmmsearch(domtblout: &Path) -> Result<Vec<PfamHit>> {
    // hmmsearch --domtblout has the same format as hmmscan --domtblout
    // but columns 0-2 and 3-5 are swapped:
    //   hmmscan:   target=domain, query=protein
    //   hmmsearch: target=protein, query=domain
    let hits = parse_domtblout_hmmsearch(domtblout)?;
    println!("  hmmsearch finished — {} domain hits", hits.len());

    Ok(hits)
}

/// Run one hmmsearch process over `hmm` (a full database or one shard).
fn run_one_hmmsearch(
    hmmsearch: &Path,
    hmm: &Path,
    query: &str,
    domtblout: &Path,
    cpu: usize,
    evalue: &str,
) -> Result<()> {
    let hmm_s = hmm.to_string_lossy();
    let dom_s = domtblout.to_string_lossy();
    let cpu_s = cpu.to_string();
    let status = Command::new(hmmsearch)
        .args([
            "--domtblout",
            dom_s.as_ref(),
            "--cpu",
            &cpu_s,
            "-E",
            evalue,
            "--domE",
            evalue,
            "--noali", // skip alignment output (saves time + disk)
            hmm_s.as_ref(),
            query,
        ])
        .stdout(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;
    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "hmmsearch failed on {}. Check that the Pfam-A.hmm file is intact.",
            hmm.display()
        )));
    }
    Ok(())
}

/// Below this many threads one hmmsearch process with `--cpu` is about as
/// fast as sharding, so keep the simple path.
const MIN_THREADS_TO_SHARD: usize = 4;
/// Upper bound on parallel shards (each re-reads its slice of Pfam-A).
const MAX_SHARDS: usize = 32;

/// How to spread `threads` over hmmsearch processes: (shards, --cpu each).
fn shard_plan(threads: usize) -> (usize, usize) {
    let threads = threads.max(1);
    if threads < MIN_THREADS_TO_SHARD {
        return (1, threads);
    }
    let shards = threads.min(MAX_SHARDS);
    (shards, (threads / shards).max(1))
}

/// Split an HMMER3 text database into `n` shards, round-robin by profile so
/// each shard gets a similar mix of profile sizes. Shards are cached beside
/// the database (`<db>.shards<n>/`) and reused while the database's size and
/// mtime are unchanged; if that directory is not writable they go under
/// `fallback_dir` instead.
fn shard_hmm_db(hmm: &Path, n: usize, fallback_dir: &Path) -> Result<Vec<PathBuf>> {
    let meta = std::fs::metadata(hmm).map_err(MycoNoteError::Io)?;
    let mtime = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let stamp = format!("{} {}\n", meta.len(), mtime);
    let file_name = hmm
        .file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_else(|| "pfam.hmm".to_string());
    let dir_name = format!("{file_name}.shards{n}");

    let beside = hmm.parent().map(|p| p.join(&dir_name));
    let candidates: Vec<PathBuf> = beside
        .into_iter()
        .chain(std::iter::once(fallback_dir.join(&dir_name)))
        .collect();

    let shard_paths = |dir: &Path| -> Vec<PathBuf> {
        (0..n).map(|i| dir.join(format!("shard{i}.hmm"))).collect()
    };

    // Reuse a valid cache from any candidate location.
    for dir in &candidates {
        let paths = shard_paths(dir);
        let stamp_ok = std::fs::read_to_string(dir.join("STAMP"))
            .map(|s| s == stamp)
            .unwrap_or(false);
        if stamp_ok && paths.iter().all(|p| p.is_file()) {
            return Ok(paths);
        }
    }

    let mut last_err = None;
    for dir in &candidates {
        match write_shards(hmm, dir, &shard_paths(dir), &stamp) {
            Ok(paths) => return Ok(paths),
            Err(e) => last_err = Some(e),
        }
    }
    Err(last_err
        .unwrap_or_else(|| MycoNoteError::InvalidFormat("could not write Pfam shards".to_string())))
}

fn write_shards(hmm: &Path, dir: &Path, paths: &[PathBuf], stamp: &str) -> Result<Vec<PathBuf>> {
    use std::io::{BufRead, BufReader, BufWriter};
    std::fs::create_dir_all(dir).map_err(MycoNoteError::Io)?;
    // Invalidate first so a half-written cache is never reused.
    let _ = std::fs::remove_file(dir.join("STAMP"));
    println!(
        "  Splitting {} into {} shards (one-time)…",
        hmm.display(),
        paths.len()
    );

    let mut writers = paths
        .iter()
        .map(|p| std::fs::File::create(p).map(BufWriter::new))
        .collect::<std::io::Result<Vec<_>>>()
        .map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(std::fs::File::open(hmm).map_err(MycoNoteError::Io)?);
    let mut current = 0usize;
    let mut profiles = 0usize;
    for line in reader.lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        writeln!(writers[current], "{line}").map_err(MycoNoteError::Io)?;
        // Each HMMER3 profile record ends with a line that is exactly "//".
        if line.trim_end() == "//" {
            profiles += 1;
            current = (current + 1) % writers.len();
        }
    }
    for w in &mut writers {
        w.flush().map_err(MycoNoteError::Io)?;
    }
    if profiles < paths.len() {
        return Err(MycoNoteError::InvalidFormat(format!(
            "{} has only {} profiles; too few to split into {} shards",
            hmm.display(),
            profiles,
            paths.len()
        )));
    }
    std::fs::write(dir.join("STAMP"), stamp).map_err(MycoNoteError::Io)?;
    Ok(paths.to_vec())
}

fn concat_files(parts: &[PathBuf], out: &Path) -> Result<()> {
    let mut f = std::fs::File::create(out).map_err(MycoNoteError::Io)?;
    for p in parts {
        let mut src = std::fs::File::open(p).map_err(MycoNoteError::Io)?;
        std::io::copy(&mut src, &mut f).map_err(MycoNoteError::Io)?;
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// hmmscan (fallback — slower but works everywhere)
// ─────────────────────────────────────────────────────────────────────────────

fn run_hmmscan(
    query_fa: &Path,
    pfam_hmm: &Path,
    domtblout: &Path,
    config: &AnnotateConfig,
) -> Result<Vec<PfamHit>> {
    let hmmscan = which::which("hmmscan")
        .map_err(|_| MycoNoteError::UnsupportedFormat("hmmscan not found".to_string()))?;

    println!("  Running hmmscan (Pfam-A) — fallback mode");

    let pfam_s = pfam_hmm.to_string_lossy();
    let query_s = query_fa.to_string_lossy();
    let dom_s = domtblout.to_string_lossy();

    let status = Command::new(&hmmscan)
        .args([
            "--domtblout",
            dom_s.as_ref(),
            "--cpu",
            &config.threads.to_string(),
            "-E",
            &config.evalue.to_string(),
            "--domE",
            &config.evalue.to_string(),
            "--noali",
            pfam_s.as_ref(),
            query_s.as_ref(),
        ])
        .stdout(std::process::Stdio::null())
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat(
            "hmmscan failed. Check that Pfam-A.hmm is properly pressed (hmmpress).".to_string(),
        ));
    }

    let hits = parse_domtblout_hmmscan(domtblout)?;
    println!("  hmmscan finished — {} domain hits", hits.len());

    Ok(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// HMM pressing helper
// ─────────────────────────────────────────────────────────────────────────────

fn press_pfam(pfam_hmm: &Path) -> Result<()> {
    let hmmpress = which::which("hmmpress").map_err(|_| {
        MycoNoteError::UnsupportedFormat("hmmpress not found. Install hmmer.".to_string())
    })?;

    let pfam_s = pfam_hmm.to_string_lossy();
    println!("  Pressing Pfam-A.hmm (one-time setup)...");
    let status = Command::new(&hmmpress)
        .args(["-f", pfam_s.as_ref()])
        .status()
        .map_err(MycoNoteError::Io)?;

    if !status.success() {
        return Err(MycoNoteError::InvalidFormat("hmmpress failed.".to_string()));
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Parser for hmmsearch --domtblout format
// ─────────────────────────────────────────────────────────────────────────────
// hmmsearch domtblout columns (space-delimited):
//   0: target name (PROTEIN)    1: target accession  2: target length
//   3: query name (DOMAIN)      4: query accession   5: query length
//   6: full E-value   7: full score   8: full bias
//   9: dom#  10: ndom
//  11: dom c-Evalue  12: dom i-Evalue  13: dom score  14: dom bias
//  15: hmm from  16: hmm to  17: ali from  18: ali to  19: env from  20: env to
//  21: acc  22+: target description

fn parse_domtblout_hmmsearch(path: &Path) -> Result<Vec<PfamHit>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut hits: Vec<PfamHit> = Vec::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 22 {
            continue;
        }

        // In hmmsearch: col 0 = protein (target), col 3 = domain (query)
        let query_id = f[0].to_string(); // protein name
        let domain_name = f[3].to_string(); // domain name
        let domain_acc = f[4].to_string(); // domain accession (PFxxxxx.xx)
        let evalue: f64 = f[12].parse().unwrap_or(f64::MAX); // domain i-Evalue
        let bitscore: f64 = f[13].parse().unwrap_or(0.0); // domain score
        let ali_from: u64 = f[17].parse().unwrap_or(0);
        let ali_to: u64 = f[18].parse().unwrap_or(0);

        // Normalize domain ID to PFxxxxx.xx format
        let domain_id = normalize_pfam_acc(&domain_acc, &domain_name);

        hits.push(PfamHit {
            query_id,
            domain_id,
            domain_name,
            evalue,
            bitscore,
            ali_from,
            ali_to,
        });
    }

    dedup_hits(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// Parser for hmmscan --domtblout format
// ─────────────────────────────────────────────────────────────────────────────
// hmmscan domtblout columns (space-delimited):
//   0: target name (DOMAIN)     1: target accession  2: target length
//   3: query name (PROTEIN)     4: query accession   5: query length
//   ...same structure after that

fn parse_domtblout_hmmscan(path: &Path) -> Result<Vec<PfamHit>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut hits: Vec<PfamHit> = Vec::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 22 {
            continue;
        }

        // In hmmscan: col 0 = domain (target), col 3 = protein (query)
        let domain_name = f[0].to_string();
        let domain_acc = f[1].to_string();
        let query_id = f[3].to_string();
        let evalue: f64 = f[12].parse().unwrap_or(f64::MAX);
        let bitscore: f64 = f[13].parse().unwrap_or(0.0);
        let ali_from: u64 = f[17].parse().unwrap_or(0);
        let ali_to: u64 = f[18].parse().unwrap_or(0);

        let domain_id = normalize_pfam_acc(&domain_acc, &domain_name);

        hits.push(PfamHit {
            query_id,
            domain_id,
            domain_name,
            evalue,
            bitscore,
            ali_from,
            ali_to,
        });
    }

    dedup_hits(hits)
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Normalize a Pfam accession to "PFxxxxx.xx" format.
fn normalize_pfam_acc(acc: &str, name: &str) -> String {
    if acc.starts_with("PF") {
        acc.to_string()
    } else if name.starts_with("PF") {
        name.to_string()
    } else {
        // Use the accession as-is (might be a Pfam-B or clan ID)
        acc.to_string()
    }
}

/// Deduplicate: keep best hit per (query, domain) pair by lowest e-value.
fn dedup_hits(mut hits: Vec<PfamHit>) -> Result<Vec<PfamHit>> {
    hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut seen = std::collections::HashSet::new();
    hits.retain(|h| seen.insert((h.query_id.clone(), h.domain_id.clone())));
    Ok(hits)
}

fn write_pfam_tsv(path: &Path, hits: &[PfamHit]) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(
        f,
        "query_id\tdomain_id\tdomain_name\tevalue\tbitscore\tali_from\tali_to"
    )
    .map_err(MycoNoteError::Io)?;
    for h in hits {
        writeln!(
            f,
            "{}\t{}\t{}\t{:.2e}\t{:.1}\t{}\t{}",
            h.query_id, h.domain_id, h.domain_name, h.evalue, h.bitscore, h.ali_from, h.ali_to
        )
        .map_err(MycoNoteError::Io)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fake_hmm(n: usize) -> String {
        (0..n)
            .map(|i| format!("HMMER3/f [3.1b2 | February 2015]\nNAME  dom{i}\nACC   PF{i:05}.1\nHMM ...\n//\n"))
            .collect()
    }

    #[test]
    fn shard_plan_keeps_small_runs_single() {
        assert_eq!(shard_plan(0), (1, 1));
        assert_eq!(shard_plan(2), (1, 2));
        assert_eq!(shard_plan(8), (8, 1));
        assert_eq!(shard_plan(64), (32, 2));
    }

    #[test]
    fn shards_partition_every_profile_round_robin() {
        let dir = tempfile::tempdir().unwrap();
        let db = dir.path().join("Pfam-A.hmm");
        std::fs::write(&db, fake_hmm(10)).unwrap();
        let shards = shard_hmm_db(&db, 4, dir.path()).unwrap();
        assert_eq!(shards.len(), 4);

        let mut names: Vec<String> = Vec::new();
        let mut counts = Vec::new();
        for s in &shards {
            let text = std::fs::read_to_string(s).unwrap();
            counts.push(text.matches("\n//\n").count());
            names.extend(
                text.lines()
                    .filter_map(|l| l.strip_prefix("NAME  ").map(String::from)),
            );
            assert!(text.starts_with("HMMER3/f"));
        }
        assert_eq!(counts, vec![3, 3, 2, 2]);
        names.sort_by_key(|n| n[3..].parse::<usize>().unwrap());
        let expected: Vec<String> = (0..10).map(|i| format!("dom{i}")).collect();
        assert_eq!(names, expected);
    }

    #[test]
    fn shard_cache_is_reused_then_rebuilt_when_db_changes() {
        let dir = tempfile::tempdir().unwrap();
        let db = dir.path().join("Pfam-A.hmm");
        std::fs::write(&db, fake_hmm(6)).unwrap();
        let first = shard_hmm_db(&db, 2, dir.path()).unwrap();
        // Tamper with a shard: a reused cache keeps the tampered content.
        std::fs::write(&first[0], "tampered").unwrap();
        let again = shard_hmm_db(&db, 2, dir.path()).unwrap();
        assert_eq!(std::fs::read_to_string(&again[0]).unwrap(), "tampered");
        // A different-sized database invalidates the stamp and rebuilds.
        std::fs::write(&db, fake_hmm(8)).unwrap();
        let rebuilt = shard_hmm_db(&db, 2, dir.path()).unwrap();
        assert!(std::fs::read_to_string(&rebuilt[0])
            .unwrap()
            .contains("NAME  dom0"));
    }

    #[test]
    fn too_few_profiles_to_shard_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let db = dir.path().join("tiny.hmm");
        std::fs::write(&db, fake_hmm(2)).unwrap();
        assert!(shard_hmm_db(&db, 4, dir.path()).is_err());
    }
}
