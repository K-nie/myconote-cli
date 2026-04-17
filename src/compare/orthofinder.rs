//! OrthoFinder wrapper — runs the `orthofinder` binary on a staged directory
//! of per-species protein FASTAs and parses the relevant outputs back into
//! MycoNote-native data structures.
//!
//! We pick OrthoFinder over rolling our own ortholog clustering because it's
//! the research-grade standard (Emms & Kelly 2019, Genome Biology), handles
//! paralog separation via gene-tree inference, and its output shape is well
//! documented and stable. The trade-off is adding a conda/pip dependency —
//! but since the MycoNote pipeline already requires BUSCO, hmmscan, IQ-TREE,
//! and minimap2, one more external tool is a fair price for trustworthy
//! ortholog output.

use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

/// One row of the OrthoFinder `Orthogroups.tsv` file after parsing: a cluster
/// identifier plus the list of gene IDs contributed by each input genome
/// (keyed by genome name).
#[derive(Debug, Clone, Default)]
pub struct Orthogroup {
    pub id: String,
    /// genome_name → gene IDs present in this orthogroup for that genome
    pub members: HashMap<String, Vec<String>>,
}

impl Orthogroup {
    /// Number of genomes contributing at least one gene to this orthogroup.
    pub fn n_genomes(&self) -> usize {
        self.members.values().filter(|v| !v.is_empty()).count()
    }

    /// True iff every genome contributes exactly one gene (the set used for
    /// species-tree inference).
    pub fn is_single_copy(&self, expected_genomes: usize) -> bool {
        self.members.len() == expected_genomes
            && self.members.values().all(|v| v.len() == 1)
    }
}

/// Aggregate pan-genome shape using Tettelin-style bins.
#[derive(Debug, Default, Clone)]
pub struct PangenomeSummary {
    pub total_clusters: usize,
    pub core: usize,
    pub soft_core: usize,
    pub shell: usize,
    pub cloud: usize,
    pub singletons: usize,
    pub single_copy_core: usize,
}

/// Check whether OrthoFinder is on PATH.
pub fn orthofinder_available() -> bool {
    which::which("orthofinder").is_ok()
}

/// Invoke OrthoFinder on `proteins_dir`. Returns the path to the
/// `Results_<date>` directory it creates inside `proteins_dir`.
///
/// OrthoFinder writes its results into a sibling directory of the input and
/// tags it with the current date, so we scan for it after the run and hand
/// back the newest match rather than trying to predict the exact name.
pub fn run_orthofinder(
    proteins_dir: &Path,
    threads: usize,
    sensitive: bool,
    msa: bool,
) -> Result<PathBuf> {
    if !orthofinder_available() {
        return Err(MycoNoteError::ExternalTool(
            "orthofinder not found in PATH. Install with:\n  \
             conda install -c bioconda orthofinder\n\
             See https://github.com/davidemms/OrthoFinder for details."
                .to_string(),
        ));
    }

    let mut cmd = Command::new("orthofinder");
    cmd.arg("-f").arg(proteins_dir);
    cmd.arg("-t").arg(threads.to_string());
    cmd.arg("-a").arg(threads.to_string()); // analysis threads
    cmd.arg("-S").arg(if sensitive {
        "diamond_ultra_sens"
    } else {
        "diamond"
    });
    // Tree inference method. OrthoFinder's default is `-M msa` which
    // requires FAMSA/MAFFT + FastTree. FAMSA v2.4.1 segfaults on macOS
    // arm64 from the bioconda channel, so default to `dendroblast` (which
    // uses DIAMOND distances directly — no MSA step). `dendroblast` is
    // slightly less accurate but substantially faster and avoids the
    // platform-specific crash. Pass `--msa` to opt back into the MSA path
    // when FAMSA/MAFFT are known-working.
    if msa {
        cmd.arg("-M").arg("msa");
    } else {
        cmd.arg("-M").arg("dendroblast");
    }
    // We only want the standard orthogroup outputs; skipping per-species gene
    // trees (`-og` flag is implicit when we don't pass `-y`) would speed
    // things up further but we keep tree inference on since the species tree
    // is a free byproduct users typically want.

    println!("  Invoking OrthoFinder (this may take a while)…");
    let status = cmd
        .status()
        .map_err(|e| MycoNoteError::ExternalTool(format!("orthofinder: {}", e)))?;

    if !status.success() {
        return Err(MycoNoteError::ExternalTool(
            "OrthoFinder exited with non-zero status".to_string(),
        ));
    }

    // OrthoFinder writes proteins_dir/OrthoFinder/Results_<DateTag>/ — pick
    // the newest match so repeated runs don't pick up a stale earlier one.
    let of_parent = proteins_dir.join("OrthoFinder");
    if !of_parent.is_dir() {
        return Err(MycoNoteError::ExternalTool(format!(
            "OrthoFinder finished but no output directory found at {}",
            of_parent.display()
        )));
    }
    let mut newest: Option<(PathBuf, std::time::SystemTime)> = None;
    for entry in std::fs::read_dir(&of_parent).map_err(MycoNoteError::Io)? {
        let entry = entry.map_err(MycoNoteError::Io)?;
        let path = entry.path();
        if !path.is_dir() { continue; }
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if !name.starts_with("Results_") { continue; }
        let mtime = entry
            .metadata()
            .and_then(|m| m.modified())
            .unwrap_or(std::time::SystemTime::UNIX_EPOCH);
        match newest {
            Some((_, t)) if t >= mtime => {}
            _ => newest = Some((path, mtime)),
        }
    }
    newest
        .map(|(p, _)| p)
        .ok_or_else(|| {
            MycoNoteError::ExternalTool(format!(
                "OrthoFinder output directory not found under {}",
                of_parent.display()
            ))
        })
}

/// Parse `Orthogroups/Orthogroups.tsv`. Header row lists genome names; each
/// subsequent row is `orthogroup_id\tgenome1_gene_list\tgenome2_gene_list...`
/// where each gene list is comma-separated (or empty when absent).
pub fn parse_orthogroups_tsv(tsv: &Path) -> Result<Vec<Orthogroup>> {
    let f = std::fs::File::open(tsv).map_err(MycoNoteError::Io)?;
    let reader = BufReader::new(f);
    let mut lines = reader.lines();

    let header = match lines.next() {
        Some(Ok(h)) => h,
        _ => return Err(MycoNoteError::InvalidFormat(
            "Empty Orthogroups.tsv".to_string(),
        )),
    };
    let cols: Vec<String> = header.split('\t').map(|s| s.to_string()).collect();
    if cols.len() < 2 {
        return Err(MycoNoteError::InvalidFormat(
            "Orthogroups.tsv header too short".to_string(),
        ));
    }
    let genome_names: Vec<String> = cols[1..].to_vec();

    let mut orthogroups = Vec::new();
    for line in lines {
        let line = line.map_err(MycoNoteError::Io)?;
        if line.trim().is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() { continue; }

        let mut og = Orthogroup {
            id: fields[0].to_string(),
            members: HashMap::new(),
        };
        for (i, col) in fields.iter().enumerate().skip(1) {
            if i - 1 >= genome_names.len() { break; }
            let genome = &genome_names[i - 1];
            let genes: Vec<String> = col
                .split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();
            og.members.insert(genome.clone(), genes);
        }
        orthogroups.push(og);
    }
    Ok(orthogroups)
}

/// Summarise orthogroups into Tettelin-style core/soft-core/shell/cloud bins.
/// `soft_core_frac` sets the threshold (default 0.95) for the soft-core cut;
/// `cloud_frac` (default 0.15) bounds the cloud from the shell.
pub fn summarise_pangenome(
    orthogroups: &[Orthogroup],
    total_genomes: usize,
    soft_core_frac: f64,
    cloud_frac: f64,
) -> PangenomeSummary {
    let mut s = PangenomeSummary::default();
    s.total_clusters = orthogroups.len();
    let soft_core_n =
        ((total_genomes as f64) * soft_core_frac).ceil() as usize;
    let cloud_n = ((total_genomes as f64) * cloud_frac).ceil() as usize;

    for og in orthogroups {
        let n = og.n_genomes();
        if n == total_genomes {
            s.core += 1;
            if og.is_single_copy(total_genomes) {
                s.single_copy_core += 1;
            }
        } else if n >= soft_core_n {
            s.soft_core += 1;
        } else if n <= cloud_n {
            if n == 1 {
                s.singletons += 1;
            }
            s.cloud += 1;
        } else {
            s.shell += 1;
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_tmp(label: &str, body: &str) -> PathBuf {
        use std::io::Write as _;
        let mut p = std::env::temp_dir();
        p.push(format!("myconote_of_test_{}_{}.tsv", std::process::id(), label));
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        p
    }

    #[test]
    fn parse_typical_orthogroups_tsv() {
        // Simulated OrthoFinder output: 3 genomes, 4 orthogroups covering
        // the usual category spread (core, shell, singleton, paralog pair).
        let body = "\
Orthogroup\tgenomeA\tgenomeB\tgenomeC
OG0000000\tgA_1\tgB_1\tgC_1
OG0000001\tgA_2\tgB_2, gB_2b\t
OG0000002\tgA_3\t\t
OG0000003\t\tgB_3\tgC_3
";
        let p = write_tmp("ok", body);
        let ogs = parse_orthogroups_tsv(&p).unwrap();
        assert_eq!(ogs.len(), 4);

        // OG0000000 = core single-copy
        assert_eq!(ogs[0].members.get("genomeA").unwrap().len(), 1);
        assert_eq!(ogs[0].n_genomes(), 3);
        assert!(ogs[0].is_single_copy(3));

        // OG0000001 = 2 genomes, one has a paralog pair (size 2)
        assert_eq!(ogs[1].n_genomes(), 2);
        assert_eq!(ogs[1].members.get("genomeB").unwrap().len(), 2);
        assert!(!ogs[1].is_single_copy(3));

        // OG0000002 = singleton (just genomeA)
        assert_eq!(ogs[2].n_genomes(), 1);

        // OG0000003 = two genomes, neither is genomeA
        assert_eq!(ogs[3].n_genomes(), 2);
        assert!(ogs[3].members.get("genomeA").unwrap().is_empty());

        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn pangenome_summary_classifies_correctly() {
        // Build 8 orthogroups across 10 synthetic genomes to hit each bin.
        let mut ogs = Vec::new();
        for (idx, n) in [10, 10, 10, 10, 9, 5, 2, 1].iter().enumerate() {
            let mut og = Orthogroup {
                id: format!("OG{:07}", idx),
                members: HashMap::new(),
            };
            for i in 0..*n {
                og.members
                    .insert(format!("genome{}", i), vec![format!("g{}_{}", i, idx)]);
            }
            ogs.push(og);
        }
        let s = summarise_pangenome(&ogs, 10, 0.95, 0.15);
        assert_eq!(s.total_clusters, 8);
        // 4 orthogroups at full prevalence → core.
        assert_eq!(s.core, 4);
        assert_eq!(s.single_copy_core, 4);
        // 1 at 9/10 = 90% (below 95% soft-core threshold) → shell, not soft-core.
        // Soft-core needs ceil(10 * 0.95) = 10 genomes, so 9/10 falls to shell.
        assert_eq!(s.soft_core, 0);
        assert_eq!(s.shell, 2); // 9-genome and 5-genome orthogroups
        // Cloud bound = ceil(10 * 0.15) = 2. 2/10 and 1/10 are cloud.
        assert_eq!(s.cloud, 2);
        assert_eq!(s.singletons, 1);
    }
}
