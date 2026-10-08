//! Structural pre-validation and repair of a GFF3 feature graph before it is
//! handed to `table2asn`.
//!
//! `table2asn` (and NCBI's downstream validators) reject GFF3s whose feature
//! hierarchy is broken — the classic failures being a CDS or exon whose
//! `Parent` points at an mRNA that was never written (funannotate issue #290,
//! "CDS not in mRNA"), and the related "N features reference missing parents"
//! that MycoNote's own benchmark `submit` hit (5 such features). Those errors
//! surface only after a full table2asn run, with a cryptic message.
//!
//! This module walks the feature graph first and:
//!   * reports every structural problem with its feature type and location;
//!   * optionally (`--fix-structure`) repairs the safe, unambiguous cases by
//!     synthesising the missing `gene` / `mRNA` parents that the orphaned
//!     children need — **without ever changing a child's coordinates**.
//!
//! The conservative default is report-and-continue: a problem is printed but
//! the pipeline proceeds (table2asn then gives its own verdict). `--fix-structure`
//! switches on the repair pass and writes a corrected GFF3 that is used for the
//! rest of the submission.

use crate::parser::gff::GFFRecord;
use std::collections::{HashMap, HashSet};

/// A single structural problem found in the feature graph.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuralIssue {
    pub kind: IssueKind,
    /// GFF3 feature type the problem is attached to (e.g. "CDS", "mRNA").
    pub feature_type: String,
    /// `seqid:start-end` of the offending feature.
    pub location: String,
    /// Human-readable detail, including the IDs involved.
    pub detail: String,
    /// Whether `--fix-structure` can repair this case automatically.
    pub repairable: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IssueKind {
    /// A `gene` / `mRNA` feature carries no `ID` attribute.
    MissingId,
    /// A feature's `Parent` references an ID that does not exist.
    OrphanParent,
    /// A CDS / exon / mRNA that should have a `Parent` has none.
    MissingParent,
    /// An mRNA's `Parent` resolves to a feature that is not a `gene`.
    MrnaParentNotGene,
    /// A CDS interval falls outside the bounds of its parent mRNA.
    CdsOutsideMrna,
}

impl IssueKind {
    fn label(&self) -> &'static str {
        match self {
            IssueKind::MissingId => "missing ID",
            IssueKind::OrphanParent => "orphan Parent (references missing ID)",
            IssueKind::MissingParent => "missing Parent",
            IssueKind::MrnaParentNotGene => "mRNA Parent is not a gene",
            IssueKind::CdsOutsideMrna => "CDS outside mRNA bounds",
        }
    }
}

#[derive(Debug, Default)]
pub struct StructuralReport {
    pub issues: Vec<StructuralIssue>,
}

impl StructuralReport {
    pub fn is_clean(&self) -> bool {
        self.issues.is_empty()
    }

    /// Number of issues `--fix-structure` could repair automatically.
    pub fn repairable_count(&self) -> usize {
        self.issues.iter().filter(|i| i.repairable).count()
    }

    pub fn print_summary(&self) {
        if self.is_clean() {
            println!("  Structure: OK (feature hierarchy is valid)");
            return;
        }
        println!(
            "  Structure: {} problem(s) found ({} auto-repairable with --fix-structure)",
            self.issues.len(),
            self.repairable_count()
        );
        for issue in self.issues.iter().take(20) {
            println!(
                "    {} [{}] {} — {}",
                issue.feature_type,
                issue.location,
                issue.kind.label(),
                issue.detail
            );
        }
        if self.issues.len() > 20 {
            println!("    ... and {} more", self.issues.len() - 20);
        }
    }
}

fn loc(rec: &GFFRecord) -> String {
    format!("{}:{}-{}", rec.seqid, rec.start, rec.end)
}

/// First Parent ID (a GFF3 `Parent` may be a comma-separated list).
fn first_parent(rec: &GFFRecord) -> Option<&str> {
    rec.attributes
        .get("Parent")
        .map(|p| p.split(',').next().unwrap_or(p.as_str()))
}

fn is_child_type(t: &str) -> bool {
    matches!(t, "CDS" | "exon" | "five_prime_UTR" | "three_prime_UTR")
}

/// Walk the feature graph and report every structural problem, without
/// changing anything.
pub fn check_structure(records: &[GFFRecord]) -> StructuralReport {
    let mut report = StructuralReport::default();

    // All declared IDs, and ID → feature_type.
    let mut id_type: HashMap<&str, &str> = HashMap::new();
    for rec in records {
        if let Some(id) = rec.attributes.get("ID") {
            id_type.insert(id.as_str(), rec.feature_type.as_str());
        }
    }
    // mRNA ID → (start, end) for the CDS-containment check.
    let mut mrna_bounds: HashMap<&str, (u64, u64)> = HashMap::new();
    for rec in records {
        if rec.feature_type == "mRNA" {
            if let Some(id) = rec.attributes.get("ID") {
                mrna_bounds.insert(id.as_str(), (rec.start, rec.end));
            }
        }
    }

    for rec in records {
        let t = rec.feature_type.as_str();

        // Required IDs on gene / mRNA.
        if matches!(t, "gene" | "mRNA") && !rec.attributes.contains_key("ID") {
            report.issues.push(StructuralIssue {
                kind: IssueKind::MissingId,
                feature_type: t.to_string(),
                location: loc(rec),
                detail: format!("{t} feature has no ID attribute"),
                repairable: true,
            });
        }

        match first_parent(rec) {
            Some(pid) => {
                // Parent must resolve to an existing ID.
                if !id_type.contains_key(pid) {
                    report.issues.push(StructuralIssue {
                        kind: IssueKind::OrphanParent,
                        feature_type: t.to_string(),
                        location: loc(rec),
                        detail: format!("Parent '{pid}' does not exist in this GFF3"),
                        // CDS/exon/mRNA orphans can be re-parented safely.
                        repairable: is_child_type(t) || t == "mRNA",
                    });
                } else if t == "mRNA" && id_type.get(pid) != Some(&"gene") {
                    // mRNA's parent exists but is not a gene.
                    report.issues.push(StructuralIssue {
                        kind: IssueKind::MrnaParentNotGene,
                        feature_type: t.to_string(),
                        location: loc(rec),
                        detail: format!(
                            "Parent '{pid}' is a '{}', expected a gene",
                            id_type.get(pid).copied().unwrap_or("?")
                        ),
                        repairable: false,
                    });
                }

                // CDS must lie within its parent mRNA.
                if t == "CDS" {
                    if let Some(&(ms, me)) = mrna_bounds.get(pid) {
                        if rec.start < ms || rec.end > me {
                            report.issues.push(StructuralIssue {
                                kind: IssueKind::CdsOutsideMrna,
                                feature_type: t.to_string(),
                                location: loc(rec),
                                detail: format!(
                                    "CDS {}-{} is outside mRNA '{pid}' bounds {ms}-{me}",
                                    rec.start, rec.end
                                ),
                                // Never auto-"fix" by moving coordinates.
                                repairable: false,
                            });
                        }
                    }
                }
            }
            None => {
                // Children and mRNAs are expected to carry a Parent.
                if is_child_type(t) || t == "mRNA" {
                    report.issues.push(StructuralIssue {
                        kind: IssueKind::MissingParent,
                        feature_type: t.to_string(),
                        location: loc(rec),
                        detail: format!("{t} feature has no Parent attribute"),
                        repairable: true,
                    });
                }
            }
        }
    }

    report
}

/// Outcome of a repair pass: the rewritten record set plus a log of what was
/// synthesised or changed (one line per action, for the report and README).
#[derive(Debug, Default)]
pub struct RepairLog {
    pub actions: Vec<String>,
}

/// Repair the safe, unambiguous structural problems and return the corrected
/// record set. Synthesised `gene`/`mRNA` parents span exactly their children,
/// so **no child coordinate is ever altered**. Cases that cannot be repaired
/// without guessing (CDS outside its mRNA, an mRNA parented to a non-gene) are
/// left untouched and remain in the returned `check_structure` report for the
/// user to resolve by hand.
///
/// The synthesised parents are prepended (genes, then mRNAs) so the output is
/// parent-before-child, which keeps IGV/JBrowse and tabix-style readers happy.
pub fn repair_structure(records: Vec<GFFRecord>) -> (Vec<GFFRecord>, RepairLog) {
    let mut log = RepairLog::default();
    let mut records = records;

    // ── 1. Give any ID-less gene/mRNA a deterministic synthetic ID ───────────
    let mut synth_id = 0usize;
    for rec in records.iter_mut() {
        if matches!(rec.feature_type.as_str(), "gene" | "mRNA")
            && !rec.attributes.contains_key("ID")
        {
            synth_id += 1;
            let id = format!(
                "myco_fix_{}_{:04}",
                rec.feature_type.to_lowercase(),
                synth_id
            );
            log.actions.push(format!(
                "assigned ID '{id}' to {} at {}",
                rec.feature_type,
                loc(rec)
            ));
            rec.attributes.insert("ID".to_string(), id);
        }
    }

    // Current ID set and ID → feature_type.
    let id_set: HashSet<String> = records
        .iter()
        .filter_map(|r| r.attributes.get("ID").cloned())
        .collect();

    // ── 2. Assign a Parent to parentless CDS/exon and parentless mRNA ────────
    // A parentless child cannot be grouped with siblings reliably, so each is
    // given its own synthesised parent chain. This is lossless for coordinates
    // (the new parents span exactly the child) and only over-splits the rare
    // genuinely multi-exon gene that arrived with no Parent links at all.
    let mut parentless_counter = 0usize;
    for rec in records.iter_mut() {
        let t = rec.feature_type.clone();
        let needs_parent = (is_child_type(&t) || t == "mRNA")
            && rec.attributes.get("Parent").map_or(true, |p| p.is_empty());
        if needs_parent {
            parentless_counter += 1;
            let synth_parent = format!("myco_fix_orphan_{parentless_counter:04}");
            log.actions.push(format!(
                "parented {} at {} to synthesised '{synth_parent}'",
                t,
                loc(rec)
            ));
            rec.attributes.insert("Parent".to_string(), synth_parent);
        }
    }

    // ── 3. Synthesise the missing parents that children now reference ────────
    // Recompute after step 2 so the newly-invented Parent IDs are included.
    let mut children_by_missing_parent: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, rec) in records.iter().enumerate() {
        if let Some(pid) = rec.attributes.get("Parent") {
            let pid = pid.split(',').next().unwrap_or(pid).to_string();
            if !id_set.contains(&pid) {
                children_by_missing_parent.entry(pid).or_default().push(idx);
            }
        }
    }

    // Deterministic order for reproducible output.
    let mut missing_parents: Vec<String> = children_by_missing_parent.keys().cloned().collect();
    missing_parents.sort();

    let mut new_genes: Vec<GFFRecord> = Vec::new();
    let mut new_mrnas: Vec<GFFRecord> = Vec::new();
    let mut gene_counter = 0usize;

    for pid in missing_parents {
        let child_idxs = &children_by_missing_parent[&pid];
        // Span = min start .. max end over the children; template fields taken
        // from the first child so seqid / source / strand match.
        let first = &records[child_idxs[0]];
        let seqid = first.seqid.clone();
        let source = first.source.clone();
        let strand = first.strand;
        let span_start = child_idxs.iter().map(|&i| records[i].start).min().unwrap();
        let span_end = child_idxs.iter().map(|&i| records[i].end).max().unwrap();

        // Are the children mRNAs (→ synthesise a gene) or CDS/exon (→ an mRNA)?
        let children_are_transcripts = child_idxs
            .iter()
            .any(|&i| records[i].feature_type == "mRNA");

        if children_are_transcripts {
            // The missing parent is a gene.
            let mut gene =
                make_feature(&seqid, &source, "gene", span_start, span_end, strand, &pid);
            gene.attributes.insert(
                "Note".to_string(),
                "synthesised_by_myconote_fix_structure".to_string(),
            );
            log.actions.push(format!(
                "synthesised gene '{pid}' ({seqid}:{span_start}-{span_end}) for {} mRNA child(ren)",
                child_idxs.len()
            ));
            new_genes.push(gene);
        } else {
            // The missing parent is an mRNA; it also needs a gene above it.
            gene_counter += 1;
            let gene_id = format!("{pid}.gene");
            let gene_id = if id_set.contains(&gene_id) {
                format!("myco_fix_gene_{gene_counter:04}")
            } else {
                gene_id
            };
            let mut gene = make_feature(
                &seqid, &source, "gene", span_start, span_end, strand, &gene_id,
            );
            gene.attributes.insert(
                "Note".to_string(),
                "synthesised_by_myconote_fix_structure".to_string(),
            );

            let mut mrna =
                make_feature(&seqid, &source, "mRNA", span_start, span_end, strand, &pid);
            mrna.attributes
                .insert("Parent".to_string(), gene_id.clone());
            mrna.attributes.insert(
                "Note".to_string(),
                "synthesised_by_myconote_fix_structure".to_string(),
            );

            log.actions.push(format!(
                "synthesised gene '{gene_id}' + mRNA '{pid}' ({seqid}:{span_start}-{span_end}) for {} CDS/exon child(ren)",
                child_idxs.len()
            ));
            new_genes.push(gene);
            new_mrnas.push(mrna);
        }
    }

    // Parent-before-child ordering: synthesised genes, synthesised mRNAs, then
    // the (possibly Parent-rewritten) original records in their original order.
    let mut out = Vec::with_capacity(new_genes.len() + new_mrnas.len() + records.len());
    out.append(&mut new_genes);
    out.append(&mut new_mrnas);
    out.extend(records);

    (out, log)
}

/// Build a minimal GFF3 feature with the given coordinates and ID.
fn make_feature(
    seqid: &str,
    source: &str,
    feature_type: &str,
    start: u64,
    end: u64,
    strand: char,
    id: &str,
) -> GFFRecord {
    let mut attributes = HashMap::new();
    attributes.insert("ID".to_string(), id.to_string());
    GFFRecord {
        seqid: seqid.to_string(),
        source: if source.is_empty() {
            "myconote".to_string()
        } else {
            source.to_string()
        },
        feature_type: feature_type.to_string(),
        start,
        end,
        score: None,
        strand,
        phase: None,
        attributes,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::gff::GFFRecord;

    fn rec(line: &str) -> GFFRecord {
        GFFRecord::from_line(line, 1).unwrap()
    }

    #[test]
    fn clean_gff_passes_unchanged() {
        let records = vec![
            rec("c1\tmyco\tgene\t100\t900\t.\t+\t.\tID=g1"),
            rec("c1\tmyco\tmRNA\t100\t900\t.\t+\t.\tID=g1.t1;Parent=g1"),
            rec("c1\tmyco\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=g1.t1"),
            rec("c1\tmyco\tCDS\t500\t900\t.\t+\t0\tID=cds2;Parent=g1.t1"),
        ];
        let report = check_structure(&records);
        assert!(report.is_clean(), "unexpected issues: {:?}", report.issues);

        // Repair on a clean set must be a no-op (same count, same IDs present).
        let n = records.len();
        let (fixed, log) = repair_structure(records);
        assert_eq!(fixed.len(), n);
        assert!(
            log.actions.is_empty(),
            "no repairs expected: {:?}",
            log.actions
        );
    }

    #[test]
    fn orphan_cds_missing_parent_is_detected() {
        // CDS points at mRNA 'g1.t1' that was never written.
        let records = vec![
            rec("c1\tmyco\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=g1.t1"),
            rec("c1\tmyco\tCDS\t500\t900\t.\t+\t0\tID=cds2;Parent=g1.t1"),
        ];
        let report = check_structure(&records);
        assert!(!report.is_clean());
        assert_eq!(report.issues.len(), 2);
        assert!(report
            .issues
            .iter()
            .all(|i| i.kind == IssueKind::OrphanParent && i.repairable));
    }

    #[test]
    fn orphan_cds_is_repaired_into_valid_parented_structure() {
        let cds1 = "c1\tmyco\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=g1.t1";
        let cds2 = "c1\tmyco\tCDS\t500\t900\t.\t+\t0\tID=cds2;Parent=g1.t1";
        let records = vec![rec(cds1), rec(cds2)];

        let (fixed, log) = repair_structure(records);
        // A gene and an mRNA were synthesised.
        assert_eq!(
            log.actions.len(),
            1,
            "one synthesis action: {:?}",
            log.actions
        );
        let genes: Vec<_> = fixed.iter().filter(|r| r.feature_type == "gene").collect();
        let mrnas: Vec<_> = fixed.iter().filter(|r| r.feature_type == "mRNA").collect();
        assert_eq!(genes.len(), 1);
        assert_eq!(mrnas.len(), 1);

        // The mRNA has the id the CDS referenced, and spans both CDS exactly.
        let mrna = mrnas[0];
        assert_eq!(mrna.attributes.get("ID").unwrap(), "g1.t1");
        assert_eq!((mrna.start, mrna.end), (100, 900));
        // The mRNA is parented to the synthesised gene.
        let gene = genes[0];
        assert_eq!(
            mrna.attributes.get("Parent").unwrap(),
            gene.attributes.get("ID").unwrap()
        );
        // Gene spans the same extent; CDS coordinates are untouched.
        assert_eq!((gene.start, gene.end), (100, 900));
        let cds_coords: Vec<(u64, u64)> = fixed
            .iter()
            .filter(|r| r.feature_type == "CDS")
            .map(|r| (r.start, r.end))
            .collect();
        assert_eq!(cds_coords, vec![(100, 300), (500, 900)]);

        // The repaired set now validates clean.
        let report = check_structure(&fixed);
        assert!(
            report.is_clean(),
            "still broken after repair: {:?}",
            report.issues
        );
    }

    #[test]
    fn parentless_cds_gets_a_chain() {
        let records = vec![rec("c1\tmyco\tCDS\t10\t99\t.\t-\t0\tID=loneCDS")];
        let report = check_structure(&records);
        assert_eq!(report.issues.len(), 1);
        assert_eq!(report.issues[0].kind, IssueKind::MissingParent);

        let (fixed, _log) = repair_structure(records);
        assert!(check_structure(&fixed).is_clean());
        assert_eq!(fixed.iter().filter(|r| r.feature_type == "gene").count(), 1);
        assert_eq!(fixed.iter().filter(|r| r.feature_type == "mRNA").count(), 1);
    }

    #[test]
    fn cds_outside_mrna_is_flagged_but_not_repaired() {
        let records = vec![
            rec("c1\tmyco\tgene\t100\t400\t.\t+\t.\tID=g1"),
            rec("c1\tmyco\tmRNA\t100\t400\t.\t+\t.\tID=g1.t1;Parent=g1"),
            // CDS runs past the mRNA end (400).
            rec("c1\tmyco\tCDS\t100\t900\t.\t+\t0\tID=cds1;Parent=g1.t1"),
        ];
        let report = check_structure(&records);
        let outside: Vec<_> = report
            .issues
            .iter()
            .filter(|i| i.kind == IssueKind::CdsOutsideMrna)
            .collect();
        assert_eq!(outside.len(), 1);
        assert!(!outside[0].repairable);

        // Repair must not touch the CDS coordinates.
        let (fixed, _log) = repair_structure(records);
        let cds = fixed.iter().find(|r| r.feature_type == "CDS").unwrap();
        assert_eq!((cds.start, cds.end), (100, 900));
    }
}
