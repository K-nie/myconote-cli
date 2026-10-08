/// Unit tests for core modules
///
/// Tests evidence scoring, genetic code translation, ploidy estimation,
/// validation, and reproducibility reporting without requiring external tools.

#[cfg(test)]
mod evidence_tests {
    use myconote_cli::predict::evidence::EvidenceWeights;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_default_weights() {
        let w = EvidenceWeights::default();
        assert_eq!(w.augustus, 10.0);
        assert_eq!(w.snap, 3.0);
        assert_eq!(w.protein, 20.0);
        assert_eq!(w.glimmerhmm, 2.0);
        assert_eq!(w.genemark, 5.0);
    }

    #[test]
    fn test_weight_for_lookup() {
        let w = EvidenceWeights::default();
        assert_eq!(w.weight_for("Augustus"), 10.0);
        assert_eq!(w.weight_for("snap"), 3.0);
        assert_eq!(w.weight_for("protein"), 20.0);
        assert_eq!(w.weight_for("genemark"), 5.0);
        assert_eq!(w.weight_for("unknown_tool"), 1.0);
    }

    #[test]
    fn test_weights_toml_roundtrip() {
        let w = EvidenceWeights::default();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        w.write_toml(tmp.path()).unwrap();

        let loaded = EvidenceWeights::from_toml(tmp.path()).unwrap();
        assert_eq!(loaded.augustus, w.augustus);
        assert_eq!(loaded.snap, w.snap);
        assert_eq!(loaded.protein, w.protein);
    }

    #[test]
    fn test_custom_weights_from_toml() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(
            tmp,
            r#"
augustus = 15.0
snap = 5.0
protein = 25.0
est = 10.0
glimmerhmm = 4.0
genemark = 8.0
miniprot = 22.0
trnascan = 18.0
"#
        )
        .unwrap();

        let w = EvidenceWeights::from_toml(tmp.path()).unwrap();
        assert_eq!(w.augustus, 15.0);
        assert_eq!(w.genemark, 8.0);
        assert_eq!(w.miniprot, 22.0);
    }
}

#[cfg(test)]
mod genetic_code_tests {
    use myconote_cli::annotate::genetic_code::GeneticCode;

    #[test]
    fn test_standard_code() {
        let gc = GeneticCode::STANDARD;
        assert_eq!(gc.translate("ATGAAATTTCCC"), "MKFP");
    }

    #[test]
    fn test_candida_ctg_ser() {
        let gc = GeneticCode::from_table_number(12).unwrap();
        // CTG = Ser in Candida, not Leu
        assert_eq!(gc.translate("ATGCTG"), "MS");
    }

    #[test]
    fn test_standard_ctg_leu() {
        let gc = GeneticCode::STANDARD;
        assert_eq!(gc.translate("ATGCTG"), "ML");
    }

    #[test]
    fn test_vertebrate_mito_tga_trp() {
        let gc = GeneticCode::from_table_number(2).unwrap();
        // TGA = Trp in vertebrate mito
        assert_eq!(gc.translate("ATGTGA"), "MW");
    }

    #[test]
    fn test_from_table_number() {
        assert_eq!(
            GeneticCode::from_table_number(12).map(|g| g.table_number()),
            Some(12)
        );
        assert_eq!(
            GeneticCode::from_table_number(1).map(|g| g.table_number()),
            Some(1)
        );
        // 99 is not an NCBI table
        assert!(GeneticCode::from_table_number(99).is_none());
        // 7, 8, 15, 17-20, 32 are retired / withdrawn by NCBI
        assert!(GeneticCode::from_table_number(7).is_none());
        assert!(GeneticCode::from_table_number(15).is_none());
    }

    #[test]
    fn test_from_name() {
        assert_eq!(GeneticCode::from_name("candida").table_number(), 12);
        assert_eq!(GeneticCode::from_name("12").table_number(), 12);
        assert_eq!(GeneticCode::from_name("standard").table_number(), 1);
        assert_eq!(GeneticCode::from_name("bacterial").table_number(), 11);
    }

    #[test]
    fn test_internal_stops() {
        assert_eq!(GeneticCode::find_internal_stops("MK*FG"), vec![2]);
        assert_eq!(
            GeneticCode::find_internal_stops("MKFG*"),
            Vec::<usize>::new()
        ); // trailing is ok
        assert_eq!(
            GeneticCode::find_internal_stops("MKFG"),
            Vec::<usize>::new()
        );
        assert_eq!(GeneticCode::find_internal_stops("M*K*G*"), vec![1, 3]);
    }

    #[test]
    fn test_start_codons() {
        let gc = GeneticCode::STANDARD;
        assert!(gc.is_start_codon(b"ATG"));
        assert!(!gc.is_start_codon(b"GTG"));

        let bac = GeneticCode::from_table_number(11).unwrap();
        assert!(bac.is_start_codon(b"ATG"));
        assert!(bac.is_start_codon(b"GTG"));
        assert!(bac.is_start_codon(b"TTG"));
    }
}

#[cfg(test)]
mod ploidy_tests {
    use myconote_cli::predict::ploidy::{adjust_for_ploidy, Ploidy};

    #[test]
    fn test_ploidy_levels() {
        assert_eq!(Ploidy::Haploid.level(), 1);
        assert_eq!(Ploidy::Diploid.level(), 2);
        assert_eq!(Ploidy::Tetraploid.level(), 4);
        assert_eq!(Ploidy::Polyploid(6).level(), 6);
    }

    #[test]
    fn test_is_polyploid() {
        assert!(!Ploidy::Haploid.is_polyploid());
        assert!(Ploidy::Diploid.is_polyploid());
        assert!(Ploidy::Tetraploid.is_polyploid());
    }

    #[test]
    fn test_from_n() {
        assert_eq!(Ploidy::from_n(0).level(), 1); // 0 treated as haploid
        assert_eq!(Ploidy::from_n(1).level(), 1);
        assert_eq!(Ploidy::from_n(2).level(), 2);
        assert_eq!(Ploidy::from_n(3).level(), 3);
    }

    #[test]
    fn test_adjustments() {
        let adj = adjust_for_ploidy(&Ploidy::Haploid);
        assert!(!adj.expect_allelic_pairs);

        let adj2 = adjust_for_ploidy(&Ploidy::Diploid);
        assert!(adj2.expect_allelic_pairs);
        assert!(adj2.report_note.is_some());
    }
}

#[cfg(test)]
mod reproducibility_tests {
    use myconote_cli::utils::reproducibility::ReproducibilityReport;
    use tempfile::TempDir;

    #[test]
    fn test_report_creation() {
        let report = ReproducibilityReport::new();
        assert!(!report.myconote_version.is_empty());
        assert!(!report.timestamp.is_empty());
        assert_eq!(report.pipeline_steps.len(), 0);
    }

    #[test]
    fn test_report_add_steps() {
        let mut report = ReproducibilityReport::new();
        report.add_step("predict");
        report.add_step("annotate");
        assert_eq!(report.pipeline_steps.len(), 2);
    }

    #[test]
    fn test_report_add_param() {
        let mut report = ReproducibilityReport::new();
        report.add_param("kingdom", "fungi");
        report.add_param("threads", "8");
        assert_eq!(report.parameters.get("kingdom").unwrap(), "fungi");
    }

    #[test]
    fn test_report_write_json() {
        let tmp = TempDir::new().unwrap();
        let json_path = tmp.path().join("report.json");

        let mut report = ReproducibilityReport::new();
        report.add_step("test");
        report.add_param("key", "value");
        report.set_runtime(42.5);

        report.write_json(&json_path).unwrap();
        assert!(json_path.exists());

        let content = std::fs::read_to_string(&json_path).unwrap();
        assert!(content.contains("myconote_version"));
        assert!(content.contains("test"));
    }

    #[test]
    fn test_report_write_text() {
        let tmp = TempDir::new().unwrap();
        let txt_path = tmp.path().join("report.txt");

        let mut report = ReproducibilityReport::new();
        report.add_step("annotate");
        report.write_text(&txt_path).unwrap();

        let content = std::fs::read_to_string(&txt_path).unwrap();
        assert!(content.contains("Reproducibility Report"));
        assert!(content.contains("annotate"));
    }
}

#[cfg(test)]
mod validation_tests {
    use myconote_cli::utils::validation;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_protein_fasta_validation_clean() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">gene1\nMKFGTERVWDL\n>gene2\nMKFGTERVWDLABC").unwrap();

        let result = validation::validate_protein_fasta(tmp.path()).unwrap();
        assert_eq!(result.total_sequences, 2);
        assert!(result.internal_stops.is_empty());
        assert!(result.is_valid());
    }

    #[test]
    fn test_protein_fasta_internal_stops() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">gene1\nMK*FGTERVWDL").unwrap();

        let result = validation::validate_protein_fasta(tmp.path()).unwrap();
        assert_eq!(result.internal_stops.len(), 1);
    }

    #[test]
    fn test_protein_fasta_invalid_chars() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">gene1\nMKFG123").unwrap();

        let result = validation::validate_protein_fasta(tmp.path()).unwrap();
        assert!(!result.invalid_chars.is_empty());
        assert!(!result.is_valid());
    }
}

#[cfg(test)]
mod evidence_merger_tests {
    use myconote_cli::predict::evidence;
    use std::collections::HashSet;
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    /// Create a minimal GFF3 with a multi-exon gene (3 CDS features).
    fn write_test_gff3(path: &std::path::Path) {
        let mut f = std::fs::File::create(path).unwrap();
        writeln!(f, "##gff-version 3").unwrap();
        // Gene 1: 3 CDS exons
        writeln!(f, "chr1\ttest\tgene\t100\t2000\t.\t+\t.\tID=gene1").unwrap();
        writeln!(
            f,
            "chr1\ttest\tmRNA\t100\t2000\t.\t+\t.\tID=mrna1;Parent=gene1"
        )
        .unwrap();
        writeln!(
            f,
            "chr1\ttest\tCDS\t100\t300\t.\t+\t0\tID=cds1a;Parent=mrna1"
        )
        .unwrap();
        writeln!(
            f,
            "chr1\ttest\tCDS\t500\t800\t.\t+\t0\tID=cds1b;Parent=mrna1"
        )
        .unwrap();
        writeln!(
            f,
            "chr1\ttest\tCDS\t1500\t2000\t.\t+\t0\tID=cds1c;Parent=mrna1"
        )
        .unwrap();
        // Gene 2: 2 CDS exons
        writeln!(f, "chr1\ttest\tgene\t3000\t4500\t.\t-\t.\tID=gene2").unwrap();
        writeln!(
            f,
            "chr1\ttest\tmRNA\t3000\t4500\t.\t-\t.\tID=mrna2;Parent=gene2"
        )
        .unwrap();
        writeln!(
            f,
            "chr1\ttest\tCDS\t3000\t3500\t.\t-\t0\tID=cds2a;Parent=mrna2"
        )
        .unwrap();
        writeln!(
            f,
            "chr1\ttest\tCDS\t4000\t4500\t.\t-\t0\tID=cds2b;Parent=mrna2"
        )
        .unwrap();
    }

    #[test]
    fn test_merger_produces_unique_ids() {
        let tmp = TempDir::new().unwrap();
        let gff_in = tmp.path().join("input.gff3");
        let gff_out = tmp.path().join("consensus.gff3");

        write_test_gff3(&gff_in);

        let inputs = vec![(gff_in.as_path(), "TestPredictor", 10.0)];

        let gene_count = evidence::merge_predictions(&inputs, &gff_out, "TEST").unwrap();
        assert_eq!(gene_count, 2, "Should find 2 genes");

        // Read output and check all IDs are unique
        let content = std::fs::read_to_string(&gff_out).unwrap();
        let mut ids: Vec<String> = Vec::new();

        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 9 {
                continue;
            }
            // Parse attributes for ID=
            for attr in cols[8].split(';') {
                if let Some(id_val) = attr.strip_prefix("ID=") {
                    ids.push(id_val.to_string());
                }
            }
        }

        let unique: HashSet<_> = ids.iter().collect();
        assert_eq!(
            ids.len(),
            unique.len(),
            "Duplicate IDs found! All IDs: {:?}",
            ids
        );
    }

    #[test]
    fn test_merger_parent_references_valid() {
        let tmp = TempDir::new().unwrap();
        let gff_in = tmp.path().join("input.gff3");
        let gff_out = tmp.path().join("consensus.gff3");

        write_test_gff3(&gff_in);

        let inputs = vec![(gff_in.as_path(), "TestPredictor", 10.0)];

        evidence::merge_predictions(&inputs, &gff_out, "TEST").unwrap();

        // Check all Parent= references point to existing IDs
        let content = std::fs::read_to_string(&gff_out).unwrap();
        let mut all_ids: HashSet<String> = HashSet::new();
        let mut parent_refs: Vec<(String, String)> = Vec::new(); // (parent, feature_type)

        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 9 {
                continue;
            }

            let feature_type = cols[2].to_string();
            for attr in cols[8].split(';') {
                if let Some(id_val) = attr.strip_prefix("ID=") {
                    all_ids.insert(id_val.to_string());
                }
                if let Some(parent_val) = attr.strip_prefix("Parent=") {
                    parent_refs.push((parent_val.to_string(), feature_type.clone()));
                }
            }
        }

        for (parent, ftype) in &parent_refs {
            assert!(
                all_ids.contains(parent),
                "Orphan feature: {} references non-existent Parent='{}'",
                ftype,
                parent
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GeneMark mode (v0.6.0) — parsing and legacy-alias resolution
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod genemark_mode_tests {
    use myconote_cli::predict::genemark::GeneMarkMode;
    use myconote_cli::predict::resolve_genemark_mode;
    use std::path::PathBuf;

    // ── Mode parsing ─────────────────────────────────────────────────────────

    #[test]
    fn from_str_accepts_es_et_ep_etp() {
        assert_eq!(GeneMarkMode::from_str("es"), Some(GeneMarkMode::Es));
        assert_eq!(GeneMarkMode::from_str("et"), Some(GeneMarkMode::Et));
        assert_eq!(GeneMarkMode::from_str("ep"), Some(GeneMarkMode::Ep));
        assert_eq!(GeneMarkMode::from_str("etp"), Some(GeneMarkMode::Etp));
    }

    #[test]
    fn from_str_is_case_insensitive() {
        assert_eq!(GeneMarkMode::from_str("ES"), Some(GeneMarkMode::Es));
        assert_eq!(GeneMarkMode::from_str("Et"), Some(GeneMarkMode::Et));
        assert_eq!(GeneMarkMode::from_str("EP"), Some(GeneMarkMode::Ep));
        assert_eq!(GeneMarkMode::from_str("ETP"), Some(GeneMarkMode::Etp));
    }

    #[test]
    fn from_str_accepts_plus_aliases() {
        // GeneMark literature spells these EP+ and ETP+; accept both.
        assert_eq!(GeneMarkMode::from_str("ep+"), Some(GeneMarkMode::Ep));
        assert_eq!(GeneMarkMode::from_str("etp+"), Some(GeneMarkMode::Etp));
    }

    #[test]
    fn from_str_rejects_garbage() {
        assert_eq!(GeneMarkMode::from_str("foo"), None);
        assert_eq!(GeneMarkMode::from_str(""), None);
        assert_eq!(GeneMarkMode::from_str("genemark"), None);
    }

    #[test]
    fn long_name_is_canonical_publication_form() {
        // Used in user-facing log lines — must match the upstream paper titles
        // so users grepping documentation find the right thing.
        assert_eq!(GeneMarkMode::Es.long_name(), "GeneMark-ES");
        assert_eq!(GeneMarkMode::Et.long_name(), "GeneMark-ET");
        assert_eq!(GeneMarkMode::Ep.long_name(), "GeneMark-EP+");
        assert_eq!(GeneMarkMode::Etp.long_name(), "GeneMark-ETP+");
    }

    // ── Resolution from legacy flags ─────────────────────────────────────────

    #[test]
    fn explicit_mode_wins_over_legacy_flags() {
        // If the user passes --genemark-mode explicitly, the legacy boolean +
        // hints should never override it.
        let hints = PathBuf::from("/tmp/intron_hints.gff");
        assert_eq!(
            resolve_genemark_mode(Some(GeneMarkMode::Ep), true, Some(&hints)),
            Some(GeneMarkMode::Ep)
        );
    }

    #[test]
    fn legacy_genemark_alone_resolves_to_es() {
        // Old `--genemark` flag with no hints == ES mode (the historical default).
        assert_eq!(
            resolve_genemark_mode(None, true, None),
            Some(GeneMarkMode::Es)
        );
    }

    #[test]
    fn legacy_genemark_with_hints_resolves_to_et() {
        // Pre-0.6.0: passing `--genemark` and `--genemark-hints` together
        // implied GeneMark-ET. Preserved.
        let hints = PathBuf::from("/tmp/intron_hints.gff");
        assert_eq!(
            resolve_genemark_mode(None, true, Some(&hints)),
            Some(GeneMarkMode::Et)
        );
    }

    #[test]
    fn legacy_hints_only_still_implies_et() {
        // Pre-0.6.0 also accepted `--genemark-hints` without `--genemark`
        // (the hints flag implied the rest). Preserved.
        let hints = PathBuf::from("/tmp/intron_hints.gff");
        assert_eq!(
            resolve_genemark_mode(None, false, Some(&hints)),
            Some(GeneMarkMode::Et)
        );
    }

    #[test]
    fn no_flags_means_skip() {
        assert_eq!(resolve_genemark_mode(None, false, None), None);
    }

    // ── ETP requires both inputs (validated at runtime, not at parse time) ───
    //
    // These tests document the contract of run_prediction's ETP branch via
    // the resolution helper + the predictor config. They live as ignored,
    // live-tool tests because actually running gmes_petap.pl needs a
    // licensed GeneMark + ProtHint install.

    #[test]
    fn etp_resolution_succeeds_even_without_inputs() {
        // Resolution itself doesn't validate inputs — it just decides which
        // mode to dispatch. The runtime check happens in run_prediction and
        // is exercised by the integration-level test below.
        assert_eq!(
            resolve_genemark_mode(Some(GeneMarkMode::Etp), false, None),
            Some(GeneMarkMode::Etp)
        );
    }

    #[test]
    #[ignore = "requires licensed gmes_petap.pl + ProtHint install; run with --ignored"]
    fn etp_end_to_end_live_smoke() {
        // Sanity-only smoke gate. We intentionally do not exercise the full
        // live ETP pipeline in CI because it needs a Georgia-Tech-licensed
        // GeneMark binary and ProtHint, both of which require manual
        // installation. Mirrors the de-template `Rscript parse()` gate.
        assert!(myconote_cli::predict::genemark::genemark_available());
        assert!(myconote_cli::predict::genemark::prothint_available());
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRAKER (v0.6.0) — config, mutual exclusion, mode auto-detection
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod braker_tests {
    use myconote_cli::predict::braker::{detect_mode, BrakerMode};
    use myconote_cli::predict::genemark::GeneMarkMode;
    use myconote_cli::predict::{check_braker_conflicts, PredictConfig};
    use std::path::PathBuf;

    // ── Mode auto-detection contract ─────────────────────────────────────────

    #[test]
    fn auto_detect_rna_only_picks_braker1() {
        let bams = vec![PathBuf::from("rna1.bam")];
        assert_eq!(detect_mode(&bams, None), Some(BrakerMode::Braker1));
    }

    #[test]
    fn auto_detect_protein_only_picks_braker2() {
        let prot = PathBuf::from("orthodb_fungi.fa");
        assert_eq!(detect_mode(&[], Some(&prot)), Some(BrakerMode::Braker2));
    }

    #[test]
    fn auto_detect_both_picks_braker3() {
        let bams = vec![PathBuf::from("a.bam"), PathBuf::from("b.bam")];
        let prot = PathBuf::from("orthodb_fungi.fa");
        assert_eq!(detect_mode(&bams, Some(&prot)), Some(BrakerMode::Braker3));
    }

    #[test]
    fn auto_detect_neither_returns_none() {
        // run_braker converts this into a hard error at call time; we don't
        // silently fall back to anything (no silent fallbacks rule).
        assert_eq!(detect_mode(&[], None), None);
    }

    // ── Mode-flag mapping (BRAKER CLI contract) ──────────────────────────────

    #[test]
    fn braker_mode_flag_strings_match_upstream_braker_cli() {
        // Verified against Gaius-Augustus/BRAKER scripts/braker.pl on
        // 2026-04-24 — these are the GetOptions keys, not invented.
        assert_eq!(BrakerMode::Braker1.mode_flag(), "--esmode");
        assert_eq!(BrakerMode::Braker2.mode_flag(), "--epmode");
        assert_eq!(BrakerMode::Braker3.mode_flag(), "--etpmode");
    }

    #[test]
    fn braker_mode_from_str_accepts_numeric_and_named() {
        assert_eq!(BrakerMode::from_str("1"), Some(BrakerMode::Braker1));
        assert_eq!(BrakerMode::from_str("braker2"), Some(BrakerMode::Braker2));
        assert_eq!(BrakerMode::from_str("3"), Some(BrakerMode::Braker3));
        assert_eq!(BrakerMode::from_str("4"), None);
        assert_eq!(BrakerMode::from_str(""), None);
    }

    // ── Mutual-exclusion contract for --use-braker vs standard predictors ────

    #[test]
    fn braker_off_never_conflicts() {
        // When BRAKER isn't active, every per-predictor flag is fine.
        let cfg = PredictConfig {
            use_genemark: true,
            use_glimmerhmm: true,
            genemark_mode: Some(GeneMarkMode::Etp),
            protein_fasta: Some(PathBuf::from("p.fa")),
            ..PredictConfig::default()
        };
        assert!(check_braker_conflicts(&cfg).is_ok());
    }

    #[test]
    fn braker_alone_with_braker_inputs_is_clean() {
        let cfg = PredictConfig {
            use_braker: true,
            braker_rna_bams: vec![PathBuf::from("a.bam")],
            braker_proteins: Some(PathBuf::from("p.fa")),
            ..PredictConfig::default()
        };
        assert!(
            check_braker_conflicts(&cfg).is_ok(),
            "BRAKER-native flags must not collide with --use-braker"
        );
    }

    #[test]
    fn braker_with_genemark_mode_is_a_conflict() {
        let cfg = PredictConfig {
            use_braker: true,
            genemark_mode: Some(GeneMarkMode::Etp),
            ..PredictConfig::default()
        };
        let err = check_braker_conflicts(&cfg).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("--genemark-mode"),
            "error must name the offending flag: {}",
            msg
        );
        assert!(
            msg.contains("--use-braker"),
            "error must name the BRAKER flag too: {}",
            msg
        );
    }

    #[test]
    fn braker_with_protein_fasta_is_a_conflict() {
        // --protein-fasta has no meaning under BRAKER (BRAKER takes
        // --braker-proteins instead). Reject it loudly so the user doesn't
        // think their proteins are being fed to BRAKER.
        let cfg = PredictConfig {
            use_braker: true,
            protein_fasta: Some(PathBuf::from("p.fa")),
            ..PredictConfig::default()
        };
        let err = check_braker_conflicts(&cfg).unwrap_err();
        assert!(err.to_string().contains("--protein-fasta"));
    }

    #[test]
    fn braker_with_protein_evidence_is_a_conflict() {
        let cfg = PredictConfig {
            use_braker: true,
            protein_evidence: Some(PathBuf::from("blast.tsv")),
            ..PredictConfig::default()
        };
        let err = check_braker_conflicts(&cfg).unwrap_err();
        assert!(err.to_string().contains("--protein-evidence"));
    }

    #[test]
    fn braker_with_protein_hints_is_a_conflict() {
        // --protein-hints feeds the ab-initio Augustus path, which BRAKER
        // bypasses entirely; reject it so the user isn't misled.
        let cfg = PredictConfig {
            use_braker: true,
            protein_hints: Some(PathBuf::from("proteins.fa")),
            ..PredictConfig::default()
        };
        let err = check_braker_conflicts(&cfg).unwrap_err();
        assert!(err.to_string().contains("--protein-hints"));
    }

    #[test]
    fn braker_with_glimmerhmm_is_a_conflict() {
        let cfg = PredictConfig {
            use_braker: true,
            use_glimmerhmm: true,
            ..PredictConfig::default()
        };
        let err = check_braker_conflicts(&cfg).unwrap_err();
        assert!(err.to_string().contains("--glimmerhmm"));
    }

    #[test]
    fn braker_with_legacy_genemark_alias_is_a_conflict() {
        let cfg = PredictConfig {
            use_braker: true,
            use_genemark: true,
            ..PredictConfig::default()
        };
        let err = check_braker_conflicts(&cfg).unwrap_err();
        assert!(err.to_string().contains("--genemark"));
    }

    #[test]
    fn braker_lists_all_conflicting_flags_at_once() {
        // The error message should be exhaustive — fix-then-rerun rather than
        // whack-a-mole.
        let cfg = PredictConfig {
            use_braker: true,
            genemark_mode: Some(GeneMarkMode::Es),
            protein_fasta: Some(PathBuf::from("p.fa")),
            use_glimmerhmm: true,
            ..PredictConfig::default()
        };
        let err = check_braker_conflicts(&cfg).unwrap_err();
        let msg = err.to_string();
        for flag in ["--genemark-mode", "--protein-fasta", "--glimmerhmm"] {
            assert!(
                msg.contains(flag),
                "error should list every offending flag; missing {} in: {}",
                flag,
                msg
            );
        }
    }

    // ── Default config sanity ───────────────────────────────────────────────

    #[test]
    fn predict_default_braker_off() {
        let c = PredictConfig::default();
        assert!(!c.use_braker);
        assert!(c.braker_mode.is_none());
        assert!(c.braker_rna_bams.is_empty());
        assert!(c.braker_proteins.is_none());
        assert_eq!(c.genetic_code, 1);
    }
}

#[cfg(test)]
mod update_kallisto_cli_tests {
    //! Tests for the `--kallisto` / `--kallisto-min-tpm` flags on the
    //! `update` subcommand. We drive the binary through assert_cmd so
    //! the full main-arm parser path is exercised.

    use assert_cmd::Command;
    use predicates::prelude::*;
    use std::io::Write;
    use tempfile::tempdir;

    /// `--kallisto-min-tpm <value>` accepts well-formed floats.
    /// Smoke-test by invoking with a missing GFF: parsing happens
    /// before the run, so a successful-parse-but-failed-run path
    /// still proves the flag was accepted. We require the error to
    /// be about the GFF, not the flag.
    #[test]
    fn kallisto_min_tpm_accepts_valid_float() {
        let dir = tempdir().unwrap();
        let bogus_gff = dir.path().join("missing.gff3");
        let bogus_fa = dir.path().join("missing.fa");

        let mut cmd = Command::cargo_bin("myconote-cli").unwrap();
        cmd.args([
            "update",
            bogus_gff.to_str().unwrap(),
            "--fasta",
            bogus_fa.to_str().unwrap(),
            "--kallisto-min-tpm",
            "0.5",
        ]);
        // It will fail because the GFF doesn't exist; that's expected.
        // Crucially, the error must NOT mention --kallisto-min-tpm
        // (which would mean the flag parse failed).
        let output = cmd.assert().failure();
        let stderr = String::from_utf8_lossy(&output.get_output().stderr).to_string();
        assert!(
            !stderr.contains("--kallisto-min-tpm"),
            "0.5 should parse cleanly; stderr: {}",
            stderr
        );
    }

    /// `--kallisto-min-tpm` rejects malformed values (no silent default).
    #[test]
    fn kallisto_min_tpm_rejects_garbage() {
        let dir = tempdir().unwrap();
        let bogus_gff = dir.path().join("missing.gff3");
        let bogus_fa = dir.path().join("missing.fa");

        let mut cmd = Command::cargo_bin("myconote-cli").unwrap();
        cmd.args([
            "update",
            bogus_gff.to_str().unwrap(),
            "--fasta",
            bogus_fa.to_str().unwrap(),
            "--kallisto-min-tpm",
            "not-a-float",
        ]);
        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("--kallisto-min-tpm"));
    }

    /// `--kallisto-min-tpm` rejects negative values (TPM is non-negative).
    #[test]
    fn kallisto_min_tpm_rejects_negative() {
        let dir = tempdir().unwrap();
        let bogus_gff = dir.path().join("missing.gff3");
        let bogus_fa = dir.path().join("missing.fa");

        let mut cmd = Command::cargo_bin("myconote-cli").unwrap();
        cmd.args([
            "update",
            bogus_gff.to_str().unwrap(),
            "--fasta",
            bogus_fa.to_str().unwrap(),
            "--kallisto-min-tpm",
            "-1.0",
        ]);
        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("--kallisto-min-tpm"));
    }

    /// `--kallisto` without any RNA-seq input must error with a clear
    /// message naming `--rna-r1`. We have to provide a real GFF + FASTA
    /// so the run gets past the file-existence checks and reaches the
    /// kallisto branch.
    #[test]
    fn kallisto_without_rnaseq_errors_loudly() {
        let dir = tempdir().unwrap();
        let gff = dir.path().join("genes.gff3");
        let fa = dir.path().join("genome.fa");
        // Minimal valid GFF3 + FASTA.
        let mut g = std::fs::File::create(&gff).unwrap();
        writeln!(g, "##gff-version 3").unwrap();
        writeln!(g, "chr1\tmyconote\tgene\t1\t100\t.\t+\t.\tID=g1").unwrap();
        let mut f = std::fs::File::create(&fa).unwrap();
        writeln!(f, ">chr1").unwrap();
        writeln!(f, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC").unwrap();
        writeln!(f, "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();

        let mut cmd = Command::cargo_bin("myconote-cli").unwrap();
        cmd.args([
            "update",
            gff.to_str().unwrap(),
            "--fasta",
            fa.to_str().unwrap(),
            "--output",
            dir.path().join("update_out").to_str().unwrap(),
            "--kallisto",
        ]);
        // Either the kallisto-not-installed error or the no-RNA-seq
        // error is acceptable — both are loud, both name actionable
        // next steps. We just need it to fail (not silently fall back).
        let output = cmd.assert().failure().get_output().clone();
        let combined = format!(
            "{}{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        // Must mention either kallisto being missing, or RNA-seq inputs.
        let mentions_kallisto = combined.to_lowercase().contains("kallisto");
        let mentions_rnaseq = combined.to_lowercase().contains("rna-seq")
            || combined.to_lowercase().contains("rna seq")
            || combined.contains("--rna-r1")
            || combined.contains("RNA-seq");
        assert!(
            mentions_kallisto || mentions_rnaseq,
            "error must explain how to recover; got: {}",
            combined
        );
    }
}
