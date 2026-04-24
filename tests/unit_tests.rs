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
