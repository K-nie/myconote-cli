# myconote-cli: a high-performance, end-to-end genome annotation pipeline for eukaryotic genomes

**Benjamin Narh-Madey**^1^

^1^ Hittinger Lab, Laboratory of Genetics, University of Wisconsin-Madison, Madison, WI 53706, USA

**Correspondence:** narhmadey@wisc.edu

---

## Abstract

**Motivation:** Genome annotation remains a critical bottleneck in eukaryotic genomics. Existing pipelines such as funannotate are widely used but constrained by Python's runtime overhead, limited annotation source integration, and lack of built-in output validation. As genome sequencing costs decline and long-read assemblies become routine, the community needs annotation tools that are faster, more comprehensive, and produce submission-ready outputs.

**Results:** We present myconote-cli, a genome annotation pipeline written in Rust that takes a eukaryotic genome assembly from raw contigs to NCBI-ready submission in a single tool. myconote-cli integrates 15 annotation sources (MMseqs2, Pfam, InterProScan, EggNOG, BUSCO, CAZyme, MEROPS, tRNAscan-SE, antiSMASH, and others), supports five eukaryotic kingdoms with tuned defaults, implements 18 NCBI genetic code tables for accurate translation of non-standard organisms, and provides built-in output validation, reproducibility reporting, and NCBI submission preparation. On a 12.9 Mb *Brettanomyces bruxellensis* genome, myconote-cli annotated 5,218 genes with 98.3% receiving functional descriptions, completing the full pipeline (prediction through annotation) in under 15 minutes. Pfam domain search using hmmsearch achieved a 6-fold speedup over the conventional hmmscan approach while maintaining equivalent sensitivity. The tool is distributed as a single 6.1 MB binary with Docker and Singularity containers, includes 89 automated tests, and provides an interactive tutorial system for new users.

**Availability:** myconote-cli is freely available under the MIT licence at https://github.com/K-nie/myconote-cli. Docker images, documentation, and a Singularity definition file are provided for reproducible deployment.

**Keywords:** genome annotation, gene prediction, functional annotation, Rust, fungi, eukaryotes

---

## 1. Introduction

Genome annotation --- the process of identifying genes and assigning them biological function --- is a foundational step in genomics research. For eukaryotic organisms, annotation typically requires a multi-stage pipeline: repeat masking, ab initio gene prediction, evidence-based model refinement, and functional characterization against reference databases. Each stage depends on specialized external tools, and the orchestration of these tools into a coherent, reproducible workflow remains a significant practical challenge.

The most widely used pipeline for fungal genomes is funannotate (Palmer and Stajich, 2020), which wraps Augustus, SNAP, GeneMark, and EvidenceModeler in a Python framework. While funannotate has been cited in over 500 publications, it has several recognized limitations: (i) Python's interpreted nature imposes runtime overhead, particularly for I/O-intensive parsing of large GFF3 and FASTA files; (ii) annotation sources are limited --- it lacks tRNA prediction, carbohydrate-active enzyme (CAZyme) annotation, protease family classification, and genetic code support for non-standard organisms such as the *Candida* CTG clade; (iii) output validation is minimal, with no automated checks for GFF3 compliance prior to NCBI submission; and (iv) reproducibility is not tracked, with no automated logging of tool versions, database dates, or parameter settings.

Here we present myconote-cli, a genome annotation pipeline written in Rust that addresses these limitations while maintaining compatibility with the established tool ecosystem. myconote-cli provides 21 CLI commands spanning the full annotation lifecycle, integrates 15 distinct annotation sources, supports five eukaryotic kingdoms, and produces NCBI-submission-ready outputs with automated validation. The tool is designed to be both a research instrument and a teaching platform, with an integrated interactive tutorial system modelled after R's swirl package.

---

## 2. Design and Implementation

### 2.1 Architecture

myconote-cli is implemented in Rust (2021 edition) and compiles to a single statically-linked binary of 6.1 MB. The architecture follows a modular design with 27 source directories, each encapsulating a distinct capability (Table 1). External bioinformatics tools (Augustus, SNAP, HMMER, MMseqs2, etc.) are invoked as subprocesses, allowing myconote-cli to benefit from their established algorithms while providing a unified interface.

Performance-critical operations leverage Rust's zero-cost abstractions: data-parallel iteration via rayon, memory-mapped file I/O via memmap2, zero-copy parsing via nom, and release-mode link-time optimization (LTO). The binary is compiled with `opt-level=3`, `codegen-units=1`, and symbol stripping for maximum throughput.

### 2.2 Pipeline stages

The annotation pipeline consists of seven sequential stages, each implemented as an independent CLI subcommand:

1. **Sort** (`myconote-cli sort`): Renames and length-sorts contigs, optionally filtering scaffolds below a minimum size threshold. Clean, sequential identifiers (scaffold_001, scaffold_002, ...) prevent naming conflicts in downstream tools.

2. **Mask** (`myconote-cli mask`): Identifies and soft-masks repetitive elements. Five masking engines are supported: minimap2 self-alignment (no database required), RepeatMasker with species-specific libraries, RepeatModeler for de novo repeat library construction, and combined modes for maximum sensitivity.

3. **Train** (`myconote-cli train`): Trains organism-specific gene prediction models using RNA-seq data. Reads are assembled with Trinity, aligned to the genome with minimap2, loaded into a PASA transcript database, and used to train Augustus and SNAP species models.

4. **Predict** (`myconote-cli predict`): Calls genes using multiple ab initio predictors (Augustus, SNAP, GlimmerHMM, GeneMark-ES) and optional protein-to-genome alignment evidence (miniprot or Exonerate). An Evidence Modeler-style consensus algorithm merges overlapping predictions using configurable per-source weights specified in TOML format.

5. **Update** (`myconote-cli update`): Refines gene models using transcript evidence. When PASA is available, full isoform-aware updates including UTR extension and alternative splicing are performed. A lightweight fallback uses minimap2 coverage to extend UTR boundaries when PASA is not installed.

6. **Annotate** (`myconote-cli annotate`): Assigns functional annotations from 15 sources (Section 2.3). All sources are optional --- the pipeline degrades gracefully when individual tools or databases are unavailable.

7. **Submit** (`myconote-cli submit`): Validates the annotated GFF3 for NCBI compliance (ID uniqueness, parent-child consistency, coordinate ordering), generates an NCBI feature table (.tbl), and optionally runs table2asn to produce a Sequin submission file (.sqn).

### 2.3 Annotation sources

myconote-cli integrates 15 annotation sources, significantly exceeding the scope of existing pipelines (Table 2):

| Source | Tool | Evidence type |
|--------|------|---------------|
| Swiss-Prot homology | MMseqs2 | Product names, UniProt accessions |
| Pfam domains | hmmsearch | Protein domain architecture |
| InterProScan | EBI REST API | InterPro, TIGRFAM, Gene3D, SMART, Superfamily |
| GO terms | UniProt + InterProScan | Gene Ontology functional categories |
| BUSCO | BUSCO 5 | Genome/proteome completeness |
| EggNOG | eggNOG-mapper | COG/NOG categories, KEGG pathways |
| CAZymes | dbCAN (DIAMOND + HMMER) | Carbohydrate-active enzyme families |
| Secretome | SignalP/DeepSig + DeepTMHMM | Signal peptides, transmembrane topology |
| BGC clusters | antiSMASH | Secondary metabolite gene clusters |
| Proteases | MEROPS (DIAMOND) | Peptidase families and clans |
| tRNA genes | tRNAscan-SE | tRNA prediction (eukaryotic, mitochondrial) |
| Protein evidence | miniprot/Exonerate | Protein-to-genome spliced alignment |
| Genetic codes | NCBI tables (built-in) | 18 translation tables |
| Output validation | Built-in | GFF3/FASTA/GenBank compliance |
| Reproducibility | Built-in | Workflow reports (JSON + text) |

### 2.4 Pfam acceleration with hmmsearch

A key performance optimization in myconote-cli is the use of `hmmsearch` rather than `hmmscan` for Pfam domain annotation. The conventional approach (hmmscan) searches each protein sequence against the full Pfam-A HMM database, incurring per-sequence overhead for loading and preprocessing the ~20,000 profile HMMs. In contrast, hmmsearch reverses the operation: each HMM profile is searched against the full protein database at once, amortizing the sequence loading cost across all profiles. As documented by the HMMER authors (Eddy, 2011), hmmsearch is asymptotically faster when the number of query profiles is large relative to the target sequences. For a typical fungal proteome (5,000 proteins) searched against Pfam-A (20,795 profiles), we observe a ~6-fold wall-clock speedup (7 minutes vs. 45 minutes on identical hardware) with equivalent sensitivity.

### 2.5 Genetic code support

myconote-cli implements 18 NCBI translation tables, enabling accurate protein translation for organisms with non-standard genetic codes. This is critical for the *Candida* CTG clade (Table 12: CTG encodes serine rather than leucine), yeast mitochondrial genomes (Table 3), and other deviations from the standard code. To our knowledge, myconote-cli is the first eukaryotic annotation pipeline to support alternative genetic codes natively in the translation step.

### 2.6 Ploidy awareness

For polyploid genomes, myconote-cli estimates ploidy from the assembly-to-expected-size ratio and adjusts overlap tolerance in the evidence merger accordingly. Allelic duplicates are detected by protein self-alignment (MMseqs2 or DIAMOND) and reported in a dedicated table, allowing users to collapse or retain allelic pairs based on their analysis goals.

### 2.7 Quality assurance

Three quality assurance mechanisms are built into the pipeline:

- **Output validation**: Every GFF3 output is checked for ID uniqueness, parent-child reference integrity, coordinate ordering, and feature hierarchy. Protein FASTA outputs are scanned for internal stop codons, invalid amino acid characters, and minimum length.

- **NCBI compliance**: The `submit` command performs pre-flight validation against NCBI GenBank requirements before generating submission files, catching errors that would otherwise result in rejection.

- **Reproducibility reports**: Each pipeline run generates a JSON and human-readable report documenting the myconote-cli version, external tool versions, database download dates, input file checksums, command-line parameters, system information, and runtime duration.

### 2.8 Interactive tutorial

myconote-cli includes an interactive, self-paced tutorial system (`myconote-cli learn`) modelled after R's swirl package. Eight lessons cover the full pipeline from basic concepts to NCBI submission, using five question types (multiple choice, free text with fuzzy matching, true/false, fill-in-the-blank, step ordering). Progress is persisted across sessions, and users can resume, skip, or request hints at any point. This feature is designed to lower the barrier to entry for graduate students and researchers new to genome annotation.

---

## 3. Results

### 3.1 Benchmark: *Brettanomyces bruxellensis*

We validated myconote-cli on a *Brettanomyces bruxellensis* genome assembly (12.9 Mb, 30 contigs) using the full pipeline on a MacBook Pro (Apple M3 Pro, 18 GB RAM).

**Gene prediction:** Augustus with the *S. cerevisiae* S288C species model predicted 5,218 protein-coding genes (median length 1,161 bp, N50 1,854 bp). The evidence merger produced a clean consensus GFF3 with zero duplicate feature IDs and zero orphan parent references, passing NCBI validation with no errors.

**Functional annotation** (Table 3):

| Annotation source | Genes annotated | Coverage | Wall-clock time |
|-------------------|----------------|----------|-----------------|
| MMseqs2 (Swiss-Prot) | 3,915 | 75.0% | 1.5 min |
| Pfam (hmmsearch) | 4,464 | 85.6% | 7 min |
| InterProScan (EBI API) | 5,129 | 98.3% | 5 min (cached) |
| GO terms (UniProt + IPR) | 4,856 | 93.1% | included |
| **Combined** | **5,129** | **98.3%** | **~15 min total** |

The 98.3% annotation rate reflects the combined contribution of InterProScan (which covers dozens of member databases) supplemented by MMseqs2 Swiss-Prot hits for product names.

**Pfam performance comparison:**

| Method | Domain hits | Wall-clock time | Speedup |
|--------|------------|-----------------|---------|
| hmmscan (conventional) | 10,608 | 45 min | 1x |
| hmmsearch (myconote-cli) | 10,686 | 7 min | 6.4x |

The slight increase in hit count with hmmsearch (+0.7%) is consistent with the HMMER documentation noting minor sensitivity differences due to different E-value calibration contexts.

### 3.2 Comparison with funannotate

Table 4 summarizes the feature comparison between myconote-cli v0.1.0 and funannotate v1.8:

| Capability | funannotate | myconote-cli |
|------------|-------------|--------------|
| CLI commands | 8 | 21 |
| Annotation sources | 8 | 15 |
| Gene predictors | 4 | 6 (+ protein evidence) |
| Kingdoms supported | 1 (fungi-focused) | 5 (per-kingdom defaults) |
| Genetic code tables | 1 (standard only) | 18 |
| tRNA prediction | No | Yes (tRNAscan-SE) |
| CAZyme annotation | No | Yes (dbCAN) |
| Protease annotation | No | Yes (MEROPS) |
| Pfam method | hmmscan | hmmsearch (6x faster) |
| NCBI submission prep | tbl2asn wrapper | Validation + .tbl + table2asn |
| Output validation | No | GFF3 + FASTA compliance |
| Reproducibility reports | No | JSON + text |
| Configurable EVM weights | No | TOML file |
| Ploidy awareness | No | Detection + allelic filtering |
| Protein evidence | Exonerate (limited) | miniprot/Exonerate (fully wired) |
| Phylogenetics | No | IQ-TREE 2 (built-in) |
| Synteny visualization | No | Ribbon diagrams |
| Genome browser | No | JBrowse2 / UCSC |
| Format conversions | Limited | 15+ formats |
| Interactive tutorial | No | 8 lessons (swirl-style) |
| Containers | Docker | Docker + Singularity |
| CI/CD | No | GitHub Actions |
| Tests | Community-tested | 89 automated tests |
| Binary size | ~50 MB (Python) | 6.1 MB (Rust) |
| Language | Python 3 | Rust 2021 |

### 3.3 Scalability

myconote-cli's Rust implementation provides consistent performance advantages on I/O-bound operations common in genome annotation. GFF3 parsing, protein extraction, and result merging all benefit from memory-mapped I/O, zero-copy parsing, and data-parallel iteration. The compiled binary eliminates Python startup overhead (~0.5s per invocation), which compounds across the hundreds of subprocess calls in a typical pipeline run.

---

## 4. Discussion

myconote-cli represents a ground-up reimplementation of the eukaryotic genome annotation pipeline in a systems programming language, designed to address the practical limitations encountered by researchers using existing tools. Three design decisions merit discussion.

**Breadth vs. depth of annotation.** By integrating 15 annotation sources into a single pipeline, myconote-cli reduces the manual effort required to achieve comprehensive functional characterization. The 98.3% annotation rate achieved on *B. bruxellensis* with InterProScan demonstrates that the annotation ceiling is set by database coverage, not pipeline capability. For organisms with less representation in reference databases, the graceful degradation model ensures that partial results are always produced.

**Genetic code correctness.** The silent mistranslation of proteins in organisms with non-standard genetic codes is an underappreciated source of error in published annotations. For the ~400 known species in the *Candida* CTG clade, every CTG codon is translated as leucine (standard code) rather than serine (Table 12), systematically corrupting protein sequences and downstream functional assignments. myconote-cli's native genetic code support eliminates this class of error.

**Reproducibility as a first-class concern.** The absence of reproducibility tracking in existing pipelines makes it difficult to determine, months after an analysis, which versions of which tools and databases produced a given result. myconote-cli's automated reproducibility reports provide a complete provenance record for every pipeline run, addressing a growing concern in the genomics community (Gruning et al., 2018).

**Limitations.** myconote-cli is a new tool (v0.1.0) without the community validation that comes from widespread use. While our automated test suite covers 89 test cases and the pipeline has been validated on real genomes, additional benchmarking across diverse taxa and genome sizes is needed. We encourage the community to test myconote-cli alongside established tools and report results.

---

## 5. Conclusion

myconote-cli provides a fast, comprehensive, and reproducible genome annotation pipeline for eukaryotic genomes. Its Rust implementation delivers significant performance improvements over Python-based alternatives, while its broader annotation source integration, genetic code support, output validation, and NCBI submission preparation close practical gaps that have limited existing tools. The integrated tutorial system lowers the barrier to entry for new users, and Docker/Singularity containers ensure reproducible deployment across computing environments.

---

## Data Availability

myconote-cli is freely available under the MIT licence at https://github.com/K-nie/myconote-cli. The *Brettanomyces bruxellensis* test dataset, including sorted genome, masked genome, and gene predictions, is included in the repository under `brettanomyces_test/`. Docker images are available at `ghcr.io/k-nie/myconote-cli`.

---

## Funding

This work was supported by the Hittinger Lab, Laboratory of Genetics, University of Wisconsin-Madison.

---

## References

Eddy, S.R. (2011) Accelerated profile HMM searches. *PLoS Computational Biology*, 7(10), e1002195.

Gruning, B. et al. (2018) Practical computational reproducibility in the life sciences. *Cell Systems*, 6(6), 631--636.

Holt, C. and Yandell, M. (2011) MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. *BMC Bioinformatics*, 12, 491.

Korf, I. (2004) Gene finding in novel genomes. *BMC Bioinformatics*, 5, 59.

Lomsadze, A. et al. (2005) Gene identification in novel eukaryotic genomes by self-training algorithm. *Nucleic Acids Research*, 33(20), 6494--6506.

Lowe, T.M. and Chan, P.P. (2016) tRNAscan-SE On-line: integrating search and context for analysis of transfer RNA genes. *Nucleic Acids Research*, 44(W1), W54--W57.

Manni, M. et al. (2021) BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. *Molecular Biology and Evolution*, 38(10), 4647--4654.

Mirdita, M. et al. (2019) MMseqs2 desktop and local web server app for fast, interactive sequence searches. *Bioinformatics*, 35(16), 2856--2858.

Palmer, J.M. and Stajich, J.E. (2020) Funannotate v1.8: eukaryotic genome annotation. *Zenodo*. https://doi.org/10.5281/zenodo.4054262.

Stanke, M. et al. (2006) Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources. *BMC Bioinformatics*, 7, 62.

---

## Supplementary Tables

**Table S1.** Complete list of 30+ external tools integrated by myconote-cli, with version requirements, installation methods, and pipeline stages in which each tool is used.

**Table S2.** Full NCBI genetic code table support in myconote-cli, showing codon reassignments for each of the 18 implemented translation tables.

**Table S3.** Evidence Modeler default weights and their biological rationale.

---

## Figure Legends

**Figure 1.** myconote-cli pipeline architecture. Seven sequential stages (sort, mask, train, predict, update, annotate, submit) process a genome assembly from raw contigs to NCBI-ready submission. Arrows indicate data flow; optional stages are shown with dashed borders. The annotate stage integrates 15 sources (right panel) with graceful degradation when individual tools are unavailable.

**Figure 2.** Annotation completeness on *Brettanomyces bruxellensis*. (A) Fraction of 5,218 predicted genes receiving functional descriptions from each annotation source. (B) Overlap between MMseqs2 Swiss-Prot hits, Pfam domain annotations, and InterProScan results. (C) GO term coverage by evidence source.

**Figure 3.** Pfam search performance comparison. Wall-clock time for searching 5,218 *B. bruxellensis* proteins against Pfam-A (20,795 profiles) using hmmscan (conventional) vs. hmmsearch (myconote-cli). Both methods were run with 4 threads on identical hardware (Apple M3 Pro, 18 GB RAM). Error bars show standard deviation across three runs.

**Figure 4.** Feature comparison between myconote-cli and funannotate. Radar plot showing normalized scores across eight capability dimensions: annotation sources, prediction tools, output formats, kingdoms supported, performance, validation, reproducibility, and user experience.
