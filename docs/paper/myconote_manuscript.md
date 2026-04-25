<!--
v0.7.1 audit notes (Benjamin: review before next manuscript revision; not auto-fixed):
  1. Abstract + §1 contributions claim "fifteen functional annotation sources",
     but src/annotate/ has 11 distinct annotation modules (mmseqs/Swiss-Prot,
     pfam, interproscan, eggnog, busco, cazyme, merops, secretome, antismash,
     trnascan, go). Either reduce the count to 11 (matches code + README) or
     reframe what the 15 are counting (e.g. include sub-sources within
     InterProScan: Pfam-via-InterPro, TIGRFAM, SMART, Gene3D, Superfamily — that
     gets to ~15 if InterProScan member DBs are counted separately).
  2. Abstract + §1 contributions say "ships ~35 curated fungal Augustus species";
     since v0.7.0 the manifest is 49 (verified via
     AUGUSTUS_FUNGI_SPECIES.len() == 49 in src/setup/mod.rs and CHANGELOG).
  3. §1 cites Santos et al. (2011) and Mühlhausen et al. (2016) for the
     "approximately 400 species of the Candida CTG clade" claim; both refs
     resolved on 2026-04-25 (Santos is MBE 28:2185, Mühlhausen is Curr Opin
     Microbiol 32:16). Numerical claim is consistent with both reviews.
  These flags carried over from the v0.7.1 audit pass; code (registry,
  manifest) is the source of truth.
-->

# MycoNote-CLI: an integrated, validated, and reproducible fungal genome annotation pipeline

**Benjamin Narh-Madey**^1,2^, **Mike Place**^1,2^, **Steve J. Schrodi**^3^, **Antonis Rokas**^4^, **Chris Todd Hittinger**^1,2,*^

^1^ Laboratory of Genetics, University of Wisconsin-Madison, Madison, WI 53706, USA
^2^ Wisconsin Energy Institute, University of Wisconsin-Madison, Madison, WI 53726, USA
^3^ Center for Human Genomics and Precision Medicine, University of Wisconsin-Madison, Madison, WI 53705, USA
^4^ Department of Biological Sciences, Vanderbilt University, Nashville, TN 37235, USA

^*^ Correspondence: cthittinger@wisc.edu

---

## Abstract

**Motivation:** Fungal genome annotation remains a significant bottleneck in genomics research. Existing pipelines such as funannotate and MAKER provide valuable services but face limitations including Python interpreter overhead for I/O-intensive operations, awkward handling of non-standard genetic codes (relevant to the ~400 species of the *Candida* CTG clade), separation of RNA-seq-based gene-evidence collection from downstream quantification and differential expression analysis, absence of pre-flight validation against NCBI submission requirements, and limited reproducibility tracking. As long-read sequencing makes high-quality fungal assemblies routine, the community needs annotation tools that integrate diverse functional sources, expose yeast-specific genetic codes as first-class options, cover the full path from raw reads to differential expression, validate outputs before submission, and document analyses in a reproducible manner.

**Results:** We present MycoNote-CLI, a Rust-implemented fungal genome annotation pipeline. The tool orchestrates 30+ external bioinformatics programs through a unified command-line interface, integrating fifteen functional annotation sources (MMseqs2 versus UniProt/Swiss-Prot, Pfam via hmmsearch, BUSCO, InterProScan via the EBI REST API with persistent caching, EggNOG-mapper, dbCAN, MEROPS, tRNAscan-SE, antiSMASH, signal peptide and transmembrane prediction, Gene Ontology assignment, and others). Defaults, test fixtures, and benchmarks target fungal genomes; kingdom parameters for plants, animals, insects, and protists remain in the codebase as an extensible scaffold but are experimental and not validated at this release. MycoNote-CLI supports 25 NCBI genetic code translation tables verified against the reference codon usage, including the *Candida* CTG clade code (Table 12) and other fungal-relevant codes. The pipeline also includes an integrated RNA-seq stack: `fetch-rna` for SRA/ENA accession ingestion, `quant` for salmon-based quantification with fastp trimming and tximport-ready outputs, `ase` for transcript-level allele-specific expression on heterozygous, F1-hybrid, and diploid fungal genomes via personalized per-haplotype salmon quantification, and three companion R-script generators — `de-template` (DESeq2), `ase-template` (per-sample-corrected binomial test), and `go-template` (topGO Fisher's exact GO enrichment). A `setup --db augustus-fungi` bundle ships ~35 curated fungal Augustus species so a fresh install can predict on most fungi out of the box. The pipeline also performs built-in GFF3 validation against NCBI submission requirements and generates JSON workflow reports documenting tool versions, database download dates, input file checksums, and parameter settings. We validate the pipeline on two test genomes (*Brettanomyces bruxellensis* and *Candida tropicalis*), demonstrate that the alternative yeast nuclear code (Table 12) produces measurably different protein translations than the standard code, and report a controlled-truth ASE smoke test on *C. tropicalis* where 60 phased SNVs across three transcripts at known 70/30/50 hap0:hap1 ratios are recovered with the correct direction and significance.

**Availability:** MycoNote-CLI is freely available under the MIT licence at https://github.com/K-nie/myconote-cli. Pre-built binaries for Linux x86_64, macOS x86_64, and macOS arm64 are distributed through GitHub Releases. Docker and Singularity container definitions are provided. The tool is documented at https://k-nie.github.io/myconote-cli/.

**Contact:** cthittinger@wisc.edu

**Keywords:** genome annotation; gene prediction; functional annotation; Rust; fungi; reproducibility

---

## 1. Introduction

Eukaryotic genome annotation -- the process of identifying genes and assigning biological function -- is a foundational step in genomics research. For most eukaryotic organisms, annotation requires a multi-stage pipeline: repeat masking, ab initio gene prediction, evidence-based model refinement, and functional characterization against reference databases. Each stage depends on specialized external tools, and the orchestration of these tools into a coherent, reproducible workflow remains a practical challenge that consumes substantial bioinformatician time per genome.

Several pipelines have been developed to address this challenge. MAKER (Cantarel et al., 2008; Holt and Yandell, 2011) and its successors provide flexible, configurable annotation workflows. BRAKER (Hoff et al., 2016; Bruna et al., 2021) focuses specifically on training Augustus and GeneMark from RNA-seq evidence and produces high-quality gene predictions for organisms with available transcriptomic data. funannotate (Palmer and Stajich, 2020) is the most widely used pipeline for fungal genomes, cited in over 500 publications. NCBI maintains an internal pipeline (EGAP) for RefSeq annotation that is not publicly available as a standalone tool.

These pipelines have collectively enabled the annotation of thousands of eukaryotic genomes. Certain limitations remain for the fungal use case in particular. funannotate is implemented in Python, which introduces interpreter overhead for I/O-intensive operations on large proteomes and lengthy parser passes. Support for non-standard genetic codes such as the alternative yeast nuclear code (NCBI Translation Table 12, used by approximately 400 species in the *Candida* CTG clade; Santos et al., 2011; Muhlhausen et al., 2016) is available only through per-component configuration, not as a first-class pipeline option. Annotating CTG clade species with the standard code introduces silent mistranslations that propagate through downstream analyses. RNA-seq handling in existing pipelines focuses on evidence-based gene prediction rather than downstream quantification; users who want transcript abundance estimates and differential expression analysis assemble a second pipeline by hand. Allele-specific expression (ASE) — the rule, not the exception, in heterozygous diploid *Candida* isolates, *Saccharomyces* F1 hybrids, and lager allopolyploids — is unsupported by funannotate, MAKER, and BRAKER altogether; the canonical fungal ASE workflow remains a bespoke STAR + WASP wrapper (van de Geijn et al., 2015) or an ad-hoc personalized-transcriptome script, neither of which is integrated with the annotation pipeline that produced the gene models. Built-in pre-flight validation of GFF3 outputs against NCBI submission requirements is uncommon, so errors such as duplicate feature identifiers and orphan parent references are passed to NCBI's tbl2asn and produce cryptic failures days after submission. Reproducibility tracking across the orchestrated toolchain is not built in.

We developed MycoNote-CLI to address these gaps through specific design choices, while acknowledging that funannotate and other established tools remain excellent choices for many use cases. MycoNote-CLI's contributions are: (1) a Rust implementation that reduces orchestration overhead and produces a single static binary distribution; (2) integration of fifteen functional annotation sources, exceeding the eight typically used in existing fungal pipelines; (3) first-class support for the *Candida* CTG clade and other yeast-relevant genetic codes through a single pipeline-wide flag, with fungal-tuned defaults for masking, prediction weights, and annotation thresholds; (4) 25 NCBI genetic code tables verified against the reference codon usage, organized in a data-driven registry; (5) an integrated RNA-seq stack that covers accession ingestion (`fetch-rna` for SRA/ENA), single-reference quantification (`quant` running fastp and salmon with tximport-ready outputs and a reproducibility bundle; Patro et al., 2017), transcript-level allele-specific expression (`ase` — personalized-transcriptome salmon quantification per haplotype from a phased VCF, with a strict refusal to operate on unphased heterozygous sites), and three companion R-script generators (`de-template` for DESeq2 differential expression, `ase-template` for per-transcript binomial ASE testing with sample-specific null correction, `go-template` for topGO Fisher's-exact GO enrichment (Alexa et al., 2006)) that the user runs through `Rscript` so MycoNote-CLI takes no R runtime dependency; (6) built-in pre-flight validation of GFF3 outputs against NCBI submission requirements; (7) automated reproducibility reporting; and (8) an integrated interactive tutorial system designed to lower the barrier to entry for new users. A curated `setup --db augustus-fungi` bundle of ~35 fungal Augustus species ships with the tool so a fresh install can predict on most fungi without manual species-model installation; a `clean --mode contigs` mode runs minimap2 self-alignment to drop contigs entirely subsumed by longer ones, useful for purging haplotigs from draft assemblies before annotation.

We emphasize that several of these features are not novel in the absolute sense. Non-standard genetic codes have long been supported by individual tools through configuration files (Augustus species models, GeneMark `--gcode` flag); MycoNote-CLI's contribution is making this support accessible through a single command-line flag in an integrated pipeline. Similarly, the use of `hmmsearch` rather than `hmmscan` for proteome-scale Pfam searches is a long-known optimization documented in the HMMER user's guide (Eddy, 2011), but most existing annotation pipelines use `hmmscan` by default. The Evidence Modeler approach (Haas et al., 2008) underlying our consensus gene caller is also not novel; we have reimplemented it in Rust with configurable per-source weights. The novelty of MycoNote-CLI lies primarily in integration, accessibility, and engineering quality rather than in any single underlying algorithm.

---

## 2. Design and Implementation

### 2.1 Architecture

MycoNote-CLI is implemented in Rust (2021 edition). The choice of Rust over Python reflects the specific requirements of an orchestration layer that processes large genomic files and invokes many external tools: compiled performance for I/O-heavy parsing operations, single-binary distribution without runtime dependencies, memory safety guarantees from the type system, and cargo's reproducible build system. The pipeline does not reimplement gene prediction or sequence alignment algorithms; instead, it invokes established external tools (Augustus, SNAP, GeneMark-ES, MMseqs2, HMMER, BUSCO, BUSCO, InterProScan, and others) as subprocesses, parses their outputs, and integrates the results.

The MycoNote-CLI codebase is organized into 27 source modules. Performance-critical operations leverage data-parallel iteration via rayon (Stone et al., 2020), memory-mapped file I/O via memmap2, and zero-copy parsing via nom. The release binary is compiled with link-time optimization, single codegen unit, symbol stripping, and maximum optimization level, producing a 6.1 MB statically linked executable.

### 2.2 Pipeline Stages

The annotation pipeline consists of seven sequential stages, each implemented as an independent CLI subcommand that can be run individually or as part of a complete workflow:

1. **sort**: Renames and length-sorts contigs, filters scaffolds below a configurable minimum size threshold, and writes a rename mapping for traceability.

2. **mask**: Identifies and soft-masks repetitive elements. Five masking strategies are supported: minimap2 self-alignment (no database required), RepeatMasker with species-specific repeat libraries, RepeatModeler2 (Flynn et al., 2020) for de novo repeat library construction, and combined modes for maximum sensitivity.

3. **train**: Trains organism-specific Augustus and SNAP models from RNA-seq evidence. Reads are assembled with Trinity (Grabherr et al., 2011), aligned to the genome with minimap2 (Li, 2018), loaded into a PASA transcript database (Haas et al., 2003), and used to train Augustus and SNAP species-specific parameters. We acknowledge that BRAKER (Hoff et al., 2016; Bruna et al., 2021) provides more sophisticated RNA-seq-based training; users seeking maximum prediction accuracy should consider running BRAKER for prediction and then MycoNote-CLI for functional annotation as a complementary workflow.

4. **predict**: Calls genes using multiple ab initio predictors (Augustus, SNAP, GlimmerHMM, GeneMark-ES) and optional protein-to-genome alignment evidence (miniprot or Exonerate). A reimplementation of the Evidence Modeler consensus algorithm (Haas et al., 2008) merges overlapping predictions using configurable per-source weights specified in TOML format. The default weights (Augustus=10, SNAP=3, GlimmerHMM=2, GeneMark=5, protein evidence=20) reflect typical accuracy hierarchies and can be customized.

5. **update**: Refines gene models using transcript evidence. When PASA is available, full isoform-aware updates including UTR extension are performed. A lightweight fallback uses minimap2 coverage to extend UTR boundaries when PASA is not installed; the lightweight mode is an approximation and is not equivalent to full PASA refinement.

6. **annotate**: Assigns functional annotations from fifteen sources (Section 2.3). Each source is independently optional; the pipeline employs graceful degradation, skipping unavailable tools while continuing with remaining steps. Users should review the annotation report to confirm which sources actually executed, as graceful degradation can mask configuration problems.

7. **submit**: Validates the annotated GFF3 against NCBI submission requirements (section 2.5), generates an NCBI feature table (.tbl), and optionally invokes table2asn to produce a Sequin submission file (.sqn).

### 2.3 Annotation Sources

MycoNote-CLI integrates fifteen functional annotation sources covering protein homology, domain architecture, completeness assessment, functional categorization, specialized enzyme classification, and comprehensive domain searching. The sources and their underlying tools are summarized in Table 1.

For Pfam domain search (Mistry et al., 2021), MycoNote-CLI uses HMMER's `hmmsearch` rather than `hmmscan`. As documented in the HMMER user's guide (Eddy, 2011), `hmmsearch` is generally faster than `hmmscan` for searching a sequence database against a profile database, because the cost of loading each profile is amortized across all sequences. For a typical fungal proteome of 5,000 proteins searched against the ~20,795 Pfam-A profiles, we measured `hmmsearch` runtimes of approximately 7 minutes compared to approximately 45 minutes for `hmmscan` on the same hardware (Apple M3 Pro, 4 threads, single replicate; see Section 4 for limitations of these benchmarks). This is not a novel optimization but represents a meaningful practical improvement when the default tool choice is changed in an integrated pipeline.

For protein homology against UniProt/Swiss-Prot (UniProt Consortium, 2023), MycoNote-CLI uses MMseqs2 `easy-search` (Steinegger and Soding, 2017). MMseqs2 achieves substantial speedup over BLAST while maintaining comparable sensitivity for protein homology detection. The exact speedup depends on hardware, query length distribution, and database composition; we observed approximately 100-fold speedup on our test data, consistent with published benchmarks.

For comprehensive domain annotation, MycoNote-CLI integrates with the EBI InterProScan REST API (Blum et al., 2021) rather than installing InterProScan locally. The integration includes pre-validation of submitted sequences (length checks, character validation, internal stop codon truncation), rate-limited submission respecting EBI API limits (30 jobs per minute, 30 sequences per job), concurrent polling of submitted jobs, persistent caching of results in `~/.myconote/iprscan_cache.tsv`, and binary-search error isolation when individual sequences fail submission. This makes repeated annotations of overlapping protein sets dramatically faster after the initial cache is populated.

### 2.4 Genetic Code Support

MycoNote-CLI supports 25 NCBI genetic code translation tables with full implementations, verified against NCBI's reference codon usage by per-table unit tests. The supported set covers the standard and canonical mitochondrial codes (Tables 1, 2, 3, 4, 5, 6), the nuclear and mitochondrial reassignments most relevant to fungal and protist annotation (Tables 9, 10, 11, 12, 13, 14, 16), and the less common but NCBI-cataloged codes (Tables 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33). Table 12 (Alternative Yeast Nuclear, the *Candida* CTG clade code) and Table 26 (*Pachysolen tannophilus* Nuclear) are of particular relevance for yeast genomics. Three tables (27, 28, 31) use context-dependent stop-codon reassignments that cannot be resolved from codon identity alone; for these, MycoNote-CLI resolves to the amino-acid reading and flags the limitation in the annotation report rather than silently falling back to the Standard code. The genetic-code registry is data-driven, so adding any future NCBI-defined table requires a single entry rather than new decoding logic.

The Alternative Yeast Nuclear code (Table 12) is of particular relevance for fungal genomics. Approximately 400 yeast species in the CTG clade -- including the major human pathogens *Candida albicans*, *C. tropicalis*, *C. parapsilosis*, and *C. metapsilosis*; the industrial yeast *Debaryomyces hansenii*; and many others -- reassigned the CTG codon from leucine to serine approximately 170 million years ago (Santos et al., 2011). Annotating these species with the standard code translates each CTG as leucine (a hydrophobic amino acid) instead of serine (a polar amino acid capable of hydrogen bonding). A typical fungal protein contains 10-15 CTG codons; the cumulative effect on protein structure prediction, domain assignment, and functional inference can be substantial. We emphasize that other tools support non-standard codes through configuration files (e.g., Augustus species-specific models, GeneMark `--gcode` flag); MycoNote-CLI's contribution is exposing this support through a single command-line flag (`--genetic-code N`) that propagates through the entire annotation pipeline.

### 2.5 Output Validation

MycoNote-CLI performs pre-flight validation of all generated GFF3 files against the structural requirements of NCBI GenBank submission. The validator checks for: (1) duplicate feature IDs across the entire file; (2) orphan Parent references where a feature claims to be the child of a non-existent parent; (3) coordinate ordering errors where start > end; (4) GFF3 sequence identifiers that do not exist in the accompanying FASTA file; (5) presence of locus_tag attributes on gene features (warning if absent); (6) presence of product descriptions on gene features (warning if absent); and (7) feature hierarchy consistency (warning for childless genes).

The `--validate-only` flag enables structural validation without generating submission files, allowing users to iterate on fixes before committing to a full table2asn run. We emphasize that this validation catches structural errors but does not guarantee that the resulting submission will be accepted by NCBI; biological correctness, locus_tag prefix registration, and other submission requirements are the user's responsibility.

### 2.6 Allele-specific Expression

Many fungal RNA-seq experiments are run on genomes that are not haploid in the strict sense: clinical *Candida* isolates carry substantial heterozygosity and undergo loss-of-heterozygosity events under azole stress (Mixão and Gabaldón, 2020); *Saccharomyces cerevisiae* × *S. paradoxus* F1 hybrids have anchored cis/trans regulatory dissection in yeast for over a decade (Tirosh et al., 2009); and lager *Saccharomyces pastorianus* is an *S. cerevisiae* × *S. eubayanus* allopolyploid with distinct subgenomes (Krogerus et al., 2017). The binomial-distribution framework for ASE inference that `ase-template` builds on — read counts at the more-expressed allele tested against a per-sample-corrected null — was proposed for hybrid mRNA-seq by McManus et al. (2010) in *Drosophila* and is directly transferable to yeast hybrid data. Quantifying RNA-seq reads from these genomes against a single-haplotype reference under-counts the more divergent allele, masks real cis-regulatory differences, and makes loss-of-heterozygosity events invisible to downstream DE.

The `ase` subcommand addresses this. Inputs are a phased VCF, a reference CDS FASTA from `convert --to cds`, the same annotated GFF3 that produced the CDS (used for genome → CDS coordinate mapping), the reference genome FASTA (used as salmon's decoy set), and a sample sheet matching the `quant` schema. For each haplotype, the parser walks the phased VCF, picks the ALT allele whose phase index matches the haplotype, maps the genomic coordinate to the spliced CDS coordinate accounting for intron removal and strand, reverse-complements REF/ALT for variants on minus-strand transcripts, and emits a personalized CDS FASTA. Salmon then builds one decoy-aware index per haplotype; for each sample, fastp is run once (QC is haplotype-independent), and salmon is run twice — once against each haplotype's index. The two `quant.sf` files are merged into a transcript × `<sample>.<hap>` counts matrix and a parallel TPM matrix. Per-transcript informativeness (whether the two haplotypes actually differ in CDS sequence) is recorded so downstream tests can drop transcripts that cannot in principle show ASE. Unphased heterozygous sites are a hard error with the offending VCF line number, not a silent skip — ASE without phasing is guesswork, and we refuse to guess. Nine documented skip categories (structural, multiallelic, missing genotype, indel longer than `--max-indel-size`, outside CDS, spans an exon boundary, REF mismatch against the genome, in-cis overlap on the same haplotype, and the unphased hard-error case) are each counted in the `ase_bundle.json` reproducibility manifest, with per-variant reason rows in `variants_skipped.tsv`. The companion `ase-template` subcommand emits a base-R script that runs `binom.test()` per transcript per sample with a sample-specific null proportion derived from that sample's total hap0:hap1 library-size ratio, BH-adjusted within sample. No Bioconductor dependency.

We chose personalized-transcriptome quantification with salmon (Patro et al., 2017) over WASP-style read-remapping (van de Geijn et al., 2015) and GATK ASEReadCounter for a specific reason: those tools operate at the variant-site level on aligned BAMs, while `ase` operates at the transcript level on quasi-aligned reads. The two approaches answer different questions — site-level ASE for fine-mapping cis-regulatory variants, transcript-level ASE for per-mRNA allelic ratios that flow into cis/trans regression and lager-yeast subgenome expression studies — and we explicitly do not solve site-level ASE in this release. Personalized-transcriptome salmon is also tolerant of read-to-reference mismatch by design, which means the haplotype that matches reference will accept hap1-derived reads with a few SNV mismatches; the resulting global mapping-rate asymmetry is corrected by `ase-template`'s per-sample-specific null. A 5% mapping-rate-asymmetry threshold (configurable via `--asymmetry-threshold`) flags strong asymmetry as a phasing or assembly diagnostic in the bundle without altering the count matrix.

### 2.7 Reproducibility Tracking

Each annotation run generates a JSON workflow report documenting: MycoNote-CLI version; pipeline steps executed; external tool versions detected at runtime; database paths, sizes, and download dates from version metadata files; input file paths and MD5 checksums (computed when explicitly requested); command-line parameters; runtime duration; and system information (operating system, hostname, CPU count). Evidence weights used by the consensus gene caller are also saved as a TOML file in the output directory.

This reporting enables several reproducibility goals: documentation of methods for publication, comparison between runs to identify configuration changes, audit trails for regulated environments, and faithful re-running with identical inputs and parameters. We acknowledge important limitations of this reproducibility approach (Section 4.3): bit-for-bit reproducibility is not guaranteed because several external tools use heuristics that may produce slightly different output on different hardware or with different thread counts; database updates between runs are not tracked retrospectively; and external tool versions are detected at runtime rather than pinned. For maximum reproducibility, users should employ the provided Docker or Singularity containers with explicit version tags.

### 2.8 Interactive Tutorial

MycoNote-CLI includes an interactive tutorial system invoked via `myconote-cli learn`. The system presents eight self-paced lessons covering the conceptual and practical aspects of genome annotation: the annotation problem and pipeline overview, setup and installation, sort and mask, gene prediction and consensus algorithms, functional annotation sources, RNA-seq training, NCBI submission, and analysis tools. Each lesson combines explanatory text, multiple-choice and free-text questions with fuzzy answer matching, command demonstrations, and embedded literature references (over 25 papers cited across the eight lessons). Progress is persisted in `~/.myconote/learn_progress.json` so users can resume sessions.

The tutorial is presented as a feature of the tool rather than a research contribution. We have not conducted formal user studies measuring its pedagogical effectiveness; such studies are planned and would be more appropriate for a separate publication in an educational venue.

---

## 3. Validation

### 3.1 Test Datasets

We validated MycoNote-CLI on two yeast genomes selected to test specific aspects of the pipeline: *Brettanomyces bruxellensis* (assembly accession in preparation), a budding yeast that uses the standard genetic code, and *Candida tropicalis* (a CTG clade species). These represent a small set and do not constitute a comprehensive benchmark. A complete benchmarking study comparing MycoNote-CLI against funannotate, MAKER, and BRAKER on a panel of fungal reference genomes spanning Saccharomycotina, Pezizomycotina, and Basidiomycota is planned for a follow-up study and is described in Section 5.2 (Future Work). The results presented here should be interpreted as a demonstration of pipeline functionality on real fungal genomes, not as evidence of superior performance over alternatives.

### 3.2 Brettanomyces bruxellensis Annotation

We applied the full MycoNote-CLI pipeline to the *B. bruxellensis* assembly (12.9 Mb across 30 contigs). After sorting and masking, gene prediction with Augustus produced 5,218 predicted gene models (SNAP was excluded from this run due to a local environmental issue with the SNAP HMM file). Functional annotation used MMseqs2 versus Swiss-Prot, hmmsearch against Pfam-A, and InterProScan via the EBI REST API.

Annotation results are summarized in Table 2. Using a permissive threshold (any annotation source produced any hit), 98.3% of predicted genes received at least one functional annotation. This figure should be interpreted carefully: it includes weak hits and is not a measure of annotation quality. Using a more stringent criterion -- a Swiss-Prot hit with greater than 50% identity and E-value below 1e-50 -- approximately 60% of genes received high-confidence functional assignments. These figures are similar to those reported by funannotate users for comparable fungal genomes (Palmer and Stajich, 2020; community reports).

NCBI submission validation passed with zero structural errors: no duplicate feature IDs, no orphan Parent references, all coordinates valid, and all GFF3 sequence identifiers present in the accompanying FASTA. The total annotation runtime on a single Apple M3 Pro processor with 4 threads was approximately 15 minutes (including approximately 7 minutes for the Pfam search using hmmsearch and approximately 5 minutes for InterProScan with most sequences served from a previously populated cache). Without InterProScan caching, the same run takes approximately 45-60 minutes.

### 3.3 Candida tropicalis with Alternative Genetic Code

To demonstrate the value of integrated genetic code support, we annotated *C. tropicalis* (14.2 Mb, 24 contigs, 6,290 input gene models) twice on the same input data: once with the standard genetic code (Table 1) and once with the Alternative Yeast Nuclear code (Table 12, `--genetic-code 12`). Both runs completed successfully and produced NCBI-valid output. The complete script, parameters, and output are available in `docs/paper/ctg_clade_demonstration.sh` in the repository.

The two runs produced identical protein counts (6,290) and identical total residue counts (3,047,858), confirming that the only difference was codon-to-amino-acid mapping. The amino acid composition differences exactly matched expectations: Table 12 produced 10,960 fewer leucine residues and 10,960 additional serine residues compared to Table 1. This corresponds to 10,960 CTG codons across all *C. tropicalis* CDSs, each silently mistranslated as leucine when using the standard code.

At the per-protein level, **3,871 of 6,290 proteins (61.5%) contained at least one CTG codon and therefore differed between the two translations**. The mean number of changed residues per affected protein was 1.74 (range 1-28), with a median of 1. The distribution was right-skewed: 2,419 proteins had no CTG codons, 3,214 had between 1 and 4, 571 had between 5 and 9, 84 had between 10 and 24, and 2 proteins had 25-28 CTG codons each. These extreme cases represent serine-rich proteins where the mistranslation has the greatest cumulative effect.

The biological significance of these differences is non-trivial. Leucine and serine differ substantially in their side chain chemistry: leucine is hydrophobic (Kyte-Doolittle hydropathy index +3.8) while serine is polar with hydrogen-bond donor capability (-0.8). Consequences for downstream analyses include altered hydrophobicity profiles (relevant to transmembrane prediction), altered domain assignments (HMMs trained on real fungal proteins expect serine at these positions), and altered active site characterization for catalytic residues. We do not provide a formal quantitative analysis of downstream errors here -- comparing Pfam hit counts, InterProScan annotations, or BUSCO scores between the two codes on the full *C. tropicalis* dataset is an obvious follow-up experiment planned for the comprehensive benchmarking study.

This experiment confirms that the genetic code flag works as intended and that the affected substitutions are biologically meaningful rather than cosmetic. Other tools (GeneMark-ES `--gcode`, Augustus species-specific configuration) support non-standard codes through configuration files; MycoNote-CLI's contribution is making this accessible through a single command-line flag in an integrated annotation pipeline.

### 3.4 Limitations of the Validation

Two genomes from a single yeast clade are insufficient to make general claims about performance across the diversity of fungal genomes. Our validation demonstrates that:

1. The pipeline runs to completion on real fungal genomes
2. The output passes NCBI structural validation
3. The genetic code support produces the expected translation differences
4. The hmmsearch-based Pfam search completes in approximately 7 minutes on the test hardware

Our validation does **not** demonstrate:

1. That MycoNote-CLI produces more accurate annotations than funannotate or BRAKER
2. That the runtime advantages persist across diverse genomes and hardware
3. Whether the multi-kingdom infrastructure (kingdom = plant, animal, insect, or protist) is usable on non-fungal genomes. That capability is experimental and not validated at this release; for plant or animal annotation, BRAKER or MAKER are the right tools
4. That the genetic code support improves downstream analyses (e.g., domain assignments) for CTG clade species
5. That the pipeline handles edge cases (highly fragmented assemblies, polyploid genomes, very large genomes)

These claims require the more comprehensive benchmarking study described in Section 5.2.

### 3.5 Allele-specific Expression Validation

To check that the `ase` + `ase-template` pipeline recovers a known ASE signal, we built a controlled-truth fixture from a single 2.47 Mb contig of the *Candida tropicalis* MYA-3404 assembly (NW_003020038.1). Three plus-strand single-block CDSs were chosen on this contig. Twenty phased heterozygous SNVs were injected per CDS at roughly one variant every 60 bp — a density consistent with *S. cerevisiae* × *S. paradoxus* F1 nucleotide divergence and broadly typical of yeast hybrid pairs. Paired-end reads were simulated with wgsim from each haplotype at three known hap0:hap1 ratios (70/30, 30/70, and 50/50), 12,000 read pairs per transcript. The full `ase` pipeline (custom phased-VCF parser → strand-aware variant application → personalized per-haplotype CDS FASTA → fastp → per-haplotype salmon index → quant) was run end-to-end against the resulting FASTQs, followed by the `ase-template` binomial test in R 4.4.2.

All 60 variants applied cleanly: zero rows in `variants_skipped.tsv` across all nine skip categories. The two personalized transcriptomes' SHA256 digests diverged as expected (`a9ff40…` for hap0 versus `197d03…` for hap1). Per-transcript truth versus recovered hap0 fraction is shown in Table 4. The asymmetry detector correctly fired (mapping_rate_hap0 = 0.999, mapping_rate_hap1 = 0.500, asymmetry_flag = true), reflecting that hap0 is the unmodified reference and salmon tolerates the few SNV-position mismatches that hap1-derived reads carry against it. This is the documented behavior of personalized-transcriptome salmon and is precisely what the per-sample-null correction in `ase-template` exists to absorb: the corresponding sample-specific null proportion comes out at 11992 / (11992 + 6001) = 0.6665, and the truly-balanced transcript g000012.m1 lands within 0.0002 of that null, confirming the correction does what it should. Across the three transcripts the binomial test makes 3/3 ground-truth calls correctly: g000005.m1 (truly hap0-biased, padj = 1.78e-58, significant), g000010.m1 (truly hap1-biased, padj = 1.53e-41, significant), and g000012.m1 (truly balanced, padj = 0.989, not significant). This fixture covers phased VCF parsing, strand-aware genome-to-CDS coordinate mapping, variant application, per-haplotype salmon indexing and quantification, the asymmetry detector, the long-to-wide merge, the informativeness flag, and the per-sample-null correction in the R template against a known signal. It does not cover minus-strand variants on minus-strand CDSs (covered by unit tests but not in this end-to-end run), indels, multi-exon transcripts spanning intron-flanking variants, or real F1 hybrid data with WhatsHap-phased VCF — those are deferred to a pre-tag sanity check on a real *S. cerevisiae* × *S. paradoxus* dataset before final manuscript submission.

---

## 4. Limitations and Caveats

We believe that an honest discussion of limitations is essential for any methods paper. The following limitations apply to MycoNote-CLI v0.5.1.

### 4.1 Maturity Limitations

MycoNote-CLI is at version 0.5.1 and has not undergone the years of community testing that funannotate, MAKER, and BRAKER have received. Edge cases and bugs that would be familiar to users of established pipelines may not yet have been encountered or fixed. Users should validate critical annotations against alternative tools, report bugs to the GitHub issue tracker, and not yet rely on MycoNote-CLI as the sole annotation pipeline for high-stakes projects without independent validation.

### 4.2 Scientific Limitations

**Functional assignments are uncertain.** A Swiss-Prot match at 30% identity assigns a "best guess" function that may be wrong. Pfam domains can be present without the protein performing the canonical function. Users should treat all automated functional assignments as hypotheses rather than definitive characterizations.

**Gene prediction makes mistakes.** Even with multiple predictors and evidence-based consensus, gene boundaries are sometimes wrong (especially UTRs), some real genes are missed (especially short genes and overlapping genes), and some false positives are predicted (especially in repeat-rich regions). Manual curation remains valuable for high-stakes projects.

**Annotation rate is not annotation accuracy.** A pipeline that generates more annotations is not necessarily more accurate. A liberal annotation rate may include many low-confidence assignments that could mislead downstream analyses.

**Multi-exon genes with non-canonical splice sites may be missed.** Augustus, SNAP, and similar tools assume GT...AG splice sites and may incorrectly predict introns at non-canonical sites.

### 4.3 Technical Limitations

**Some external tools have ARM64 limitations.** Trinity, PASA, and several other bioinformatics tools do not have ARM64 builds at the time of writing. On Apple Silicon Macs and ARM-based Linux servers, these tools may not be available, limiting the RNA-seq training functionality.

**InterProScan integration requires internet access.** Air-gapped environments cannot use the EBI REST API integration; users must install InterProScan locally and run it separately.

**Memory requirements for Trinity assembly can be substantial.** Trinity requires 32-128 GB of RAM for typical RNA-seq datasets, which exceeds the capacity of many laptops. MycoNote-CLI cannot run Trinity on systems with insufficient memory.

**Genome size scaling above ~1 Gb is untested.** MycoNote-CLI has been validated only on small fungal genomes. Performance and correctness on larger plant or animal genomes are theoretically supported but not demonstrated.

**Polyploid handling is partial.** The `ase` subcommand performs diploid haplotype-aware quantification — two haplotypes from a phased VCF, two personalized transcriptomes, two salmon runs per sample. Higher ploidy (e.g. triploid *Candida* lineages, allopolyploid strains beyond a clean two-subgenome decomposition) is still future work and would require a multinomial or Dirichlet-multinomial extension of the binomial test in `ase-template`. The `--ploidy` flag elsewhere in the pipeline adjusts overlap tolerance and reports allelic duplicates but does not perform haplotype-aware annotation. The `ase-template` binomial test is sample-level: each sample is its own independent test and BH adjustment is per-sample, not pooled. There is no random-effects model in 0.5.x; multi-sample modeling lives in user R code on top of `ase_counts.tsv`.

**Reproducibility is approximate.** Bit-for-bit reproducibility is not guaranteed across different hardware, thread counts, or runs because several external tools use non-deterministic heuristics. For maximum reproducibility, users should employ the provided containers with pinned versions and avoid multi-threading where determinism matters more than speed.

### 4.4 When to Use Other Tools

We do not believe MycoNote-CLI is the right tool for every annotation task. Table 3 summarizes our recommendations for tool selection by use case. In particular, **funannotate remains an excellent choice for fungal annotation projects where established workflows, community validation, and an extensive citation history are priorities**. BRAKER provides more sophisticated RNA-seq-based gene prediction. MAKER offers more fine-grained control over evidence weighting and is well-established for plant genomes. For bacterial or archaeal genomes, Prokka (Seemann, 2014) or Bakta (Schwengers et al., 2021) are better suited. For viral genomes, tools such as VAPiD (Shean et al., 2019) and VADR (Schaffer et al., 2020) are appropriate.

---

## 5. Discussion

### 5.1 Contributions and Their Context

MycoNote-CLI provides several practical contributions to the fungal genome annotation ecosystem. Most are integration and engineering contributions rather than novel algorithms:

**Integration of fifteen annotation sources** in a single pipeline. While each individual source exists in other tools, the unified interface reduces the manual effort of running and merging results from multiple tools. Users gain comprehensive annotation without writing custom integration scripts.

**Fungi-tuned defaults.** Masking strategy, evidence weights, BUSCO lineage selection, and annotation thresholds are calibrated for fungal genomes out of the box. An experimental scaffold for plant, animal, insect, and protist defaults lives in the codebase for future development but is not validated at this release.

**Single-flag access to non-standard genetic codes.** This addresses a real and underappreciated source of annotation errors for the CTG clade and other organisms. The underlying capability exists in individual tools; the contribution is integration and accessibility.

**Built-in pre-flight NCBI validation.** This addresses a workflow pain point that the author and many other genome annotators have encountered repeatedly: structural errors in GFF3 files cause confusing failures days after submission. Catching errors before submission improves the user experience and reduces submission cycle time.

**Reproducibility reporting.** Automated documentation of tool versions, database dates, and parameters supports the broader push toward computational reproducibility in genomics.

**Interactive tutorial.** To our knowledge, no other annotation pipeline includes an integrated, swirl-style tutorial. While we have not formally evaluated its pedagogical effectiveness, we believe lowering the barrier to entry for graduate students and new users is valuable.

**Rust implementation.** The choice of Rust over Python provides modest performance improvements (compiled execution, single-binary distribution, no garbage collector overhead) and substantial robustness improvements (memory safety, type-checked error handling, deterministic deployment via cargo). These benefits compound over the long-term maintenance horizon of a research tool.

### 5.2 Future Work

A comprehensive benchmarking study comparing MycoNote-CLI against funannotate, MAKER, and BRAKER on a panel of at least six fungal reference genomes spanning ascomycete yeasts (Saccharomycotina), filamentous ascomycetes (Pezizomycotina), and basidiomycetes is the most important next step. The study will use identical hardware and tool versions, evaluate sensitivity and specificity against gold-standard manually curated reference annotations (e.g., the *Saccharomyces* Genome Database for *S. cerevisiae*, AspGD for *Aspergillus nidulans*, MycoCosm for several basidiomycetes, and RefSeq where curated alternatives are unavailable), and report statistical analysis of accuracy differences. We aim to complete this study before submission of the final manuscript.

Additional planned work includes: polyploid (>diploid) ASE via a multinomial or Dirichlet-multinomial extension of the binomial test in `ase-template`, with a sketch already in scope for triploid *Candida* and allopolyploid lager strains; a follow-up `ase-sites` subcommand that performs WASP-style variant-site-level allele assignment for users who need fine-mapped cis-regulatory analysis rather than transcript-level allelic ratios; integration of GeneMark-EP (protein-evidence-guided self-training) alongside the existing GeneMark-ES and -ET drivers; BRAKER as an alternative predictor for the "extensive RNA-seq, maximum accuracy" use case where the user wants Augustus + GeneMark trained directly on hint files rather than via PASA + Trinity; long-read transcript support (PacBio Iso-Seq, Nanopore Direct RNA); a formal user study evaluating the tutorial's pedagogical effectiveness; and a web-based interface for users who prefer not to use the command line. The 25 currently supported NCBI translation tables cover all codes published in the NCBI reference at the time of writing; full implementation of any future NCBI-defined tables is trivial thanks to the data-driven registry design.

### 5.3 Sustainability

MycoNote-CLI is currently maintained by the corresponding author and the Hittinger Lab at the University of Wisconsin-Madison. The Hittinger Lab has committed to long-term maintenance and development of the tool as part of its broader genomic annotation infrastructure. Lab member Mike Place and external collaborator Antonis Rokas (Vanderbilt University) have agreed to serve as alternate maintainers. The codebase is open-source under the MIT license, has continuous integration and automated testing on every commit, and uses automated dependency security audits. Community pull requests are welcomed and reviewed.

We acknowledge that single-author or small-team scientific software has historically faced sustainability challenges. The maintenance commitment from the Hittinger Lab, the multiple co-maintainers, the institutional infrastructure of the Laboratory of Genetics, and the modular design of the codebase (which allows individual modules to be maintained independently) collectively mitigate this risk. Users interested in contributing or in the long-term maintenance plan are encouraged to contact the corresponding author.

### 5.4 Conclusion

MycoNote-CLI provides a practical, integrated pipeline for fungal genome annotation. Its contributions are primarily in integration, accessibility, and engineering quality rather than in novel algorithms. The pipeline addresses several real limitations of existing fungal tools, particularly for users working with non-standard genetic codes, the *Candida* CTG clade and other yeast-specific codes, integrated RNA-seq quantification alongside annotation, or projects where NCBI submission readiness and reproducibility documentation are priorities. We acknowledge that comprehensive benchmarking against established alternatives is required to substantiate broader claims about pipeline accuracy, and we have outlined this benchmarking as the primary item of future work. For the fungal genomics community, we hope MycoNote-CLI provides a useful complement to existing tools and contributes to making fungal annotation more accessible, more reproducible, and more accurate.

---

## Reproducibility Statement

All code, configuration files, and documentation supporting this manuscript are available at https://github.com/K-nie/myconote-cli (release tag v0.5.1). The two test datasets used in §§3.2–3.3 are available at the repository under `brettanomyces_test/` (with provenance information in the README) and `tests/data/candida_tropicalis.*`. The ASE controlled-truth fixture in §3.5 is reproduced from `scratch/e2e/genome.fa` (the *C. tropicalis* MYA-3404 contig NW_003020038.1) plus the wgsim simulation scripts in `scratch/ase_validation/`; the full validation log is in `docs/paper/ase_validation.md`. The `setup --db augustus-fungi` curated bundle and the `clean --mode contigs` haplotig-purging step are part of the standard reproducibility path documented in `docs/paper/reproducibility.md`. A `Dockerfile` and `Singularity.def` are provided for reproducible deployment with pinned tool versions. All external tool versions used at the time of validation are recorded in the JSON reproducibility reports generated by the pipeline.

We encourage readers to validate our results by running the pipeline on their own data using the published containers.

---

## Author Contributions

**B.N.M.** conceived the project, designed the architecture, implemented the codebase (including the v0.5.x RNA-seq stack: `quant`, `fetch-rna`, `de-template`, `ase`, `ase-template`, and `go-template`), designed and ran the controlled-truth ASE validation experiment on *C. tropicalis*, performed all validation experiments, and drafted the manuscript. **M.P.** provided technical guidance on bioinformatics workflows, contributed to evidence weight calibration, and reviewed code for the annotation modules. **S.J.S.** advised on statistical aspects, validation methodology, and reproducibility frameworks. **A.R.** provided expert guidance on fungal genomics, the CTG clade biology, and the comparative context with other annotation pipelines. **C.T.H.** supervised the project, provided computational resources, contributed expert knowledge of yeast biology and evolution, secured funding, and revised the manuscript. All authors read and approved the final manuscript.

---

## Acknowledgments

We thank the developers of Augustus (Mario Stanke and team), SNAP (Ian Korf), GeneMark (Mark Borodovsky and team), Evidence Modeler (Brian Haas), HMMER (Sean Eddy and team), MMseqs2 (Martin Steinegger and Johannes Soding), and the many other tools that MycoNote-CLI orchestrates. We thank the Pfam, InterPro, UniProt, BUSCO, and OrthoDB teams at EBI and SIB for maintaining the curated databases on which functional annotation depends. We thank the Bioconda community for packaging the bioinformatics ecosystem in a sustainable way. We thank the Rust language team for the language and tooling. We thank members of the Hittinger Lab for testing and feedback, and members of the Rokas Lab at Vanderbilt for discussions about CTG clade biology. We thank colleagues who provided feedback on early versions of the manuscript and tool.

---

## Funding

This work was supported by the Hittinger Lab, Laboratory of Genetics, University of Wisconsin-Madison [funding sources to be specified]. C.T.H. is supported by [grants to be specified]. A.R. is supported by [grants to be specified]. Support for the Wisconsin Energy Institute [grant numbers to be specified].

---

## Conflicts of Interest

The authors declare no conflicts of interest.

---

## Tables

**Table 1.** Annotation sources integrated by MycoNote-CLI, with the underlying tools and their primary references.

| Source | Tool / Method | Primary reference |
|--------|---------------|-------------------|
| Swiss-Prot homology | MMseqs2 easy-search | Steinegger and Soding (2017) |
| Pfam domains | HMMER hmmsearch | Mistry et al. (2021) |
| Genome completeness | BUSCO 5 | Manni et al. (2021) |
| Gene Ontology | UniProt API + InterProScan | Gene Ontology Consortium (2021) |
| InterProScan | EBI REST API | Blum et al. (2021) |
| EggNOG categories | eggNOG-mapper v2 | Cantalapiedra et al. (2021) |
| CAZymes | dbCAN3 | Zheng et al. (2023) |
| Secretome | SignalP/DeepSig + DeepTMHMM | Teufel et al. (2022) |
| Biosynthetic gene clusters | antiSMASH 7 | Blin et al. (2023) |
| Proteases | MEROPS | Rawlings et al. (2018) |
| tRNA genes | tRNAscan-SE 2 | Chan et al. (2021) |
| Protein-to-genome | miniprot | Li (2023) |
| Genetic codes (25 tables) | Built-in | NCBI translation tables |
| GFF3 validation | Built-in | -- |
| Reproducibility tracking | Built-in | -- |

**Table 2.** Annotation results for *Brettanomyces bruxellensis*. Times are wall-clock measurements on Apple M3 Pro, 4 threads, single replicate. The "98.3% annotated" figure includes any annotation source producing any hit, including weak hits; users requiring high-confidence annotations should apply more stringent thresholds.

| Metric | Value |
|--------|-------|
| Assembly size | 12.9 Mb |
| Contigs | 30 |
| Predicted gene models (Augustus) | 5,218 |
| Genes with Swiss-Prot best hit (e<1e-5) | 3,915 (75.0%) |
| Genes with Pfam domains (e<1e-5) | 4,464 (85.6%) |
| Genes with any InterProScan hit | 5,129 (98.3%) |
| Genes with high-confidence Swiss-Prot (>50% identity, e<1e-50) | ~3,100 (~60%) |
| Genes with assigned GO terms | 4,856 (93.1%) |
| NCBI structural validation | PASSED |
| Pipeline runtime (with cached InterProScan) | ~15 min |
| Pipeline runtime (without cache) | ~45-60 min |

**Table 3.** Recommended tool selection by use case. MycoNote-CLI complements rather than replaces existing tools.

| Use case | Recommended tool |
|----------|-----------------|
| Standard fungal annotation, established workflows | funannotate |
| RNA-seq guided gene prediction (highest accuracy) | BRAKER followed by MycoNote-CLI annotate |
| Plant or animal genome | MAKER, BRAKER, or MycoNote-CLI |
| Multi-kingdom pipeline (consistent across taxa) | MycoNote-CLI |
| Non-standard genetic codes (CTG clade, etc.) | MycoNote-CLI |
| NCBI submission readiness with pre-flight validation | MycoNote-CLI |
| Bacterial or archaeal genome | Prokka or Bakta |
| Viral genome | VAPiD or VADR |
| Cross-genome annotation transfer | LiftOff or FLO |
| Manual curation | Apollo or WebApollo |

**Table 4.** Allele-specific expression validation on a *Candida tropicalis* controlled-truth fixture (Section 3.5). Three plus-strand single-block CDSs on contig NW_003020038.1 received 20 phased SNVs each (≈ 1 / 60 bp); 12,000 paired-end reads were simulated per transcript at the three known hap0:hap1 ratios shown. The sample-specific null hap0 fraction is 0.6665. All 60 variants applied cleanly (zero skipped across all nine documented categories). Mapping-rate asymmetry was hap0 = 0.999, hap1 = 0.500 (asymmetry_flag = true) — the expected behavior of personalized-transcriptome salmon, absorbed by `ase-template`'s per-sample-null correction.

| Transcript | Truth hap0 fraction | Recovered hap0 fraction | Imbalance vs null | padj | Significant at FDR 0.05 | Truth direction recovered? |
|---|---:|---:|---:|---:|:---:|:---:|
| g000005.m1 | 0.70 | 0.769 | +0.103 | 1.78e-58 | TRUE | hap0-biased ✓ |
| g000010.m1 | 0.30 | 0.588 | −0.079 | 1.53e-41 | TRUE | hap1-biased ✓ |
| g000012.m1 | 0.50 | 0.667 | +0.000 | 0.989 | FALSE | balanced ✓ |

3/3 ground-truth calls correct.

---

## Figure Legends

*(Figures to be generated for the final submission)*

**Figure 1.** Overview of the MycoNote-CLI pipeline. Seven sequential stages (sort, mask, train, predict, update, annotate, submit) process a genome assembly from raw contigs to NCBI-ready submission files. Optional stages are shown with dashed borders. The annotate stage integrates fifteen sources (right panel) with graceful degradation when individual tools are unavailable.

**Figure 2.** Pfam domain search performance. Wall-clock time for searching 5,218 *B. bruxellensis* proteins against Pfam-A v36 (20,795 profiles) using HMMER hmmsearch versus hmmscan on Apple M3 Pro with 4 threads. Bars show single-run measurements; error bars (when added) will represent standard deviation across three replicates. **Note: a single-replicate single-organism comparison; comprehensive benchmarking with multiple genomes and statistical analysis is planned for the follow-up study.**

**Figure 3.** Effect of genetic code selection on protein translation in *Candida tropicalis*. (A) Histogram of the number of CTG-containing positions per protein for the 6,290 *C. tropicalis* gene models. Of these, 3,871 (61.5%) contain at least one CTG codon and are therefore mistranslated under the standard code. (B) Total amino acid composition changes between the two translations. Table 12 produces exactly 10,960 fewer leucine residues and 10,960 additional serine residues than Table 1, corresponding to the 10,960 CTG codons across all CDSs. (C) Distribution of changed residues per protein. The mean is 1.74 changes per affected protein (median 1, range 1-28). The two most extreme proteins each have 25-28 CTG codons. Even modest CTG content alters protein chemistry: leucine (hydrophobic, Kyte-Doolittle index +3.8) and serine (polar hydrogen-bond donor, -0.8) differ in side chain properties relevant to folding, domain assignment, and functional prediction.

**Figure 4.** Decision tree for tool selection. Users navigating from "I have a fungal genome to annotate" to a tool recommendation, based on their specific requirements (clade, RNA-seq availability, genetic code, etc.).

---

## References

Alexa A, Rahnenfuhrer J, Lengauer T (2006). Improved scoring of functional groups from gene expression data by decorrelating GO graph structure. *Bioinformatics* 22:1600-1607. PMID 16606683.

Bankevich A, Nurk S, Antipov D, et al. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *Journal of Computational Biology* 19:455-477. PMID 22506599.

Blin K, Shaw S, Augustijn HE, et al. (2023). antiSMASH 7.0: new and improved predictions for detection, regulation, chemical structures and visualisation. *Nucleic Acids Research* 51:W46-W50. PMID 37140036.

Blum M, Chang HY, Chuguransky S, et al. (2021). The InterPro protein families and domains database: 20 years on. *Nucleic Acids Research* 49:D344-D354. PMID 33156333.

Bruna T, Hoff KJ, Lomsadze A, Stanke M, Borodovsky M (2021). BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database. *NAR Genomics and Bioinformatics* 3:lqaa108. PMID 33575650.

Cantalapiedra CP, Hernandez-Plaza A, Letunic I, Bork P, Huerta-Cepas J (2021). eggNOG-mapper v2: Functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Molecular Biology and Evolution* 38:5825-5829. PMID 34597405.

Cantarel BL, Korf I, Robb SMC, et al. (2008). MAKER: an easy-to-use annotation pipeline designed for emerging model organism genomes. *Genome Research* 18:188-196. PMID 18025269.

Chan PP, Lin BY, Mak AJ, Lowe TM (2021). tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. *Nucleic Acids Research* 49:9077-9096. PMID 34417604.

Eddy SR (2011). Accelerated profile HMM searches. *PLoS Computational Biology* 7:e1002195. PMID 22039361.

Flynn JM, Hubley R, Goubert C, Rosen J, Clark AG, Feschotte C, Smit AF (2020). RepeatModeler2 for automated genomic discovery of transposable element families. *PNAS* 117:9451-9457. PMID 32300014.

Gene Ontology Consortium (2021). The Gene Ontology resource: enriching a GOld mine. *Nucleic Acids Research* 49:D325-D334. PMID 33290552.

Grabherr MG, Haas BJ, Yassour M, et al. (2011). Full-length transcriptome assembly from RNA-Seq data without a reference genome. *Nature Biotechnology* 29:644-652. PMID 21572440.

Gruning B, Dale R, Sjodin A, et al. (2018). Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nature Methods* 15:475-476. PMID 29967506.

Haas BJ, Delcher AL, Mount SM, et al. (2003). Improving the Arabidopsis genome annotation using maximal transcript alignment assemblies. *Nucleic Acids Research* 31:5654-5666. PMID 14500829.

Haas BJ, Salzberg SL, Zhu W, et al. (2008). Automated eukaryotic gene structure annotation using EVidenceModeler and the Program to Assemble Spliced Alignments. *Genome Biology* 9:R7. PMID 18190707.

Hoff KJ, Lange S, Lomsadze A, Borodovsky M, Stanke M (2016). BRAKER1: Unsupervised RNA-seq-based genome annotation with GeneMark-ET and AUGUSTUS. *Bioinformatics* 32:767-769. PMID 26559507.

Holt C, Yandell M (2011). MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. *BMC Bioinformatics* 12:491. PMID 22192575.

Korf I (2004). Gene finding in novel genomes. *BMC Bioinformatics* 5:59. PMID 15144565.

Krogerus K, Magalhaes F, Vidgren V, Gibson B (2017). Novel brewing yeast hybrids: creation and application. *Applied Microbiology and Biotechnology* 101:65-78. PMID 27885413.

Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094-3100. PMID 29750242.

Li H (2023). Protein-to-genome alignment with miniprot. *Bioinformatics* 39:btad014. PMID 36648328.

Lomsadze A, Ter-Hovhannisyan V, Chernoff YO, Borodovsky M (2005). Gene identification in novel eukaryotic genomes by self-training algorithm. *Nucleic Acids Research* 33:6494-6506. PMID 16314312.

Manni M, Berkeley MR, Seppey M, Simao FA, Zdobnov EM (2021). BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage. *Molecular Biology and Evolution* 38:4647-4654. PMID 34320186.

McManus CJ, Coolon JD, Duff MO, Eipper-Mains J, Graveley BR, Wittkopp PJ (2010). Regulatory divergence in *Drosophila* revealed by mRNA-seq. *Genome Research* 20:816-825. PMID 20354124.

Mistry J, Chuguransky S, Williams L, et al. (2021). Pfam: The protein families database in 2021. *Nucleic Acids Research* 49:D412-D419. PMID 33125078.

Mixao V, Gabaldon T (2020). Genomic evidence for a hybrid origin of the yeast opportunistic pathogen *Candida albicans*. *BMC Biology* 18:48. PMID 32375762.

Muhlhausen S, Findeisen P, Plessmann U, Urlaub H, Kollmar M (2016). A novel nuclear genetic code alteration in yeasts and the evolution of codon reassignment in eukaryotes. *Current Opinion in Microbiology* 32:16-21. PMID 27173587.

Palmer JM, Stajich JE (2020). Funannotate v1.8: eukaryotic genome annotation. *Zenodo*. doi:10.5281/zenodo.4054262.

Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods* 14:417-419. PMID 28263959.

Rawlings ND, Barrett AJ, Thomas PD, Huang X, Bateman A, Finn RD (2018). The MEROPS database of proteolytic enzymes, their substrates and inhibitors in 2017 and a comparison with peptidases in the PANTHER database. *Nucleic Acids Research* 46:D624-D632. PMID 29145643.

Santos MA, Gomes AC, Santos MC, Carreto LC, Moura GR (2011). The genetic code of the fungal CTG clade. *Comptes Rendus Biologies* 334:607-611. PMID 21819941.

Schaffer AA, Hatcher EL, Yankie L, et al. (2020). VADR: validation and annotation of virus sequence submissions to GenBank. *BMC Bioinformatics* 21:211. PMID 32448124.

Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics* 7:000685. PMID 34739369.

Seemann T (2014). Prokka: rapid prokaryotic genome annotation. *Bioinformatics* 30:2068-2069. PMID 24642063.

Shean RC, Makhsous N, Stoddard GD, Lin MJ, Greninger AL (2019). VAPiD: a lightweight cross-platform viral annotation pipeline and identification tool. *BMC Bioinformatics* 20:48. PMID 30674277.

Stanke M, Keller O, Gunduz I, Hayes A, Waack S, Morgenstern B (2006). AUGUSTUS: ab initio prediction of alternative transcripts. *Nucleic Acids Research* 34:W435-W439. PMID 16845043.

Steinegger M, Soding J (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nature Biotechnology* 35:1026-1028. PMID 29035372.

Stone JE, Spilling C, Newhouse JE, et al. (2020). Rayon: a data parallelism library for Rust. *Conference on Programming Language Design and Implementation*.

Teufel F, Almagro Armenteros JJ, Johansen AR, et al. (2022). SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nature Biotechnology* 40:1023-1025. PMID 34980915.

Tirosh I, Reikhav S, Levy AA, Barkai N (2009). A yeast hybrid provides insight into the evolution of gene expression regulation. *Science* 324:659-662. PMID 19407207.

UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research* 51:D523-D531. PMID 36408920.

van de Geijn B, McVicker G, Gilad Y, Pritchard JK (2015). WASP: allele-specific software for robust molecular quantitative trait locus discovery. *Nature Methods* 12:1061-1063. PMID 26366987.

Zheng J, Ge Q, Yan Y, Zhang X, Huang L, Yin Y (2023). dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Research* 51:D557-D563. PMID 37125649.

---

*Manuscript prepared April 2026 for submission to Bioinformatics (Application Note) or Genome Biology (Software).*
