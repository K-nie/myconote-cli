# myconote ase

Quantify **allele-specific expression (ASE)** in heterozygous, hybrid, or polyploid fungal genomes. Builds personalized transcriptomes per haplotype from a phased VCF, then runs salmon independently against each haplotype so each read is assigned to the allele it matches better.

ASE is the rule, not the exception, for many fungi: diploid *Candida* isolates, *Saccharomyces* F1 hybrids, interspecies crosses, and polyploids all need this. Quantifying against a single reference under-counts the allele farther from reference — reference bias — and makes real cis-regulatory differences invisible.

## What this does

1. Parses a phased VCF. **Unphased heterozygous sites error out with a line number** — ASE requires phasing; we refuse to guess.
2. For each haplotype, applies its variants to the reference CDS to produce a personalized transcriptome FASTA.
3. Runs fastp on each sample once (QC is haplotype-independent), then salmon against each haplotype's index. Two salmon runs per sample.
4. Merges per-sample, per-haplotype `quant.sf` into a `transcript × <sample>.<hap>` counts matrix.
5. Computes a per-transcript `informative` flag (haplotypes actually differ in CDS sequence), per-haplotype variant counts, and a max mapping-rate asymmetry across samples.
6. Writes a reproducibility bundle (`ase_bundle.json`) with SHA256s of every input, tool versions, variant-application counts, and per-sample mapping rates.

## Prerequisites — install before running

You need **salmon** and **fastp** on your PATH:

```bash
conda install -c bioconda salmon fastp
```

Then prepare inputs:

1. **Reference CDS FASTA** (`cds.fa`) — from `myconote-cli convert --to cds genome.fa annotated.gff3`.
2. **Phased VCF** (`phased.vcf.gz`) — produced by your own variant-calling pipeline (GATK HaplotypeCaller → WhatsHap, or HapCUT2, or `bcftools call` with a phased genotype caller upstream). Unphased input is rejected.
3. **Annotated GFF3** (`annotated.gff3`) — same GFF3 you passed to `convert --to cds`. Needed to map genome coordinates in the VCF onto CDS positions. Must be sort/pad/clean: the same file that produced the CDS, no later edits.
4. **Sample sheet** (`samples.tsv`) — same schema as `quant`: `sample_id`, `fastq_1`, `fastq_2` (or blank for SE), optional `condition`, `batch`, and arbitrary extra columns.
5. **Reference genome FASTA** (`genome.fa`) — used as the decoy set in salmon's index.

## Quickstart

```bash
myconote-cli ase cds.fa \
    --vcf phased.vcf.gz \
    --gff3 annotated.gff3 \
    --samples samples.tsv \
    --genome genome.fa \
    -o ase_out
```

Then run the downstream binomial test with `ase-template`:

```bash
myconote-cli ase-template --ase-dir ase_out -o ase_analysis.R
Rscript ase_analysis.R
```

## Inputs

| flag | meaning |
|---|---|
| `<cds.fa>` | Reference CDS FASTA (positional, required). |
| `--vcf <phased.vcf[.gz]>` | Phased VCF. Unphased heterozygous sites cause a hard error with the offending line number. |
| `--gff3 <annotated.gff3>` | Same GFF3 that produced `<cds.fa>`. Required for genome → CDS coordinate mapping. |
| `--samples <sheet.tsv>` | Sample sheet (same schema as `quant`). |
| `--genome <genome.fa>` | Reference genome FASTA for decoy-aware indexing. |
| `--output <dir>` / `-o` | Output directory (default `ase_out`). |
| `--haplotype-names N1,N2` | Custom haplotype labels (default `hap0,hap1`). Common alternative: `paternal,maternal` or `Scer,Spar`. |
| `-k <n>` | salmon k-mer length (default `31`). |
| `--threads <n>` / `-t` | Threads per salmon / fastp run. |
| `--tmpdir <dir>` | Temp directory root for fastp-trimmed FASTQs (default `$TMPDIR`). |
| `--index-cache <dir>` | Index cache root (default `$MYCONOTE_INDEX_CACHE` or XDG). |
| `--keep-trimmed <dir>` | Persist fastp-trimmed FASTQs to this directory (off by default — deleted after quant). |
| `--max-indel-size <n>` | Skip indels longer than N bp, marked in `variants_skipped.tsv` (default `50`). |
| `--asymmetry-threshold <f>` | Flag samples where `|hap0 rate − hap1 rate| > f` (default `0.05`, i.e. 5 %). |
| `--seed <n>` | Deterministic seed (default `42`). |
| `--fastp <path>` | Override fastp binary. |
| `--salmon <path>` | Override salmon binary. |

## What gets skipped, and why

The VCF → personalized-CDS step is deliberately conservative — variants that could silently produce garbage are skipped and logged to `variants_skipped.tsv` with a reason:

- **unphased** — heterozygous genotype without the `|` phase separator. Hard error, not a skip.
- **structural** — non-SNV/MNP/indel variants (symbolic alleles, breakends).
- **multiallelic** — sites with multiple ALT alleles; not represented in the current scheme.
- **missing_genotype** — `.|.` or `./.` at this site.
- **indel_too_large** — indel length exceeds `--max-indel-size` (default 50 bp).
- **outside_cds** — variant position doesn't fall inside any CDS exon.
- **spans_exon_boundary** — indel spans an exon/intron boundary; introns aren't in the CDS.
- **ref_mismatch** — VCF's REF allele doesn't match the genome FASTA at that position (VCF / genome version skew).
- **in_cis_overlap** — multiple variants on the *same* haplotype overlap in CDS coordinates.

Counts for every category appear in `ase_bundle.json` under `variants` so downstream tools can summarize without re-deriving the list.

## Mapping-rate asymmetry

After both salmon runs finish, `|mapping_rate_hap0 − mapping_rate_hap1|` is computed for every sample. Anything over `--asymmetry-threshold` (default 5 %) emits a warning to stderr and sets `asymmetry_flag = true` in `ase_bundle.json`.

Strong asymmetry usually signals:

- A phasing error — e.g. large haplotype blocks with the wrong parent-of-origin.
- An assembly issue — one haplotype has many more accurate regions than the other.
- A sample swap — the sample sheet claims one parent as hap0 when it's actually the other.

Not necessarily biology. Chase down the warning before interpreting the counts.

## Outputs (under `<output-dir>`)

| file | contents |
|---|---|
| `cds_<hap>.fa` | One personalized CDS FASTA per haplotype (e.g. `cds_hap0.fa`, `cds_hap1.fa`). |
| `salmon/<sample>.<hap>/` | One salmon output directory per `sample × haplotype`. Each contains the usual `quant.sf`, `aux_info/`, `logs/`. |
| `ase_counts.tsv` | Wide counts matrix: rows = transcripts, columns = `<sample>.<hap>`. |
| `ase_tpm.tsv` | Same shape as `ase_counts.tsv`, TPM-normalized. |
| `ase_summary.tsv` | Per-transcript metadata: `informative` flag, `n_variants_hap0`, `n_variants_hap1`, `max_asymmetry`. |
| `variants_applied.tsv` | Per-variant audit: which haplotype, which transcript, CDS position, strand. |
| `variants_skipped.tsv` | Per-variant audit: skip reason per category above. |
| `ase_bundle.json` | Reproducibility manifest (SHA256 of every input + tool versions + every skip count). |
| `fastp/<sample>.json` | Raw fastp JSON report per sample. |

## What `ase` does NOT do

- **Doesn't phase.** Run a dedicated phaser (WhatsHap / HapCUT2) upstream. Unphased het sites are a hard error.
- **Doesn't variant-call.** Bring your own VCF.
- **Doesn't do per-site ASE.** This is transcript-level ASE (WASP / GATK ASEReadCounter operate at the site level — different problem shape).
- **Doesn't run the statistical test.** See `myconote-cli ase-template` for the R script that does the binomial analysis on `ase_counts.tsv`.
- **Doesn't support >2 haplotypes in 0.5.0.** Diploid-only for now; the code path is written to generalize to higher ploidy in 0.5.x.

## Troubleshooting

- **"Unphased heterozygous genotype at VCF line N"** — your VCF has `0/1` where ASE needs `0|1`. Run WhatsHap or HapCUT2 first, or use a caller that emits phased output (DeepVariant + read-backed phasing, GATK + `VariantFiltration` → `WhatsHap phase`).
- **"REF mismatch at chrX:12345"** — VCF was called against a different genome build than `--genome`. Re-lift the VCF, or regenerate it against the same FASTA.
- **"spans_exon_boundary" count is very high** — your indels fall on intron boundaries, which are outside the CDS. Not usually a bug; it's the reality of intronic variation in the VCF.
- **Low mapping rate on one haplotype** — check the asymmetry warning in the bundle. If it's a phasing artifact, re-phase with longer reads or more stringent thresholds. If it's a real assembly gap, document it and proceed with the informative transcripts only.
- **"salmon: index not found"** — `--index-cache` is not writable, or the cache slot got stale. Delete the cache directory or point to a fresh one.

## Related

- `myconote-cli convert --to cds` — prerequisite; produces `cds.fa`.
- `myconote-cli ase-template` — emits the binomial-test R script for `ase_counts.tsv`.
- `myconote-cli quant` — single-reference quantification (the non-ASE counterpart).
- salmon documentation: <https://salmon.readthedocs.io/>
- WhatsHap phaser: <https://whatshap.readthedocs.io/>
