# Quick Start

Annotate a eukaryotic genome from start to finish. Estimated time: 30-60 minutes (excluding database downloads).

**New to myconote-cli?** Run `myconote-cli learn` for an interactive tutorial that teaches each concept step by step.

---

## 1. Verify your installation

```bash
myconote-cli --version       # check binary is installed
myconote-cli check           # verify external tools
myconote-cli setup --check   # verify databases
```

If tools or databases are missing:
```bash
myconote-cli install --yes   # install all missing tools via conda
myconote-cli setup           # download all databases (~2.5 GB)
```

---

## 2. Sort and clean scaffolds

```bash
myconote-cli sort assembly.fa --min-length 500
```

Renames scaffolds to clean IDs (scaffold_001, scaffold_002...), sorts by length, and filters short contigs.

---

## 3. Mask repeats

```bash
myconote-cli mask assembly_sorted.fa --engine repeatmodeler --threads 8
```

Builds a de novo repeat library with RepeatModeler2, then soft-masks the genome with RepeatMasker. Soft-masked bases are lowercased so gene finders can skip them.

---

## 4. Train gene predictors (optional)

If you have RNA-seq data:

```bash
myconote-cli train assembly_masked.fa \
  --left R1.fastq.gz --right R2.fastq.gz \
  --species my_organism --threads 8
```

This assembles transcripts with Trinity, aligns them with minimap2, builds a PASA database, and trains Augustus + SNAP on your organism's gene structures.

!!! tip
    If you don't have RNA-seq, skip this step. Augustus will use its pre-trained model for your kingdom.

---

## 5. Predict gene models

```bash
myconote-cli predict assembly_masked.fa \
  --kingdom fungi \
  --locus-prefix MYORG \
  --threads 8
```

Runs Augustus + SNAP, merges predictions with the Evidence Modeler consensus, and produces a clean GFF3 with sequential locus tags.

**Advanced options:**
```bash
# Add protein evidence for better predictions
myconote-cli predict assembly_masked.fa --kingdom fungi \
  --protein-fasta swissprot.fasta --locus-prefix MYORG

# Use custom evidence weights
myconote-cli predict assembly_masked.fa --weights weights.toml

# Handle diploid/polyploid genomes
myconote-cli predict assembly_masked.fa --ploidy 2
```

---

## 6. Functional annotation

```bash
myconote-cli annotate predict_out/consensus.gff3 \
  --fasta assembly_masked.fa \
  --kingdom fungi \
  --trnascan \
  --threads 8
```

Runs MMseqs2 (Swiss-Prot homology), Pfam (domain search via hmmsearch), BUSCO (completeness), GO terms, and tRNAscan-SE (tRNA genes).

**Enable additional annotation sources:**
```bash
myconote-cli annotate predict_out/consensus.gff3 \
  --fasta assembly_masked.fa \
  --eggnog --cazyme --secretome --antismash --merops \
  --interproscan --email you@email.edu \
  --genetic-code 12    # for Candida CTG clade
```

---

## 7. Validate and prepare NCBI submission

```bash
# Check for errors first
myconote-cli submit annotate_out/annotated.gff3 \
  --fasta assembly_masked.fa \
  --organism "Genus species" --validate-only

# Generate submission files
myconote-cli submit annotate_out/annotated.gff3 \
  --fasta assembly_masked.fa \
  --organism "Genus species" \
  --strain "CBS 123" \
  --locus-prefix MYORG \
  --bioproject PRJNA123456
```

---

## 8. Explore your results

```bash
# Summary statistics with taxonomic benchmarking
myconote-cli stats annotate_out/annotated.gff3 --taxon fungi

# Generate a genome map
myconote-cli plot annotate_out/annotated.gff3 --type circular --output genome_map.png

# Interactive genome browser
myconote-cli view annotate_out/annotated.gff3 --fasta assembly_masked.fa

# Convert formats
myconote-cli convert annotate_out/annotated.gff3 --to genbank --fasta assembly_masked.fa
```

---

## Output files

| File | Description |
|------|-------------|
| `predict_out/consensus.gff3` | Predicted gene models |
| `annotate_out/annotated.gff3` | Gene models with functional annotations |
| `annotate_out/proteins.fa` | Predicted protein sequences |
| `annotate_out/annotations.tsv` | Full annotation table (for R/Python) |
| `annotate_out/annotation_report.txt` | Summary with gene counts and coverage |
| `submit_out/annotation.tbl` | NCBI feature table |

---

## Next steps

- Learn interactively: `myconote-cli learn`
- Build a phylogenetic tree: [phylogeny](analysis/phylogeny.md)
- Compare two genomes: `myconote-cli synteny a.gff3 b.gff3 --fasta1 a.fa --fasta2 b.fa`
- Full workshop tutorial: [Workshop Lesson](lesson.md)
