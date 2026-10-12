# myconote-cli annotate

```bash
myconote-cli annotate --help
```

See the [Pipeline Overview](overview.md) for how `annotate` fits into the full workflow.

## Usage

`annotate` takes the gene GFF3 as a positional argument and the genome via
`--fasta`:

```bash
myconote-cli annotate predict_out/consensus.gff3 \
    --fasta genome_masked.fa \
    --genetic-code 1 \
    --threads 8 \
    --output annotate_out/
```

Output is written as `annotated.gff3` (plus per-source tables and a BUSCO
summary) in the output directory.

## Functional sources

The **core sources run by default**: MMseqs2 (Swiss-Prot homology), Pfam domains
(`--no-pfam` to skip), and BUSCO completeness (`--no-busco` to skip). The rest
are **opt-in** — add the flag to enable each:

| Source | Flag | Database | Provides |
|--------|------|----------|----------|
| MMseqs2 | default on | Swiss-Prot | Gene names, best-hit species |
| Pfam | default on (`--no-pfam`) | Pfam-A (hmmsearch) | Protein domain families |
| BUSCO | default on (`--no-busco`) | fungi_odb10 | Completeness assessment |
| EggNOG | `--eggnog` | EggNOG / COG | GO terms, KEGG, COG categories |
| CAZyme | `--cazyme` | dbCAN (DIAMOND) | Carbohydrate-active enzymes |
| MEROPS | `--merops` | MEROPS (DIAMOND) | Protease families |
| antiSMASH | `--antismash` | antiSMASH | Secondary-metabolite clusters |
| InterProScan | `--interproscan` | EBI REST | Domain + GO cross-refs |
| secretome | `--secretome` | DeepSig + DeepTMHMM | Secreted-protein prediction |
| tRNAscan-SE | `--trnascan` | — | tRNA gene prediction |

```bash
# Core (default): MMseqs2 + Pfam + BUSCO
myconote-cli annotate predict_out/consensus.gff3 --fasta genome_masked.fa --output annotate_out/

# Full: add the opt-in sources
myconote-cli annotate predict_out/consensus.gff3 --fasta genome_masked.fa \
    --eggnog --cazyme --merops --antismash --secretome --trnascan --interproscan \
    --threads 8 --output annotate_out/
```

## EggNOG speed (v0.7.8)

The EggNOG step runs ~4× faster by default: `emapper.py` is called with
`--dmnd_iterate no` (measured ~46 min vs ~3.4 h on a 10,655-protein proteome, at
identical recall). To restore the slower, maximum-recall DIAMOND iterate path,
pass `--eggnog-iterate`.
