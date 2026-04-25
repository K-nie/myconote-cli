# myconote update

```bash
myconote-cli update --help
```

`update` refines an initial gene annotation produced by `predict` using
RNA-seq evidence. The default path runs PASA (with a lightweight
minimap2-coverage fallback when PASA is unavailable) to add UTRs,
correct intron/exon boundaries, and add alternative isoforms. See the
[Pipeline Overview](overview.md) for how `update` fits into the full
workflow.

## Inputs

At least one of the following RNA-seq inputs is required:

- `--rna-r1` (and `--rna-r2` for paired-end) — raw FASTQ reads. Trinity
  is invoked internally to assemble transcripts.
- `--transcripts <file>` — pre-assembled transcript FASTA (Trinity
  output, IsoSeq, or any cDNA FASTA). Skips the Trinity step.
- `--rna-bam <file>` — pre-aligned BAM. Used by the lightweight UTR
  fallback only; PASA re-aligns internally.

## Using Kallisto for transcript filtering

`update --kallisto` adds an abundance pre-filter step that mirrors
funannotate's update flow: every transcript is assigned a TPM by
`kallisto quant`, transcripts below `--kallisto-min-tpm` (default
**1.0** TPM) are dropped, and only the surviving IDs are written to a
`pasa_kallisto_filter.txt` sidecar passed to PASA's update step. This
removes spurious Trinity assemblies and assembly artifacts that would
otherwise feed PASA noise and inflate isoform counts.

```bash
myconote-cli update predict_out/consensus.gff3 \
    --fasta genome.fa \
    --rna-r1 reads_R1.fastq.gz \
    --rna-r2 reads_R2.fastq.gz \
    --kallisto \
    --kallisto-min-tpm 1.0 \
    --threads 16
```

### When to use it

- **Use `--kallisto`** when you have extensive RNA-seq coverage and
  expect a substantial fraction of Trinity assemblies to be assembly
  artifacts. The 1.0-TPM cut is funannotate's tested-default and
  works well at typical Illumina depths (≥30M paired reads).
- **Skip `--kallisto`** when RNA-seq coverage is low (single-replicate
  or fewer than ~10M reads) — the 1.0-TPM threshold will over-filter
  and drop legitimate low-expressed transcripts. Either run without
  `--kallisto` or drop the cut substantially with
  `--kallisto-min-tpm 0.1`.
- **Skip `--kallisto`** when your only RNA-seq input is a pre-aligned
  BAM (`--rna-bam`). Kallisto requires raw reads to pseudoalign; a BAM
  is rejected with a clear error.

### Outputs added under `--kallisto`

| Path | Contents |
|------|----------|
| `update_out/kallisto/transcripts.idx` | Kallisto index over the transcript FASTA |
| `update_out/kallisto/sample_NN/abundance.tsv` | Per-sample kallisto quant output |
| `update_out/kallisto/passing_transcripts.txt` | Transcript IDs that cleared the TPM filter |
| `update_out/pasa_kallisto_filter.txt` | The same passing-ID list at the top level |

### Dependencies

`--kallisto` requires `kallisto` on PATH. Install with:

```bash
myconote-cli install update          # installs kallisto + the rest of update's deps
# or directly:
conda install -c bioconda kallisto
```

If `--kallisto` is set and `kallisto` is missing, the run aborts with
a clear error and a pointer at `myconote-cli install update`. There is
no silent fallback.
