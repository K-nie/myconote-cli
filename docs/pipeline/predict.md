# myconote predict

```bash
myconote predict --help
```

See the [Pipeline Overview](overview.md) for how `predict` fits into the full workflow.

## Choosing a GeneMark mode

`gmes_petap.pl` ships four self-training algorithms. Pick one with
`--genemark-mode <mode>`; if you don't pass the flag, GeneMark is skipped.

| Mode | Flag                       | Extra inputs required           | Quality | When to use |
|------|----------------------------|---------------------------------|---------|-------------|
| ES   | `--genemark-mode es`       | none                            | Lowest  | No RNA-seq, no protein database. Fast smoke run. |
| ET   | `--genemark-mode et`       | `--genemark-hints <introns.gff>`| Medium  | RNA-seq available; HISAT2 / STAR splice-junctions converted to GFF. |
| EP+  | `--genemark-mode ep`       | `--protein-fasta <proteins.fa>` | High    | No RNA-seq but a curated protein DB exists. **Recommended for novel CTG-clade fungi** — the alt-yeast nuclear code makes published predictors unreliable, and protein evidence anchors the training. |
| ETP+ | `--genemark-mode etp`      | both `--genemark-hints` *and* `--protein-fasta` | Highest | Both data types in hand. Same accuracy ceiling as BRAKER3 but driven by GeneMark alone. |

EP+ and ETP+ run **ProtHint** internally to convert the protein FASTA into
splice-site hints (`prothint.gff` + `evidence.gff`) before invoking
`gmes_petap.pl --EP` / `--ETP`. ProtHint is bundled with the GeneMark-ES
installer tarball under `<install>/ProtHint/bin/` and with the bioconda
`braker3` package; there is no standalone bioconda recipe.

For the protein database, OrthoDB fungi is the standard choice:

```bash
wget https://data.orthodb.org/v11/odb11v0_proteins_fungi.fa.gz
gunzip odb11v0_proteins_fungi.fa.gz
myconote-cli predict genome_masked.fa \
    --kingdom fungi \
    --genemark-mode ep \
    --protein-fasta odb11v0_proteins_fungi.fa \
    --threads 16
```

### Backward compatibility

The pre-0.6.0 flags still work and are translated automatically:

- `--genemark` (alone)        → `--genemark-mode es`
- `--genemark` + `--genemark-hints <gff>` → `--genemark-mode et`
- `--genemark-hints <gff>` (alone) → `--genemark-mode et`

`--genemark-mode` always wins when passed alongside the legacy flags.

### License note

GeneMark itself requires a free academic license from Georgia Tech
(<http://topaz.gatech.edu/GeneMark/>). `myconote-cli install predict` will
flag both `gmes_petap.pl` and `prothint.py` under "manual download" until the
licensed tarball is on `PATH`.

## Using BRAKER for maximum accuracy

When you have **both** RNA-seq alignments and a fungal protein database,
`--use-braker` is the recommended path. BRAKER is a complete predictor that
internally cross-trains Augustus and GeneMark with hints from both evidence
types; BRAKER3 is widely held to be the current gold standard for fungal gene
prediction quality.

```bash
myconote-cli predict genome_masked.fa \
    --kingdom fungi \
    --use-braker \
    --braker-rna-bam rnaseq_sample1.bam --braker-rna-bam rnaseq_sample2.bam \
    --braker-proteins odb11v0_proteins_fungi.fa \
    --threads 16
```

`--use-braker` short-circuits the standard ab-initio + EVM consensus stack;
MycoNote-CLI does **not** run Augustus / SNAP / GlimmerHMM / GeneMark / miniprot
when BRAKER is the engine, since BRAKER produces a complete `braker.gff3`
that flows unchanged into the downstream `update`, `annotate`, and `submit`
stages.

### Mode auto-detection

By default the BRAKER mode is selected from the inputs:

| Inputs                              | Detected mode |
|-------------------------------------|---------------|
| `--braker-rna-bam` only             | BRAKER1 (`--esmode` plus RNA hints) |
| `--braker-proteins` only            | BRAKER2 (`--epmode`) |
| both `--braker-rna-bam` + `--braker-proteins` | BRAKER3 (`--etpmode`) |
| neither                             | error — at least one evidence input is required |

Pass `--braker-mode <1|2|3>` to override the auto-detection (e.g. force a
BRAKER1 run even if you happen to have a protein FASTA on disk but want to
exclude it).

### Genetic code

`--genetic-code <n>` is forwarded to BRAKER as `--translation_table=<n>`
unchanged. Set this to **12** for Candida-clade (CTG) fungi or **6** for
ciliate-style codes; default is 1 (standard).

### Dependencies

`--use-braker` requires:

- `braker.pl` (bioconda `braker3` provides this; the package name is
  `braker3` even though the binary is `braker.pl`)
- Augustus + Augustus-config writeable (`AUGUSTUS_CONFIG_PATH`); use
  `myconote-cli setup --db augustus-fungi` to populate
  `~/.myconote/augustus_config/` and export the variable
- A licensed GeneMark install (EP+/ETP+ for BRAKER2/3)
- ProtHint (bundled with the BRAKER3 conda package)

Run `myconote-cli check predict` to see which of these are missing.
