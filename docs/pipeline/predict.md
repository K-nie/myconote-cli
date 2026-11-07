# myconote-cli predict

```bash
myconote-cli predict --help
```

See the [Pipeline Overview](overview.md) for how `predict` fits into the full workflow.

## Choosing the Augustus species model (`--clade`)

Augustus needs a species model. By default a fungal run uses
`saccharomyces_cerevisiae_S288C`, but that is a Saccharomycotina model and a
poor choice off that clade — so instead of hard-coding it, pass `--clade` to get
a clade-appropriate default:

```bash
myconote-cli predict genome_masked.fa --kingdom fungi --clade pezizomycotina
```

| `--clade` | Default model |
|-----------|---------------|
| `saccharomycotina` | `saccharomyces_cerevisiae_S288C` |
| `pezizomycotina` | `aspergillus_nidulans` |
| `basidiomycota` | `cryptococcus_neoformans_neoformans_JEC21` |

An explicit `--species <model>` always wins. A fungal run with neither
`--species` nor `--clade` still uses S288C but prints a warning (BUSCO
completeness can fall below 15 % off-Saccharomycotina).

## Protein hints (`--protein-hints`)

Add protein evidence from related species as Augustus hints via miniprot. The
protein DB is prefiltered with `diamond blastx` so a whole-kingdom OrthoDB set
does not stall miniprot:

```bash
myconote-cli predict genome_masked.fa --kingdom fungi --species aspergillus_nidulans \
    --protein-hints orthodb_fungi.fa --protein-hints-min-identity 70
```

| Flag | Default | Description |
|------|---------|-------------|
| `--protein-hints <proteins.fa>` | off | Protein FASTA → miniprot → Augustus hints |
| `--protein-hints-evalue <E>` | `1e-5` | diamond e-value for the prefilter |
| `--protein-hints-max <N>` | `50000` | Max proteins fed to miniprot (`0` = unbounded) |
| `--protein-hints-min-identity <P>` | `0` | Keep only prefilter hits ≥ P% identity — raise it to drop distant, low-identity homologs (close relatives beat remote proteins) |

This path is non-fatal: if diamond/miniprot are absent or nothing survives the
prefilter, prediction proceeds ab-initio. (Distinct from `--protein-fasta`,
which feeds GeneMark EP+/ETP+ via ProtHint — see below.)

## Self-training Augustus (`--augustus-training self`)

`--augustus-training stock` (default) uses a curated pre-trained species.
`--augustus-training self` instead trains a genome-specific model on the fly
(first-pass prediction → `new_species.pl` → `etraining`). It needs the Augustus
companion scripts and the `species/generic/` template from
`setup --db augustus-fungi`; if either is missing, `predict` prints a loud
banner and falls back to the stock species (it never silently pretends to
train). On fungi with a curated model, `stock` usually wins — self-train earns
its keep on a novel fungus with no close stock species.

`--augustus-training busco` trains on the genome's complete single-copy BUSCO
genes instead. It uses the `--kingdom` lineage (`fungi_odb10`) by default;
pass `--busco-lineage <name|dir>` to pick a narrower dataset such as
`dikarya_odb10`, or a path to a lineage directory.

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

`--use-braker` short-circuits the standard ab-initio + agreement-weighted consensus stack;
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
