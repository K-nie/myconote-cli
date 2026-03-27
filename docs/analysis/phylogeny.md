# myconote phylogeny

Infer a maximum-likelihood phylogenetic tree from a multiple sequence alignment using IQ-TREE2.

```bash
myconote phylogeny --alignment proteins.aln --out phylogeny/
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--alignment` | required | Multiple sequence alignment (FASTA or PHYLIP) |
| `--out` | required | Output directory |
| `--model` | `MFP` | Substitution model (`MFP` = ModelFinder Plus auto-select) |
| `--bootstrap` | 1000 | UFBoot2 replicates |
| `--alrt` | 1000 | SH-aLRT branch support replicates |
| `--threads` | 4 | CPU threads |
| `--prefix` | `tree` | Output file prefix |
| `--partition` | — | Partition file for multi-locus analysis |

## Output files

| File | Description |
|------|-------------|
| `tree.treefile` | Best ML tree in Newick format |
| `tree.iqtree` | Full IQ-TREE2 log with model selection |
| `tree.contree` | Consensus tree |
| `tree.log` | Run log |

## Example: multi-gene phylogeny

```bash
# Align proteins first
cat 06_annotation/proteins.faa | mafft --auto - > proteins.aln

# Infer tree
myconote phylogeny \
  --alignment proteins.aln \
  --bootstrap 1000 \
  --threads 8 \
  --out phylogeny/
```

## Notes

!!! info "Model selection"
    The default `MFP` model runs ModelFinder Plus to automatically select the best-fit substitution model. This adds ~10–30% to run time but produces better-supported topologies.

IQ-TREE2 must be installed and on your PATH. Verify with `myconote check`.
