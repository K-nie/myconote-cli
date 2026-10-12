# myconote-cli stats

Report annotation statistics from a GFF3 (gene/exon/transcript counts, lengths,
structure). The GFF3 is a positional argument.

```bash
myconote-cli stats genes.gff3 --taxon fungi
myconote-cli stats genes.gff3 --format json > stats.json
myconote-cli stats genes.gff3 --chromosome scaffold_1
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `<annotation.gff3>` | required | Input annotation GFF3 (positional) |
| `--format <json\|csv\|human>` | `human` | Output format |
| `--taxon <group>` | — | Taxonomic group for benchmarking context |
| `--chromosome <chr>` | — | Focus on one or more specific sequences |
| `--region <chr:start-end>` | — | Focus on a specific genomic region |
| `--exclude <chr>` | — | Exclude a sequence from the stats |
| `--primary-only` | off | Collapse isoforms to the primary transcript only |

!!! note
    `stats` works on the annotation GFF3. For assembly-level metrics (contig
    N50, GC, BUSCO completeness), use the `sort` report and `annotate`'s BUSCO
    output.
