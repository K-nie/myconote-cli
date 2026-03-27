# myconote stats

Report assembly and annotation statistics.

```bash
myconote stats --genome genome.fas --gff annotation.gff3
```

Reports include scaffold N50, gene count, average gene length, exon count, BUSCO completeness score, and GC content.

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--genome` | required | Genome FASTA |
| `--gff` | — | Annotation GFF3 (optional) |
| `--busco-lineage` | `fungi_odb10` | BUSCO lineage dataset |
| `--out` | stdout | Write report to file |
| `--format` | `text` | `text`, `json`, or `csv` |
