# myconote sort

Standardises scaffold names and removes sequences below a minimum length threshold.

```bash
myconote sort --genome mygenome.fas --out 01_sorted/
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--genome` | required | Input genome FASTA |
| `--out` | required | Output directory |
| `--min-length` | 500 | Minimum scaffold length (bp) to retain |
| `--prefix` | `scaffold` | Scaffold name prefix |
| `--threads` | 4 | CPU threads |

## Output files

| File | Description |
|------|-------------|
| `genome.fas` | Renamed, filtered genome FASTA |
| `sort_report.txt` | Scaffold count, total length, N50 before/after |
| `name_map.tsv` | Mapping of original → new scaffold names |

## Notes

!!! tip
    Always run `sort` first. Downstream tools require clean, consistently named scaffolds.

Scaffolds are renamed to `scaffold_0001`, `scaffold_0002`, etc. by default. Use `--prefix` to customise (e.g., `--prefix chr` for chromosome-scale assemblies).
