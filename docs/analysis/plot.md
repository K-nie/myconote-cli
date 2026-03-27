# myconote plot

Generate genome visualisations — circular genome maps or linear track plots.

```bash
myconote plot circular --genome genome.fas --gff annotation.gff3 --out plots/
myconote plot linear   --genome genome.fas --gff annotation.gff3 --out plots/
```

## Circular plots

Radial plots showing chromosome arcs, gene density rings, and strand tracks.

| Flag | Default | Description |
|------|---------|-------------|
| `--genome` | required | Genome FASTA |
| `--gff` | required | Annotation GFF3 |
| `--out` | required | Output directory |
| `--format` | `png` | Output format (`png`, `svg`, `pdf`) |
| `--dpi` | 300 | Resolution for raster outputs |

## Linear plots

Track-based linear visualisations, similar to IGV or UCSC Genome Browser tracks.

```bash
myconote plot linear \
  --genome genome.fas \
  --gff annotation.gff3 \
  --region scaffold_0001:1-50000 \
  --out plots/
```

## Output files

| File | Description |
|------|-------------|
| `circular_plot.png` | Full-genome circular map |
| `linear_plot.png` | Linear track view |
