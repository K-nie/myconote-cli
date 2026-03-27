# myconote blast

Run BLAST searches against local or remote databases.

```bash
myconote blast --query proteins.faa --db swissprot --out blast/
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--query` | required | Query FASTA (nucleotide or protein) |
| `--db` | `swissprot` | Database name or path |
| `--out` | required | Output directory |
| `--evalue` | `1e-5` | E-value threshold |
| `--threads` | 4 | CPU threads |
| `--remote` | false | Use NCBI remote BLAST (requires `NCBI_EMAIL`) |

## Output files

| File | Description |
|------|-------------|
| `blast_results.tsv` | Tabular BLAST output (fmt 6) |
| `blast_summary.txt` | Hit count, coverage statistics |
