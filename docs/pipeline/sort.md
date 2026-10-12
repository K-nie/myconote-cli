# myconote-cli sort

Reorders contigs longest-first and filters out sequences below a minimum length.
By **default it is lossless**: original contig seqids are preserved end-to-end.
Renaming to NCBI-style short IDs is opt-in via `--ncbi-clean`.

```bash
# Default (lossless): reorder + drop anything under 500 bp, keep original names
myconote-cli sort assembly.fa --min-length 500

# NCBI-ready: also rename to scaffold_N and save the lift-over table
myconote-cli sort assembly.fa --min-length 500 --ncbi-clean --rename-table rename.tsv
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `<genome.fa>` | required | Input genome FASTA (positional) |
| `--output`, `-o` | `<input>_sorted.fa` | Output FASTA file |
| `--min-length <bp>` | `0` | Discard contigs shorter than this |
| `--ncbi-clean` (alias `--rename-contigs`) | off | Rename contigs to `<prefix>_N` for NCBI's short-name requirement |
| `--prefix <str>` | `scaffold` | ID prefix when `--ncbi-clean` is set |
| `--rename-table <file>` | — | Write the old→new ID table (TSV); use with `--ncbi-clean` |
| `--keep-desc` | off (strip) | Keep the original description after the ID |
| `--sort-by-name` | off (by length) | Sort alphabetically instead of by length |

## Output files

| File | Description |
|------|-------------|
| `<input>_sorted.fa` | Reordered, length-filtered genome FASTA (original seqids unless `--ncbi-clean`) |
| `<rename-table>.tsv` | old→new ID mapping — only when `--ncbi-clean` + `--rename-table` are used |

## Notes

!!! tip
    Run `sort` first to drop junk-length contigs and get a deterministic order.
    The lossless default means your downstream GFF3 keeps the accessions you
    know; only add `--ncbi-clean` when preparing an NCBI submission.

With `--ncbi-clean`, contigs are renamed to `scaffold_1`, `scaffold_2`, … (not
zero-padded). Use `--prefix` to customise (e.g. `--prefix chr` for
chromosome-scale assemblies), and `--rename-table` to keep the lift-over map so
`submit` (and your own records) can trace back to the originals.
