# myconote synteny

Two-genome ribbon diagram. Runs minimap2 asm-to-asm, chains colinear hits into consolidated syntenic blocks, and renders a self-contained HTML viewer with D3.js — no server or external hosting.

```bash
myconote-cli synteny <genome_a.gff3> <genome_b.gff3> [options]
```

## Quick examples

```bash
# Chromosome layout only (no FASTA → no ribbons, just the tracks)
myconote-cli synteny a.gff3 b.gff3 --output a_vs_b.html

# Full ribbon diagram with minimap2 alignment
myconote-cli synteny a.gff3 b.gff3 \
  --fasta1 a.fa --fasta2 b.fa \
  --label1 "A. niger" --label2 "A. fumigatus" \
  --output aniger_vs_afumigatus.html

# Label ribbons using a local gene-name TSV
myconote-cli synteny a.gff3 b.gff3 \
  --fasta1 a.fa --fasta2 b.fa \
  --names curated_names.tsv

# Auto-fetch gene names from NCBI for the given taxon
myconote-cli synteny a.gff3 b.gff3 \
  --fasta1 a.fa --fasta2 b.fa \
  --fetch-names --taxon 5061

# Keep the intermediate PAF for debugging or external dotplot tools
myconote-cli synteny a.gff3 b.gff3 --fasta1 a.fa --fasta2 b.fa --keep-paf
```

## What it does

1. **Parse chromosome sizes** from both GFF3 files — you always get a to-scale chromosome layout, even without FASTA.
2. **Align** the two genomes with `minimap2 -cx asm5 --cs` (asm-to-asm preset, ≥5 % divergence tolerated). Requires `minimap2` on `PATH`.
3. **Parse PAF** natively in Rust: query/target coords, strand, residue matches, block length, identity.
4. **Chain** adjacent colinear PAF hits into consolidated synteny blocks so the viewer shows real syntenic regions rather than the raw asm5 fragments.
5. **Extract gene intervals** from both GFF3s for the on-track gene overlay and for block-label resolution.
6. **Render** a single self-contained HTML file — D3.js from CDN, graceful offline error message if that fails.

## Requirements

| Requirement | When needed |
|---|---|
| Two GFF3 files | Always |
| `minimap2` in `PATH` | Only for alignment (pass `--fasta1`/`--fasta2`) |
| Two FASTA files | Only for alignment — without them, the output shows just chromosome tracks |
| Internet | Only to load D3.js at viewer open (offline fallback: red error banner, rest of page still works minus the visualisation) |

Install minimap2 with:

```bash
conda install -c bioconda minimap2
# or
brew install minimap2
```

## Options

| Flag | Default | Description |
|---|---|---|
| `--fasta1 <file>` | — | FASTA for genome A (enables alignment ribbons) |
| `--fasta2 <file>` | — | FASTA for genome B |
| `--output` / `-o <file>` | `synteny.html` | Output file |
| `--label1 <name>` | `Genome A` | Label for genome A track |
| `--label2 <name>` | `Genome B` | Label for genome B track |
| `--min-block <bp>` | `1000` | Minimum PAF hit length to keep |
| `--chain-gap <bp>` | `100000` | Merge adjacent colinear hits within this gap (0 = off) |
| `--threads` / `-t <n>` | `4` | Threads for minimap2 |
| `--names <file.tsv>` | — | `gene_id<TAB>display_name` mapping |
| `--fetch-names` | — | Auto-fetch gene names from NCBI/UniProt/FungiDB |
| `--taxon <id>` | — | NCBI taxon ID for online name lookup |
| `--keep-paf` | — | Keep the intermediate PAF (normally deleted after rendering) |

## Output: the HTML viewer

The emitted HTML is a single file — open it in any modern browser, no server needed. The embedded JS imports D3.js v7 from a CDN; if the CDN is unreachable a visible error banner explains the issue.

### Layout

- **Top track** — chromosomes of genome A, drawn to scale
- **Middle area** — syntenic ribbons (Bézier curves) between the two genomes
- **Bottom track** — chromosomes of genome B

Chromosomes are sorted by decreasing length on each side, gapped by 4 px.

### Toolbar controls

| Control | Effect |
|---|---|
| **Min identity** slider | Hide ribbons below this identity threshold (default 70 %) |
| **Colour by** selector | `Strand` (blue=+, red=−), `Identity` (YlOrRd gradient), `Chromosome` (Tableau palette keyed by genome-A contig) |
| **Show gene names** checkbox | On-ribbon labels for large blocks (≥ 40 px wide) |
| **Gene overlay** checkbox | Tick marks on the chromosome bars for every annotated gene — named genes highlighted in the accent colour |
| **⬇ SVG** button | Download the current view as a standalone SVG (CSS custom properties resolved to hex at export time so the file renders correctly in any viewer) |
| **⬇ PNG** button | Rasterise the current view to PNG at 2× the on-screen resolution |

### Ribbon tooltip

Hover a ribbon to see:

- Gene name (genome A ↔ genome B, resolved via overlap on the named-gene map)
- Query contig and coordinates
- Target contig and coordinates
- Strand
- Identity (percentage, with chained blocks using a length-weighted mean)
- Length (kb)

## How gene labels work

Labels come from gene-ID → display-name pairs. There are two ways to populate the map:

1. **Local TSV** via `--names file.tsv`:

   ```text
   g001559	CtAlphaBetaHydrolase
   g001560	CtAlphaHydrolase2
   ```

2. **Online fetch** via `--fetch-names [--taxon <id>]` — queries NCBI / UniProt / FungiDB for gene symbols and product names. `--taxon` narrows the search; without it the fetcher uses any species match.

At render time the tool extracts **every gene** from both GFF3s and serialises `(contig, start, end, gene_id)` into the viewer. For each PAF block, the viewer finds the largest-overlap gene on that contig and displays its name (or falls back to the ID if unnamed). This is why labels are correct even though PAF blocks key on contig names rather than gene IDs.

## How block chaining works

minimap2 `asm5` fragments a single colinear region into many short PAF hits whenever small indels, repeats, or heterogeneity break the alignment. Rendering each raw hit as its own ribbon makes the viewer look like alignment confetti.

The chainer groups hits by `(query_contig, target_contig, strand)`, sorts by query start, and merges adjacent neighbours when:

- Query-side gap ≤ `--chain-gap` (default 100 kb)
- Target-side gap ≤ `--chain-gap`
- Target coordinates remain monotonic (increasing on `+` strand, decreasing on `−`)

Merged identity is a **length-weighted mean** (`sum(residue_matches) / sum(block_len)`), so a tiny low-identity fragment joining a long high-identity stretch won't pull the score down disproportionately.

Pass `--chain-gap 0` to disable chaining and inspect the raw PAF hits.

## Troubleshooting

**"minimap2 not found in PATH"** — install it (`conda install -c bioconda minimap2` or `brew install minimap2`) and ensure the binary is on `PATH` when launching the terminal.

**Ribbons look fragmented** — increase `--chain-gap` (e.g. 250000 for divergent genomes with larger indels).

**No ribbons at all** — verify both `--fasta1` and `--fasta2` were passed, that both are valid FASTA, and that you weren't below `--min-block`. Run with `--keep-paf` and inspect the retained `.paf` file to confirm minimap2 actually produced hits.

**Viewer shows red "Could not load D3.js"** — the machine can't reach the CDN. Copy `https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js` to the same directory as the output HTML and edit the `<script src=...>` tag to point to it.

**Too many gene ticks** — leave the **Gene overlay** checkbox off. It's disabled by default for a reason: a gene-rich fungal genome pair pushes ~12,000 tick marks into the DOM.

## Related commands

- [`compare`](../pipeline/overview.md) — N-genome protein-level comparison (BLAST / MMseqs2 / MUMmer). `compare --synteny` now delegates here since a faithful 2-genome ribbon needs whole-genome FASTAs.
- [`phylogeny`](phylogeny.md) — tree-based comparison across many genomes.
- [`plot`](plot.md) — single-genome feature plots.
