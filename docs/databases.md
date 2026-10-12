# Databases

MycoNote-CLI uses several external databases for functional annotation. Use
`myconote-cli setup` to download and index them. Databases live in
`~/.myconote/dbs` by default.

```bash
myconote-cli setup --list                 # show all databases + download status
myconote-cli setup --dbs swiss-prot pfam  # download a specific set
myconote-cli setup --db augustus-fungi    # curated fungal Augustus species (~50) + the
                                          # species/generic/ self-training template
myconote-cli setup --dbs eggnog           # ~50 GB — only if you need COG/KEGG
```

---

## Database summary

| Database | Size | Used by | Notes |
|----------|------|---------|-------|
| UniProtKB/Swiss-Prot | ~250 MB | MMseqs2 homology | Curated protein annotations |
| dbCAN | ~50 MB | CAZyme annotation | Carbohydrate-active enzymes |
| MEROPS | ~80 MB | Protease annotation | Peptidase families |
| Pfam | ~350 MB | Domain annotation (hmmsearch) | Protein family HMMs |
| GO ontology (`go.obo`) | ~35 MB | GO term mapping | Gene Ontology |
| eggNOG | ~50 GB | eggNOG-mapper | Orthology + COG + KEGG. Download only if needed |
| RepeatMasker libraries | ~1 GB | Repeat masking | Downloaded by RepeatModeler2 during `mask` |
| augustus-fungi | ~small | `predict` / `train` | ~50 curated fungal species + the `species/generic/` template needed for `--augustus-training self` |
| Ollama + LLM model | 4-48 GB | `explain` (AI interpreter) | Local LLM, auto-selected by RAM from pulled models. Setup: `myconote-cli setup --dbs ollama` |
| Chat paper corpus | ~50 MB | `explain` (citations) | Q1 open-access papers (CC-BY). Setup: `myconote-cli setup --dbs chat-corpus` |

---

## Storage layout

By default databases are stored in `~/.myconote/dbs`. Override per-run with
`--db-dir <dir>`.

```bash
myconote-cli setup --dbs swiss-prot pfam --db-dir /data/myconote_databases
```

---

## Checking database status

```bash
myconote-cli setup --list    # what is available and whether it's downloaded
myconote-cli setup --check   # verify existing databases are intact
```

Re-download anything with `--force`, or preview a download with `--dry-run`.

!!! warning "eggNOG database"
    The eggNOG database is ~50 GB compressed and needs ~50 GB more during
    indexing. Download it only if you need COG/KEGG assignments beyond what the
    Swiss-Prot (MMseqs2) + dbCAN annotations provide.
