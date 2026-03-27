# Databases

Myconote_CLI uses several external databases for functional annotation. Run `myconote db setup` to download and index them.

```bash
myconote db setup --kingdom fungi        # minimal set (~5 GB)
myconote db setup --kingdom fungi --all  # full set including eggNOG (~75 GB)
```

---

## Database summary

| Database | Size | Used by | Notes |
|----------|------|---------|-------|
| UniProtKB/Swiss-Prot | ~250 MB | BLAST | Curated protein annotations |
| dbCAN | ~50 MB | CAZyme annotation | Carbohydrate-active enzymes |
| MEROPS | ~80 MB | Protease annotation | Peptidase families |
| Pfam | ~350 MB | Domain annotation | Protein family HMMs |
| GO ontology (`go.obo`) | ~35 MB | GO term mapping | Gene Ontology |
| eggNOG | ~50 GB | eggNOG-mapper | Orthology + COG + KEGG. Download only if needed |
| RepeatMasker libraries | ~1 GB | Repeat masking | Downloaded by RepeatModeler2 during `mask` |

---

## Storage layout

By default databases are stored in `~/.myconote/db/`. Override with `--db-dir` or the `MYCONOTE_DB` environment variable.

```bash
export MYCONOTE_DB=/data/myconote_databases
myconote db setup --kingdom fungi
```

---

## Checking database status

```bash
myconote db status
```

Reports which databases are present, their versions, and when they were last updated.

!!! warning "eggNOG database"
    The eggNOG database is ~50 GB compressed and requires ~50 GB additional space during indexing. Only download it if you need COG/KEGG assignments beyond what BLAST + dbCAN provide.
