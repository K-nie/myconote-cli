# CAZyme annotation fix — dbCAN HMMER path (benchmark myconote arm)

**Status: RESOLVED (v0.7.15).** CAZyme annotation in the MycoNote benchmark arm
now runs an `hmmsearch` of the predicted proteins against the dbCAN HMM database,
the same method the funannotate arm uses, so the two pipelines are directly
comparable. Verified on *S. cerevisiae* (see Verification below).

## Why this exists (the question it answers)

The functional-annotation comparison claims MycoNote and funannotate both
annotate CAZymes (carbohydrate-active enzymes). A pre-submission check found the
MycoNote arm was emitting **no** CAZyme calls — `annotate --cazyme` was silently
doing nothing — which would have made the CAZyme comparison an unfair asymmetry
(funannotate reports CAZymes; MycoNote reported none).

## Root causes (all verified with real commands on the cluster)

1. **Corrupt dbCAN database.** `~/.myconote/dbs/dbcan/CAZyDB.fa` and
   `dbCAN-HMMdb-V12.txt` were both 18,618-byte HTML error pages (`head` →
   `<!DOCTYPE html`; `grep -c "^>"` = 0) saved by a failed 2025-07-09 download.
2. **No reachable CAZyme backend.** `run_myconote.sh` passed `--cazyme` with no
   database flag. MycoNote's `cazyme::dbcan_available()` probes `which
   run_dbcan.py`; the installed dbCAN entrypoint is `run_dbcan` (v3.0.7, no
   `.py`), so the probe failed and annotate hit the "No dbCAN backend" branch.
3. **dbCAN2 → dbCAN3 incompatibility.** Even exposed, dbCAN 3.0.7's `overview.txt`
   columns (`Gene ID | EC# | HMMER | dbCAN_sub | DIAMOND`) differ from the
   dbCAN2 layout (`Gene ID | HMMER | Hotpep | DIAMOND`) MycoNote's parser read,
   so driving run_dbcan 3.0.7 unchanged would have misparsed families.

An interim DIAMOND fallback (blastp vs the full CAZyDB sequence set) was tried
and rejected: it **over-called CAZymes ~16×** — 2,174 of 5,773 *S. cerevisiae*
genes (37.7 %) — because sequence similarity to a full-length multidomain CAZy
protein is not family-domain membership. DIAMOND-blastp is not a valid CAZyme
comparator here.

## The fix (v0.7.15)

A native HMMER path in `src/annotate/cazyme.rs`:

- `run_hmmer_cazyme()` runs `hmmsearch --domtblout` of the proteins against the
  dbCAN HMM database, reusing the same engine/approach as the Pfam step.
- `parse_dbcan_domtblout()` applies the canonical dbCAN cutoffs — HMM-profile
  coverage > 0.35, and domain i-E-value < 1e-15 (tightened to < 1e-18 when the
  aligned HMM span ≤ 80 aa) — and keeps the best (lowest-E) hit per
  (gene, family). Family = HMM name with a trailing `.hmm` stripped (`GH5.hmm`
  → `GH5`).
- New flag `--cazyme-hmm <dbCAN.hmm>`. Backend precedence:
  `--cazyme-hmm` (hmmsearch) → `run_dbcan` → `--cazyme-db` (DIAMOND fallback).

This is funannotate's method: both arms `hmmsearch` the proteins against the
dbCAN HMM profiles.

## Inputs / provenance (versions pinned)

- **MycoNote** v0.7.15 (source rsync'd to the cluster; `myconote-cli --version`
  → v0.7.15; build job cluster 146776, `cargo install --locked --offline`,
  synced `~/.cargo/bin` → `~/.local/bin`, md5 `f547aad0822eb8259dbc29f7fa5aa103`).
- **dbCAN HMM database (primary):** `~/.myconote/dbs/funannotate/dbCAN.hmm`
  (100,916,453 bytes) + pressed indices `.h3{f,i,m,p}`; dbCAN HMMdb **14.0**,
  CAZyDB 2025-07-10 (provisioned earlier by `funannotate setup -d`). This is the
  identical HMM database funannotate's arm uses.
- **hmmsearch / hmmpress:** HMMER from the `myconote_dbcan` conda env (prepended
  to PATH in the run scripts).
- **DIAMOND fallback DB (retained, not used when `--cazyme-hmm` resolves):**
  `~/.myconote/dbs/dbcan/CAZyDB.dmnd` (dbCAN CAZyDB 2025-07-24, dbCAN v5-2 bundle;
  built from `CAZyDB.07242025.fa`, src md5 `47f4a68c3b8a7f9d4d1be897bed43d14`,
  dmnd md5 `7efb6928ffc20ac3f7c4875624ff0966`; build job via
  `jobs/build_cazyme_db.sh`). The corrupt originals are preserved at
  `~/.myconote/dbs/dbcan/corrupt_html_backup/`; provenance stamped in
  `~/.myconote/dbs/dbcan/dbcan.version`.

## Reproduce (exact commands)

Benchmark arm (`jobs/run_myconote.sh`) now passes:

```
myconote-cli annotate <predict>/consensus.gff3 --fasta <masked.fa> \
  --output <out>/annotate_out --kingdom fungi --locus-prefix <PFX> \
  --genetic-code 1 --db-dir ~/.myconote/dbs \
  --eggnog --eggnog-db ~/.myconote/dbs/eggnog \
  --cazyme --cazyme-hmm ~/.myconote/dbs/funannotate/dbCAN.hmm \
  --cazyme-db ~/.myconote/dbs/dbcan/CAZyDB.dmnd \
  --merops --threads 8
```

Fast CAZyme-only re-verification: `jobs/verify_cazyme_hmm.sh` (annotate on an
existing `sce` predict set), submitted via `verify_cazyme_hmm.sub` (compute node
only, never the submit node).

## Verification (sce, cluster 146777, v0.7.15)

Re-ran annotate on the existing *S. cerevisiae* predict output (5,773 predicted
genes):

| metric | DIAMOND fallback (rejected) | **HMMER (shipped)** |
|---|---|---|
| CAZyme annotations (rows) | 2,174 | **333** |
| distinct genes called | 2,174 (37.7 %) | **136 (2.36 %)** |
| rows with junk family label | 992 (46 %) | **0** |
| distinct families | 830 (inflated) | 104 |
| cazyme step wall time | ~44 s | ~41 s |

Class breakdown (HMMER): GH 237, GT 67, AA 20, CBM 6, CE 3 — GH-dominated, as
expected for a fermentative yeast with a reduced CAZyme repertoire. 136 CAZyme
genes is consistent with published dbCAN/CAZy counts for *S. cerevisiae*.

## Outputs

- `annotate_out/cazyme_hits.tsv` (`gene_id`, `cazyme_family`, `evalue`, `method`;
  method = `HMMER`).
- `annotate_out/cazyme/dbcan_hmmer.domtblout` (raw hmmsearch table).
- Per-source timing in `annotate_out/annotate_timings.tsv` (`cazyme` row).

## Decisions / caveats

- **HMMER over DIAMOND** because DIAMOND-blastp over-calls ~16× and is not the
  dbCAN-recommended method. The DIAMOND path remains only as a fallback for
  users without an HMM database and is documented as approximate.
- **dbCAN HMMdb version (14.0, CAZyDB 2025-07-10)** matches the funannotate arm
  exactly by design — both consume the same `dbCAN.hmm`.
- **antiSMASH / InterProScan remain excluded** from both arms (require network on
  execute nodes); unchanged by this fix.
- Verified on `sce` only so far; the per-genome CAZyme counts for all panel
  genomes land with the production run.
