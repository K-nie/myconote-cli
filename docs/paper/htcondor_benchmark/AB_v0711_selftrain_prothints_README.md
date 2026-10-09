# A/B: Augustus self-training and protein-hints lift (v0.7.11)

Reproducibility + methods record for the corrected A/B that tests whether two
MycoNote accuracy levers — Augustus self-training and miniprot protein hints —
actually raise gene-prediction F1 once two provisioning/filter defects are
fixed. Written against commit `5ea962a` (v0.7.11).

## Question

Earlier A/B rounds (v0.7.9–v0.7.10) showed neither lever lifting loose/strict
gene-level F1. Diagnosis found the levers never got a fair test:

1. **Self-training never ran.** Augustus `new_species.pl` aborted with
   `Could not open .../augustus_config//species/generic/generic_parameters.cfg`
   and MycoNote silently fell back to a stock species. Root cause: MycoNote's
   `setup --dbs augustus-fungi` fetched only the curated per-species dirs from
   the Augustus GitHub mirror and **omitted `species/generic/`**, the template
   `new_species.pl` copies to seed a new species.
2. **Protein hints fed whole-kingdom OrthoDB with no identity filter.** The A/B
   fed `~/.myconote/dbs/orthodb/Fungi.fa` (all-fungi, ~3.9 GB) straight to the
   DIAMOND prefilter with only an e-value cutoff. On the one completed replicate
   (cne) loose was flat and strict dropped (0.526 → 0.49).

This A/B asks: with `species/generic/` provisioned and with a percent-identity
gate on the protein prefilter, do the levers lift?

## Literature basis for the fixes

- **Self-training / config:** the community fix for the missing-generic /
  non-writable `AUGUSTUS_CONFIG_PATH` error is to copy the **whole** Augustus
  `config/` tree (generic included) to a writable location (BRAKER issue #23;
  funannotate issue #118, PR #1222 `find_augustus_base`; TAMU HPRC guide).
  funannotate `rsync`s the entire config, so it never hits the missing-generic
  bug. v0.7.11 provisions `species/generic/` the same way, within MycoNote's
  curated-download model.
- **Protein hints / identity gate:** GALBA (Brůna et al., miniprot + Augustus)
  recommends proteins of the **closest relatives**, not broad remote databases,
  and warns miniprot "may produce spurious alignments if the reference proteins
  originate from distantly related species," adding a miniprothint scorer to
  filter them. BRAKER3 (Gabriel et al.) reports higher F1 with a species-excluded
  close-relative DB than with a large remote one. funannotate's default protein
  evidence is UniProtKB/Swiss-Prot (curated/compact), mapped by exonerate at a
  default **80 % identity**. v0.7.11's `--protein-hints-min-identity` is
  MycoNote's equivalent identity gate on the DIAMOND prefilter.

## Inputs

- Genomes (6-genome benchmark subset used here): `ani` (*A. nidulans*),
  `ncr` (*N. crassa*), `cne` (*C. neoformans* JEC21). Masked FASTA, Augustus
  species, genetic code and the seqid rename-table are reused verbatim from each
  genome's baseline rep1 (`results/myconote/<id>/rep1/`) so only the one tested
  flag varies vs baseline.
- Reference annotations: `data/<id>/reference.gff3`.
- Protein DB (protein-hints arms): `~/.myconote/dbs/orthodb/Fungi.fa` (OrthoDB
  fungal partition) — deliberately the same DB as the no-gate baseline arm, so
  the identity gate is the only changed variable.
- `species/generic/` template: 8 `generic_*` files copied into
  `~/.myconote/augustus_config/species/` (stopgap identical to what the fixed
  `setup --dbs augustus-fungi` now produces).

## Commands

Harness: `~/ab_v077/run_ab.sh <genome_id> <feature> <benchmark_dir>` (non-repo),
one Condor job per (feature × genome). Per job:

```
myconote-cli predict <masked.fa> --kingdom fungi --species <SP> \
    --locus-prefix <ID> --output <out>/predict --threads 8 <FLAG>
myconote-cli annotate <out>/predict/consensus.gff3 --fasta <masked.fa> \
    --output <out>/annotate --kingdom fungi --locus-prefix <ID> \
    --genetic-code <code> --threads 8
python3 scripts/compare_annotations.py <out>/annotate/annotated.gff3 \
    <reference.gff3> --rename-table <rename_table.tsv> --output <out>/metrics.json
```

`<FLAG>` per arm:

| arm          | FLAG                                                                 |
|--------------|---------------------------------------------------------------------|
| baseline     | (none — stock species / ab-initio, from each genome's rep1)         |
| selftrain    | `--augustus-training self`                                          |
| prothints    | `--protein-hints Fungi.fa` (no-gate control)                        |
| prothints70  | `--protein-hints Fungi.fa --protein-hints-min-identity 70`          |

Submits: `condor_submit ab_selftrain.sub` (ncr, cne), `ab_prothints70.sub`
(ani, ncr, cne). Baseline protein-hints (no-gate) is clusters 145951 (ani, ncr,
cne). Build: cluster 145978 (v0.7.11, `myco_build` conda env, `cargo install
--locked --offline`), then `~/.cargo/bin` synced to `~/.local/bin` (run scripts
prepend `~/.local/bin`; a stale copy there has shadowed rebuilds before).

## Scoring

`compare_annotations.py` — Eilbeck et al. (2009) gene/exon/nucleotide F1.
Strict gene match = identical ordered CDS structure (CDS-to-CDS); loose = locus
overlap. Metric reported: `gene_level.{loose,strict}.F1`.

## Baselines (gene-level loose / strict F1)

- self-train: ncr 0.897 / 0.564, cne 0.921 / 0.526
- protein-hints: ani 0.897 / 0.445, ncr 0.897 / 0.564, cne 0.921 / 0.526

## Tool / software versions

- MycoNote-CLI v0.7.11 (`5ea962a`); Augustus 3.3.3 (not AVX2 — runs full node
  pool); SNAP; GeneMark-ES/EP (`gmes_linux_64_4`); DIAMOND (prefilter); miniprot
  (protein-to-genome); Python 3 (scorer). Built with the `myco_build` conda env.

## Outputs

- Per arm: `results/ab/<feature>/<genome>/{predict,annotate,metrics.json,*.log}`.
- Engagement evidence lives in `predict.log`: self-train prints `Trained
  species: …` on success or the loud `SELF-TRAINING DID NOT RUN` banner on
  fallback (now naming which gate failed); the gated prefilter logs
  `evalue ≤ …, identity ≥ 70%` and the kept-protein count.

## Results

PENDING — orchestration watcher `bivhy3h0d` reports loose/strict vs baseline
for all arms on completion. Fill this section with the table and the
progress-vs-lag decision once the jobs drain.

## Key decisions / caveats

- The protein-hints gate is tested on the **same** all-fungi OrthoDB as the
  no-gate control, so the comparison isolates the identity filter. It does NOT
  test the orthogonal lever the literature also favours — swapping the DB for a
  Swiss-Prot / close-relative set. If the 70 % gate alone does not lift, a DB
  swap is the next test, not another gate tweak.
- `--protein-hints-min-identity` default is 0 (off): v0.7.11 does not silently
  change existing behaviour; the gate is opt-in and the A/B picks the operating
  point.
- Self-training provisioning used a `cp` stopgap of `species/generic/` from the
  `myconote_augustus` conda env; functionally identical to the fixed `setup`
  path, chosen only to avoid a re-download before the A/B.
