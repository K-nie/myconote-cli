# myconote fetch-rna

Download public RNA-seq FASTQs by SRA/ENA accession. ENA's REST API is the default backend — it requires no credentials, no `vdb-config`, and resolves study / project / sample / experiment accessions to their constituent runs transparently. sra-toolkit (`prefetch` + `fasterq-dump`) is the fallback for runs ENA hasn't mirrored.

Output is a directory of gzipped FASTQs plus a `samples.tsv` pre-populated for `myconote-cli quant`, so the typical usage is a single hop from accession list to DE-ready count matrix.

## Quickstart

```bash
# One run
myconote-cli fetch-rna SRR12345678 -o rna/

# A whole study — ENA expands it into the constituent runs
myconote-cli fetch-rna PRJNA123456 -o rna/

# An accession list (one per line, # comments allowed)
myconote-cli fetch-rna accessions.txt -o rna/ --threads 4

# Immediately quantify against your transcriptome
myconote-cli quant cds.fa --samples rna/samples.tsv --genome genome.fa
```

## Accession kinds

| prefix | kind | notes |
|---|---|---|
| `SRR` / `ERR` / `DRR` | run | direct — one row per accession |
| `SRP` / `ERP` / `DRP` | study | expanded to all runs |
| `PRJNA` / `PRJEB` / `PRJDB` | BioProject | expanded to all runs |
| `SRS` / `ERS` / `DRS` | sample | expanded to all runs |
| `SRX` / `ERX` / `DRX` | experiment | expanded to all runs |

Any of the above can appear in an accession-list file passed on the command line.

## Options

| flag | default | meaning |
|---|---|---|
| `--output <dir>` / `-o` | `rna/` | where FASTQs and `samples.tsv` land |
| `--threads <n>` / `-t` | `1` | parallel downloads across accessions |
| `--retries <n>` | `3` | attempts per FASTQ before giving up |
| `--backend ena\|sra\|auto` | `auto` | backend selection; `auto` tries ENA first, falls back to sra-toolkit per run |
| `--no-verify-md5` | off | skip MD5 checksum verification (not recommended — ENA ships the expected MD5 in the metadata) |
| `--dry-run` | off | resolve accessions and print the URLs that would be fetched, without downloading |

## Outputs

Under `--output <dir>`:

- `{run_accession}.fastq.gz` — single-end
- `{run_accession}_R1.fastq.gz`, `{run_accession}_R2.fastq.gz` — paired-end
- `samples.tsv` — a sample sheet with `sample_id / fastq_r1 / fastq_r2 / condition / strandedness`, ready to feed `myconote-cli quant`. `sample_id` mirrors the run accession; `condition` is left blank for the user to fill in.

## Why ENA first

ENA (the European Nucleotide Archive at EBI) mirrors SRA, serves FASTQ over plain HTTPS, and publishes per-file MD5 checksums. No `vdb-config` dialog, no geographic surprises, no NCBI credentials. For the majority of public fungal RNA-seq submitted in the last five years, ENA has the FASTQs available within hours of submission.

The sra-toolkit fallback exists for the rare run that ENA hasn't mirrored yet — we try it only when ENA returns no URLs or when the user explicitly forces `--backend sra`. On macOS, `sra-toolkit` needs `vdb-config --interactive` run once to accept the EULA; the error message surfaces that when we hit the condition.

## MD5 verification

Every FASTQ is streamed to disk while its MD5 is computed on the fly. On mismatch the partial file is deleted and the download retries up to `--retries` times. The checksum ENA publishes is the same one NCBI attaches to the SRA record, so a verified download matches the original submission byte-for-byte.

`--no-verify-md5` exists for environments where MD5 is unnecessary (read-only NFS mounts, for example), but the cost of verification is modest — MD5 runs at disk speed — so we leave it on by default.

## Typical errors

- **"ENA returned 400 for accession X"** — accession may be a typo, or it may only exist in SRA. Retry with `--backend sra` (requires `sra-toolkit`).
- **"md5 mismatch for ..."** — partial download or corrupted transfer. `fetch-rna` auto-retries up to `--retries`; if it keeps failing, network MTU or corporate proxy may be truncating responses.
- **"prefetch not on PATH"** — ENA didn't have the accession and we fell through to sra-toolkit, which isn't installed. Either install it (`conda install -c bioconda sra-tools`) or confirm the accession is retrievable from ENA via `curl 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=...&result=read_run&fields=fastq_ftp&format=tsv'`.

## Related

- `myconote-cli quant` — feeds on the `samples.tsv` we emit here
- `myconote-cli check` — verifies any optional `sra-tools` install is reachable
