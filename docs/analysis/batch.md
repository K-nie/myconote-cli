# myconote batch

Annotate multiple genomes in one command. Runs the full pipeline (sort, mask, predict, annotate, submit) per genome with progress tracking, resume support, and HTCondor integration.

```bash
myconote-cli batch <genomes_dir | sample_sheet.tsv> [options]
```

## Quick examples

```bash
# Annotate all FASTAs in a directory
myconote-cli batch genomes/ --kingdom fungi --threads 8

# Use a sample sheet for per-genome settings
myconote-cli batch samples.tsv --parallel 4

# Run only specific stages
myconote-cli batch genomes/ --stages sort,mask,predict

# Resume after interruption
myconote-cli batch --resume batch_out/

# Generate HTCondor submit files
myconote-cli batch genomes/ --condor --condor-mem 64G --condor-cpus 16
```

## Input modes

### Directory of FASTAs

Point batch at a directory containing `.fa`, `.fasta`, `.fas`, or `.fna` files:

```bash
myconote-cli batch genomes/
```

All genomes use the same settings (kingdom, genetic code, etc.) from CLI flags.

### Sample sheet (TSV)

For per-genome settings, provide a tab-separated file with a header row:

```
name        fasta                kingdom    species          genetic_code    locus_prefix
isolate_A   /data/isolate_A.fa   fungi      saccharomyces    1               ISOA
isolate_B   /data/isolate_B.fa   fungi      auto             12              ISOB
isolate_C   /data/isolate_C.fa   plant      arabidopsis      1               ISOC
```

Only the `fasta` column is required (also accepts `path` or `genome` as column names). All other columns are optional and use defaults when omitted.

| Column | Required | Default | Description |
|--------|----------|---------|-------------|
| `name` | No | FASTA filename stem | Display name for this genome |
| `fasta` | Yes | -- | Path to input FASTA |
| `kingdom` | No | `fungi` | Kingdom for prediction |
| `species` | No | auto | Augustus species override |
| `genetic_code` | No | `1` | NCBI translation table |
| `locus_prefix` | No | `GENE` | Locus tag prefix |

## Options

### Pipeline options

| Flag | Default | Description |
|------|---------|-------------|
| `--output <dir>` | `batch_out` | Batch output directory |
| `--stages <list>` | `sort,mask,predict,annotate,submit` | Comma-separated stages to run |
| `--kingdom <k>` | `fungi` | Default kingdom |
| `--threads <n>` | `4` | Threads per genome |
| `--parallel <n>` | `2` | Max genomes in parallel |
| `--min-length <bp>` | `500` | Minimum contig length for sort |
| `--mask-engine <e>` | `repeatmodeler` | Masking engine |
| `--genetic-code <n>` | `1` | Default translation table |
| `--locus-prefix <s>` | `GENE` | Default locus prefix |
| `--no-compare` | -- | Skip auto-compare suggestion |
| `--resume <dir>` | -- | Resume a previous batch run |

### HTCondor options

| Flag | Default | Description |
|------|---------|-------------|
| `--condor` | -- | Generate HTCondor submit files (don't run locally) |
| `--condor-cpus <n>` | `8` | CPUs per job |
| `--condor-mem <size>` | `32G` | Memory per job |
| `--condor-disk <size>` | `50G` | Disk per job |
| `--condor-queue <name>` | -- | HTCondor accounting group |
| `--condor-extra <file>` | -- | Extra submit directives to append |

## Output structure

```
batch_out/
  status.json                    # batch state (for resume)
  isolate_A/
    isolate_A_sorted.fa
    isolate_A_masked.fa
    predict_out/
      consensus.gff3
    annotate_out/
      annotated.gff3
      proteins.fa
      annotations.tsv
    submit_out/
      annotation.tbl
  isolate_B/
    ...
```

## Dashboard

The dashboard auto-detects your environment:

### Interactive terminal (laptop, SSH session)

Live multi-progress bars showing per-genome stage progress:

```
  myconote batch -- 3 / 8 genomes complete

  v isolate_A    sort > mask > predict > annotate > submit   42m
  v isolate_B    sort > mask > predict > annotate > submit   38m
  * isolate_C    sort > mask > predict > update              27m
  * isolate_D    sort > mask > train                         14m
  . isolate_E    queued
  . isolate_F    queued

  [============---------------------------------------] 3/8 genomes (ETA ~2h15m)
```

### HPC batch job (no TTY)

Timestamped log lines, one per event:

```
[2026-04-16 14:30:22] isolate_A                      | sort       | started
[2026-04-16 14:31:05] isolate_A                      | sort       | done
[2026-04-16 14:31:05] isolate_A                      | mask       | started
...
[2026-04-16 16:45:33] isolate_E                      | predict    | FAILED -- Augustus exit code 1
...
[2026-04-16 18:12:01] isolate_A                      | ALL        | done (3720s)
```

Monitor with `tail -f slurm-12345.out` or `condor_tail -f job_0`.

## Resume

A `status.json` file is written after every stage completion. If the process is interrupted:

```bash
myconote-cli batch --resume batch_out/
```

This reads the state file, skips completed stages, and picks up where it left off. Per-genome, per-stage granularity -- if genome A completed sort+mask but failed on predict, resume reruns predict for genome A and continues with the remaining genomes.

## HTCondor integration

Generate submit files instead of running locally:

```bash
myconote-cli batch genomes/ --condor --condor-mem 64G --condor-cpus 16
```

This creates:

| File | Description |
|------|-------------|
| `condor.sub` | HTCondor submit file |
| `run_genome.sh` | Wrapper script executed by each job |
| `condor_genomes.txt` | Per-genome argument list (tab-separated) |
| `condor_logs/` | Directory for per-job log/output/error files |

Submit to the cluster:

```bash
condor_submit batch_out/condor.sub
```

Monitor:

```bash
condor_q                                      # job status
tail -f batch_out/condor_logs/job_0.out       # first job output
condor_tail -f <cluster_id>.0                 # live output
```

Each genome runs as a separate job. The wrapper script chains the pipeline stages sequentially within each job. Jobs use the shared filesystem -- no file transfer needed.

### Custom HTCondor directives

Append extra submit directives from a file:

```bash
myconote-cli batch genomes/ --condor --condor-extra extra.sub
```

Example `extra.sub`:

```
+WantGPU = true
requirements = (OpSysAndVer == "Ubuntu22")
notification = Complete
notify_user = you@wisc.edu
```

## Auto-compare

When 2 or more genomes complete successfully, batch prints a suggested compare command:

```
  Next: myconote-cli compare batch_out/**/annotate_out/*.gff3
```

Disable with `--no-compare`.

## Resource guidelines

| Organism type | RAM per genome | Disk per genome | Time per genome |
|---------------|---------------|-----------------|-----------------|
| Fungi (~30-50 MB) | 8-16 GB | 2-5 GB | 1-3 hours |
| Plant (~500 MB-1 GB) | 32-64 GB | 20-50 GB | 12-48 hours |
| Insect (~200-500 MB) | 16-32 GB | 10-20 GB | 4-12 hours |
| Small animal (~1-3 GB) | 64+ GB | 50-100 GB | 24-72 hours |

For HTCondor, set `--condor-mem` and `--condor-cpus` accordingly.
