# submit -- NCBI GenBank Submission

Validates your annotation and prepares files for NCBI GenBank submission.

## Usage

```bash
myconote-cli submit <annotated.gff3> --fasta <genome.fa> --organism "Genus species" [options]
```

## What it does

1. **Validates** your GFF3 for NCBI compliance (duplicate IDs, orphan features, coordinates)
2. **Generates** an NCBI feature table (`.tbl` format)
3. **Runs** `table2asn` to produce a Sequin file (`.sqn`) if available
4. **Creates** submission metadata files

## Options

| Flag | Description |
|------|-------------|
| `--fasta <file>` | Genome FASTA (required) |
| `--output <dir>` | Output directory (default: `submit_out`) |
| `--organism <name>` | Organism name (required) |
| `--strain <name>` | Strain name |
| `--bioproject <acc>` | BioProject accession (e.g. PRJNA123456) |
| `--biosample <acc>` | BioSample accession (e.g. SAMN12345678) |
| `--locus-prefix <str>` | Locus tag prefix (default: MYCO) |
| `--genetic-code <n>` | Translation table (default: 1) |
| `--email <address>` | Contact email (written into `template.sbt`) |
| `--contact-first <name>` | Submitter first name (`template.sbt`) |
| `--contact-last <name>` | Submitter last name (`template.sbt`) |
| `--institution <name>` | Submitter institution (`template.sbt`) |
| `--fix-structure` | Repair safe structural problems before `table2asn` (synthesise missing gene/mRNA parents for orphan CDS); child coordinates are never altered |
| `--validate-only` | Only validate (including the structure check), do not generate files |

!!! note "Complete the `.sbt` or table2asn won't load it"
    `table2asn` requires a complete submission template. Provide
    `--contact-first`, `--contact-last`, `--institution` and `--email`, or it
    rejects the stub with `Error loading template file`. MycoNote warns when any
    of these is missing.

## Validation checks

| Check | Type |
|-------|------|
| Duplicate feature IDs | Error |
| Orphan Parent references | Error |
| Coordinate ordering (start > end) | Error |
| GFF3 seqids missing from FASTA | Error |
| Missing locus_tag attributes | Warning |
| Missing product descriptions | Warning |

## Examples

```bash
# Validate only (no files generated)
myconote-cli submit genes.gff3 --fasta genome.fa --organism "Aspergillus niger" --validate-only

# Full submission prep
myconote-cli submit annotated.gff3 --fasta genome.fa \
  --organism "Candida albicans" \
  --strain SC5314 \
  --genetic-code 12 \
  --bioproject PRJNA123456 \
  --locus-prefix CALB \
  --email narhmadey@wisc.edu

# Output files
ls submit_out/
#  annotation.tbl      -- NCBI feature table
#  annotation.fsa      -- Genome FASTA copy
#  annotation.sqn      -- Sequin file (if table2asn installed)
#  template.sbt        -- Submission template
#  structure_fixed.gff3 -- repaired GFF3 (only when --fix-structure applied repairs)
```

## Before submitting to NCBI

1. Register a [BioProject](https://submit.ncbi.nlm.nih.gov/) for your genome
2. Register a [BioSample](https://submit.ncbi.nlm.nih.gov/) for your organism/strain
3. Request a [locus_tag prefix](https://www.ncbi.nlm.nih.gov/genbank/genome_locus_tag/) from NCBI
4. Run `myconote-cli submit --validate-only` to check for errors
5. Fix any errors, then run without `--validate-only` to generate files
6. Upload the `.sqn` file to [NCBI Submission Portal](https://submit.ncbi.nlm.nih.gov/)
