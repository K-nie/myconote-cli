# Output Formats

All Myconote_CLI outputs use widely supported bioinformatics formats.

## File formats produced

| Extension | Format | Compatible with |
|-----------|--------|-----------------|
| `.gff3` | GFF3 gene annotation | IGV, JBrowse2, Geneious, UCSC |
| `.faa` | Protein FASTA | BLAST, eggNOG-mapper, InterProScan |
| `.fna` | Nucleotide FASTA | BLAST, Trinity, PASA |
| `.gbk` | GenBank flat file | Geneious, NCBI submission |
| `.vcf` | Variant Call Format | GATK, bcftools |
| `.treefile` | Newick tree | FigTree, iTOL, dendroscope |
| `.tsv` | Tab-separated values | Excel, R, Python |
| `.json` | JSON | Web viewers, custom scripts |
| `.png` / `.svg` | Raster / vector plots | Publications, presentations |

## GFF3 attribute conventions

Myconote_CLI follows the [Sequence Ontology GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md). Gene features use these attributes:

| Attribute | Example | Description |
|-----------|---------|-------------|
| `ID` | `gene_0001` | Unique feature identifier |
| `Name` | `YAL001C` | Gene name (if available) |
| `product` | `ATP synthase subunit alpha` | Functional description |
| `go_terms` | `GO:0005737` | GO term assignments |
| `eggnog` | `COG0055` | eggNOG COG assignment |
| `cazyme` | `GH18` | CAZyme family (if applicable) |
