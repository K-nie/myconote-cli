/// Format conversion utilities for myconote-cli
///
/// Sub-modules:
///   gff3     — GFF3 → GTF / BED6 / BED12 / BEDGraph / TSV / Protein (.faa)
///   sequence — FASTA ↔ FASTQ, FASTA+QUAL, FASTA table, alignment formats
///   vcf      — VCF → BED / TSV / consensus FASTA / ANNOVAR / MAF

pub mod gff3;
pub mod sequence;
pub mod vcf;

pub use gff3::{
    gff3_to_gtf,
    gff3_to_bed,
    gff3_to_bed12,
    gff3_to_bedgraph,
    gff3_to_table,
    gff3_to_protein,
};

pub use sequence::{
    fasta_to_fastq,
    fastq_to_fasta,
    fasta_qual_to_fastq,
    fastq_to_fasta_qual,
    fasta_to_table,
    convert_alignment,
    Alignment,
};

pub use vcf::{
    vcf_to_bed,
    vcf_to_table,
    vcf_to_consensus,
    vcf_to_annovar,
    vcf_to_maf,
};
