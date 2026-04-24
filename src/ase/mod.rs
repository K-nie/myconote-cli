//! Allele-specific expression (ASE) for heterozygous / hybrid /
//! polyploid fungal genomes.
//!
//! Given a phased VCF + a reference CDS FASTA + a sample sheet,
//! build personalized transcriptomes per haplotype and quantify
//! expression separately against each. Produces allele-count
//! matrices suitable for downstream cis/trans regression or
//! binomial ASE tests in R.
//!
//! Implementation state: foundation landing incrementally. See
//! `scratch/ase_spec.md` for the full design; all 10 open decisions
//! were locked on 2026-04-24 before code started.

pub mod vcf;

use crate::utils::error::{MycoNoteError, Result};

/// Placeholder dispatcher — implementation arrives module by module
/// (see `scratch/ase_spec.md` release plan).
pub fn run_ase(_args: &[String]) -> Result<()> {
    Err(MycoNoteError::QuantSheet(
        "`ase` is under construction for 0.5.0. See scratch/ase_spec.md \
         for the design and current progress."
            .to_string(),
    ))
}
