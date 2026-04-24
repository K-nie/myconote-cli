//! RNA-seq expression quantification for fungal genomes.
//!
//! Full pipeline: parse a sample sheet → build (or reuse) a decoy-aware
//! salmon index → per sample, fastp QC/trim then `salmon quant` → merge
//! the per-sample `quant.sf` files into a wide count matrix → write a
//! reproducibility manifest. DE analysis stays in R.
//!
//! The pipeline is additive to the existing annotation flow: users run
//! `predict → annotate`, then `convert --to cds` to derive the
//! transcript FASTA salmon indexes, then `quant` here.
//!
//! See `scratch/rnaseq_spec.md` and `scratch/rnaseq_spec_decisions.md`
//! for the design that drives this module.

pub mod bundle;
pub mod merge;
pub mod sample_sheet;

use crate::utils::error::{MycoNoteError, Result};

/// Entry point called from `src/main.rs` when the user runs
/// `myconote-cli quant …`. Current status: stubbed during 0.3.0
/// staging — returns an informative error. Replaced with the real
/// dispatcher as sample_sheet / index / fastp / salmon / merge /
/// bundle modules land.
pub fn run_quant(_args: &[String]) -> Result<()> {
    Err(MycoNoteError::InvalidFormat(
        "`quant` is under construction for 0.3.0. Track progress at \
         scratch/rnaseq_spec.md. For now, use external salmon directly."
            .to_string(),
    ))
}
