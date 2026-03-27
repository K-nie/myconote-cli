pub mod gff;
pub mod region;
pub mod genbank;
pub mod fasta;
pub mod autodetect;

pub use gff::{GFFRecord, GFFReader};
pub use region::{Region, RegionSelector};
pub use fasta::{FastaRecord, FastaReader, read_fasta, read_fasta_index, reverse_complement};
pub use autodetect::{FileFormat, detect_format};
