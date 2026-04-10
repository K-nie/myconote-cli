pub mod autodetect;
pub mod fasta;
pub mod genbank;
pub mod gff;
pub mod region;

pub use autodetect::{detect_format, FileFormat};
pub use fasta::{read_fasta, read_fasta_index, reverse_complement, FastaReader, FastaRecord};
pub use gff::{GFFReader, GFFRecord};
pub use region::{Region, RegionSelector};
