//! Y1000+ reference bundle integration.
//!
//! Makes the Opulente et al. (2024) dataset available to every myconote
//! subcommand via a single normalised local cache at
//! `~/.myconote/y1000plus/` (override with `MYCONOTE_Y1000_CACHE`).
//!
//! Citation
//! --------
//! Opulente, D. A., LaBella, A. L., Harrison, M.-C. et al. (2024).
//! "Genomic factors shape carbon and nitrogen metabolic niche breadth across
//! Saccharomycotina yeasts." *Science* **384**(6694): eadj4503.
//! <https://doi.org/10.1126/science.adj4503>
//!
//! Figshare collection: <https://plus.figshare.com/collections/_/6714042>.

pub mod benchmark;
pub mod commands;
pub mod download;
pub mod extract;
pub mod install;
pub mod manifest;
pub mod place;
pub mod presets;
pub mod subsets;

pub use manifest::{cache_root, load as load_manifest, save as save_manifest, Manifest};
pub use presets::Preset;
pub use subsets::{format_bytes, ExtractKind, FileSpec, Subset};

/// One-stop citation string surfaced by `--list` and stamped into the
/// manifest so downstream papers can cite the source.
pub const CITATION: &str =
    "Opulente DA et al. (2024). Genomic factors shape carbon and nitrogen metabolic \
     niche breadth across Saccharomycotina yeasts. Science 384(6694): eadj4503. \
     doi:10.1126/science.adj4503";
