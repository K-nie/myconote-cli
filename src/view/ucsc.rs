/// UCSC Genome Browser integration
///
/// Two outputs:
///
/// 1. **Custom-track file** – prepends the UCSC `track` header to the
///    original GFF3 so it can be drag-dropped onto the UCSC browser or
///    uploaded via "My Data → Custom Tracks".
///
/// 2. **Custom-track URL** – builds a URL of the form
///    `https://genome.ucsc.edu/cgi-bin/hgTracks?db=<db>&hgt.customText=<url>`
///    that opens UCSC at the requested position with the data pre-loaded.
///    (The data must be hosted at a publicly reachable URL for UCSC to
///    fetch it; for private genomes use the custom-track file instead.)

use crate::utils::error::Result;
use super::ViewConfig;
use std::io::{BufRead, BufReader};
use std::fs::File;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Minimal percent-encoder for URL query-string values (RFC 3986 unreserved +
/// a few safe extras).  We avoid pulling in an external crate.
fn url_encode(s: &str) -> String {
    let mut out = String::with_capacity(s.len() * 2);
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9'
            | b'-' | b'_' | b'.' | b'~' | b':' | b'/' => {
                out.push(b as char);
            }
            _ => {
                out.push('%');
                out.push_str(&format!("{:02X}", b));
            }
        }
    }
    out
}

/// Sniff the most-used UCSC assembly DB from feature seqids in the GFF3.
/// Falls back to `"sacCer3"` (yeast) because that is the most common
/// species in fungal genomics work.
fn guess_ucsc_db(gff_path: &std::path::Path) -> String {
    // Check the first few sequence names
    let seqids = sniff_seqids(gff_path, 50);
    for s in &seqids {
        let sl = s.to_lowercase();
        if sl.starts_with("chr") && sl.len() <= 6 {
            // Looks like a chromosomal organism – could be many things.
            // Without more info, default to sacCer3.
            return "sacCer3".to_string();
        }
        if sl.contains("scaffold") || sl.contains("node") || sl.contains("contig") {
            // De-novo assembled fungal genome – no public UCSC db exists;
            // the user will need to load it as a custom assembly.
            return "sacCer3".to_string();
        }
    }
    "sacCer3".to_string()
}

/// Read at most `limit` seqid values from a GFF3.
fn sniff_seqids(gff_path: &std::path::Path, limit: usize) -> Vec<String> {
    let Ok(file) = File::open(gff_path) else { return vec![] };
    let reader = BufReader::new(file);
    let mut out = Vec::new();
    for line in reader.lines().flatten() {
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        if let Some(seqid) = line.split('\t').next() {
            if !out.contains(&seqid.to_string()) {
                out.push(seqid.to_string());
            }
        }
        if out.len() >= limit { break; }
    }
    out
}

/// Choose a sensible UCSC position string:
/// - use the region override when present
/// - otherwise land on the first chromosome with a reasonable window
fn ucsc_position(config: &ViewConfig) -> String {
    if let Some(ref r) = config.region {
        return r.clone();
    }
    // Default: first seqid, first 100 kb
    let seqids = sniff_seqids(&config.gff_path, 1);
    let chr = seqids.into_iter().next().unwrap_or_else(|| "chr1".to_string());
    format!("{}:1-100000", chr)
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Write a UCSC-ready custom-track GFF3 file (adds the `track` header line).
pub fn generate_ucsc_track_file(config: &ViewConfig) -> Result<()> {
    let gff_name = config
        .gff_path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "annotations".to_string());

    let description = config
        .title
        .clone()
        .unwrap_or_else(|| format!("Gene annotations – {}", gff_name));

    let track_header = format!(
        "track name=\"{}\" description=\"{}\" visibility=2 gffTags=on colorByStrand=\"255,128,0 0,0,255\"\n",
        gff_name, description
    );

    // Read original GFF3
    let original = std::fs::read_to_string(&config.gff_path)?;

    // Write new file
    let out_path = config.output.with_extension("ucsc.gff3");
    std::fs::write(&out_path, format!("{}{}", track_header, original))?;

    println!("✓ UCSC custom-track file: {}", out_path.display());
    println!("  Upload at: https://genome.ucsc.edu/cgi-bin/hgCustom");
    println!();

    Ok(())
}

/// Print the UCSC custom-track URL to stdout (data must be publicly hosted).
pub fn print_ucsc_url(config: &ViewConfig) -> Result<()> {
    let db  = guess_ucsc_db(&config.gff_path);
    let pos = ucsc_position(config);

    println!("🔗 UCSC Genome Browser URL");
    println!("──────────────────────────────────────────────────────────────────");
    println!("  Reference DB: {db}  (override with --ucsc-db if needed)");
    println!("  Position:     {pos}");
    println!();

    // URL when the GFF3 is hosted at a public URL
    println!("  ① If your GFF3 is at a public URL, open:");
    println!(
        "    https://genome.ucsc.edu/cgi-bin/hgTracks?db={db}&position={pos_enc}&hgt.customText=<YOUR_GFF3_URL>",
        db      = db,
        pos_enc = url_encode(&pos),
    );
    println!();

    // Track-hub stub for novel assemblies
    println!("  ② For a novel / private fungal assembly, use a UCSC Track Hub:");
    println!("    https://genome.ucsc.edu/goldenPath/help/trackHub.html");
    println!();

    // JBrowse2 online sharing (better for private genomes)
    println!("  ③ Or share interactively with JBrowse2 (no upload needed):");
    println!("    python3 -m http.server 8080");
    println!("    → http://localhost:8080/{}", config.output.display());
    println!("──────────────────────────────────────────────────────────────────");

    Ok(())
}
