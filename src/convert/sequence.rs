/// Sequence format conversions
///
/// FASTA ↔ FASTQ, FASTA + QUAL → FASTQ, FASTQ → FASTA + QUAL,
/// FASTA → TSV table, alignment format conversion
/// (FASTA alignment ↔ PHYLIP ↔ NEXUS ↔ CLUSTAL)
use crate::utils::error::{MycoNoteError, Result};
use std::io::{BufRead, Write};
use std::path::Path;

// ─────────────────────────────────────────────────────────────────────────────
// FASTA → FASTQ
// ─────────────────────────────────────────────────────────────────────────────

/// Convert FASTA to FASTQ by assigning a uniform quality score to every base.
/// Default quality character is 'I' (Phred 40).
pub fn fasta_to_fastq(input: &Path, output: &Path, default_qual: char) -> Result<usize> {
    let file = std::fs::File::open(input).map_err(MycoNoteError::Io)?;
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut count = 0;

    let mut current_id = String::new();
    let mut current_seq = String::new();

    let flush = |out: &mut std::fs::File, id: &str, seq: &str| -> Result<()> {
        if id.is_empty() || seq.is_empty() {
            return Ok(());
        }
        let qual: String = std::iter::repeat(default_qual).take(seq.len()).collect();
        writeln!(out, "@{}\n{}\n+\n{}", id, seq, qual).map_err(MycoNoteError::Io)
    };

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let line = line.trim_end();
        if line.starts_with('>') {
            flush(&mut out, &current_id, &current_seq)?;
            if !current_id.is_empty() {
                count += 1;
            }
            current_id = line[1..].to_string();
            current_seq = String::new();
        } else if !current_id.is_empty() {
            current_seq.push_str(line.trim());
        }
    }
    flush(&mut out, &current_id, &current_seq)?;
    if !current_id.is_empty() {
        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTQ → FASTA
// ─────────────────────────────────────────────────────────────────────────────

/// Strip quality scores from FASTQ, output FASTA.
pub fn fastq_to_fasta(input: &Path, output: &Path) -> Result<usize> {
    let file = std::fs::File::open(input).map_err(MycoNoteError::Io)?;
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut lines = std::io::BufReader::new(file).lines();
    let mut count = 0;

    loop {
        // FASTQ record = 4 lines: @header, seq, +, qual
        let header = match lines.next() {
            Some(Ok(l)) => l,
            Some(Err(e)) => return Err(MycoNoteError::Io(e)),
            None => break,
        };
        let seq = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };
        let _plus = lines.next(); // skip '+'
        let _qual = lines.next(); // skip quality

        if header.starts_with('@') {
            writeln!(out, ">{}", &header[1..]).map_err(MycoNoteError::Io)?;
            // Write sequence in 60-char lines
            for chunk in seq.as_bytes().chunks(60) {
                writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                    .map_err(MycoNoteError::Io)?;
            }
            count += 1;
        }
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTA + QUAL → FASTQ
// ─────────────────────────────────────────────────────────────────────────────

/// Combine a FASTA file and its companion QUAL file into a FASTQ.
/// QUAL file format: >id followed by space-separated Phred score integers.
pub fn fasta_qual_to_fastq(fasta: &Path, qual: &Path, output: &Path) -> Result<usize> {
    // Load quality scores indexed by sequence ID
    let qual_map = load_qual_file(qual)?;

    let file = std::fs::File::open(fasta).map_err(MycoNoteError::Io)?;
    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    let mut count = 0;

    let mut current_id = String::new();
    let mut current_seq = String::new();

    let flush = |out: &mut std::fs::File, id: &str, seq: &str| -> Result<()> {
        if id.is_empty() {
            return Ok(());
        }
        let bare_id = id.split_whitespace().next().unwrap_or(id);
        let quals = qual_map.get(bare_id).cloned().unwrap_or_else(|| {
            // Default to Phred 40 if no QUAL entry
            vec![40u8; seq.len()]
        });
        // Convert Phred integers to ASCII+33 characters
        let qual_str: String = quals.iter().map(|&q| (q.min(93) + 33) as char).collect();
        writeln!(out, "@{}\n{}\n+\n{}", id, seq, qual_str).map_err(MycoNoteError::Io)
    };

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let line = line.trim_end();
        if line.starts_with('>') {
            flush(&mut out, &current_id, &current_seq)?;
            if !current_id.is_empty() {
                count += 1;
            }
            current_id = line[1..].to_string();
            current_seq = String::new();
        } else if !current_id.is_empty() {
            current_seq.push_str(line.trim());
        }
    }
    flush(&mut out, &current_id, &current_seq)?;
    if !current_id.is_empty() {
        count += 1;
    }

    Ok(count)
}

fn load_qual_file(path: &Path) -> Result<std::collections::HashMap<String, Vec<u8>>> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut map = std::collections::HashMap::new();
    let mut current_id = String::new();
    let mut current_quals: Vec<u8> = Vec::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let line = line.trim();
        if line.starts_with('>') {
            if !current_id.is_empty() {
                map.insert(current_id.clone(), current_quals.clone());
            }
            current_id = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_quals = Vec::new();
        } else {
            for token in line.split_whitespace() {
                if let Ok(q) = token.parse::<u8>() {
                    current_quals.push(q);
                }
            }
        }
    }
    if !current_id.is_empty() {
        map.insert(current_id, current_quals);
    }

    Ok(map)
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTQ → FASTA + QUAL
// ─────────────────────────────────────────────────────────────────────────────

/// Split a FASTQ file into a FASTA file and a companion QUAL file.
pub fn fastq_to_fasta_qual(input: &Path, out_fasta: &Path, out_qual: &Path) -> Result<usize> {
    let file = std::fs::File::open(input).map_err(MycoNoteError::Io)?;
    let mut fa = std::fs::File::create(out_fasta).map_err(MycoNoteError::Io)?;
    let mut qf = std::fs::File::create(out_qual).map_err(MycoNoteError::Io)?;
    let mut lines = std::io::BufReader::new(file).lines();
    let mut count = 0;

    loop {
        let header = match lines.next() {
            Some(Ok(l)) if l.starts_with('@') => l[1..].to_string(),
            Some(Ok(_)) => continue,
            Some(Err(e)) => return Err(MycoNoteError::Io(e)),
            None => break,
        };
        let seq = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };
        let _plus = lines.next();
        let qual = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };

        // Write FASTA
        writeln!(fa, ">{}", header).map_err(MycoNoteError::Io)?;
        for chunk in seq.as_bytes().chunks(60) {
            writeln!(fa, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                .map_err(MycoNoteError::Io)?;
        }

        // Write QUAL (Phred scores as space-separated integers)
        let phred: Vec<String> = qual
            .bytes()
            .map(|b| (b.saturating_sub(33)).to_string())
            .collect();
        writeln!(qf, ">{}", header).map_err(MycoNoteError::Io)?;
        // wrap at 60 scores per line
        for chunk in phred.chunks(60) {
            writeln!(qf, "{}", chunk.join(" ")).map_err(MycoNoteError::Io)?;
        }

        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTA → TSV table
// ─────────────────────────────────────────────────────────────────────────────

/// Write a two-column TSV: sequence_id \t sequence.
/// Useful for loading sequences into databases or spreadsheets.
pub fn fasta_to_table(input: &Path, output: &Path) -> Result<usize> {
    use crate::parser::fasta::FastaReader;

    let mut out = std::fs::File::create(output).map_err(MycoNoteError::Io)?;
    writeln!(out, "id\tsequence").map_err(MycoNoteError::Io)?;
    let mut count = 0;

    for rec_res in FastaReader::from_path(input)? {
        let rec = rec_res?;
        writeln!(out, "{}\t{}", rec.id, rec.sequence).map_err(MycoNoteError::Io)?;
        count += 1;
    }

    Ok(count)
}

// ─────────────────────────────────────────────────────────────────────────────
// Alignment format conversion
// ─────────────────────────────────────────────────────────────────────────────

/// In-memory representation of a multiple sequence alignment.
#[derive(Debug, Default)]
pub struct Alignment {
    pub sequences: Vec<(String, String)>, // (id, aligned_sequence)
}

impl Alignment {
    pub fn len(&self) -> usize {
        self.sequences.len()
    }
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// Alignment length (all seqs should be equal length)
    pub fn aln_len(&self) -> usize {
        self.sequences.first().map(|(_, s)| s.len()).unwrap_or(0)
    }
}

/// Convert an alignment file between FASTA, PHYLIP, NEXUS, or CLUSTAL formats.
pub fn convert_alignment(
    input: &Path,
    output: &Path,
    in_format: &str,
    out_format: &str,
) -> Result<usize> {
    let aln = read_alignment(input, in_format)?;
    let n = aln.len();
    write_alignment(output, &aln, out_format)?;
    Ok(n)
}

fn read_alignment(path: &Path, format: &str) -> Result<Alignment> {
    match format.to_lowercase().as_str() {
        "fasta" | "fa" | "aln" => read_fasta_alignment(path),
        "phylip" | "phy" => read_phylip(path),
        "nexus" | "nex" => read_nexus_alignment(path),
        "clustal" | "clustalw" => read_clustal(path),
        other => Err(MycoNoteError::UnsupportedFormat(format!(
            "Unknown alignment input format '{}'. Use: fasta, phylip, nexus, clustal",
            other
        ))),
    }
}

fn write_alignment(path: &Path, aln: &Alignment, format: &str) -> Result<()> {
    match format.to_lowercase().as_str() {
        "fasta" | "fa" => write_fasta_alignment(path, aln),
        "phylip" | "phy" => write_phylip(path, aln),
        "nexus" | "nex" => write_nexus_alignment(path, aln),
        "clustal" | "clustalw" => write_clustal(path, aln),
        other => Err(MycoNoteError::UnsupportedFormat(format!(
            "Unknown alignment output format '{}'. Use: fasta, phylip, nexus, clustal",
            other
        ))),
    }
}

// ── FASTA alignment reader/writer ────────────────────────────────────────────

fn read_fasta_alignment(path: &Path) -> Result<Alignment> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut aln = Alignment::default();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let line = line.trim_end();
        if line.starts_with('>') {
            if !current_id.is_empty() {
                aln.sequences
                    .push((current_id.clone(), current_seq.clone()));
            }
            current_id = line[1..].to_string();
            current_seq = String::new();
        } else if !current_id.is_empty() {
            current_seq.push_str(line.trim());
        }
    }
    if !current_id.is_empty() {
        aln.sequences.push((current_id, current_seq));
    }
    Ok(aln)
}

fn write_fasta_alignment(path: &Path, aln: &Alignment) -> Result<()> {
    let mut out = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    for (id, seq) in &aln.sequences {
        writeln!(out, ">{}", id).map_err(MycoNoteError::Io)?;
        for chunk in seq.as_bytes().chunks(60) {
            writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                .map_err(MycoNoteError::Io)?;
        }
    }
    Ok(())
}

// ── PHYLIP reader/writer ──────────────────────────────────────────────────────

fn read_phylip(path: &Path) -> Result<Alignment> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut lines = std::io::BufReader::new(file).lines();
    let mut aln = Alignment::default();

    // First line: ntaxa nchars
    let header = lines
        .next()
        .ok_or_else(|| MycoNoteError::InvalidFormat("Empty PHYLIP file".into()))??;
    let parts: Vec<&str> = header.split_whitespace().collect();
    if parts.len() < 2 {
        return Err(MycoNoteError::InvalidFormat(
            "PHYLIP header must have ntaxa and nchars".into(),
        ));
    }

    for line in lines {
        let line = line.map_err(MycoNoteError::Io)?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        // PHYLIP interleaved or sequential — try to split at first whitespace
        let (id, seq) = if line.len() > 10 {
            // strict PHYLIP: first 10 chars = name
            let (id_part, seq_part) = line.split_at(10.min(line.len()));
            (id_part.trim().to_string(), seq_part.replace(' ', ""))
        } else {
            let mut parts = line.splitn(2, char::is_whitespace);
            (
                parts.next().unwrap_or("").to_string(),
                parts.next().unwrap_or("").replace(' ', ""),
            )
        };

        if !id.is_empty() && !seq.is_empty() {
            aln.sequences.push((id, seq));
        }
    }

    Ok(aln)
}

fn write_phylip(path: &Path, aln: &Alignment) -> Result<()> {
    let mut out = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(out, " {} {}", aln.len(), aln.aln_len()).map_err(MycoNoteError::Io)?;
    for (id, seq) in &aln.sequences {
        // Pad or truncate ID to 10 chars (strict PHYLIP)
        let padded_id = format!("{:<10}", &id[..id.len().min(10)]);
        writeln!(out, "{}{}", padded_id, seq).map_err(MycoNoteError::Io)?;
    }
    Ok(())
}

// ── NEXUS alignment reader/writer ─────────────────────────────────────────────

fn read_nexus_alignment(path: &Path) -> Result<Alignment> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut aln = Alignment::default();
    let mut in_matrix = false;

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let line = line.trim();

        if line.to_uppercase().starts_with("MATRIX") {
            in_matrix = true;
            continue;
        }
        if line == ";" && in_matrix {
            in_matrix = false;
            continue;
        }

        if in_matrix && !line.is_empty() && !line.starts_with('[') {
            let mut parts = line.splitn(2, char::is_whitespace);
            let id = parts.next().unwrap_or("").to_string();
            let seq = parts.next().unwrap_or("").replace(' ', "");
            if !id.is_empty() && !seq.is_empty() {
                // If the ID already exists (interleaved format), append sequence
                if let Some(entry) = aln.sequences.iter_mut().find(|(i, _)| i == &id) {
                    entry.1.push_str(&seq);
                } else {
                    aln.sequences.push((id, seq));
                }
            }
        }
    }
    Ok(aln)
}

fn write_nexus_alignment(path: &Path, aln: &Alignment) -> Result<()> {
    let mut out = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(out, "#NEXUS\nBegin data;").map_err(MycoNoteError::Io)?;
    writeln!(
        out,
        "  Dimensions ntax={} nchar={};",
        aln.len(),
        aln.aln_len()
    )
    .map_err(MycoNoteError::Io)?;
    writeln!(out, "  Format datatype=DNA missing=? gap=-;").map_err(MycoNoteError::Io)?;
    writeln!(out, "  Matrix").map_err(MycoNoteError::Io)?;
    for (id, seq) in &aln.sequences {
        writeln!(out, "    {} {}", id, seq).map_err(MycoNoteError::Io)?;
    }
    writeln!(out, "  ;\nEnd;").map_err(MycoNoteError::Io)
}

// ── CLUSTAL reader/writer ─────────────────────────────────────────────────────

fn read_clustal(path: &Path) -> Result<Alignment> {
    let file = std::fs::File::open(path).map_err(MycoNoteError::Io)?;
    let mut order: Vec<String> = Vec::new();
    let mut seq_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();

    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(MycoNoteError::Io)?;
        let line = line.trim_end();

        // Skip CLUSTAL header line, blank lines, and conservation lines
        if line.starts_with("CLUSTAL") || line.is_empty() {
            continue;
        }
        if line.starts_with(' ')
            || line.starts_with('*')
            || line.starts_with(':')
            || line.starts_with('.')
        {
            continue;
        }

        let mut parts = line.splitn(2, char::is_whitespace);
        let id = parts.next().unwrap_or("").to_string();
        let seq = parts.next().unwrap_or("").trim().replace(' ', "");

        if !id.is_empty() && !seq.is_empty() {
            if !seq_map.contains_key(&id) {
                order.push(id.clone());
            }
            seq_map.entry(id).or_default().push_str(&seq);
        }
    }

    let sequences = order
        .into_iter()
        .filter_map(|id| seq_map.remove(&id).map(|seq| (id, seq)))
        .collect();

    Ok(Alignment { sequences })
}

fn write_clustal(path: &Path, aln: &Alignment) -> Result<()> {
    let mut out = std::fs::File::create(path).map_err(MycoNoteError::Io)?;
    writeln!(out, "CLUSTAL W (1.83) multiple sequence alignment").map_err(MycoNoteError::Io)?;
    writeln!(out).map_err(MycoNoteError::Io)?;

    let block_size = 60usize;
    let aln_len = aln.aln_len();
    let name_width = aln
        .sequences
        .iter()
        .map(|(id, _)| id.len())
        .max()
        .unwrap_or(10)
        + 2;

    let mut pos = 0;
    while pos < aln_len {
        let end = (pos + block_size).min(aln_len);
        for (id, seq) in &aln.sequences {
            let block = &seq[pos..end.min(seq.len())];
            writeln!(out, "{:<width$}{}", id, block, width = name_width)
                .map_err(MycoNoteError::Io)?;
        }
        writeln!(out).map_err(MycoNoteError::Io)?;
        pos += block_size;
    }
    Ok(())
}
