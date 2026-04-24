//! sra-toolkit fallback — `prefetch` + `fasterq-dump` wrapper.
//!
//! Invoked only when (a) the user forces `--backend sra`, or
//! (b) ENA returned no URLs for an accession (rare but real for
//! runs still being mirrored, or for geographically restricted data).
//!
//! This path is strictly worse than ENA — sra-toolkit requires
//! `vdb-config`, hits NCBI's SRA store instead of ENA's mirror, is
//! slower on average, and the binaries on macOS are a known install
//! wart. We surface those pains in error messages rather than silently
//! eating them.

use crate::utils::error::{MycoNoteError, Result};
use std::path::{Path, PathBuf};

/// Check that both `prefetch` and `fasterq-dump` resolve on PATH. The
/// caller can call this once before looping over accessions to fail
/// fast with a clear "install sra-toolkit" message instead of
/// discovering missing binaries one accession at a time.
pub fn ensure_sratoolkit() -> Result<()> {
    for bin in ["prefetch", "fasterq-dump"] {
        if which::which(bin).is_err() {
            return Err(MycoNoteError::QuantTool {
                tool: "sra-toolkit".to_string(),
                message: format!(
                    "`{bin}` not found on PATH. Install sra-toolkit (`conda install \
                     -c bioconda sra-tools`), then run `vdb-config --interactive` \
                     once to accept the EULA. On macOS you may also need to set \
                     `NCBI_SETTINGS` to the generated config location."
                ),
            });
        }
    }
    Ok(())
}

/// Output paths produced by `fetch_sra` for one run.
#[derive(Debug)]
pub struct SraRunOutput {
    pub run_accession: String,
    pub fastq_r1: PathBuf,
    pub fastq_r2: Option<PathBuf>,
}

/// Pull one run via `prefetch` + `fasterq-dump`. Output files land
/// directly under `outdir/` with `fasterq-dump`'s default naming:
///   single-end: `{run}.fastq.gz`
///   paired-end: `{run}_1.fastq.gz`, `{run}_2.fastq.gz`
///
/// We rename the paired outputs to the project's `_R1` / `_R2`
/// convention so downstream `samples.tsv` is uniform regardless of
/// which backend was used.
pub fn fetch_sra(run_accession: &str, outdir: &Path) -> Result<SraRunOutput> {
    ensure_sratoolkit()?;
    std::fs::create_dir_all(outdir)?;

    // `prefetch <acc>` downloads the compressed .sra archive into a
    // sub-directory under the current dir. We funnel the work into
    // outdir so later paths are predictable.
    let prefetch_status = duct::cmd!("prefetch", run_accession, "-O", outdir)
        .stderr_to_stdout()
        .unchecked()
        .run()
        .map_err(|e| MycoNoteError::QuantTool {
            tool: "prefetch".to_string(),
            message: format!("{run_accession}: spawn failed: {e}"),
        })?;
    if !prefetch_status.status.success() {
        return Err(MycoNoteError::QuantTool {
            tool: "prefetch".to_string(),
            message: format!(
                "{run_accession}: prefetch exited {}:\n{}",
                prefetch_status.status,
                String::from_utf8_lossy(&prefetch_status.stdout)
            ),
        });
    }

    // `fasterq-dump --split-files` emits separate R1/R2 files for
    // paired-end runs (using `_1` / `_2` suffixes). `--gzip` isn't a
    // standard fasterq-dump flag; we gzip after the fact if needed.
    let fqd_status = duct::cmd!(
        "fasterq-dump",
        run_accession,
        "--split-files",
        "--skip-technical",
        "-O",
        outdir
    )
    .stderr_to_stdout()
    .unchecked()
    .run()
    .map_err(|e| MycoNoteError::QuantTool {
        tool: "fasterq-dump".to_string(),
        message: format!("{run_accession}: spawn failed: {e}"),
    })?;
    if !fqd_status.status.success() {
        return Err(MycoNoteError::QuantTool {
            tool: "fasterq-dump".to_string(),
            message: format!(
                "{run_accession}: fasterq-dump exited {}:\n{}",
                fqd_status.status,
                String::from_utf8_lossy(&fqd_status.stdout)
            ),
        });
    }

    // Normalize output naming. fasterq-dump produces uncompressed
    // .fastq; we compress and rename to the `_R1` / `_R2` convention.
    let paired_1 = outdir.join(format!("{run_accession}_1.fastq"));
    let paired_2 = outdir.join(format!("{run_accession}_2.fastq"));
    let single = outdir.join(format!("{run_accession}.fastq"));

    if paired_1.exists() && paired_2.exists() {
        let r1 = outdir.join(format!("{run_accession}_R1.fastq.gz"));
        let r2 = outdir.join(format!("{run_accession}_R2.fastq.gz"));
        gzip_rename(&paired_1, &r1)?;
        gzip_rename(&paired_2, &r2)?;
        Ok(SraRunOutput {
            run_accession: run_accession.to_string(),
            fastq_r1: r1,
            fastq_r2: Some(r2),
        })
    } else if single.exists() {
        let out = outdir.join(format!("{run_accession}.fastq.gz"));
        gzip_rename(&single, &out)?;
        Ok(SraRunOutput {
            run_accession: run_accession.to_string(),
            fastq_r1: out,
            fastq_r2: None,
        })
    } else {
        Err(MycoNoteError::QuantTool {
            tool: "fasterq-dump".to_string(),
            message: format!(
                "{run_accession}: fasterq-dump exited successfully but produced \
                 no recognized output under {}",
                outdir.display()
            ),
        })
    }
}

/// Gzip-compress `src` into `dest` (streaming), then remove the
/// uncompressed source. Avoids shelling out to `gzip` for portability
/// and uses the already-vendored `flate2`.
fn gzip_rename(src: &Path, dest: &Path) -> Result<()> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::{Read, Write};

    let mut input = std::fs::File::open(src)?;
    let out = std::fs::File::create(dest)?;
    let mut enc = GzEncoder::new(out, Compression::default());

    let mut buf = [0u8; 65_536];
    loop {
        let n = input.read(&mut buf)?;
        if n == 0 {
            break;
        }
        enc.write_all(&buf[..n])?;
    }
    enc.finish()?;
    drop(input);
    std::fs::remove_file(src)?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn ensure_sratoolkit_complains_when_missing() {
        // On CI and most dev boxes neither binary is installed; the
        // function's job is to say so clearly rather than silently
        // fail later inside prefetch.
        if which::which("prefetch").is_err() {
            let err = ensure_sratoolkit().unwrap_err();
            let msg = format!("{err}");
            assert!(msg.contains("prefetch") || msg.contains("fasterq-dump"));
            assert!(
                msg.contains("conda install"),
                "should hint at install: {msg}"
            );
        }
    }

    #[test]
    fn gzip_rename_roundtrips_content() {
        use std::io::Read;
        let tmp = TempDir::new().unwrap();
        let src = tmp.path().join("raw.fastq");
        let dst = tmp.path().join("raw.fastq.gz");
        std::fs::write(&src, b"@r\nACGT\n+\nIIII\n").unwrap();
        gzip_rename(&src, &dst).unwrap();
        assert!(!src.exists(), "source should be removed");
        assert!(dst.exists(), "gzip target missing");

        let f = std::fs::File::open(&dst).unwrap();
        let mut dec = flate2::read::GzDecoder::new(f);
        let mut s = String::new();
        dec.read_to_string(&mut s).unwrap();
        assert_eq!(s, "@r\nACGT\n+\nIIII\n");
    }
}
