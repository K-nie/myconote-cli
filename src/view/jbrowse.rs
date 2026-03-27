/// JBrowse2 visualisation output
///
/// Generates a self-contained HTML file that loads JBrowse2 from the
/// official CDN and embeds the GFF3 annotation features as inline JSON
/// (using JBrowse2's `FromConfigAdapter`).  No local server or separate
/// index files are required for most fungal-scale annotation files.
///
/// For assemblies without a FASTA file, chromosome sizes are derived
/// from the maximum coordinate seen in the GFF3.  When a FASTA path is
/// supplied the config references it via an `IndexedFastaAdapter` (the
/// `.fai` index must already exist alongside the FASTA).

use crate::parser::gff::GFFReader;
use crate::utils::error::Result;
use super::ViewConfig;
use serde_json::{json, Value};
use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// Internal data types
// ─────────────────────────────────────────────────────────────────────────────

struct ParsedGFF {
    /// Top-level features (genes, etc.) with nested subfeatures
    features: Vec<Value>,
    /// Chromosome → max observed coordinate (used when no FASTA is given)
    chrom_sizes: HashMap<String, u64>,
    /// Total number of features parsed (for size-warning logic)
    total_count: usize,
}

// ─────────────────────────────────────────────────────────────────────────────
// GFF3 parser → JBrowse2 feature tree
// ─────────────────────────────────────────────────────────────────────────────

fn percent_encode_chrom_sizes(sizes: &HashMap<String, u64>) -> String {
    let mut sorted: Vec<(&String, &u64)> = sizes.iter().collect();
    sorted.sort_by_key(|(n, _)| n.as_str());
    let raw: String = sorted
        .iter()
        .map(|(n, s)| format!("{}\t{}\n", n, s))
        .collect();
    // Encode only the characters that break data: URIs
    raw.chars()
        .map(|c| match c {
            '\t' => "%09".to_string(),
            '\n' => "%0A".to_string(),
            ' '  => "%20".to_string(),
            '#'  => "%23".to_string(),
            _    => c.to_string(),
        })
        .collect()
}

fn strand_to_int(c: char) -> i8 {
    match c {
        '+' => 1,
        '-' => -1,
        _   => 0,
    }
}

/// Parse the GFF3 file and return feature hierarchy + chromosome sizes.
/// `names` is an optional ID → display-name map loaded from `--names` or the API.
fn parse_gff(
    gff_path: &std::path::Path,
    names: &std::collections::HashMap<String, String>,
) -> Result<ParsedGFF> {
    // flat map: id → raw JBrowse2 feature value (subfeatures empty initially)
    let mut flat: HashMap<String, Value> = HashMap::new();
    let mut children_of: HashMap<String, Vec<String>> = HashMap::new();
    let mut top_level_ids: Vec<String> = Vec::new();
    let mut chrom_sizes: HashMap<String, u64> = HashMap::new();
    let mut id_counter: usize = 0;

    let reader = GFFReader::from_path(gff_path)?;

    for result in reader {
        let record = match result {
            Ok(r)  => r,
            Err(_) => continue,   // skip malformed lines
        };

        // Track chromosome max position
        let max = chrom_sizes.entry(record.seqid.clone()).or_insert(0);
        if record.end > *max { *max = record.end; }

        // Determine a stable unique ID
        let id = record.id()
            .cloned()
            .unwrap_or_else(|| {
                id_counter += 1;
                format!("feat_{}_{}", record.feature_type, id_counter)
            });

        // Human-readable display name — resolution priority:
        //   1. names map (from --names file or API fetch)
        //   2. GFF3 Name / gene / gene_name attribute
        //   3. Feature ID (fallback)
        let name: Option<String> = names.get(&id)
            .cloned()
            .or_else(|| record.attributes.get("Name").cloned())
            .or_else(|| record.attributes.get("gene").cloned())
            .or_else(|| record.attributes.get("gene_name").cloned())
            .or_else(|| record.id().cloned());

        let feature = json!({
            "uniqueId":  id,
            "refName":   record.seqid,
            "start":     record.start.saturating_sub(1),  // GFF3 → 0-based
            "end":       record.end,
            "type":      record.feature_type,
            "name":      name,
            "strand":    strand_to_int(record.strand),
            "subfeatures": []
        });

        let parent = record.parent().cloned();
        flat.insert(id.clone(), feature);

        match parent {
            Some(p) => children_of.entry(p).or_default().push(id),
            None    => top_level_ids.push(id),
        }
    }

    let total_count = flat.len();

    // ── Build the hierarchy ──────────────────────────────────────────────────
    fn attach_children(
        feat: &mut Value,
        own_id: &str,
        children_of: &HashMap<String, Vec<String>>,
        flat: &mut HashMap<String, Value>,
    ) {
        if let Some(child_ids) = children_of.get(own_id) {
            let mut subs: Vec<Value> = Vec::new();
            for cid in child_ids {
                if let Some(mut child) = flat.remove(cid) {
                    attach_children(&mut child, cid, children_of, flat);
                    subs.push(child);
                }
            }
            if let Some(obj) = feat.as_object_mut() {
                obj.insert("subfeatures".to_string(), Value::Array(subs));
            }
        }
    }

    let mut features: Vec<Value> = Vec::new();
    for tid in &top_level_ids {
        if let Some(mut feat) = flat.remove(tid) {
            attach_children(&mut feat, tid, &children_of, &mut flat);
            features.push(feat);
        }
    }
    // Any remaining orphaned features (parent not found)
    for (_, feat) in flat {
        features.push(feat);
    }

    Ok(ParsedGFF { features, chrom_sizes, total_count })
}

// ─────────────────────────────────────────────────────────────────────────────
// JBrowse2 config builders
// ─────────────────────────────────────────────────────────────────────────────

fn build_assembly(
    name: &str,
    chrom_sizes: &HashMap<String, u64>,
    fasta_path: Option<&std::path::Path>,
) -> Value {
    let track_id = format!("{}-ReferenceSequenceTrack", name);

    let adapter = if let Some(fa) = fasta_path {
        let fa_str  = fa.to_string_lossy();
        let fai_str = format!("{}.fai", fa_str);
        json!({
            "type": "IndexedFastaAdapter",
            "fastaLocation": { "uri": fa_str },
            "faiLocation":   { "uri": fai_str }
        })
    } else {
        // Encode chrom sizes into a data: URI so the HTML is fully self-contained
        let encoded = percent_encode_chrom_sizes(chrom_sizes);
        json!({
            "type": "ChromSizesAdapter",
            "chromSizesLocation": {
                "uri": format!("data:text/plain,{}", encoded)
            }
        })
    };

    json!({
        "name": name,
        "sequence": {
            "type":    "ReferenceSequenceTrack",
            "trackId": track_id,
            "adapter": adapter
        }
    })
}

fn build_annotation_track(assembly_name: &str, features: Vec<Value>) -> Value {
    let track_id   = format!("{}-annotations", assembly_name);
    let display_id = format!("{}-LinearBasicDisplay", track_id);
    json!({
        "type":          "FeatureTrack",
        "trackId":       track_id,
        "name":          "Gene Annotations",
        "assemblyNames": [assembly_name],
        "adapter": {
            "type":     "FromConfigAdapter",
            "features": features
        },
        "displays": [{
            "type":      "LinearBasicDisplay",
            "displayId": display_id
        }]
    })
}

fn build_default_session(
    assembly_name: &str,
    chrom_sizes: &HashMap<String, u64>,
    region_override: Option<&str>,
) -> Value {
    // Choose the first (largest) chromosome for the default view
    let mut sorted: Vec<(&String, &u64)> = chrom_sizes.iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(a.1));

    let (default_ref, default_end) = sorted
        .first()
        .map(|(n, s)| (n.as_str(), **s))
        .unwrap_or(("chr1", 1_000_000));

    let (ref_name, start, end) = if let Some(r) = region_override {
        parse_region_str(r).unwrap_or((default_ref.to_string(), 0, default_end))
    } else {
        (default_ref.to_string(), 0, default_end)
    };

    let track_id   = format!("{}-annotations", assembly_name);
    let display_id = format!("{}-LinearBasicDisplay", track_id);

    json!({
        "name": "MycoNote session",
        "view": {
            "id":   "LinearGenomeView-1",
            "type": "LinearGenomeView",
            "displayedRegions": [{
                "refName":      ref_name,
                "start":        start,
                "end":          end,
                "assemblyName": assembly_name
            }],
            "tracks": [{
                "id":            format!("{}-track", track_id),
                "type":          "FeatureTrack",
                "configuration": track_id,
                "displays": [{
                    "id":     display_id,
                    "type":   "LinearBasicDisplay",
                    "height": 300
                }]
            }]
        }
    })
}

fn parse_region_str(r: &str) -> Option<(String, u64, u64)> {
    // Format: "seqid:start-end"  (1-based, like samtools)
    let (seq, coords) = r.split_once(':')?;
    let (s, e) = coords.split_once('-')?;
    let start: u64 = s.replace(',', "").parse().ok()?;
    let end:   u64 = e.replace(',', "").parse().ok()?;
    Some((seq.to_string(), start.saturating_sub(1), end))
}

// ─────────────────────────────────────────────────────────────────────────────
// HTML template
// ─────────────────────────────────────────────────────────────────────────────

fn render_html(title: &str, config_json: &str) -> String {
    // Escape </script> inside the JSON blob to avoid breaking the HTML parser
    let safe_json = config_json.replace("</script>", "<\\/script>");
    format!(
r#"<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{title}</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      background: #0d1117;
      color: #e6edf3;
      height: 100vh;
      display: flex;
      flex-direction: column;
    }}
    header {{
      background: linear-gradient(135deg, #161b22 0%, #0d1117 100%);
      border-bottom: 1px solid #30363d;
      padding: 10px 20px;
      display: flex;
      align-items: center;
      gap: 14px;
      flex-shrink: 0;
    }}
    .logo {{ font-size: 24px; }}
    .brand {{ font-size: 17px; font-weight: 700; color: #58a6ff; letter-spacing: -0.3px; }}
    .subtitle {{ font-size: 12px; color: #8b949e; margin-top: 1px; }}
    .badge {{
      margin-left: auto;
      background: #21262d;
      border: 1px solid #30363d;
      border-radius: 6px;
      padding: 3px 10px;
      font-size: 11px;
      color: #8b949e;
    }}
    #jbrowse-container {{
      flex: 1;
      background: #ffffff;
      overflow: hidden;
    }}
    .splash {{
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      height: 100%;
      background: #ffffff;
      color: #555;
      gap: 12px;
    }}
    .spinner {{
      width: 40px; height: 40px;
      border: 3px solid #e0e0e0;
      border-top-color: #58a6ff;
      border-radius: 50%;
      animation: spin 0.8s linear infinite;
    }}
    @keyframes spin {{ to {{ transform: rotate(360deg); }} }}
    .error-box {{
      background: #fff8f8;
      border: 1px solid #ffcdd2;
      border-radius: 8px;
      padding: 20px 24px;
      max-width: 600px;
      text-align: left;
    }}
    .error-box h3 {{ color: #c62828; margin-bottom: 8px; }}
    .error-box pre {{ background: #f5f5f5; padding: 10px; border-radius: 4px; font-size: 12px; overflow: auto; }}
    .error-box .tip {{ margin-top: 12px; font-size: 13px; color: #555; }}
  </style>
</head>
<body>
  <header>
    <span class="logo">🍄</span>
    <div>
      <div class="brand">MycoNote</div>
      <div class="subtitle">{title}</div>
    </div>
    <span class="badge">JBrowse2 · powered by MycoNote CLI</span>
  </header>

  <div id="jbrowse-container">
    <div class="splash" id="splash">
      <div class="spinner"></div>
      <span>Loading JBrowse2…</span>
    </div>
  </div>

  <!-- JBrowse2 React Linear Genome View (CDN) -->
  <script src="https://unpkg.com/react@17/umd/react.production.min.js" crossorigin></script>
  <script src="https://unpkg.com/react-dom@17/umd/react-dom.production.min.js" crossorigin></script>
  <script src="https://unpkg.com/@jbrowse/react-linear-genome-view@2/dist/react-linear-genome-view.umd.production.min.js" crossorigin></script>

  <script>
    (function () {{
      var config = {safe_json};

      function init() {{
        var splash = document.getElementById('splash');
        var container = document.getElementById('jbrowse-container');

        try {{
          var JB = JBrowseReactLinearGenomeView;
          var state = JB.createViewState(config);

          splash.remove();

          ReactDOM.render(
            React.createElement(JB.JBrowseLinearGenomeView, {{ viewState: state }}),
            container
          );
        }} catch (err) {{
          container.innerHTML = [
            '<div class="splash"><div class="error-box">',
            '<h3>⚠ JBrowse2 failed to initialise</h3>',
            '<p>This usually means the CDN scripts could not be loaded.</p>',
            '<pre>' + String(err) + '</pre>',
            '<p class="tip">Make sure you have an internet connection, then reload.<br>',
            'Alternatively, serve this file with:<br>',
            '<code>python3 -m http.server 8080</code></p>',
            '</div></div>'
          ].join('');
          console.error('[MycoNote] JBrowse2 init error:', err);
        }}
      }}

      if (document.readyState === 'loading') {{
        document.addEventListener('DOMContentLoaded', init);
      }} else {{
        init();
      }}
    }})();
  </script>
</body>
</html>"#,
        title = title,
        safe_json = safe_json,
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Generate a self-contained JBrowse2 HTML viewer for the given GFF3 file.
pub fn generate_jbrowse_html(config: &ViewConfig) -> Result<()> {
    println!("🧬 Generating JBrowse2 viewer…");

    let parsed = parse_gff(&config.gff_path, &config.names)?;

    // Warn if the annotation is very large (inline JSON will be heavy)
    if parsed.total_count > 80_000 {
        eprintln!(
            "⚠  Large annotation ({} features). The HTML file will be large (~{} MB).",
            parsed.total_count,
            parsed.total_count / 5000
        );
        eprintln!(
            "   For better performance, consider using --region to limit the view,\n   \
             or install JBrowse2 Desktop: https://jbrowse.org/jb2/download/"
        );
    }

    let assembly_name = config
        .assembly_name
        .clone()
        .unwrap_or_else(|| {
            config.gff_path
                .file_stem()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "genome".to_string())
        });

    let title = config
        .title
        .clone()
        .unwrap_or_else(|| format!("{} — Gene Annotations", assembly_name));

    // Build JBrowse2 config object
    let assembly = build_assembly(
        &assembly_name,
        &parsed.chrom_sizes,
        config.fasta_path.as_deref(),
    );
    let track   = build_annotation_track(&assembly_name, parsed.features);
    let session = build_default_session(
        &assembly_name,
        &parsed.chrom_sizes,
        config.region.as_deref(),
    );

    let jbrowse_cfg = json!({
        "assembly":       assembly,
        "tracks":         [track],
        "defaultSession": session
    });

    let config_json = serde_json::to_string(&jbrowse_cfg)
        .map_err(|e| crate::utils::error::MycoNoteError::InvalidFormat(e.to_string()))?;

    // Render and write HTML
    let html = render_html(&title, &config_json);
    std::fs::write(&config.output, html)?;

    // User-facing summary
    let abs = config.output.canonicalize().unwrap_or_else(|_| config.output.clone());
    println!("✓ JBrowse2 HTML saved:  {}", abs.display());
    println!("  Open in browser:      file://{}", abs.display());
    println!();
    println!("  {} chromosomes/contigs  •  {} top-level features",
        parsed.chrom_sizes.len(),
        parsed.total_count,
    );
    if let Some(ref reg) = config.region {
        println!("  Focused on region:    {}", reg);
    }
    if config.fasta_path.is_none() {
        println!();
        println!("💡 Tip: supply --fasta <genome.fa> to enable the reference sequence track.");
        println!("   (The .fai index must exist alongside the FASTA file.)");
    }
    println!();
    println!("📡 For interactive team sharing, serve locally:");
    println!("   python3 -m http.server 8080");
    println!("   → http://localhost:8080/{}", config.output.display());

    Ok(())
}
