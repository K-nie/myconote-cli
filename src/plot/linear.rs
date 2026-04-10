use crate::parser::gff::{GFFReader, GFFRecord};
use crate::parser::region::{Region, RegionSelector};
use crate::plot::PlotConfig;
use crate::utils::error::Result;
use plotters::prelude::*;
use std::path::Path;

pub fn plot<P: AsRef<Path>>(gff_path: P, config: &PlotConfig) -> Result<()> {
    println!("Generating linear genome plot: {}", config.output_path);

    // Parse the GFF file
    let mut features: Vec<GFFRecord> = Vec::new();
    let reader = GFFReader::from_path(gff_path)?;

    // Determine region to plot
    let selector = if let Some(region_str) = &config.region {
        let region = Region::parse(region_str)?;
        let mut sel = RegionSelector::new();
        sel.include_regions.push(region);
        sel.include_all = false;
        sel
    } else {
        RegionSelector::new()
    };

    // Collect features
    let mut min_pos = u64::MAX;
    let mut max_pos = 0;
    let mut chromosome = String::new();

    for result in reader {
        let record = result?;

        if selector.should_include_feature(&record.seqid, record.start, record.end) {
            if chromosome.is_empty() {
                chromosome = record.seqid.clone();
            }

            min_pos = min_pos.min(record.start);
            max_pos = max_pos.max(record.end);
            features.push(record);
        }
    }

    if features.is_empty() {
        println!("No features found in the specified region");
        return Ok(());
    }

    // Add padding (10% on each side)
    let range = (max_pos - min_pos) as f64;
    let padding = (range * 0.1) as u64;
    let plot_min = min_pos.saturating_sub(padding);
    let plot_max = max_pos + padding;

    // Create the plot
    let root =
        BitMapBackend::new(&config.output_path, (config.width, config.height)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            config
                .title
                .as_deref()
                .unwrap_or(&format!("Genome Map - {}", chromosome)),
            ("sans-serif", 30).into_font(),
        )
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(plot_min as f64..plot_max as f64, 0.0..10.0)?;

    chart
        .configure_mesh()
        .x_desc("Position (bp)")
        .y_desc("")
        .disable_y_mesh()
        .disable_y_axis()
        .x_labels(10)
        .x_label_formatter(&|v| {
            if *v >= 1_000_000.0 {
                format!("{:.1} Mb", *v / 1_000_000.0)
            } else if *v >= 1_000.0 {
                format!("{:.1} Kb", *v / 1_000.0)
            } else {
                format!("{}", *v as u64)
            }
        })
        .draw()?;

    // Group features by type
    let mut genes: Vec<&GFFRecord> = Vec::new();
    let mut transcripts: Vec<&GFFRecord> = Vec::new();
    let mut cds_features: Vec<&GFFRecord> = Vec::new();
    let mut exons: Vec<&GFFRecord> = Vec::new();

    for feature in &features {
        match feature.feature_type.as_str() {
            "gene" => genes.push(feature),
            "mRNA" | "transcript" => transcripts.push(feature),
            "CDS" => cds_features.push(feature),
            "exon" => exons.push(feature),
            _ => {}
        }
    }

    let mut y_pos = 2.0;
    let track_height = 1.5;
    let track_spacing = 2.0;

    // Draw genes
    if config.show_genes && !genes.is_empty() {
        // Add track label
        chart.draw_series(std::iter::once(Text::new(
            "Genes",
            (plot_min as f64, y_pos + track_height / 2.0),
            ("sans-serif", 12).into_font(),
        )))?;

        chart.draw_series(genes.iter().map(|gene| {
            let y0 = y_pos;
            let y1 = y_pos + track_height;
            Rectangle::new(
                [(gene.start as f64, y0), (gene.end as f64, y1)],
                RGBColor(100, 149, 237).filled().stroke_width(2), // Cornflower blue with border
            )
        }))?;

        // Add gene label only once per gene, not repeated
        for gene in genes {
            let mid = (gene.start + gene.end) as f64 / 2.0;
            let label = gene.id().unwrap_or(&"gene".to_string()).clone();
            chart.draw_series(std::iter::once(Text::new(
                label,
                (mid, y_pos + track_height + 0.3),
                ("sans-serif", 12).into_font(),
            )))?;
        }
        y_pos += track_spacing;
    }

    // Draw transcripts
    if config.show_genes && !transcripts.is_empty() {
        chart.draw_series(std::iter::once(Text::new(
            "Transcripts",
            (plot_min as f64, y_pos + track_height / 2.0),
            ("sans-serif", 12).into_font(),
        )))?;

        chart.draw_series(transcripts.iter().map(|transcript| {
            let y0 = y_pos;
            let y1 = y_pos + track_height * 0.7;
            Rectangle::new(
                [(transcript.start as f64, y0), (transcript.end as f64, y1)],
                RGBColor(255, 165, 0).filled().stroke_width(1), // Orange
            )
        }))?;
        y_pos += track_spacing;
    }

    // Draw exons
    if config.show_exons && !exons.is_empty() {
        chart.draw_series(std::iter::once(Text::new(
            "Exons",
            (plot_min as f64, y_pos + track_height / 2.0),
            ("sans-serif", 12).into_font(),
        )))?;

        chart.draw_series(exons.iter().map(|exon| {
            let y0 = y_pos;
            let y1 = y_pos + track_height;
            Rectangle::new(
                [(exon.start as f64, y0), (exon.end as f64, y1)],
                RGBColor(255, 99, 71).filled().stroke_width(1), // Tomato
            )
        }))?;
        y_pos += track_spacing;
    }

    // Draw CDS
    if config.show_cds && !cds_features.is_empty() {
        chart.draw_series(std::iter::once(Text::new(
            "CDS",
            (plot_min as f64, y_pos + track_height / 2.0),
            ("sans-serif", 12).into_font(),
        )))?;

        chart.draw_series(cds_features.iter().map(|cds| {
            let y0 = y_pos;
            let y1 = y_pos + track_height;
            Rectangle::new(
                [(cds.start as f64, y0), (cds.end as f64, y1)],
                RGBColor(50, 205, 50).filled().stroke_width(1), // Lime green
            )
        }))?;
    }

    root.present()?;
    println!("✓ Plot saved to {}", config.output_path);

    Ok(())
}
