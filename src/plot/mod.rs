use crate::utils::error::Result;
use std::path::Path;

pub mod circular;
pub mod linear;

#[derive(Debug, Clone)]
pub struct PlotConfig {
    pub output_path: String,
    pub width: u32,
    pub height: u32,
    pub title: Option<String>,
    pub show_genes: bool,
    pub show_cds: bool,
    pub show_exons: bool,
    pub color_scheme: String,
    pub region: Option<String>,
}

impl Default for PlotConfig {
    fn default() -> Self {
        Self {
            output_path: "output.png".to_string(),
            width: 1200,
            height: 800,
            title: None,
            show_genes: true,
            show_cds: true,
            show_exons: true,
            color_scheme: "default".to_string(),
            region: None,
        }
    }
}

impl PlotConfig {
    pub fn new(output: &str) -> Self {
        Self {
            output_path: output.to_string(),
            ..Default::default()
        }
    }

    pub fn with_size(mut self, width: u32, height: u32) -> Self {
        self.width = width;
        self.height = height;
        self
    }

    pub fn with_title(mut self, title: &str) -> Self {
        self.title = Some(title.to_string());
        self
    }

    pub fn with_region(mut self, region: &str) -> Self {
        self.region = Some(region.to_string());
        self
    }

    pub fn hide_genes(mut self) -> Self {
        self.show_genes = false;
        self
    }

    pub fn hide_cds(mut self) -> Self {
        self.show_cds = false;
        self
    }

    pub fn hide_exons(mut self) -> Self {
        self.show_exons = false;
        self
    }
}

pub fn generate_linear_plot<P: AsRef<Path>>(gff_path: P, config: &PlotConfig) -> Result<()> {
    linear::plot(gff_path, config)
}

pub fn generate_circular_plot<P: AsRef<Path>>(gff_path: P, config: &PlotConfig) -> Result<()> {
    circular::plot(gff_path, config)
}
