use crate::utils::error::{MycoNoteError, Result};
use std::collections::HashSet;

#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    pub chromosome: String,
    pub start: Option<u64>,
    pub end: Option<u64>,
}

impl Region {
    pub fn new(chromosome: &str) -> Self {
        Self {
            chromosome: chromosome.to_string(),
            start: None,
            end: None,
        }
    }

    pub fn with_range(chromosome: &str, start: u64, end: u64) -> Self {
        Self {
            chromosome: chromosome.to_string(),
            start: Some(start),
            end: Some(end),
        }
    }

    pub fn contains(&self, chr: &str, pos: u64) -> bool {
        if chr != self.chromosome {
            return false;
        }

        match (self.start, self.end) {
            (Some(s), Some(e)) => pos >= s && pos <= e,
            _ => true, // No range specified means whole chromosome
        }
    }

    pub fn parse(region_str: &str) -> Result<Self> {
        // Format: "chr1" or "chr1:1000-2000"
        if region_str.contains(':') {
            let parts: Vec<&str> = region_str.split(':').collect();
            if parts.len() != 2 {
                return Err(MycoNoteError::ParseError {
                    line: 0,
                    message: format!("Invalid region format: {}. Use chr:start-end", region_str),
                });
            }

            let chromosome = parts[0].to_string();
            let range_parts: Vec<&str> = parts[1].split('-').collect();

            if range_parts.len() != 2 {
                return Err(MycoNoteError::ParseError {
                    line: 0,
                    message: format!("Invalid range format: {}. Use start-end", parts[1]),
                });
            }

            let start = range_parts[0]
                .parse::<u64>()
                .map_err(|_| MycoNoteError::ParseError {
                    line: 0,
                    message: format!("Invalid start coordinate: {}", range_parts[0]),
                })?;

            let end = range_parts[1]
                .parse::<u64>()
                .map_err(|_| MycoNoteError::ParseError {
                    line: 0,
                    message: format!("Invalid end coordinate: {}", range_parts[1]),
                })?;

            if start > end {
                return Err(MycoNoteError::ParseError {
                    line: 0,
                    message: format!("Start ({}) cannot be greater than end ({})", start, end),
                });
            }

            Ok(Region::with_range(&chromosome, start, end))
        } else {
            // Whole chromosome
            Ok(Region::new(region_str))
        }
    }
}

#[derive(Debug, Default)]
pub struct RegionSelector {
    pub include_regions: Vec<Region>,
    pub exclude_chromosomes: HashSet<String>,
    pub include_all: bool,
}

impl RegionSelector {
    pub fn new() -> Self {
        RegionSelector {
            include_regions: Vec::new(),
            exclude_chromosomes: HashSet::new(),
            include_all: true, // Change this to true by default!
        }
    }

    pub fn with_chromosomes(chromosomes: &[String]) -> Self {
        let mut selector = RegionSelector::new();
        selector.include_all = false; // Turn off include_all when specific chromosomes are given
        for chr in chromosomes {
            selector.include_regions.push(Region::new(chr));
        }
        selector
    }

    pub fn with_regions(regions: &[String]) -> Result<Self> {
        let mut selector = RegionSelector::new();
        selector.include_all = false; // Turn off include_all when specific regions are given
        for region_str in regions {
            selector.include_regions.push(Region::parse(region_str)?);
        }
        Ok(selector)
    }

    pub fn exclude_chromosome(&mut self, chromosome: &str) {
        self.exclude_chromosomes.insert(chromosome.to_string());
    }

    pub fn should_include(&self, chromosome: &str, position: u64) -> bool {
        // First check if chromosome is excluded
        if self.exclude_chromosomes.contains(chromosome) {
            return false;
        }

        // If include_all is true and no specific regions, include everything
        if self.include_all && self.include_regions.is_empty() {
            return true;
        }

        // Check if position falls within any included region
        for region in &self.include_regions {
            if region.contains(chromosome, position) {
                return true;
            }
        }

        false
    }

    pub fn should_include_feature(&self, chromosome: &str, start: u64, end: u64) -> bool {
        // First check if chromosome is excluded
        if self.exclude_chromosomes.contains(chromosome) {
            return false;
        }

        // If include_all is true and no specific regions, include everything
        if self.include_all && self.include_regions.is_empty() {
            return true;
        }

        // Check if feature overlaps with any included region
        for region in &self.include_regions {
            if region.chromosome == chromosome {
                match (region.start, region.end) {
                    (Some(rs), Some(re)) => {
                        // Check for overlap
                        if end >= rs && start <= re {
                            return true;
                        }
                    }
                    _ => return true, // Whole chromosome
                }
            }
        }

        false
    }
}
