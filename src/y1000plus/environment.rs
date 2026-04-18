//! Isolation-source ontology lookup for Y1000+ species.
//!
//! Parses `budding_yeasts_isolations.owl` from the `environment` subset into
//! a species → human-readable-niche-label map. The OWL is structured as:
//!
//! ```text
//! <owl:Class rdf:about="…#Soil_environment_type">
//!   <rdfs:label>Soil</rdfs:label>
//! </owl:Class>
//!
//! <owl:NamedIndividual rdf:about="…">
//!   <rdf:type rdf:resource="…#Soil_environment_type"/>
//!   <rdfs:label>yHMPu5000038081_candida_takata_180604</rdfs:label>
//! </owl:NamedIndividual>
//! ```
//!
//! We walk the XML once, first collecting Class URI → label, then pairing
//! each NamedIndividual's type URI with its species label. Some types use
//! opaque `webprotege.stanford.edu/R…` URIs without a friendly label; in
//! those cases we fall back to a cleaned-up tail of the URI fragment.
//!
//! Citation: Opulente DA et al. (2024). Science 384(6694): eadj4503.

use crate::utils::error::{MycoNoteError, Result};
use crate::y1000plus::metabolism::normalise_species;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use xml::reader::{EventReader, XmlEvent};

/// A species' (possibly missing) inferred ecological niche.
#[derive(Debug, Clone)]
pub struct Niche {
    pub species_pretty: String,
    pub niche_label: String,
}

#[derive(Debug, Clone, Default)]
pub struct EnvironmentIndex {
    /// normalised "genus species" → niche entry.
    pub by_species: HashMap<String, Niche>,
}

impl EnvironmentIndex {
    pub fn lookup(&self, kegg_species_stem: &str) -> Option<&Niche> {
        self.by_species.get(&normalise_species(kegg_species_stem))
    }
    pub fn len(&self) -> usize {
        self.by_species.len()
    }
}

pub fn load_index(cache_root: &Path) -> Result<EnvironmentIndex> {
    let root = cache_root.join("environment");
    let owl = root.join("budding_yeasts_isolations.owl");
    if !owl.exists() {
        return Err(MycoNoteError::UnsupportedFormat(format!(
            "No isolation ontology at {}. Install with \
             `setup --y1000plus --include environment`.",
            owl.display()
        )));
    }
    parse_owl(&owl)
}

/// Single streaming pass — collects Class URI → label as we see Classes,
/// and NamedIndividual → (type_uri, species_label) pairs, then joins.
fn parse_owl(path: &Path) -> Result<EnvironmentIndex> {
    let f = File::open(path).map_err(MycoNoteError::Io)?;
    let parser = EventReader::new(BufReader::new(f));

    #[derive(Default)]
    struct Block {
        /// rdf:about URI of the current Class or NamedIndividual.
        about: String,
        /// rdf:type rdf:resource for a NamedIndividual (URI of the class).
        type_uri: String,
        /// rdfs:label text collected inside this block.
        label: String,
    }

    // URI → friendly label, collected from Class blocks.
    let mut class_labels: HashMap<String, String> = HashMap::new();
    // (species_label, type_uri) pairs, collected from NamedIndividual blocks.
    let mut individuals: Vec<(String, String)> = Vec::new();

    let mut current_block: Option<&'static str> = None; // "class" | "individual"
    let mut block = Block::default();
    let mut in_label = false;

    for ev in parser {
        let ev = ev.map_err(|e| MycoNoteError::InvalidFormat(format!("OWL parse: {e}")))?;
        match ev {
            XmlEvent::StartElement {
                name, attributes, ..
            } => {
                let local = name.local_name.as_str();
                match local {
                    "Class" | "NamedIndividual" => {
                        block = Block::default();
                        block.about = attr_value(&attributes, "about").unwrap_or_default();
                        current_block = Some(if local == "Class" {
                            "class"
                        } else {
                            "individual"
                        });
                    }
                    "type" if current_block == Some("individual") => {
                        if let Some(uri) = attr_value(&attributes, "resource") {
                            block.type_uri = uri;
                        }
                    }
                    "label" => {
                        in_label = true;
                        block.label.clear();
                    }
                    _ => {}
                }
            }
            XmlEvent::Characters(text) | XmlEvent::CData(text) => {
                if in_label {
                    block.label.push_str(&text);
                }
            }
            XmlEvent::EndElement { name } => {
                let local = name.local_name.as_str();
                match local {
                    "label" => in_label = false,
                    "Class" if current_block == Some("class") => {
                        if !block.about.is_empty() && !block.label.is_empty() {
                            class_labels
                                .insert(block.about.clone(), block.label.trim().to_string());
                        }
                        current_block = None;
                    }
                    "NamedIndividual" if current_block == Some("individual") => {
                        if !block.type_uri.is_empty() && !block.label.is_empty() {
                            individuals
                                .push((block.label.trim().to_string(), block.type_uri.clone()));
                        }
                        current_block = None;
                    }
                    _ => {}
                }
            }
            _ => {}
        }
    }

    // Resolve each species' type URI to a friendly niche label.
    let mut index = EnvironmentIndex::default();
    for (species_label, type_uri) in individuals {
        let niche_label = class_labels
            .get(&type_uri)
            .cloned()
            .unwrap_or_else(|| uri_tail_to_label(&type_uri));
        let niche_label = humanise_niche(&niche_label);
        let key = normalise_species(&species_label);
        if key.is_empty() {
            continue;
        }
        index.by_species.insert(
            key,
            Niche {
                species_pretty: species_label,
                niche_label,
            },
        );
    }
    Ok(index)
}

fn attr_value(attrs: &[xml::attribute::OwnedAttribute], local: &str) -> Option<String> {
    attrs
        .iter()
        .find(|a| a.name.local_name == local)
        .map(|a| a.value.clone())
}

/// Fallback label when a type URI has no associated rdfs:label: take the
/// URI fragment (after `#` or final `/`) and humanise it.
fn uri_tail_to_label(uri: &str) -> String {
    let tail = uri
        .rsplit_once('#')
        .map(|(_, t)| t)
        .unwrap_or_else(|| uri.rsplit_once('/').map(|(_, t)| t).unwrap_or(uri));
    tail.to_string()
}

fn humanise_niche(raw: &str) -> String {
    // Replace underscores with spaces, collapse trailing "_environment_type",
    // and title-case the first character. We deliberately keep the raw
    // taxonomic terms (anus_of_cat, Soil_environment_type) intact so the
    // user can see the literal ontology term — only gently prettified.
    let mut t = raw.trim().to_string();
    for suffix in &["_environment_type", "_environmental_type", "_type"] {
        if t.ends_with(suffix) {
            t.truncate(t.len() - suffix.len());
            break;
        }
    }
    t = t.replace('_', " ");
    if t.is_empty() {
        return "unknown".to_string();
    }
    // Capitalise first letter only (leave rest untouched).
    let mut chars = t.chars();
    let first = chars.next().unwrap().to_ascii_uppercase();
    format!("{first}{}", chars.collect::<String>())
}
