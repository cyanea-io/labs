//! Metagenomics analysis for the Cyanea bioinformatics ecosystem.
//!
//! `cyanea-meta` provides comprehensive tools for metagenomic data analysis:
//!
//! - **Taxonomy** — [`taxonomy::TaxonomyDB`] with k-mer classification, LCA queries,
//!   and lineage retrieval for NCBI-style taxonomic trees
//! - **Profiling** — [`profile::TaxonomicProfile`] for abundance estimation and
//!   profile aggregation across samples
//! - **Diversity** — Community diversity metrics: alpha diversity (Shannon, Simpson,
//!   Chao1, ACE, Fisher's alpha), beta diversity (Bray-Curtis, Jaccard),
//!   and rarefaction curves
//! - **Composition** — Compositional data analysis: CLR/ILR transforms,
//!   differential abundance (ALDEx2-style, ANCOM)
//! - **Functional** — Functional profiling: mapping taxa to KEGG/COG functions,
//!   pathway abundance, functional diversity
//! - **Binning** — Metagenomic binning by tetranucleotide frequency and coverage,
//!   with bin quality assessment (completeness/contamination)
//! - **Assembly** — Quality-control metrics: N50, L50, N90, L90, GC content,
//!   mean/median lengths, auN
//!
//! # Example: Classification and Profiling
//!
//! ```ignore
//! use cyanea_meta::taxonomy::{TaxonomyDB, TaxonRank};
//! use cyanea_meta::profile::profile_from_classifications;
//!
//! // Build a taxonomy
//! let mut db = TaxonomyDB::new(1);
//! db.add_node(2, 1, TaxonRank::Phylum, "Bacteria").unwrap();
//!
//! // Index reference sequences
//! db.add_reference(b"ACGTACGTACGTACGTACGTACGTACGTACG", 2).unwrap();
//!
//! // Classify a metagenomic read
//! let read = b"ACGTACGTACGTACGTACGTACGTACGTACG";
//! if let Ok(Some(taxid)) = db.classify_sequence(read) {
//!     println!("Classified to taxon {}", taxid);
//! }
//!
//! // Build a taxonomic profile from classifications
//! let classifications = vec![2, 2, 3, 4];
//! let profile = profile_from_classifications(&classifications).unwrap();
//! ```
//!
//! # Example: Diversity Analysis
//!
//! ```ignore
//! use cyanea_meta::diversity::{alpha_diversity, rarefaction_curve};
//!
//! let counts = vec![100, 50, 25, 10];
//! let ad = alpha_diversity(&counts).unwrap();
//! println!("Shannon: {}, Simpson: {}", ad.shannon, ad.simpson);
//!
//! let depths = vec![10, 50, 100, 150];
//! let rare = rarefaction_curve(&counts, &depths).unwrap();
//! ```

pub mod error;
pub mod taxonomy;
pub mod profile;
pub mod diversity;
pub mod composition;
pub mod functional;
pub mod binning;
pub mod assembly;

// Re-export common types
pub use error::{MetaError, Result};
pub use taxonomy::{TaxonomyDB, TaxonNode, TaxonRank};
pub use profile::TaxonomicProfile;
pub use diversity::{AlphaDiversity, BetaDiversityMatrix};
pub use composition::CompositionTransform;
pub use functional::FunctionalProfile;
pub use binning::{Contig, Bin};
pub use assembly::AssemblyStats;
