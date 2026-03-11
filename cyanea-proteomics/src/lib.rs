//! Proteomics and mass spectrometry analysis for the Cyanea bioinformatics ecosystem.
//!
//! `cyanea-proteomics` provides comprehensive tools for proteomics data analysis:
//!
//! - **Spectrum** — [`spectrum::MassSpectrum`] with peak operations, filtering,
//!   normalization, deisotoping, and run-level statistics
//! - **Peptide** — [`peptide::Peptide`] with amino acid masses, modifications,
//!   fragment ion generation (b/y/a series), and in-silico digestion
//!   (trypsin, LysC, chymotrypsin, AspN, GluC)
//! - **MGF** — [`mgf::parse_mgf`] / [`mgf::write_mgf`] for Mascot Generic Format
//! - **mzML** — [`mzml::parse_mzml_text`] for simplified mzML parsing
//! - **Search** — [`search::score_peptide`] / [`search::search_all`] for
//!   peptide-spectrum matching with XCorr and hyperscore
//! - **Protein** — [`protein::infer_proteins`] for parsimony-based protein inference
//! - **Quantification** — [`quantification::spectral_counting`],
//!   [`quantification::intensity_quantification`], [`quantification::quantify_tmt`]
//!   for label-free and TMT/iTRAQ quantification
//! - **FDR** — [`fdr::estimate_fdr`] / [`fdr::filter_fdr`] for target-decoy
//!   FDR estimation and q-value calculation
//! - **mzTab** — [`mztab::write_mztab_proteins`] / [`mztab::write_mztab_psms`]
//!   for standardized result output
//!
//! # Example: Database Search Pipeline
//!
//! ```ignore
//! use cyanea_proteomics::mgf::parse_mgf;
//! use cyanea_proteomics::peptide::{digest, DigestConfig};
//! use cyanea_proteomics::search::{search_all, generate_decoys, SearchConfig};
//! use cyanea_proteomics::fdr::filter_fdr;
//! use cyanea_proteomics::protein::infer_proteins;
//!
//! // Parse spectra
//! let spectra = parse_mgf(mgf_text).unwrap();
//!
//! // Digest protein database
//! let config = DigestConfig::default();
//! let peptides = digest(b"MAAAKPEPTIDEKSEQUENCER", &config).unwrap();
//!
//! // Search
//! let search_config = SearchConfig::default();
//! let target_psms = search_all(&spectra, &peptides, &search_config).unwrap();
//!
//! // Generate decoys and search again
//! let decoys = generate_decoys(&peptides);
//! let decoy_psms = search_all(&spectra, &decoys, &search_config).unwrap();
//!
//! // FDR filtering at 1%
//! let filtered = filter_fdr(&target_psms, &decoy_psms, 0.01).unwrap();
//! ```

pub mod error;
pub mod spectrum;
pub mod peptide;
pub mod mgf;
pub mod mzml;
pub mod search;
pub mod protein;
pub mod quantification;
pub mod fdr;
pub mod mztab;

// Re-export common types
pub use error::{ProteomicsError, Result};
pub use spectrum::{MassSpectrum, Peak, MsLevel, Precursor};
pub use peptide::{Peptide, Protease, DigestConfig, Modification};
pub use search::Psm;
pub use protein::{ProteinGroup, ProteinEntry};
pub use quantification::{ProteinQuant, QuantMethod, TmtPlex};
pub use fdr::{FdrResult, FdrConfig};
