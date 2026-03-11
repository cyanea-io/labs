//! Bundled sample datasets for the Cyanea bioinformatics ecosystem.
//!
//! Provides small, self-contained demo datasets for every domain — genomics,
//! alignment, epigenomics, single-cell, chemistry, phylogenetics, metagenomics,
//! and structural biology. All data is generated in-memory (no external files).
//!
//! Designed for:
//! - **Notebooks**: "open and run" experiences with real-looking data
//! - **Tests**: consistent, reproducible data for integration tests
//! - **Tutorials**: learning material for each cyanea crate
//!
//! # Example
//!
//! ```
//! use cyanea_datasets::genomics;
//!
//! let (name, seq) = genomics::ecoli_16s_rrna();
//! assert!(seq.len() > 500);
//!
//! let vcf = genomics::demo_vcf_chr22();
//! assert!(vcf.contains("#CHROM"));
//! ```

pub mod alignment;
pub mod chemistry;
pub mod epigenomics;
pub mod genomics;
pub mod metagenomics;
pub mod phylogenetics;
pub mod single_cell;
pub mod structural;
