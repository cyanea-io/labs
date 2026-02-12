//! Omics data structures for the Cyanea bioinformatics ecosystem.
//!
//! This crate provides core types for working with omics data:
//!
//! - **Genomic coordinates** — [`Strand`], [`GenomicPosition`], [`GenomicInterval`]
//! - **Interval collections** — [`IntervalSet`] with overlap queries
//! - **Expression matrices** — Dense [`ExpressionMatrix`] (features × samples)
//! - **Sparse matrices** — [`SparseMatrix`] in COO format
//! - **Variants** — VCF-style [`Variant`] representation
//! - **Gene annotations** — [`Gene`], [`Transcript`], [`Exon`] hierarchy
//!
//! # Quick start
//!
//! ```
//! use cyanea_omics::ExpressionMatrix;
//! use cyanea_core::Summarizable;
//!
//! let matrix = ExpressionMatrix::new(
//!     vec![vec![1.0, 2.0], vec![3.0, 4.0]],
//!     vec!["gene1".into(), "gene2".into()],
//!     vec!["sample_a".into(), "sample_b".into()],
//! ).unwrap();
//!
//! assert_eq!(matrix.shape(), (2, 2));
//! assert_eq!(matrix.get(0, 1), Some(2.0));
//! assert_eq!(matrix.summary(), "ExpressionMatrix: 2 features \u{00d7} 2 samples");
//! ```

pub mod genomic;
pub mod interval;
pub mod expr;
pub mod sparse;
pub mod variant;
pub mod annotation;
pub mod single_cell;
#[cfg(feature = "h5ad")]
pub mod h5ad;
#[cfg(feature = "zarr")]
pub mod zarr;

pub use genomic::{GenomicInterval, GenomicPosition, Strand};
pub use interval::IntervalSet;
pub use expr::ExpressionMatrix;
pub use sparse::SparseMatrix;
pub use variant::{Variant, VariantFilter, VariantType, Zygosity};
pub use annotation::{Exon, Gene, GeneType, Transcript};
pub use single_cell::ColumnData;
#[cfg(feature = "h5ad")]
pub use h5ad::{read_h5ad, write_h5ad};
#[cfg(feature = "zarr")]
pub use zarr::{read_zarr, write_zarr};
