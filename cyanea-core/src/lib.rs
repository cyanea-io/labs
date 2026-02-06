//! Shared primitives, traits, and utilities for the Cyanea bioinformatics ecosystem.
//!
//! `cyanea-core` provides the foundation that all other Cyanea crates build on:
//!
//! - **Error types** — [`CyaneaError`] and [`Result`] for structured error handling
//! - **Traits** — Core abstractions like [`Sequence`], [`ContentAddressable`], [`Compressible`]
//! - **Hashing** — SHA-256 content addressing for data integrity
//! - **Compression** — zstd and gzip with algorithm auto-detection
//! - **Memory mapping** — Zero-copy file access (std feature only)

pub mod error;
pub mod traits;
pub mod hash;

#[cfg(feature = "std")]
pub mod compress;

#[cfg(feature = "std")]
pub mod mmap;

pub use error::{CyaneaError, Result};
pub use traits::*;
