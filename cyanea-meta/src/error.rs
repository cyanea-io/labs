//! Error types for metagenomics operations.

use thiserror::Error;

/// Errors specific to metagenomics analysis.
#[derive(Debug, Error)]
pub enum MetaError {
    /// Taxonomy database is missing or malformed.
    #[error("taxonomy error: {0}")]
    Taxonomy(String),

    /// Classification failed (no hits, ambiguous assignment).
    #[error("classification error: {0}")]
    Classification(String),

    /// Profile is invalid or incompatible.
    #[error("profile error: {0}")]
    Profile(String),

    /// Binning operation failed.
    #[error("binning error: {0}")]
    Binning(String),

    /// Functional annotation error.
    #[error("functional error: {0}")]
    Functional(String),

    /// Wrapped core error.
    #[error(transparent)]
    Core(#[from] cyanea_core::CyaneaError),
}

/// Convenience alias for metagenomics results.
pub type Result<T> = std::result::Result<T, MetaError>;
