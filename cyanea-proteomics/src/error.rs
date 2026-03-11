//! Error types for proteomics operations.

use thiserror::Error;

/// Errors specific to proteomics analysis.
#[derive(Debug, Error)]
pub enum ProteomicsError {
    /// Mass spectrum is invalid or malformed.
    #[error("spectrum error: {0}")]
    Spectrum(String),

    /// File parsing failed (mzML, mzXML, MGF, mzTab).
    #[error("parse error: {0}")]
    Parse(String),

    /// Peptide or protein sequence error.
    #[error("peptide error: {0}")]
    Peptide(String),

    /// Database search or scoring failed.
    #[error("search error: {0}")]
    Search(String),

    /// Quantification error.
    #[error("quantification error: {0}")]
    Quantification(String),

    /// FDR estimation error.
    #[error("fdr error: {0}")]
    Fdr(String),

    /// Wrapped core error.
    #[error(transparent)]
    Core(#[from] cyanea_core::CyaneaError),
}

/// Convenience alias for proteomics results.
pub type Result<T> = std::result::Result<T, ProteomicsError>;
