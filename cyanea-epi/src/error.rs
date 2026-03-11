//! Error types for epigenomics analysis.

use thiserror::Error;

/// Errors that can occur during epigenomics analysis.
#[derive(Debug, Error)]
pub enum EpiError {
    /// Invalid peak or parameter
    #[error("invalid peak: {0}")]
    InvalidPeak(String),

    /// Invalid motif data
    #[error("invalid motif: {0}")]
    InvalidMotif(String),

    /// Invalid chromatin state model
    #[error("invalid chromatin state: {0}")]
    InvalidChromatinState(String),

    /// Insufficient data for operation
    #[error("insufficient data: {0}")]
    InsufficientData(String),

    /// Convergence failure
    #[error("convergence failed: {0}")]
    ConvergenceFailed(String),

    /// Invalid input parameters
    #[error("invalid input: {0}")]
    InvalidInput(String),

    /// Underlying Cyanea error
    #[error("cyanea error: {0}")]
    CyaneaError(#[from] cyanea_core::CyaneaError),
}

/// Convenience alias for Result with EpiError.
pub type Result<T> = std::result::Result<T, EpiError>;
