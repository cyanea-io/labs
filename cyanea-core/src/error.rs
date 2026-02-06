//! Structured error types for the Cyanea ecosystem.

use thiserror::Error;

/// Unified error type for all Cyanea operations.
#[derive(Debug, Error)]
pub enum CyaneaError {
    /// I/O error (file not found, permission denied, etc.)
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Parse error (malformed input data)
    #[error("parse error: {0}")]
    Parse(String),

    /// Invalid input (bad arguments, out-of-range values)
    #[error("invalid input: {0}")]
    InvalidInput(String),

    /// Compression or decompression failure
    #[error("compression error: {0}")]
    Compression(String),

    /// Hashing failure
    #[error("hash error: {0}")]
    Hash(String),

    /// Catch-all for other errors
    #[error("{0}")]
    Other(String),
}

/// Convenience alias used throughout the Cyanea ecosystem.
pub type Result<T> = std::result::Result<T, CyaneaError>;
