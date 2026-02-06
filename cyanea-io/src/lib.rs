//! File format parsing for the Cyanea bioinformatics ecosystem.
//!
//! Currently supports:
//! - **CSV/TSV** — via the `csv` feature (enabled by default)

#[cfg(feature = "csv")]
pub mod csv;

#[cfg(feature = "csv")]
pub use csv::{csv_preview, parse_csv_info, CsvInfo};
