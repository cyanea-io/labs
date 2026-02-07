//! File format parsing for the Cyanea bioinformatics ecosystem.
//!
//! Supported formats:
//! - **CSV/TSV** — via the `csv` feature (enabled by default)
//! - **VCF** — Variant Call Format, via the `vcf` feature
//! - **BED** — Browser Extensible Data, via the `bed` feature
//! - **GFF3** — General Feature Format, via the `gff` feature

#[cfg(feature = "csv")]
pub mod csv;

#[cfg(feature = "vcf")]
pub mod vcf;

#[cfg(feature = "bed")]
pub mod bed;

#[cfg(feature = "gff")]
pub mod gff;

// Re-exports for convenience.

#[cfg(feature = "csv")]
pub use csv::{csv_preview, parse_csv_info, CsvInfo};

#[cfg(feature = "vcf")]
pub use vcf::{parse_vcf, vcf_stats, VcfStats};

#[cfg(feature = "bed")]
pub use bed::{bed_stats, parse_bed, parse_bed_intervals, BedRecord, BedStats};

#[cfg(feature = "gff")]
pub use gff::{gff3_stats, parse_gff3, GffStats};
