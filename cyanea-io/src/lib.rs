//! File format parsing for the Cyanea bioinformatics ecosystem.
//!
//! Supported formats:
//! - **CSV/TSV** — via the `csv` feature (enabled by default)
//! - **VCF** — Variant Call Format, via the `vcf` feature
//! - **BED** — Browser Extensible Data, via the `bed` feature
//! - **GFF3** — General Feature Format, via the `gff` feature
//! - **SAM** — Sequence Alignment/Map (text format), via the `sam` feature

#[cfg(feature = "csv")]
pub mod csv;

#[cfg(feature = "vcf")]
pub mod vcf;

#[cfg(feature = "bed")]
pub mod bed;

#[cfg(feature = "bed")]
pub mod bedpe;

#[cfg(feature = "gff")]
pub mod gff;

#[cfg(feature = "sam")]
pub mod sam;

#[cfg(feature = "sam")]
pub mod pileup;

#[cfg(feature = "bam")]
pub mod bam;

#[cfg(feature = "parquet")]
pub mod parquet;

#[cfg(feature = "cram")]
pub mod cram;

// Re-exports for convenience.

#[cfg(feature = "csv")]
pub use csv::{csv_preview, parse_csv_info, CsvInfo};

#[cfg(feature = "vcf")]
pub use vcf::{parse_vcf, vcf_stats, VcfStats};

#[cfg(feature = "bed")]
pub use bed::{bed_stats, parse_bed, parse_bed_intervals, BedRecord, BedStats};

#[cfg(feature = "bed")]
pub use bedpe::{bedpe_stats, parse_bedpe, BedpeRecord, BedpeStats};

#[cfg(feature = "gff")]
pub use gff::{gff3_stats, parse_gff3, GffStats};

#[cfg(feature = "sam")]
pub use sam::{
    filter_proper_pairs, pair_sam_records, paired_sam_stats, parse_sam, sam_stats,
    sam_stats_from_path, PairedSamStats, SamPair, SamRecord, SamStats,
};

#[cfg(feature = "sam")]
pub use pileup::{depth_stats, pileup, pileup_region, pileup_to_mpileup, DepthStats, Pileup, PileupColumn};

#[cfg(all(feature = "sam", feature = "vcf"))]
pub use pileup::call_snps;

#[cfg(feature = "bam")]
pub use bam::{parse_bam, bam_stats, BamReference};

#[cfg(feature = "parquet")]
pub use parquet::{
    parquet_info, parquet_interval_stats, parquet_variant_stats, read_intervals_parquet,
    read_variants_parquet, write_intervals_parquet, write_variants_parquet, ParquetInfo,
};

#[cfg(feature = "cram")]
pub use cram::{cram_stats, cram_stats_default, parse_cram, parse_cram_default, CramConfig};
