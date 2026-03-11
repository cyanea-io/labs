//! Epigenomics analysis for the Cyanea bioinformatics ecosystem.
//!
//! `cyanea-epi` provides comprehensive tools for analyzing epigenomic data from
//! ChIP-seq, ATAC-seq, and MNase-seq experiments:
//!
//! - **Peak calling** — MACS2-style narrow and broad peak calling with local background estimation
//! - **Signal analysis** — Pileup construction, normalization, smoothing, and replicate correlation
//! - **Motif discovery** — K-mer enrichment-based motif discovery and PWM scanning
//! - **Chromatin states** — ChromHMM-like state learning via EM and genomic segmentation
//! - **Differential analysis** — DESeq2-style differential binding with negative binomial test
//! - **Nucleosome positioning** — MNase-seq nucleosome calling with periodicity detection
//! - **ATAC-seq QC** — TSS enrichment, FRiP, NFR ratio, and fragment size metrics
//!
//! # Example
//!
//! ```
//! use cyanea_epi::pileup::{build_pileup, TagPileup};
//!
//! // Build a pileup from aligned reads
//! let reads = vec![
//!     ("chr1".to_string(), 100, 50),
//!     ("chr1".to_string(), 150, 50),
//!     ("chr1".to_string(), 200, 50),
//! ];
//!
//! let pileup = build_pileup(&reads, 200);
//! assert!(pileup.coverage.contains_key("chr1"));
//! ```

pub mod error;
pub mod peaks;
pub mod pileup;
pub mod motifs;
pub mod chromatin;
pub mod differential;
pub mod nucleosome;
pub mod accessibility;

// Re-export error types
pub use error::{EpiError, Result};

// Re-export peak types
pub use peaks::{Peak, PeakCallParams, PeakSet, PeakStats, call_peaks, call_broad_peaks};

// Re-export pileup types and functions
pub use pileup::{
    TagPileup, build_pileup, normalize_pileup, smooth_pileup, pileup_correlation, fingerprint,
};

// Re-export motif types and functions
pub use motifs::{
    Motif, MotifMatch, DiscoveryParams, discover_motifs, scan_sequence, parse_meme, write_meme,
    compare_motifs, motif_enrichment,
};

// Re-export chromatin types and functions
pub use chromatin::{
    ChromatinState, ChromHMMModel, ChromHMMParams, ChromatinSegmentation,
    learn_chromatin_states, segment_genome, state_enrichment,
};

// Re-export differential types and functions
pub use differential::{
    DiffResult, differential_peaks, count_reads_in_peaks, ma_plot_data,
};

// Re-export nucleosome types and functions
pub use nucleosome::{
    NucleosomePosition, NucleosomeParams, call_nucleosomes, nfr_score, periodicity,
};

// Re-export accessibility types and functions
pub use accessibility::{
    InsertSizeMetrics, AtacQcResult, tss_enrichment, fragment_size_distribution,
    insert_size_metrics, frip, atacqc,
};
