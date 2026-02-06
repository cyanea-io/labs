//! Sequence I/O and manipulation for the Cyanea bioinformatics ecosystem.
//!
//! - **FASTA/FASTQ parsing** — Streaming parser via [`fasta`] module
//! - **Sequence types** — Planned DNA/RNA/Protein types in [`sequence`]

pub mod fasta;
pub mod sequence;

pub use fasta::{parse_fasta_stats, FastaStats};
