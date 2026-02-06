//! Biological sequence types.
//!
//! This module will provide strongly-typed sequence representations:
//!
//! - `DnaSequence` — validated DNA (A, C, G, T, N)
//! - `RnaSequence` — validated RNA (A, C, G, U, N)
//! - `ProteinSequence` — amino acid sequences
//!
//! Each type will implement [`cyanea_core::Sequence`] and support
//! operations like reverse complement, translation, and k-mer iteration.
