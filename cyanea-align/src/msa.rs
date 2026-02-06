//! Multiple sequence alignment (not yet implemented).
//!
//! This module will provide algorithms for aligning three or more sequences
//! simultaneously, producing a multi-row alignment that reveals conserved regions
//! across an entire family of sequences.
//!
//! # Planned features
//!
//! - **Progressive alignment** — Build a guide tree from pairwise distances
//!   (e.g. UPGMA or neighbour-joining), then progressively align sequences/profiles
//!   bottom-up. This is the strategy used by ClustalW/ClustalOmega.
//!
//! - **Profile HMMs** — Represent a multiple alignment as a hidden Markov model
//!   with match, insert, and delete states. Enables sensitive homology search and
//!   iterative alignment refinement. Based on the HMMER approach (Eddy, 2011).
//!
//! - **Iterative refinement** — After an initial progressive alignment, repeatedly
//!   remove each sequence and re-align it to the profile of the remaining sequences
//!   until convergence. Inspired by MUSCLE (Edgar, 2004).
//!
//! # References
//!
//! - Sievers F. et al. (2011) "Fast, scalable generation of high-quality protein
//!   multiple sequence alignments using Clustal Omega." *Mol Syst Biol* 7:539.
//! - Edgar R.C. (2004) "MUSCLE: multiple sequence alignment with high accuracy and
//!   high throughput." *Nucleic Acids Res* 32(5):1792-1797.
//! - Eddy S.R. (2011) "Accelerated Profile HMM Searches." *PLoS Comput Biol*
//!   7(10):e1002195.
