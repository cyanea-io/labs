//! SIMD-accelerated alignment kernels (not yet implemented).
//!
//! This module will provide vectorised pairwise alignment using SIMD intrinsics,
//! targeting both x86-64 (SSE4.1, AVX2, AVX-512) and AArch64 (NEON) instruction sets.
//!
//! # Planned features
//!
//! - **Striped DP** — The Farrar (2007) striped Smith-Waterman using query-profile
//!   vectors and inter-lane shuffles for the dependency chain. Expected 8-16x speedup
//!   over scalar for typical protein alignments.
//!
//! - **Banded alignment** — Restrict the DP matrix to a diagonal band of width 2k+1,
//!   reducing work from O(mn) to O(m·k). Combinable with SIMD for another ~4x.
//!
//! - **Score-only mode** — When only the alignment score is needed (no traceback),
//!   the DP can operate in O(n) space using two rows and skip the traceback matrix
//!   entirely.
//!
//! - **Auto dispatch** — Runtime detection of the best available instruction set via
//!   `std::is_x86_feature_detected!` / `std::is_aarch64_feature_detected!`, with
//!   transparent fallback to scalar. Inspired by the approach in Rognes (2011)
//!   SWIPE and the parasail library.
//!
//! # References
//!
//! - Farrar M. (2007) "Striped Smith-Waterman speeds database searches six times
//!   over other SIMD implementations." *Bioinformatics* 23(2):156-161.
//! - Rognes T. (2011) "Faster Smith-Waterman database searches with inter-sequence
//!   SIMD parallelisation." *BMC Bioinformatics* 12:221.
