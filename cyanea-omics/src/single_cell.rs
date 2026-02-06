//! Stub module for single-cell omics data structures.
//!
//! This module will provide an AnnData-like in-memory representation for
//! single-cell data, inspired by:
//!
//! - **AnnData** (Python) — the standard data structure in the scverse ecosystem
//! - **scanpy** — single-cell analysis toolkit built on AnnData
//! - **Seurat** (R) — widely used single-cell analysis framework
//!
//! Planned structures:
//!
//! - `AnnData` — the central container, holding:
//!   - `X` — the primary data matrix (typically counts or normalized expression),
//!     backed by [`crate::SparseMatrix`] or [`crate::ExpressionMatrix`]
//!   - `obs` — per-cell (observation) metadata
//!   - `var` — per-gene (variable) metadata
//!   - `obsm` / `varm` — multi-dimensional annotations (e.g. PCA, UMAP embeddings)
//!   - `obsp` / `varp` — pairwise annotations (e.g. distance matrices, graphs)
//!   - `layers` — alternative data matrices aligned to the same cells and genes
//!   - `uns` — unstructured metadata
//!
//! - `CellType` — enum or ontology-backed cell type annotation
//! - `QualityMetrics` — per-cell QC metrics (n_genes, total_counts, pct_mito, etc.)
//!
//! This module is intentionally empty for now. Implementation will follow once
//! the foundational matrix types and I/O layer are stabilized.
