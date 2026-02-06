//! Dimensionality reduction (planned).
//!
//! This module will provide implementations for common dimensionality reduction
//! techniques used in bioinformatics:
//!
//! - **PCA** (Principal Component Analysis) — linear projection onto directions
//!   of maximum variance. Requires eigenvalue decomposition or SVD.
//!
//! - **t-SNE** (t-distributed Stochastic Neighbor Embedding) — non-linear
//!   embedding that preserves local structure. Good for visualization of
//!   high-dimensional data (e.g. single-cell RNA-seq).
//!
//! - **UMAP** (Uniform Manifold Approximation and Projection) — non-linear
//!   embedding that preserves both local and global structure. Faster than
//!   t-SNE for large datasets.
//!
//! These methods depend on linear algebra primitives (matrix multiplication,
//! eigendecomposition, SVD) that are not yet available in the workspace.
//! Implementation is deferred until a suitable linear algebra foundation is
//! in place.
