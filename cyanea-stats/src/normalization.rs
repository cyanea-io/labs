//! Count normalization for RNA-seq and related assays.
//!
//! All functions operate on row-major `&[f64]` slices with dimensions
//! `(n_genes, n_samples)`, matching `cyanea_omics::ExpressionMatrix` layout.
//!
//! - [`cpm`] — Counts per million
//! - [`tpm`] — Transcripts per million (requires gene lengths)
//! - [`fpkm`] — Fragments per kilobase of transcript per million mapped reads
//! - [`size_factors`] — DESeq2-style median-of-ratios normalization factors
//! - [`normalize_by_size_factors`] — Divide counts by per-sample size factors

use cyanea_core::{CyaneaError, Result};

use crate::descriptive;

// ── Helpers ──────────────────────────────────────────────────────────────────

fn validate_matrix(counts: &[f64], n_genes: usize, n_samples: usize) -> Result<()> {
    if n_genes == 0 || n_samples == 0 {
        return Err(CyaneaError::InvalidInput(
            "normalization: matrix must have at least 1 gene and 1 sample".into(),
        ));
    }
    if counts.len() != n_genes * n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "normalization: counts length ({}) != n_genes ({}) * n_samples ({})",
            counts.len(),
            n_genes,
            n_samples,
        )));
    }
    Ok(())
}

fn column_sums(counts: &[f64], n_genes: usize, n_samples: usize) -> Vec<f64> {
    let mut sums = vec![0.0; n_samples];
    for i in 0..n_genes {
        let row = &counts[i * n_samples..(i + 1) * n_samples];
        for (j, &v) in row.iter().enumerate() {
            sums[j] += v;
        }
    }
    sums
}

// ── CPM ──────────────────────────────────────────────────────────────────────

/// Counts per million: `CPM_ij = count_ij / library_size_j * 1e6`.
pub fn cpm(counts: &[f64], n_genes: usize, n_samples: usize) -> Result<Vec<f64>> {
    validate_matrix(counts, n_genes, n_samples)?;
    let lib_sizes = column_sums(counts, n_genes, n_samples);
    let mut out = vec![0.0; counts.len()];
    for i in 0..n_genes {
        for j in 0..n_samples {
            let idx = i * n_samples + j;
            out[idx] = if lib_sizes[j] > 0.0 {
                counts[idx] / lib_sizes[j] * 1e6
            } else {
                0.0
            };
        }
    }
    Ok(out)
}

// ── TPM ──────────────────────────────────────────────────────────────────────

/// Transcripts per million.
///
/// 1. Divide each count by gene length (→ reads per kilobase, RPK).
/// 2. Sum RPK per sample, then scale each sample's RPK values to sum to 1M.
///
/// `gene_lengths` must have `n_genes` elements (in bases or kilobases — units
/// are consistent because the per-kilobase step cancels in the ratio).
pub fn tpm(
    counts: &[f64],
    n_genes: usize,
    n_samples: usize,
    gene_lengths: &[f64],
) -> Result<Vec<f64>> {
    validate_matrix(counts, n_genes, n_samples)?;
    if gene_lengths.len() != n_genes {
        return Err(CyaneaError::InvalidInput(format!(
            "tpm: gene_lengths length ({}) != n_genes ({})",
            gene_lengths.len(),
            n_genes,
        )));
    }

    // RPK: count / (length / 1000)
    let mut rpk = vec![0.0; counts.len()];
    for i in 0..n_genes {
        let len_kb = gene_lengths[i] / 1000.0;
        if len_kb <= 0.0 {
            return Err(CyaneaError::InvalidInput(format!(
                "tpm: gene_lengths[{i}] must be positive",
            )));
        }
        for j in 0..n_samples {
            rpk[i * n_samples + j] = counts[i * n_samples + j] / len_kb;
        }
    }

    // Per-sample RPK sums
    let rpk_sums = column_sums(&rpk, n_genes, n_samples);

    // Scale to 1M
    let mut out = vec![0.0; counts.len()];
    for i in 0..n_genes {
        for j in 0..n_samples {
            let idx = i * n_samples + j;
            out[idx] = if rpk_sums[j] > 0.0 {
                rpk[idx] / rpk_sums[j] * 1e6
            } else {
                0.0
            };
        }
    }
    Ok(out)
}

// ── FPKM ─────────────────────────────────────────────────────────────────────

/// Fragments per kilobase of transcript per million mapped reads.
///
/// `FPKM_ij = count_ij * 1e9 / (library_size_j * length_i)`
pub fn fpkm(
    counts: &[f64],
    n_genes: usize,
    n_samples: usize,
    gene_lengths: &[f64],
) -> Result<Vec<f64>> {
    validate_matrix(counts, n_genes, n_samples)?;
    if gene_lengths.len() != n_genes {
        return Err(CyaneaError::InvalidInput(format!(
            "fpkm: gene_lengths length ({}) != n_genes ({})",
            gene_lengths.len(),
            n_genes,
        )));
    }
    let lib_sizes = column_sums(counts, n_genes, n_samples);
    let mut out = vec![0.0; counts.len()];
    for i in 0..n_genes {
        if gene_lengths[i] <= 0.0 {
            return Err(CyaneaError::InvalidInput(format!(
                "fpkm: gene_lengths[{i}] must be positive",
            )));
        }
        for j in 0..n_samples {
            let idx = i * n_samples + j;
            out[idx] = if lib_sizes[j] > 0.0 {
                counts[idx] * 1e9 / (lib_sizes[j] * gene_lengths[i])
            } else {
                0.0
            };
        }
    }
    Ok(out)
}

// ── Size factors (median-of-ratios) ──────────────────────────────────────────

/// DESeq2-style size factors via the median-of-ratios method (Anders & Huber 2010).
///
/// 1. Compute the geometric mean of each gene across samples.
/// 2. For each sample, compute the ratio `count / geometric_mean` for every gene.
/// 3. The size factor for a sample is the median of those ratios.
///
/// Genes with any zero count are excluded from the geometric mean calculation.
/// Returns one size factor per sample.
pub fn size_factors(counts: &[f64], n_genes: usize, n_samples: usize) -> Result<Vec<f64>> {
    validate_matrix(counts, n_genes, n_samples)?;

    // Compute geometric means, excluding genes with any zero.
    let mut geo_means = Vec::with_capacity(n_genes);
    let mut usable_genes = Vec::with_capacity(n_genes);

    for i in 0..n_genes {
        let row = &counts[i * n_samples..(i + 1) * n_samples];
        if row.iter().any(|&v| v <= 0.0) {
            continue;
        }
        let log_sum: f64 = row.iter().map(|v| v.ln()).sum();
        let geo_mean = (log_sum / n_samples as f64).exp();
        geo_means.push(geo_mean);
        usable_genes.push(i);
    }

    if usable_genes.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "size_factors: no genes with all non-zero counts".into(),
        ));
    }

    // For each sample, compute median of ratios.
    let mut factors = Vec::with_capacity(n_samples);
    for j in 0..n_samples {
        let ratios: Vec<f64> = usable_genes
            .iter()
            .zip(geo_means.iter())
            .map(|(&gene_i, &gm)| counts[gene_i * n_samples + j] / gm)
            .collect();
        let med = descriptive::median(&ratios)?;
        factors.push(med);
    }

    Ok(factors)
}

/// Divide each count by the corresponding sample's size factor.
pub fn normalize_by_size_factors(
    counts: &[f64],
    n_genes: usize,
    n_samples: usize,
    factors: &[f64],
) -> Result<Vec<f64>> {
    validate_matrix(counts, n_genes, n_samples)?;
    if factors.len() != n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "normalize_by_size_factors: factors length ({}) != n_samples ({})",
            factors.len(),
            n_samples,
        )));
    }
    let mut out = vec![0.0; counts.len()];
    for i in 0..n_genes {
        for j in 0..n_samples {
            let idx = i * n_samples + j;
            out[idx] = counts[idx] / factors[j];
        }
    }
    Ok(out)
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    #[test]
    fn cpm_column_sums_to_1m() {
        // 2 genes, 3 samples
        let counts = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0];
        let result = cpm(&counts, 2, 3).unwrap();
        // Column sums: [50, 70, 90]
        for j in 0..3 {
            let col_sum: f64 = (0..2).map(|i| result[i * 3 + j]).sum();
            assert!((col_sum - 1e6).abs() < 1.0, "col {j} sum={col_sum}");
        }
    }

    #[test]
    fn tpm_column_sums_to_1m() {
        let counts = [100.0, 200.0, 300.0, 400.0];
        let lengths = [1000.0, 2000.0];
        let result = tpm(&counts, 2, 2, &lengths).unwrap();
        for j in 0..2 {
            let col_sum: f64 = (0..2).map(|i| result[i * 2 + j]).sum();
            assert!((col_sum - 1e6).abs() < 1.0, "col {j} sum={col_sum}");
        }
    }

    #[test]
    fn tpm_length_normalization() {
        // Two genes with same count but different lengths:
        // shorter gene should have higher TPM
        let counts = [100.0, 100.0];
        let lengths = [500.0, 2000.0];
        let result = tpm(&counts, 2, 1, &lengths).unwrap();
        assert!(result[0] > result[1], "shorter gene should have higher TPM");
    }

    #[test]
    fn fpkm_known_values() {
        // 1 gene, 1 sample: count=100, lib_size=100, length=1000
        // FPKM = 100 * 1e9 / (100 * 1000) = 1_000_000
        let counts = [100.0];
        let lengths = [1000.0];
        let result = fpkm(&counts, 1, 1, &lengths).unwrap();
        assert!((result[0] - 1_000_000.0).abs() < TOL);
    }

    #[test]
    fn fpkm_to_tpm_relationship() {
        // TPM_i = FPKM_i / sum(FPKM_j) * 1e6
        let counts = [100.0, 200.0, 50.0, 300.0];
        let lengths = [1000.0, 2000.0];
        let fpkm_vals = fpkm(&counts, 2, 2, &lengths).unwrap();
        let tpm_vals = tpm(&counts, 2, 2, &lengths).unwrap();
        for j in 0..2 {
            let fpkm_sum: f64 = (0..2).map(|i| fpkm_vals[i * 2 + j]).sum();
            for i in 0..2 {
                let tpm_from_fpkm = fpkm_vals[i * 2 + j] / fpkm_sum * 1e6;
                assert!(
                    (tpm_from_fpkm - tpm_vals[i * 2 + j]).abs() < 1.0,
                    "gene {i} sample {j}: tpm_from_fpkm={tpm_from_fpkm}, tpm={}", tpm_vals[i * 2 + j]
                );
            }
        }
    }

    #[test]
    fn size_factors_equal_libraries() {
        // Two identical samples → size factors should be ~1
        let counts = [10.0, 10.0, 20.0, 20.0, 30.0, 30.0];
        let sf = size_factors(&counts, 3, 2).unwrap();
        assert!((sf[0] - 1.0).abs() < TOL);
        assert!((sf[1] - 1.0).abs() < TOL);
    }

    #[test]
    fn size_factors_doubled_library() {
        // Second sample has 2x counts → size factor ~2
        let counts = [10.0, 20.0, 20.0, 40.0, 30.0, 60.0];
        let sf = size_factors(&counts, 3, 2).unwrap();
        let ratio = sf[1] / sf[0];
        assert!((ratio - 2.0).abs() < TOL, "ratio={ratio}");
    }

    #[test]
    fn size_factors_skip_zeros() {
        // Gene 0 has a zero in sample 0 → excluded from calculation
        let counts = [0.0, 10.0, 20.0, 20.0, 30.0, 30.0];
        let sf = size_factors(&counts, 3, 2).unwrap();
        // Without gene 0: genes 1,2 are equal → factors ~1
        assert!((sf[0] - 1.0).abs() < TOL);
        assert!((sf[1] - 1.0).abs() < TOL);
    }

    #[test]
    fn normalize_roundtrip() {
        let counts = [10.0, 20.0, 30.0, 60.0];
        let sf = size_factors(&counts, 2, 2).unwrap();
        let normed = normalize_by_size_factors(&counts, 2, 2, &sf).unwrap();
        // After normalization, re-computing size factors should give ~1
        let sf2 = size_factors(&normed, 2, 2).unwrap();
        for &s in &sf2 {
            assert!((s - 1.0).abs() < 1e-4, "s={s}");
        }
    }

    #[test]
    fn dimension_mismatch() {
        assert!(cpm(&[1.0, 2.0], 3, 1).is_err());
        assert!(tpm(&[1.0, 2.0], 2, 1, &[100.0]).is_err()); // lengths wrong
        assert!(fpkm(&[1.0], 1, 1, &[100.0, 200.0]).is_err());
        assert!(normalize_by_size_factors(&[1.0, 2.0], 1, 2, &[1.0]).is_err());
    }

    #[test]
    fn empty_matrix() {
        assert!(cpm(&[], 0, 0).is_err());
    }
}
