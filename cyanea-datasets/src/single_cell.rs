//! Single-cell demo dataset: small PBMC-like expression matrix.

/// Demo single-cell expression matrix: 50 cells × 20 genes.
///
/// Simulates a small PBMC dataset with 3 cell types:
/// - T cells (cells 0–19): high CD3D, CD3E, IL7R
/// - B cells (cells 20–34): high CD19, MS4A1, CD79A
/// - Monocytes (cells 35–49): high CD14, LYZ, FCGR3A
///
/// Returns (gene_names, cell_barcodes, count_matrix) where matrix is genes × cells.
pub fn demo_pbmc_50() -> DemoSingleCell {
    let genes = vec![
        "CD3D", "CD3E", "IL7R", "CD8A", "CD8B",    // T cell markers
        "CD19", "MS4A1", "CD79A", "CD79B", "PAX5",   // B cell markers
        "CD14", "LYZ", "FCGR3A", "S100A8", "S100A9", // Monocyte markers
        "ACTB", "GAPDH", "B2M", "MALAT1", "TMSB4X",  // Housekeeping
    ];

    let cells: Vec<String> = (0..50).map(|i| format!("CELL_{:03}", i)).collect();

    // Generate counts: marker genes high in their cell type, low elsewhere
    let mut matrix = vec![vec![0.0f64; 50]; 20];

    // Seed-based pseudo-random for reproducibility
    let mut rng = 42u64;
    let next = |state: &mut u64| -> f64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*state >> 33) as f64) / (u32::MAX as f64)
    };

    for cell in 0..50 {
        let cell_type = if cell < 20 { 0 } else if cell < 35 { 1 } else { 2 };

        for gene in 0..20 {
            let base = next(&mut rng) * 2.0; // background noise

            let signal = match (cell_type, gene) {
                // T cell markers high in T cells
                (0, 0..=4) => 8.0 + next(&mut rng) * 6.0,
                // B cell markers high in B cells
                (1, 5..=9) => 7.0 + next(&mut rng) * 5.0,
                // Monocyte markers high in monocytes
                (2, 10..=14) => 9.0 + next(&mut rng) * 7.0,
                // Housekeeping genes expressed in all
                (_, 15..=19) => 5.0 + next(&mut rng) * 3.0,
                _ => 0.0,
            };

            matrix[gene][cell] = (base + signal).round();
        }
    }

    DemoSingleCell {
        genes: genes.iter().map(|s| s.to_string()).collect(),
        cells,
        matrix,
        cell_types: (0..50)
            .map(|i| {
                if i < 20 {
                    "T_cell".to_string()
                } else if i < 35 {
                    "B_cell".to_string()
                } else {
                    "Monocyte".to_string()
                }
            })
            .collect(),
    }
}

/// A demo single-cell dataset.
#[derive(Debug, Clone)]
pub struct DemoSingleCell {
    /// Gene names (rows).
    pub genes: Vec<String>,
    /// Cell barcodes (columns).
    pub cells: Vec<String>,
    /// Count matrix (genes × cells).
    pub matrix: Vec<Vec<f64>>,
    /// Ground-truth cell type labels.
    pub cell_types: Vec<String>,
}

impl DemoSingleCell {
    /// Number of genes.
    pub fn num_genes(&self) -> usize {
        self.genes.len()
    }

    /// Number of cells.
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }

    /// Get expression vector for a gene.
    pub fn gene_expression(&self, gene_name: &str) -> Option<&Vec<f64>> {
        self.genes
            .iter()
            .position(|g| g == gene_name)
            .map(|i| &self.matrix[i])
    }

    /// Total UMI counts per cell.
    pub fn total_counts_per_cell(&self) -> Vec<f64> {
        (0..self.num_cells())
            .map(|c| self.matrix.iter().map(|row| row[c]).sum())
            .collect()
    }

    /// Number of detected genes per cell.
    pub fn genes_per_cell(&self) -> Vec<usize> {
        (0..self.num_cells())
            .map(|c| self.matrix.iter().filter(|row| row[c] > 0.0).count())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_demo_pbmc() {
        let data = demo_pbmc_50();
        assert_eq!(data.num_genes(), 20);
        assert_eq!(data.num_cells(), 50);
        assert_eq!(data.matrix.len(), 20);
        assert!(data.matrix.iter().all(|row| row.len() == 50));
    }

    #[test]
    fn test_cell_types() {
        let data = demo_pbmc_50();
        assert_eq!(data.cell_types.len(), 50);
        let t_cells = data.cell_types.iter().filter(|t| *t == "T_cell").count();
        let b_cells = data.cell_types.iter().filter(|t| *t == "B_cell").count();
        let mono = data.cell_types.iter().filter(|t| *t == "Monocyte").count();
        assert_eq!(t_cells, 20);
        assert_eq!(b_cells, 15);
        assert_eq!(mono, 15);
    }

    #[test]
    fn test_marker_enrichment() {
        let data = demo_pbmc_50();
        let cd3d = data.gene_expression("CD3D").unwrap();
        // T cells (0–19) should have higher CD3D than B cells (20–34)
        let t_mean: f64 = cd3d[0..20].iter().sum::<f64>() / 20.0;
        let b_mean: f64 = cd3d[20..35].iter().sum::<f64>() / 15.0;
        assert!(t_mean > b_mean);
    }

    #[test]
    fn test_total_counts() {
        let data = demo_pbmc_50();
        let counts = data.total_counts_per_cell();
        assert_eq!(counts.len(), 50);
        assert!(counts.iter().all(|&c| c > 0.0));
    }

    #[test]
    fn test_genes_per_cell() {
        let data = demo_pbmc_50();
        let detected = data.genes_per_cell();
        assert_eq!(detected.len(), 50);
        assert!(detected.iter().all(|&g| g > 5));
    }

    #[test]
    fn test_gene_expression_not_found() {
        let data = demo_pbmc_50();
        assert!(data.gene_expression("NONEXISTENT").is_none());
    }
}
