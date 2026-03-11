//! Spatial transcriptomics platform data structures — Visium, MERFISH, Slide-seq.
//!
//! Provides typed containers for the major spatial transcriptomics platforms,
//! handling platform-specific coordinate systems, resolution, and metadata.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Visium (10x Genomics)
// ---------------------------------------------------------------------------

/// A 10x Visium spatial transcriptomics dataset.
///
/// Visium uses a hexagonal array of capture spots (55 µm diameter, ~100 µm
/// center-to-center) on a glass slide. Each spot captures mRNA from the tissue
/// section above it.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct VisiumData {
    /// Gene names (features).
    pub genes: Vec<String>,
    /// Spot barcodes.
    pub barcodes: Vec<String>,
    /// Expression matrix: `counts[spot][gene]`.
    pub counts: Vec<Vec<f64>>,
    /// Spot spatial coordinates (pixel space): `(x, y)`.
    pub spot_coords: Vec<(f64, f64)>,
    /// Array row/col positions on the Visium grid.
    pub array_positions: Vec<(u32, u32)>,
    /// Whether each spot is under tissue.
    pub in_tissue: Vec<bool>,
    /// Scale factors for converting between coordinate systems.
    pub scale_factors: VisiumScaleFactors,
}

/// Visium scale factors for coordinate conversion.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct VisiumScaleFactors {
    /// Spot diameter in full-resolution pixels.
    pub spot_diameter_fullres: f64,
    /// Tissue hires image scale factor.
    pub tissue_hires_scalef: f64,
    /// Tissue lowres image scale factor.
    pub tissue_lowres_scalef: f64,
    /// Fiducial diameter in full-resolution pixels.
    pub fiducial_diameter_fullres: f64,
}

impl Default for VisiumScaleFactors {
    fn default() -> Self {
        Self {
            spot_diameter_fullres: 177.0,
            tissue_hires_scalef: 0.17,
            tissue_lowres_scalef: 0.034,
            fiducial_diameter_fullres: 328.0,
        }
    }
}

impl VisiumData {
    /// Create a new Visium dataset.
    pub fn new(
        genes: Vec<String>,
        barcodes: Vec<String>,
        counts: Vec<Vec<f64>>,
        spot_coords: Vec<(f64, f64)>,
        array_positions: Vec<(u32, u32)>,
        in_tissue: Vec<bool>,
    ) -> Result<Self> {
        let n_spots = barcodes.len();
        let n_genes = genes.len();
        if counts.len() != n_spots {
            return Err(CyaneaError::InvalidInput(format!(
                "counts rows ({}) != barcodes ({})",
                counts.len(),
                n_spots
            )));
        }
        for (i, row) in counts.iter().enumerate() {
            if row.len() != n_genes {
                return Err(CyaneaError::InvalidInput(format!(
                    "counts row {} has {} cols, expected {}",
                    i,
                    row.len(),
                    n_genes
                )));
            }
        }
        if spot_coords.len() != n_spots {
            return Err(CyaneaError::InvalidInput(
                "spot_coords length mismatch".into(),
            ));
        }
        if array_positions.len() != n_spots {
            return Err(CyaneaError::InvalidInput(
                "array_positions length mismatch".into(),
            ));
        }
        if in_tissue.len() != n_spots {
            return Err(CyaneaError::InvalidInput(
                "in_tissue length mismatch".into(),
            ));
        }
        Ok(Self {
            genes,
            barcodes,
            counts,
            spot_coords,
            array_positions,
            in_tissue,
            scale_factors: VisiumScaleFactors::default(),
        })
    }

    /// Number of spots.
    pub fn n_spots(&self) -> usize {
        self.barcodes.len()
    }

    /// Number of spots under tissue.
    pub fn n_tissue_spots(&self) -> usize {
        self.in_tissue.iter().filter(|&&t| t).count()
    }

    /// Number of genes.
    pub fn n_genes(&self) -> usize {
        self.genes.len()
    }

    /// Filter to only in-tissue spots.
    pub fn filter_tissue(&self) -> Self {
        let indices: Vec<usize> = self
            .in_tissue
            .iter()
            .enumerate()
            .filter(|(_, &t)| t)
            .map(|(i, _)| i)
            .collect();
        Self {
            genes: self.genes.clone(),
            barcodes: indices.iter().map(|&i| self.barcodes[i].clone()).collect(),
            counts: indices.iter().map(|&i| self.counts[i].clone()).collect(),
            spot_coords: indices.iter().map(|&i| self.spot_coords[i]).collect(),
            array_positions: indices.iter().map(|&i| self.array_positions[i]).collect(),
            in_tissue: vec![true; indices.len()],
            scale_factors: self.scale_factors.clone(),
        }
    }

    /// Total UMI counts per spot.
    pub fn total_counts_per_spot(&self) -> Vec<f64> {
        self.counts.iter().map(|row| row.iter().sum()).collect()
    }

    /// Genes detected (count > 0) per spot.
    pub fn genes_detected_per_spot(&self) -> Vec<usize> {
        self.counts
            .iter()
            .map(|row| row.iter().filter(|&&v| v > 0.0).count())
            .collect()
    }
}

// ---------------------------------------------------------------------------
// MERFISH
// ---------------------------------------------------------------------------

/// A MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization) dataset.
///
/// MERFISH provides subcellular-resolution spatial gene expression with
/// single-molecule sensitivity. Each molecule is localized to x, y (and
/// optionally z) coordinates and assigned to a cell via segmentation.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MerfishData {
    /// Gene names in the panel.
    pub genes: Vec<String>,
    /// Cell IDs (unique identifiers for segmented cells).
    pub cell_ids: Vec<String>,
    /// Cell centroids: `(x, y)`.
    pub cell_centroids: Vec<(f64, f64)>,
    /// Cell-by-gene count matrix: `counts[cell][gene]`.
    pub counts: Vec<Vec<f64>>,
    /// Optional: cell volumes (µm³) from 3D segmentation.
    pub cell_volumes: Option<Vec<f64>>,
    /// Optional: per-cell FOV (field of view) assignment.
    pub fov_ids: Option<Vec<u32>>,
    /// Blank control barcode counts per cell (for error rate estimation).
    pub blank_counts: Option<Vec<f64>>,
}

impl MerfishData {
    /// Create a new MERFISH dataset.
    pub fn new(
        genes: Vec<String>,
        cell_ids: Vec<String>,
        cell_centroids: Vec<(f64, f64)>,
        counts: Vec<Vec<f64>>,
    ) -> Result<Self> {
        let n_cells = cell_ids.len();
        let n_genes = genes.len();
        if counts.len() != n_cells {
            return Err(CyaneaError::InvalidInput(format!(
                "counts rows ({}) != cells ({})",
                counts.len(),
                n_cells
            )));
        }
        for (i, row) in counts.iter().enumerate() {
            if row.len() != n_genes {
                return Err(CyaneaError::InvalidInput(format!(
                    "counts row {} has {} cols, expected {}",
                    i,
                    row.len(),
                    n_genes
                )));
            }
        }
        if cell_centroids.len() != n_cells {
            return Err(CyaneaError::InvalidInput(
                "cell_centroids length mismatch".into(),
            ));
        }
        Ok(Self {
            genes,
            cell_ids,
            cell_centroids,
            counts,
            cell_volumes: None,
            fov_ids: None,
            blank_counts: None,
        })
    }

    /// Number of cells.
    pub fn n_cells(&self) -> usize {
        self.cell_ids.len()
    }

    /// Number of genes in the panel.
    pub fn n_genes(&self) -> usize {
        self.genes.len()
    }

    /// Estimate per-cell false positive rate from blank barcodes.
    pub fn estimated_fpr(&self) -> Option<Vec<f64>> {
        let blanks = self.blank_counts.as_ref()?;
        let n_genes = self.n_genes() as f64;
        Some(
            blanks
                .iter()
                .zip(self.counts.iter())
                .map(|(&blank, row)| {
                    let total: f64 = row.iter().sum();
                    if total > 0.0 {
                        (blank / n_genes) / (total / n_genes)
                    } else {
                        0.0
                    }
                })
                .collect(),
        )
    }

    /// Total transcript counts per cell.
    pub fn total_counts_per_cell(&self) -> Vec<f64> {
        self.counts.iter().map(|row| row.iter().sum()).collect()
    }
}

// ---------------------------------------------------------------------------
// Slide-seq
// ---------------------------------------------------------------------------

/// A Slide-seq / Slide-seqV2 spatial transcriptomics dataset.
///
/// Slide-seq uses DNA-barcoded beads (10 µm) randomly deposited on a glass
/// surface. Bead positions are decoded via sequencing-by-synthesis, giving
/// near-cellular resolution.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SlideseqData {
    /// Gene names.
    pub genes: Vec<String>,
    /// Bead barcodes.
    pub barcodes: Vec<String>,
    /// Bead spatial coordinates (µm): `(x, y)`.
    pub bead_coords: Vec<(f64, f64)>,
    /// Expression matrix: `counts[bead][gene]`.
    pub counts: Vec<Vec<f64>>,
}

impl SlideseqData {
    /// Create a new Slide-seq dataset.
    pub fn new(
        genes: Vec<String>,
        barcodes: Vec<String>,
        bead_coords: Vec<(f64, f64)>,
        counts: Vec<Vec<f64>>,
    ) -> Result<Self> {
        let n_beads = barcodes.len();
        let n_genes = genes.len();
        if counts.len() != n_beads {
            return Err(CyaneaError::InvalidInput(format!(
                "counts rows ({}) != beads ({})",
                counts.len(),
                n_beads
            )));
        }
        for (i, row) in counts.iter().enumerate() {
            if row.len() != n_genes {
                return Err(CyaneaError::InvalidInput(format!(
                    "counts row {} has {} cols, expected {}",
                    i,
                    row.len(),
                    n_genes
                )));
            }
        }
        if bead_coords.len() != n_beads {
            return Err(CyaneaError::InvalidInput(
                "bead_coords length mismatch".into(),
            ));
        }
        Ok(Self {
            genes,
            barcodes,
            bead_coords,
            counts,
        })
    }

    /// Number of beads.
    pub fn n_beads(&self) -> usize {
        self.barcodes.len()
    }

    /// Number of genes.
    pub fn n_genes(&self) -> usize {
        self.genes.len()
    }

    /// Total UMI counts per bead.
    pub fn total_counts_per_bead(&self) -> Vec<f64> {
        self.counts.iter().map(|row| row.iter().sum()).collect()
    }

    /// Filter beads by minimum total count.
    pub fn filter_min_counts(&self, min_counts: f64) -> Self {
        let totals = self.total_counts_per_bead();
        let indices: Vec<usize> = totals
            .iter()
            .enumerate()
            .filter(|(_, &t)| t >= min_counts)
            .map(|(i, _)| i)
            .collect();
        Self {
            genes: self.genes.clone(),
            barcodes: indices.iter().map(|&i| self.barcodes[i].clone()).collect(),
            bead_coords: indices.iter().map(|&i| self.bead_coords[i]).collect(),
            counts: indices.iter().map(|&i| self.counts[i].clone()).collect(),
        }
    }
}

// ---------------------------------------------------------------------------
// Generic spatial coordinates conversion
// ---------------------------------------------------------------------------

/// Convert any platform's coordinates to `SpatialPoint` for use with spatial
/// analysis functions (Moran's I, neighbor graphs, etc.).
pub fn visium_to_spatial_points(data: &VisiumData) -> Vec<super::spatial::SpatialPoint> {
    data.spot_coords
        .iter()
        .enumerate()
        .map(|(i, &(x, y))| super::spatial::SpatialPoint { x, y, index: i })
        .collect()
}

/// Convert MERFISH cell centroids to `SpatialPoint`.
pub fn merfish_to_spatial_points(data: &MerfishData) -> Vec<super::spatial::SpatialPoint> {
    data.cell_centroids
        .iter()
        .enumerate()
        .map(|(i, &(x, y))| super::spatial::SpatialPoint { x, y, index: i })
        .collect()
}

/// Convert Slide-seq bead coordinates to `SpatialPoint`.
pub fn slideseq_to_spatial_points(data: &SlideseqData) -> Vec<super::spatial::SpatialPoint> {
    data.bead_coords
        .iter()
        .enumerate()
        .map(|(i, &(x, y))| super::spatial::SpatialPoint { x, y, index: i })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_visium(n_spots: usize, n_genes: usize) -> VisiumData {
        let genes: Vec<String> = (0..n_genes).map(|i| format!("Gene{}", i)).collect();
        let barcodes: Vec<String> = (0..n_spots).map(|i| format!("BARCODE{}", i)).collect();
        let counts = vec![vec![1.0; n_genes]; n_spots];
        let spot_coords: Vec<(f64, f64)> = (0..n_spots).map(|i| (i as f64 * 100.0, 0.0)).collect();
        let array_positions: Vec<(u32, u32)> = (0..n_spots).map(|i| (0, i as u32)).collect();
        let in_tissue = vec![true; n_spots];
        VisiumData::new(genes, barcodes, counts, spot_coords, array_positions, in_tissue).unwrap()
    }

    #[test]
    fn test_visium_new() {
        let v = make_visium(10, 5);
        assert_eq!(v.n_spots(), 10);
        assert_eq!(v.n_genes(), 5);
        assert_eq!(v.n_tissue_spots(), 10);
    }

    #[test]
    fn test_visium_filter_tissue() {
        let mut v = make_visium(10, 5);
        v.in_tissue[0] = false;
        v.in_tissue[5] = false;
        let filtered = v.filter_tissue();
        assert_eq!(filtered.n_spots(), 8);
        assert!(filtered.in_tissue.iter().all(|&t| t));
    }

    #[test]
    fn test_visium_total_counts() {
        let v = make_visium(5, 10);
        let totals = v.total_counts_per_spot();
        assert_eq!(totals.len(), 5);
        assert!((totals[0] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_visium_genes_detected() {
        let mut v = make_visium(3, 5);
        v.counts[0] = vec![1.0, 0.0, 2.0, 0.0, 3.0]; // 3 detected
        let detected = v.genes_detected_per_spot();
        assert_eq!(detected[0], 3);
        assert_eq!(detected[1], 5); // all 1.0
    }

    #[test]
    fn test_visium_validation() {
        let result = VisiumData::new(
            vec!["G1".into()],
            vec!["B1".into(), "B2".into()],
            vec![vec![1.0]], // 1 row but 2 barcodes
            vec![(0.0, 0.0), (1.0, 0.0)],
            vec![(0, 0), (0, 1)],
            vec![true, true],
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_merfish_new() {
        let m = MerfishData::new(
            vec!["Slc17a7".into(), "Gad1".into()],
            vec!["cell_0".into(), "cell_1".into(), "cell_2".into()],
            vec![(0.0, 0.0), (10.0, 5.0), (20.0, 10.0)],
            vec![vec![5.0, 0.0], vec![0.0, 8.0], vec![3.0, 2.0]],
        )
        .unwrap();
        assert_eq!(m.n_cells(), 3);
        assert_eq!(m.n_genes(), 2);
    }

    #[test]
    fn test_merfish_fpr() {
        let mut m = MerfishData::new(
            vec!["G1".into(), "G2".into()],
            vec!["c0".into(), "c1".into()],
            vec![(0.0, 0.0), (1.0, 1.0)],
            vec![vec![10.0, 10.0], vec![5.0, 5.0]],
        )
        .unwrap();
        // No blanks → None
        assert!(m.estimated_fpr().is_none());
        m.blank_counts = Some(vec![0.5, 0.2]);
        let fpr = m.estimated_fpr().unwrap();
        assert_eq!(fpr.len(), 2);
        assert!(fpr[0] > 0.0);
    }

    #[test]
    fn test_merfish_total_counts() {
        let m = MerfishData::new(
            vec!["G1".into(), "G2".into(), "G3".into()],
            vec!["c0".into()],
            vec![(0.0, 0.0)],
            vec![vec![3.0, 4.0, 5.0]],
        )
        .unwrap();
        let totals = m.total_counts_per_cell();
        assert!((totals[0] - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_slideseq_new() {
        let s = SlideseqData::new(
            vec!["G1".into(), "G2".into()],
            vec!["bead_0".into(), "bead_1".into()],
            vec![(100.0, 200.0), (300.0, 400.0)],
            vec![vec![1.0, 2.0], vec![3.0, 4.0]],
        )
        .unwrap();
        assert_eq!(s.n_beads(), 2);
        assert_eq!(s.n_genes(), 2);
    }

    #[test]
    fn test_slideseq_filter_min_counts() {
        let s = SlideseqData::new(
            vec!["G1".into(), "G2".into()],
            vec!["b0".into(), "b1".into(), "b2".into()],
            vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)],
            vec![vec![1.0, 0.0], vec![10.0, 5.0], vec![0.5, 0.0]],
        )
        .unwrap();
        let filtered = s.filter_min_counts(5.0);
        assert_eq!(filtered.n_beads(), 1);
        assert_eq!(filtered.barcodes[0], "b1");
    }

    #[test]
    fn test_coordinate_conversions() {
        let v = make_visium(5, 2);
        let pts = visium_to_spatial_points(&v);
        assert_eq!(pts.len(), 5);
        assert!((pts[1].x - 100.0).abs() < 1e-10);

        let m = MerfishData::new(
            vec!["G1".into()],
            vec!["c0".into(), "c1".into()],
            vec![(5.0, 10.0), (15.0, 20.0)],
            vec![vec![1.0], vec![2.0]],
        )
        .unwrap();
        let mpts = merfish_to_spatial_points(&m);
        assert_eq!(mpts.len(), 2);
        assert!((mpts[0].y - 10.0).abs() < 1e-10);

        let s = SlideseqData::new(
            vec!["G1".into()],
            vec!["b0".into()],
            vec![(42.0, 99.0)],
            vec![vec![1.0]],
        )
        .unwrap();
        let spts = slideseq_to_spatial_points(&s);
        assert_eq!(spts.len(), 1);
        assert!((spts[0].x - 42.0).abs() < 1e-10);
    }
}
