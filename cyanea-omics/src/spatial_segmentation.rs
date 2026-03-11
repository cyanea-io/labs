//! Cell segmentation for spatial transcriptomics — watershed, Voronoi, expansion.
//!
//! Provides cell segmentation algorithms for assigning transcripts or pixels
//! to cells based on spatial coordinates and marker signals.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A segmented cell with boundary and properties.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SegmentedCell {
    /// Unique cell ID.
    pub cell_id: usize,
    /// Centroid (x, y).
    pub centroid: (f64, f64),
    /// Area in coordinate units².
    pub area: f64,
    /// Boundary vertices as (x, y) pairs (convex hull approximation).
    pub boundary: Vec<(f64, f64)>,
    /// Number of transcripts assigned to this cell.
    pub n_transcripts: usize,
}

/// Parameters for nucleus expansion segmentation.
#[derive(Debug, Clone)]
pub struct ExpansionParams {
    /// Maximum expansion radius from nucleus centroid (µm).
    pub max_radius: f64,
    /// Minimum distance between cell boundaries (µm).
    pub min_gap: f64,
}

impl Default for ExpansionParams {
    fn default() -> Self {
        Self {
            max_radius: 15.0,
            min_gap: 1.0,
        }
    }
}

/// Result of a segmentation run.
#[derive(Debug, Clone)]
pub struct SegmentationResult {
    /// Segmented cells.
    pub cells: Vec<SegmentedCell>,
    /// Number of transcripts assigned to cells.
    pub assigned_transcripts: usize,
    /// Number of transcripts not assigned (background).
    pub unassigned_transcripts: usize,
}

// ---------------------------------------------------------------------------
// Voronoi-based segmentation
// ---------------------------------------------------------------------------

/// Segment transcripts into cells using Voronoi tessellation around seed points.
///
/// Each transcript is assigned to the nearest seed (nucleus centroid). Seeds
/// are provided as `(x, y)` coordinates. Transcripts are `(x, y)` positions.
///
/// An optional `max_radius` clips assignments beyond that distance (transcripts
/// too far from any seed are unassigned).
pub fn voronoi_segmentation(
    seeds: &[(f64, f64)],
    transcripts: &[(f64, f64)],
    max_radius: Option<f64>,
) -> Result<SegmentationResult> {
    if seeds.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "need at least 1 seed for segmentation".into(),
        ));
    }

    let n_seeds = seeds.len();
    let mut cell_transcripts: Vec<Vec<(f64, f64)>> = vec![Vec::new(); n_seeds];
    let mut unassigned = 0usize;

    for &(tx, ty) in transcripts {
        let mut best_idx = 0;
        let mut best_dist = f64::MAX;
        for (si, &(sx, sy)) in seeds.iter().enumerate() {
            let d = ((tx - sx).powi(2) + (ty - sy).powi(2)).sqrt();
            if d < best_dist {
                best_dist = d;
                best_idx = si;
            }
        }
        if let Some(max_r) = max_radius {
            if best_dist > max_r {
                unassigned += 1;
                continue;
            }
        }
        cell_transcripts[best_idx].push((tx, ty));
    }

    let mut cells = Vec::with_capacity(n_seeds);
    let assigned = transcripts.len() - unassigned;

    for (i, pts) in cell_transcripts.iter().enumerate() {
        let centroid = if pts.is_empty() {
            seeds[i]
        } else {
            let cx = pts.iter().map(|p| p.0).sum::<f64>() / pts.len() as f64;
            let cy = pts.iter().map(|p| p.1).sum::<f64>() / pts.len() as f64;
            (cx, cy)
        };

        let boundary = convex_hull(pts);
        let area = polygon_area(&boundary);

        cells.push(SegmentedCell {
            cell_id: i,
            centroid,
            area,
            boundary,
            n_transcripts: pts.len(),
        });
    }

    Ok(SegmentationResult {
        cells,
        assigned_transcripts: assigned,
        unassigned_transcripts: unassigned,
    })
}

// ---------------------------------------------------------------------------
// Nucleus expansion segmentation
// ---------------------------------------------------------------------------

/// Segment transcripts using nucleus expansion — grow circles from seed points
/// until they meet neighboring cells or reach `max_radius`.
///
/// This is a simplified version of the Baysor/Cellpose expansion approach.
/// Each transcript is assigned to the nearest seed within the effective radius,
/// where the effective radius is `min(max_radius, half_distance_to_nearest_seed)`.
pub fn expansion_segmentation(
    seeds: &[(f64, f64)],
    transcripts: &[(f64, f64)],
    params: &ExpansionParams,
) -> Result<SegmentationResult> {
    if seeds.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "need at least 1 seed for segmentation".into(),
        ));
    }

    let n_seeds = seeds.len();

    // Compute effective radius per seed: min(max_radius, half nearest-neighbor distance - min_gap/2)
    let mut effective_radii = vec![params.max_radius; n_seeds];
    for i in 0..n_seeds {
        let mut min_dist = f64::MAX;
        for j in 0..n_seeds {
            if i == j {
                continue;
            }
            let d = ((seeds[i].0 - seeds[j].0).powi(2) + (seeds[i].1 - seeds[j].1).powi(2)).sqrt();
            if d < min_dist {
                min_dist = d;
            }
        }
        let half_nn = (min_dist - params.min_gap) / 2.0;
        if half_nn > 0.0 && half_nn < effective_radii[i] {
            effective_radii[i] = half_nn;
        }
    }

    let mut cell_transcripts: Vec<Vec<(f64, f64)>> = vec![Vec::new(); n_seeds];
    let mut unassigned = 0usize;

    for &(tx, ty) in transcripts {
        let mut best_idx = None;
        let mut best_dist = f64::MAX;
        for (si, &(sx, sy)) in seeds.iter().enumerate() {
            let d = ((tx - sx).powi(2) + (ty - sy).powi(2)).sqrt();
            if d <= effective_radii[si] && d < best_dist {
                best_dist = d;
                best_idx = Some(si);
            }
        }
        match best_idx {
            Some(idx) => cell_transcripts[idx].push((tx, ty)),
            None => unassigned += 1,
        }
    }

    let assigned = transcripts.len() - unassigned;
    let mut cells = Vec::with_capacity(n_seeds);

    for (i, pts) in cell_transcripts.iter().enumerate() {
        let centroid = if pts.is_empty() {
            seeds[i]
        } else {
            let cx = pts.iter().map(|p| p.0).sum::<f64>() / pts.len() as f64;
            let cy = pts.iter().map(|p| p.1).sum::<f64>() / pts.len() as f64;
            (cx, cy)
        };

        let boundary = convex_hull(pts);
        let area = polygon_area(&boundary);

        cells.push(SegmentedCell {
            cell_id: i,
            centroid,
            area,
            boundary,
            n_transcripts: pts.len(),
        });
    }

    Ok(SegmentationResult {
        cells,
        assigned_transcripts: assigned,
        unassigned_transcripts: unassigned,
    })
}

// ---------------------------------------------------------------------------
// Watershed-style segmentation on a grid
// ---------------------------------------------------------------------------

/// Watershed segmentation on a 2D intensity grid.
///
/// `grid` is a row-major intensity image (e.g. DAPI). `seeds` are pixel
/// coordinates `(row, col)` of detected nuclei. Returns cell label for each
/// pixel (0 = background/boundary).
pub fn watershed_grid(
    grid: &[Vec<f64>],
    rows: usize,
    cols: usize,
    seeds: &[(usize, usize)],
) -> Result<Vec<Vec<usize>>> {
    if grid.len() != rows {
        return Err(CyaneaError::InvalidInput(
            "grid rows mismatch".into(),
        ));
    }
    for r in grid {
        if r.len() != cols {
            return Err(CyaneaError::InvalidInput(
                "grid cols mismatch".into(),
            ));
        }
    }
    if seeds.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "need at least 1 seed".into(),
        ));
    }

    // Label grid: 0 = unlabeled
    let mut labels = vec![vec![0usize; cols]; rows];

    // Priority queue: (priority inverted for min-heap via BTreeMap, row, col, label)
    // Use a simple BFS from seeds, processing in order of decreasing intensity.
    // Simplified: use a vec-based queue sorted by intensity.
    struct QueueItem {
        intensity: f64,
        row: usize,
        col: usize,
        label: usize,
    }

    let mut queue: Vec<QueueItem> = Vec::new();

    for (i, &(r, c)) in seeds.iter().enumerate() {
        let label = i + 1;
        if r < rows && c < cols {
            labels[r][c] = label;
            queue.push(QueueItem {
                intensity: grid[r][c],
                row: r,
                col: c,
                label,
            });
        }
    }

    // Sort by intensity descending (process high-intensity pixels first)
    queue.sort_by(|a, b| b.intensity.partial_cmp(&a.intensity).unwrap_or(std::cmp::Ordering::Equal));

    let mut head = 0;
    while head < queue.len() {
        let row = queue[head].row;
        let col = queue[head].col;
        let label = queue[head].label;
        head += 1;

        // 4-connected neighbors
        let neighbors: [(i64, i64); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];
        for (dr, dc) in &neighbors {
            let nr = row as i64 + dr;
            let nc = col as i64 + dc;
            if nr >= 0 && nr < rows as i64 && nc >= 0 && nc < cols as i64 {
                let nr = nr as usize;
                let nc = nc as usize;
                if labels[nr][nc] == 0 {
                    labels[nr][nc] = label;
                    queue.push(QueueItem {
                        intensity: grid[nr][nc],
                        row: nr,
                        col: nc,
                        label,
                    });
                }
            }
        }
    }

    Ok(labels)
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

/// Compute convex hull of 2D points (Andrew's monotone chain).
fn convex_hull(points: &[(f64, f64)]) -> Vec<(f64, f64)> {
    if points.len() < 3 {
        return points.to_vec();
    }

    let mut pts: Vec<(f64, f64)> = points.to_vec();
    pts.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });
    pts.dedup();

    if pts.len() < 3 {
        return pts;
    }

    fn cross(o: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    }

    let mut lower = Vec::new();
    for &p in &pts {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 0.0 {
            lower.pop();
        }
        lower.push(p);
    }

    let mut upper = Vec::new();
    for &p in pts.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 0.0 {
            upper.pop();
        }
        upper.push(p);
    }

    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// Area of a simple polygon via the shoelace formula.
fn polygon_area(vertices: &[(f64, f64)]) -> f64 {
    if vertices.len() < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    let n = vertices.len();
    for i in 0..n {
        let j = (i + 1) % n;
        area += vertices[i].0 * vertices[j].1;
        area -= vertices[j].0 * vertices[i].1;
    }
    (area / 2.0).abs()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_voronoi_basic() {
        let seeds = vec![(0.0, 0.0), (10.0, 0.0)];
        let transcripts = vec![
            (1.0, 0.0),
            (2.0, 1.0),
            (8.0, 0.0),
            (9.0, 1.0),
            (5.0, 0.0), // equidistant, goes to first
        ];
        let result = voronoi_segmentation(&seeds, &transcripts, None).unwrap();
        assert_eq!(result.cells.len(), 2);
        assert_eq!(result.assigned_transcripts, 5);
        assert_eq!(result.unassigned_transcripts, 0);
        // First cell should have transcripts near (0,0)
        assert!(result.cells[0].n_transcripts >= 2);
    }

    #[test]
    fn test_voronoi_max_radius() {
        let seeds = vec![(0.0, 0.0)];
        let transcripts = vec![(1.0, 0.0), (100.0, 0.0)];
        let result = voronoi_segmentation(&seeds, &transcripts, Some(5.0)).unwrap();
        assert_eq!(result.assigned_transcripts, 1);
        assert_eq!(result.unassigned_transcripts, 1);
    }

    #[test]
    fn test_expansion_basic() {
        let seeds = vec![(0.0, 0.0), (30.0, 0.0)];
        let params = ExpansionParams {
            max_radius: 15.0,
            min_gap: 1.0,
        };
        // Transcripts near seeds and far away
        let transcripts = vec![
            (2.0, 0.0),   // near seed 0
            (28.0, 0.0),  // near seed 1
            (100.0, 0.0), // far from both
        ];
        let result = expansion_segmentation(&seeds, &transcripts, &params).unwrap();
        assert_eq!(result.cells.len(), 2);
        assert_eq!(result.unassigned_transcripts, 1);
    }

    #[test]
    fn test_expansion_overlapping_seeds() {
        // Two seeds close together — effective radii should shrink
        let seeds = vec![(0.0, 0.0), (10.0, 0.0)];
        let params = ExpansionParams {
            max_radius: 20.0,
            min_gap: 2.0,
        };
        let transcripts = vec![(5.0, 0.0)]; // midpoint
        let result = expansion_segmentation(&seeds, &transcripts, &params).unwrap();
        // Effective radius = (10.0 - 2.0) / 2.0 = 4.0, midpoint at 5.0 is outside both
        assert_eq!(result.unassigned_transcripts, 1);
    }

    #[test]
    fn test_watershed_basic() {
        // 4x4 grid with two intensity peaks
        let grid = vec![
            vec![1.0, 2.0, 1.0, 1.0],
            vec![2.0, 5.0, 2.0, 1.0],
            vec![1.0, 2.0, 2.0, 5.0],
            vec![1.0, 1.0, 1.0, 2.0],
        ];
        let seeds = vec![(1, 1), (2, 3)]; // two nuclei
        let labels = watershed_grid(&grid, 4, 4, &seeds).unwrap();
        assert_eq!(labels.len(), 4);
        // Seed locations should have their labels
        assert_eq!(labels[1][1], 1);
        assert_eq!(labels[2][3], 2);
        // All pixels should be labeled
        for row in &labels {
            for &l in row {
                assert!(l > 0);
            }
        }
    }

    #[test]
    fn test_watershed_single_seed() {
        let grid = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let labels = watershed_grid(&grid, 2, 2, &[(0, 0)]).unwrap();
        // All pixels should be label 1
        assert!(labels.iter().all(|row| row.iter().all(|&l| l == 1)));
    }

    #[test]
    fn test_convex_hull() {
        let pts = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5)];
        let hull = convex_hull(&pts);
        assert_eq!(hull.len(), 4); // interior point excluded
    }

    #[test]
    fn test_polygon_area() {
        // Unit square
        let square = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let area = polygon_area(&square);
        assert!((area - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_segmentation_result_counts() {
        let seeds = vec![(0.0, 0.0), (50.0, 50.0)];
        let transcripts: Vec<(f64, f64)> = (0..100)
            .map(|i| {
                if i < 60 {
                    (i as f64 * 0.1, 0.0)
                } else {
                    (50.0 + (i - 60) as f64 * 0.1, 50.0)
                }
            })
            .collect();
        let result = voronoi_segmentation(&seeds, &transcripts, None).unwrap();
        assert_eq!(
            result.assigned_transcripts + result.unassigned_transcripts,
            100
        );
    }
}
