//! OTU/ASV table operations for metagenomics.
//!
//! Create, filter, rarefy, collapse, and merge OTU/ASV abundance tables.

use cyanea_core::{CyaneaError, Result};

/// An OTU/ASV abundance table (samples × OTUs).
#[derive(Debug, Clone)]
pub struct OtuTable {
    /// Count matrix: `counts[sample][otu]`.
    pub counts: Vec<Vec<usize>>,
    /// Sample identifiers.
    pub sample_names: Vec<String>,
    /// OTU/ASV identifiers.
    pub otu_names: Vec<String>,
}

impl OtuTable {
    /// Create a new OTU table.
    ///
    /// # Errors
    ///
    /// Returns an error if dimensions are inconsistent (number of sample names
    /// doesn't match rows, or any row has the wrong number of OTUs).
    pub fn new(
        counts: Vec<Vec<usize>>,
        sample_names: Vec<String>,
        otu_names: Vec<String>,
    ) -> Result<Self> {
        if counts.len() != sample_names.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "number of samples ({}) does not match number of count rows ({})",
                sample_names.len(),
                counts.len()
            )));
        }
        for (i, row) in counts.iter().enumerate() {
            if row.len() != otu_names.len() {
                return Err(CyaneaError::InvalidInput(format!(
                    "row {} has {} OTUs, expected {}",
                    i,
                    row.len(),
                    otu_names.len()
                )));
            }
        }
        Ok(Self {
            counts,
            sample_names,
            otu_names,
        })
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.counts.len()
    }

    /// Number of OTUs/ASVs.
    pub fn n_otus(&self) -> usize {
        self.otu_names.len()
    }

    /// Total counts per sample.
    pub fn total_counts(&self) -> Vec<usize> {
        self.counts.iter().map(|row| row.iter().sum()).collect()
    }

    /// Relative abundance per sample (each row sums to 1.0).
    pub fn relative_abundance(&self) -> Vec<Vec<f64>> {
        self.counts
            .iter()
            .map(|row| {
                let total: usize = row.iter().sum();
                if total == 0 {
                    vec![0.0; row.len()]
                } else {
                    row.iter().map(|&c| c as f64 / total as f64).collect()
                }
            })
            .collect()
    }

    /// Keep only OTUs with total count ≥ `min_count` across all samples.
    pub fn filter_min_count(&self, min_count: usize) -> Self {
        let n_otus = self.n_otus();
        let mut keep = vec![false; n_otus];
        for j in 0..n_otus {
            let total: usize = self.counts.iter().map(|row| row[j]).sum();
            keep[j] = total >= min_count;
        }
        self.filter_otus(&keep)
    }

    /// Keep only OTUs present (count > 0) in at least `min_fraction` of samples.
    pub fn filter_min_prevalence(&self, min_fraction: f64) -> Self {
        let n_otus = self.n_otus();
        let n_samples = self.n_samples();
        let mut keep = vec![false; n_otus];
        for j in 0..n_otus {
            let n_present = self.counts.iter().filter(|row| row[j] > 0).count();
            keep[j] = n_present as f64 / n_samples as f64 >= min_fraction;
        }
        self.filter_otus(&keep)
    }

    /// Rarefy (subsample) each sample to exactly `depth` reads.
    ///
    /// Samples with fewer than `depth` total reads are excluded.
    /// Uses a deterministic approach based on proportional allocation
    /// with rounding, ensuring exactly `depth` total counts per sample.
    ///
    /// # Errors
    ///
    /// Returns an error if no samples have sufficient depth.
    pub fn rarefy(&self, depth: usize) -> Result<Self> {
        let mut new_counts = Vec::new();
        let mut new_names = Vec::new();

        for (i, row) in self.counts.iter().enumerate() {
            let total: usize = row.iter().sum();
            if total < depth {
                continue;
            }

            // Proportional allocation with controlled rounding.
            let mut rarefied = vec![0usize; row.len()];
            let mut allocated = 0usize;

            // First pass: floor allocation.
            let mut remainders: Vec<(f64, usize)> = Vec::new();
            for (j, &c) in row.iter().enumerate() {
                let proportion = c as f64 / total as f64;
                let exact = proportion * depth as f64;
                let floored = exact.floor() as usize;
                rarefied[j] = floored;
                allocated += floored;
                remainders.push((exact - floored as f64, j));
            }

            // Distribute remaining counts to highest-remainder OTUs.
            remainders.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
            let remaining = depth - allocated;
            for k in 0..remaining {
                if k < remainders.len() {
                    rarefied[remainders[k].1] += 1;
                }
            }

            new_counts.push(rarefied);
            new_names.push(self.sample_names[i].clone());
        }

        if new_counts.is_empty() {
            return Err(CyaneaError::InvalidInput(format!(
                "no samples have at least {} reads",
                depth
            )));
        }

        Ok(Self {
            counts: new_counts,
            sample_names: new_names,
            otu_names: self.otu_names.clone(),
        })
    }

    /// Collapse OTUs by taxonomy at the given level.
    ///
    /// `taxonomy[i]` is the taxonomic lineage for OTU `i`, split by level.
    /// OTUs sharing the same taxonomy prefix up to `level` are summed.
    ///
    /// # Errors
    ///
    /// Returns an error if `taxonomy` length doesn't match the number of OTUs,
    /// or any lineage is shorter than `level`.
    pub fn collapse_taxonomy(
        &self,
        level: usize,
        taxonomy: &[Vec<String>],
    ) -> Result<Self> {
        if taxonomy.len() != self.n_otus() {
            return Err(CyaneaError::InvalidInput(format!(
                "taxonomy length {} does not match OTU count {}",
                taxonomy.len(),
                self.n_otus()
            )));
        }
        for (i, lineage) in taxonomy.iter().enumerate() {
            if lineage.len() < level {
                return Err(CyaneaError::InvalidInput(format!(
                    "taxonomy for OTU {} has {} levels, need at least {}",
                    i,
                    lineage.len(),
                    level
                )));
            }
        }

        // Group OTUs by taxonomy prefix.
        let mut groups: Vec<(String, Vec<usize>)> = Vec::new();
        for (j, lineage) in taxonomy.iter().enumerate() {
            let prefix = lineage[..level].join(";");
            if let Some(group) = groups.iter_mut().find(|(k, _)| k == &prefix) {
                group.1.push(j);
            } else {
                groups.push((prefix, vec![j]));
            }
        }

        let new_otu_names: Vec<String> = groups.iter().map(|(k, _)| k.clone()).collect();
        let mut new_counts = vec![vec![0usize; groups.len()]; self.n_samples()];

        for (gi, (_, otu_indices)) in groups.iter().enumerate() {
            for (si, row) in self.counts.iter().enumerate() {
                for &oi in otu_indices {
                    new_counts[si][gi] += row[oi];
                }
            }
        }

        Ok(Self {
            counts: new_counts,
            sample_names: self.sample_names.clone(),
            otu_names: new_otu_names,
        })
    }

    /// Merge two OTU tables, combining samples.
    ///
    /// Both tables must share the same OTU names (in the same order).
    ///
    /// # Errors
    ///
    /// Returns an error if OTU names differ between the tables.
    pub fn merge(&self, other: &OtuTable) -> Result<Self> {
        if self.otu_names != other.otu_names {
            return Err(CyaneaError::InvalidInput(
                "cannot merge tables with different OTU names".into(),
            ));
        }
        let mut counts = self.counts.clone();
        counts.extend(other.counts.clone());

        let mut sample_names = self.sample_names.clone();
        sample_names.extend(other.sample_names.clone());

        Ok(Self {
            counts,
            sample_names,
            otu_names: self.otu_names.clone(),
        })
    }

    fn filter_otus(&self, keep: &[bool]) -> Self {
        let new_otu_names: Vec<String> = self
            .otu_names
            .iter()
            .enumerate()
            .filter(|(j, _)| keep[*j])
            .map(|(_, name)| name.clone())
            .collect();

        let new_counts: Vec<Vec<usize>> = self
            .counts
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(j, _)| keep[*j])
                    .map(|(_, &c)| c)
                    .collect()
            })
            .collect();

        Self {
            counts: new_counts,
            sample_names: self.sample_names.clone(),
            otu_names: new_otu_names,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_table() -> OtuTable {
        OtuTable::new(
            vec![
                vec![10, 5, 0, 1],
                vec![20, 0, 3, 2],
                vec![15, 10, 1, 0],
            ],
            vec!["S1".into(), "S2".into(), "S3".into()],
            vec!["OTU1".into(), "OTU2".into(), "OTU3".into(), "OTU4".into()],
        )
        .unwrap()
    }

    #[test]
    fn otu_table_creation() {
        let table = sample_table();
        assert_eq!(table.n_samples(), 3);
        assert_eq!(table.n_otus(), 4);
    }

    #[test]
    fn relative_abundance_sums_to_one() {
        let table = sample_table();
        let rel = table.relative_abundance();
        for row in &rel {
            let sum: f64 = row.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, "row sums to {}", sum);
        }
    }

    #[test]
    fn filter_min_count_removes() {
        let table = sample_table();
        // OTU3 total: 0+3+1 = 4, OTU4 total: 1+2+0 = 3.
        let filtered = table.filter_min_count(5);
        // OTU1 (45), OTU2 (15) should remain; OTU3 (4), OTU4 (3) filtered.
        assert_eq!(filtered.n_otus(), 2);
        assert!(filtered.otu_names.contains(&"OTU1".to_string()));
        assert!(filtered.otu_names.contains(&"OTU2".to_string()));
    }

    #[test]
    fn filter_prevalence_removes() {
        let table = sample_table();
        // OTU2 is present in S1, S3 (2/3 = 0.67); OTU3 in S2, S3 (2/3); OTU4 in S1, S2 (2/3)
        // OTU1 in all 3 (1.0). With threshold 0.8, only OTU1 remains.
        let filtered = table.filter_min_prevalence(0.8);
        assert_eq!(filtered.n_otus(), 1);
        assert_eq!(filtered.otu_names[0], "OTU1");
    }

    #[test]
    fn rarefy_correct_depth() {
        let table = sample_table();
        let depth = 10;
        let rarefied = table.rarefy(depth).unwrap();
        for row in &rarefied.counts {
            let total: usize = row.iter().sum();
            assert_eq!(total, depth, "rarefied sample total should be {}", depth);
        }
    }

    #[test]
    fn merge_combines_samples() {
        let t1 = OtuTable::new(
            vec![vec![10, 5]],
            vec!["S1".into()],
            vec!["OTU1".into(), "OTU2".into()],
        )
        .unwrap();
        let t2 = OtuTable::new(
            vec![vec![20, 15]],
            vec!["S2".into()],
            vec!["OTU1".into(), "OTU2".into()],
        )
        .unwrap();
        let merged = t1.merge(&t2).unwrap();
        assert_eq!(merged.n_samples(), 2);
        assert_eq!(merged.counts[0], vec![10, 5]);
        assert_eq!(merged.counts[1], vec![20, 15]);
    }
}
