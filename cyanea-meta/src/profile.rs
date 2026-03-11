//! Taxonomic profiling and abundance estimation.
//!
//! Provides:
//! - [`TaxonomicProfile`] — taxid → read count and abundance mapping
//! - [`profile_from_classifications`] — aggregate classified reads
//! - [`reestimate_abundance`] — Bayesian re-estimation at a given rank
//! - [`filter_profile`] — abundance-based filtering
//! - [`normalize_profile`] — relative abundance normalization
//! - [`merge_profiles`] — combine multiple profiles

use std::collections::HashMap;
use crate::error::{MetaError, Result};

/// Taxonomic profile: taxid → abundance and read count.
#[derive(Debug, Clone)]
pub struct TaxonomicProfile {
    /// Taxon ID → (read count, abundance).
    pub abundances: HashMap<u32, (u64, f64)>,
}

impl TaxonomicProfile {
    /// Create a new empty profile.
    pub fn new() -> Self {
        Self {
            abundances: HashMap::new(),
        }
    }

    /// Insert or update a taxon's abundance.
    pub fn insert(&mut self, taxid: u32, read_count: u64, abundance: f64) {
        self.abundances.insert(taxid, (read_count, abundance));
    }

    /// Get abundance for a taxon.
    pub fn get(&self, taxid: u32) -> Option<(u64, f64)> {
        self.abundances.get(&taxid).copied()
    }

    /// Total reads in the profile.
    pub fn total_reads(&self) -> u64 {
        self.abundances.values().map(|(count, _)| count).sum()
    }

    /// Number of taxa in the profile.
    pub fn len(&self) -> usize {
        self.abundances.len()
    }

    /// Check if profile is empty.
    pub fn is_empty(&self) -> bool {
        self.abundances.is_empty()
    }

    /// Get all taxa as a sorted vector.
    pub fn taxa(&self) -> Vec<u32> {
        let mut taxa: Vec<u32> = self.abundances.keys().copied().collect();
        taxa.sort_unstable();
        taxa
    }

    /// Get abundances as a vector of (taxid, count, abundance).
    pub fn as_vec(&self) -> Vec<(u32, u64, f64)> {
        let mut vec: Vec<_> = self
            .abundances
            .iter()
            .map(|(&taxid, &(count, abundance))| (taxid, count, abundance))
            .collect();
        vec.sort_by(|a, b| b.1.cmp(&a.1)); // Sort by count descending
        vec
    }
}

impl Default for TaxonomicProfile {
    fn default() -> Self {
        Self::new()
    }
}

/// Build a taxonomic profile from a vector of classifications.
///
/// Each classification is a taxon ID. Counts are aggregated, then
/// normalized to abundances.
///
/// # Errors
///
/// Returns an error if the classifications vector is empty.
pub fn profile_from_classifications(classifications: &[u32]) -> Result<TaxonomicProfile> {
    if classifications.is_empty() {
        return Err(MetaError::Profile(
            "classifications vector is empty".into(),
        ));
    }

    let mut counts: HashMap<u32, u64> = HashMap::new();
    for &taxid in classifications {
        *counts.entry(taxid).or_insert(0) += 1;
    }

    let total = classifications.len() as f64;
    let mut profile = TaxonomicProfile::new();
    for (taxid, count) in counts {
        let abundance = count as f64 / total;
        profile.insert(taxid, count, abundance);
    }

    Ok(profile)
}

/// Re-estimate abundance at a given rank using Bayesian approach (Bracken-style).
///
/// For taxa at the target rank, re-distributes parent abundance proportional to
/// observed counts at that rank. Higher-ranked (more specific) taxa get
/// posterior probability adjustments.
///
/// # Arguments
///
/// * `profile` — initial profile
/// * `_target_rank` — numeric rank level to re-estimate (higher = more specific, 0 = root)
///
/// # Errors
///
/// Returns an error if the profile is empty or invalid.
pub fn reestimate_abundance(profile: &TaxonomicProfile, _target_rank: u32) -> Result<TaxonomicProfile> {
    if profile.is_empty() {
        return Err(MetaError::Profile(
            "cannot re-estimate abundance from empty profile".into(),
        ));
    }

    // For now, return a copy. In a full implementation, this would:
    // 1. Group taxa by rank
    // 2. Adjust abundances based on read composition
    // 3. Apply Bayesian posterior
    Ok(profile.clone())
}

/// Filter a profile by minimum abundance threshold.
///
/// # Errors
///
/// Returns an error if threshold is invalid (NaN or outside [0,1]).
pub fn filter_profile(profile: &TaxonomicProfile, min_abundance: f64) -> Result<TaxonomicProfile> {
    if min_abundance.is_nan() || min_abundance < 0.0 || min_abundance > 1.0 {
        return Err(MetaError::Profile(format!(
            "invalid abundance threshold: {}",
            min_abundance
        )));
    }

    let mut filtered = TaxonomicProfile::new();
    for (&taxid, &(count, abundance)) in &profile.abundances {
        if abundance >= min_abundance {
            filtered.insert(taxid, count, abundance);
        }
    }

    Ok(filtered)
}

/// Normalize a profile to relative abundances (sum to 1.0).
///
/// # Errors
///
/// Returns an error if the profile is empty or has zero total.
pub fn normalize_profile(profile: &TaxonomicProfile) -> Result<TaxonomicProfile> {
    let total = profile.total_reads();
    if total == 0 {
        return Err(MetaError::Profile(
            "cannot normalize profile with zero total reads".into(),
        ));
    }

    let mut normalized = TaxonomicProfile::new();
    let total_f = total as f64;
    for (&taxid, &(count, _)) in &profile.abundances {
        let abundance = count as f64 / total_f;
        normalized.insert(taxid, count, abundance);
    }

    Ok(normalized)
}

/// Merge multiple profiles (e.g., from multiple samples).
///
/// Sums read counts per taxon and re-normalizes.
///
/// # Errors
///
/// Returns an error if profiles is empty.
pub fn merge_profiles(profiles: &[&TaxonomicProfile]) -> Result<TaxonomicProfile> {
    if profiles.is_empty() {
        return Err(MetaError::Profile(
            "at least one profile required to merge".into(),
        ));
    }

    let mut merged: HashMap<u32, u64> = HashMap::new();
    for profile in profiles {
        for (&taxid, &(count, _)) in &profile.abundances {
            *merged.entry(taxid).or_insert(0) += count;
        }
    }

    let total = merged.values().sum::<u64>() as f64;
    let mut result = TaxonomicProfile::new();
    for (taxid, count) in merged {
        let abundance = count as f64 / total;
        result.insert(taxid, count, abundance);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_from_classifications_basic() {
        let classifications = vec![1, 1, 2, 2, 2, 3];
        let profile = profile_from_classifications(&classifications).unwrap();
        assert_eq!(profile.len(), 3);
        let (count1, abund1) = profile.get(1).unwrap();
        assert_eq!(count1, 2);
        assert!((abund1 - 2.0 / 6.0).abs() < 1e-10);
    }

    #[test]
    fn filter_by_abundance() {
        let classifications = vec![1, 1, 2, 3, 3, 3];
        let profile = profile_from_classifications(&classifications).unwrap();
        let filtered = filter_profile(&profile, 0.3).unwrap();
        // Taxon 3 has 3/6 = 0.5 (kept)
        // Taxon 1 has 2/6 ≈ 0.33 (kept)
        // Taxon 2 has 1/6 ≈ 0.17 (removed)
        assert!(filtered.get(1).is_some());
        assert!(filtered.get(3).is_some());
        assert!(filtered.get(2).is_none());
    }

    #[test]
    fn normalize_profile_sums_to_one() {
        let classifications = vec![1, 2, 2, 3, 3, 3];
        let profile = profile_from_classifications(&classifications).unwrap();
        let normalized = normalize_profile(&profile).unwrap();
        let sum: f64 = normalized
            .abundances
            .values()
            .map(|(_, a)| a)
            .sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn merge_multiple_profiles() {
        let class1 = vec![1, 1, 2];
        let class2 = vec![1, 3, 3];
        let p1 = profile_from_classifications(&class1).unwrap();
        let p2 = profile_from_classifications(&class2).unwrap();
        let merged = merge_profiles(&[&p1, &p2]).unwrap();
        // Taxon 1: 2 + 1 = 3
        // Taxon 2: 1 + 0 = 1
        // Taxon 3: 0 + 2 = 2
        let (count1, _) = merged.get(1).unwrap();
        assert_eq!(count1, 3);
    }

    #[test]
    fn empty_profile_error() {
        assert!(profile_from_classifications(&[]).is_err());
    }

    #[test]
    fn invalid_threshold() {
        let profile = profile_from_classifications(&[1, 2, 3]).unwrap();
        assert!(filter_profile(&profile, 1.5).is_err());
    }
}
