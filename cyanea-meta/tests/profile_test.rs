/// Integration tests for profile module.

use cyanea_meta::profile::{
    profile_from_classifications, filter_profile, normalize_profile, merge_profiles,
};

#[test]
fn profile_from_classifications_basic() {
    let classifications = vec![1, 1, 2, 2, 2, 3];
    let profile = profile_from_classifications(&classifications).unwrap();

    assert_eq!(profile.len(), 3);
    assert_eq!(profile.total_reads(), 6);

    let (count1, abund1) = profile.get(1).unwrap();
    assert_eq!(count1, 2);
    assert!((abund1 - 2.0 / 6.0).abs() < 1e-10);

    let (count2, abund2) = profile.get(2).unwrap();
    assert_eq!(count2, 3);
    assert!((abund2 - 3.0 / 6.0).abs() < 1e-10);
}

#[test]
fn filter_by_abundance() {
    let classifications = vec![1, 1, 2, 3, 3, 3];
    let profile = profile_from_classifications(&classifications).unwrap();

    // Taxon 3 has 3/6 = 0.5 (should be kept with threshold 0.3)
    // Taxon 1 has 2/6 ≈ 0.33 (should be kept with threshold 0.3)
    // Taxon 2 has 1/6 ≈ 0.17 (should be removed with threshold 0.3)
    let filtered = filter_profile(&profile, 0.3).unwrap();

    assert!(filtered.get(1).is_some());
    assert!(filtered.get(3).is_some());
    assert!(filtered.get(2).is_none());
    assert_eq!(filtered.len(), 2);
}

#[test]
fn filter_strict_threshold() {
    let classifications = vec![1, 2, 2, 2, 3, 3];
    let profile = profile_from_classifications(&classifications).unwrap();
    // Taxon 1: 1/6 ≈ 0.167
    // Taxon 2: 3/6 = 0.5
    // Taxon 3: 2/6 ≈ 0.333

    let filtered = filter_profile(&profile, 0.4).unwrap();
    // Only taxon 2 should remain
    assert_eq!(filtered.len(), 1);
    assert!(filtered.get(2).is_some());
}

#[test]
fn normalize_profile_sums_to_one() {
    let classifications = vec![1, 2, 2, 3, 3, 3];
    let profile = profile_from_classifications(&classifications).unwrap();

    let normalized = normalize_profile(&profile).unwrap();
    let sum: f64 = normalized.abundances.values().map(|(_, a)| a).sum();
    assert!((sum - 1.0).abs() < 1e-10);
}

#[test]
fn merge_two_profiles() {
    let class1 = vec![1, 1, 2];
    let class2 = vec![1, 3, 3];

    let p1 = profile_from_classifications(&class1).unwrap();
    let p2 = profile_from_classifications(&class2).unwrap();

    let merged = merge_profiles(&[&p1, &p2]).unwrap();

    // Taxon 1: 2 + 1 = 3 out of 6 total
    // Taxon 2: 1 + 0 = 1 out of 6 total
    // Taxon 3: 0 + 2 = 2 out of 6 total
    let (count1, abund1) = merged.get(1).unwrap();
    assert_eq!(count1, 3);
    assert!((abund1 - 0.5).abs() < 1e-10);

    let (count3, abund3) = merged.get(3).unwrap();
    assert_eq!(count3, 2);
    assert!((abund3 - (2.0 / 6.0)).abs() < 1e-10);
}

#[test]
fn merge_multiple_profiles() {
    let p1_class = vec![1, 1];
    let p2_class = vec![2, 2, 2];
    let p3_class = vec![1, 3];

    let p1 = profile_from_classifications(&p1_class).unwrap();
    let p2 = profile_from_classifications(&p2_class).unwrap();
    let p3 = profile_from_classifications(&p3_class).unwrap();

    let merged = merge_profiles(&[&p1, &p2, &p3]).unwrap();

    // Total reads: 2 + 3 + 2 = 7
    // Taxon 1: 2 + 0 + 1 = 3
    // Taxon 2: 0 + 3 + 0 = 3
    // Taxon 3: 0 + 0 + 1 = 1
    let (count1, _) = merged.get(1).unwrap();
    let (count2, _) = merged.get(2).unwrap();
    let (count3, _) = merged.get(3).unwrap();

    assert_eq!(count1, 3);
    assert_eq!(count2, 3);
    assert_eq!(count3, 1);
}

#[test]
fn as_vec_sorted_by_count() {
    let classifications = vec![1, 2, 2, 3, 3, 3];
    let profile = profile_from_classifications(&classifications).unwrap();

    let vec = profile.as_vec();
    // Should be sorted by count descending: 3 (count 3), 2 (count 2), 1 (count 1)
    assert_eq!(vec[0].0, 3); // taxid
    assert_eq!(vec[0].1, 3); // count
    assert_eq!(vec[1].0, 2);
    assert_eq!(vec[1].1, 2);
    assert_eq!(vec[2].0, 1);
    assert_eq!(vec[2].1, 1);
}

#[test]
fn empty_classifications_error() {
    assert!(profile_from_classifications(&[]).is_err());
}

#[test]
fn invalid_filter_threshold() {
    let profile = profile_from_classifications(&[1, 2, 3]).unwrap();
    assert!(filter_profile(&profile, 1.5).is_err());
    assert!(filter_profile(&profile, -0.1).is_err());
}

#[test]
fn empty_merge_error() {
    assert!(merge_profiles(&[]).is_err());
}
