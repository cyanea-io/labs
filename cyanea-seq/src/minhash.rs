//! MinHash and FracMinHash sketching for rapid genome comparison.
//!
//! Implements bottom-k MinHash and scaled FracMinHash for estimating Jaccard
//! similarity, containment, and average nucleotide identity (ANI) between
//! DNA sequences without full alignment.
//!
//! # Algorithms
//!
//! - **MinHash** — bottom-k sketch: keeps the `sketch_size` smallest hash values
//!   from the set of canonical k-mers. Jaccard similarity is estimated as the
//!   fraction of shared hashes in the union of the two sketches.
//!
//! - **FracMinHash** — scaled sketch: keeps all hash values below
//!   `u64::MAX / scale`. This allows unbiased containment estimation and works
//!   well for genomes of very different sizes.
//!
//! Both use canonical k-mers (min of forward and reverse complement hash) so
//! that a sequence and its reverse complement produce identical sketches.
//!
//! # Example
//!
//! ```
//! use cyanea_seq::minhash::{MinHash, FracMinHash};
//!
//! let seq_a = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
//! let seq_b = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
//!
//! let sketch_a = MinHash::from_sequence(seq_a, 7, 100).unwrap();
//! let sketch_b = MinHash::from_sequence(seq_b, 7, 100).unwrap();
//!
//! let jaccard = sketch_a.jaccard(&sketch_b).unwrap();
//! assert!((jaccard - 1.0).abs() < 1e-10);
//! ```

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// MurmurHash3-style hash function (inline, no external deps)
// ---------------------------------------------------------------------------

/// MurmurHash3 finalizer-based 64-bit hash.
fn murmurhash3_64(key: &[u8], seed: u64) -> u64 {
    let mut h = seed;
    for chunk in key.chunks(8) {
        let mut k = 0u64;
        for (i, &b) in chunk.iter().enumerate() {
            k |= (b as u64) << (i * 8);
        }
        k = k.wrapping_mul(0xff51afd7ed558ccd);
        k ^= k >> 33;
        k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
        h ^= k;
        h = h.wrapping_mul(5).wrapping_add(0xe6546b64);
    }
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    h
}

// ---------------------------------------------------------------------------
// DNA complement and canonical k-mer helpers
// ---------------------------------------------------------------------------

/// DNA complement for canonical k-mer computation (ACGT only).
#[inline]
fn dna_complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

/// Compute the reverse complement of a k-mer in place into a buffer.
fn reverse_complement(kmer: &[u8], buf: &mut Vec<u8>) {
    buf.clear();
    buf.extend(kmer.iter().rev().map(|&b| dna_complement(b)));
}

/// Normalize a byte to uppercase.
#[inline]
fn to_upper(b: u8) -> u8 {
    b.to_ascii_uppercase()
}

/// Check if a byte is a valid DNA base (ACGT, case-insensitive).
#[inline]
fn is_dna(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// Hash a k-mer canonically: hash both forward and reverse complement, return
/// the minimum. This ensures a sequence and its reverse complement produce the
/// same sketch.
fn canonical_kmer_hash(kmer: &[u8], rc_buf: &mut Vec<u8>, seed: u64) -> u64 {
    // Uppercase the k-mer for consistent hashing
    let upper: Vec<u8> = kmer.iter().map(|&b| to_upper(b)).collect();
    let h_fwd = murmurhash3_64(&upper, seed);
    reverse_complement(&upper, rc_buf);
    let h_rc = murmurhash3_64(rc_buf, seed);
    h_fwd.min(h_rc)
}

// ---------------------------------------------------------------------------
// MinHash — bottom-k sketch
// ---------------------------------------------------------------------------

/// Bottom-k MinHash sketch for rapid genome comparison.
///
/// Keeps the `sketch_size` smallest canonical k-mer hash values from a DNA
/// sequence. Jaccard similarity between two sketches is estimated from the
/// overlap of their bottom-k hash sets.
#[derive(Debug, Clone)]
pub struct MinHash {
    /// K-mer size.
    k: usize,
    /// Number of minimum hashes to keep.
    sketch_size: usize,
    /// Sorted bottom-k hash values (ascending).
    hashes: Vec<u64>,
}

const HASH_SEED: u64 = 42;

impl MinHash {
    /// Create a new empty MinHash sketch.
    ///
    /// # Errors
    ///
    /// Returns an error if `k == 0` or `sketch_size == 0`.
    pub fn new(k: usize, sketch_size: usize) -> Result<Self> {
        if k == 0 {
            return Err(CyaneaError::InvalidInput(
                "k-mer size must be at least 1".into(),
            ));
        }
        if sketch_size == 0 {
            return Err(CyaneaError::InvalidInput(
                "sketch size must be at least 1".into(),
            ));
        }
        Ok(Self {
            k,
            sketch_size,
            hashes: Vec::new(),
        })
    }

    /// Build a MinHash sketch from a DNA sequence.
    ///
    /// Non-ACGT bases act as k-mer break points (k-mers spanning them are
    /// skipped). Input is case-insensitive.
    ///
    /// # Errors
    ///
    /// Returns an error if `k == 0` or `sketch_size == 0`.
    pub fn from_sequence(seq: &[u8], k: usize, sketch_size: usize) -> Result<Self> {
        let mut mh = Self::new(k, sketch_size)?;
        mh.add_sequence(seq);
        Ok(mh)
    }

    /// Add k-mers from a DNA sequence to this sketch.
    ///
    /// Non-ACGT bases act as k-mer break points. The sketch is maintained as a
    /// sorted bottom-k set.
    pub fn add_sequence(&mut self, seq: &[u8]) {
        if seq.len() < self.k {
            return;
        }

        let mut rc_buf = Vec::with_capacity(self.k);

        for window in seq.windows(self.k) {
            // Skip windows containing non-DNA bases
            if !window.iter().all(|&b| is_dna(b)) {
                continue;
            }

            let h = canonical_kmer_hash(window, &mut rc_buf, HASH_SEED);
            self.insert_hash(h);
        }
    }

    /// Insert a hash value into the bottom-k set.
    fn insert_hash(&mut self, h: u64) {
        // If we haven't filled the sketch yet, insert and maintain sort
        if self.hashes.len() < self.sketch_size {
            let pos = self.hashes.binary_search(&h).unwrap_or_else(|p| p);
            // Skip duplicates
            if pos < self.hashes.len() && self.hashes[pos] == h {
                return;
            }
            self.hashes.insert(pos, h);
        } else if h < self.hashes[self.hashes.len() - 1] {
            // h is smaller than the current max in our bottom-k set
            let pos = self.hashes.binary_search(&h).unwrap_or_else(|p| p);
            // Skip duplicates
            if pos < self.hashes.len() && self.hashes[pos] == h {
                return;
            }
            self.hashes.insert(pos, h);
            self.hashes.truncate(self.sketch_size);
        }
        // else: h >= current max and sketch is full, discard
    }

    /// Estimate Jaccard similarity between this sketch and another.
    ///
    /// Uses the merge-based estimator: merge both sorted hash arrays, count
    /// how many of the bottom-k values from the union appear in both sketches.
    ///
    /// # Errors
    ///
    /// Returns an error if the sketches have different `k` values.
    pub fn jaccard(&self, other: &MinHash) -> Result<f64> {
        if self.k != other.k {
            return Err(CyaneaError::InvalidInput(format!(
                "incompatible k-mer sizes: {} vs {}",
                self.k, other.k
            )));
        }
        if self.hashes.is_empty() && other.hashes.is_empty() {
            return Ok(1.0);
        }
        if self.hashes.is_empty() || other.hashes.is_empty() {
            return Ok(0.0);
        }

        // Merge the two sorted hash lists, take the smallest `sketch_size`
        // from the union, count how many appear in both.
        let max_size = self.sketch_size.max(other.sketch_size);
        let mut i = 0;
        let mut j = 0;
        let mut union_count = 0;
        let mut intersection_count = 0;

        while union_count < max_size && i < self.hashes.len() && j < other.hashes.len() {
            let a = self.hashes[i];
            let b = other.hashes[j];
            if a < b {
                i += 1;
            } else if a > b {
                j += 1;
            } else {
                // Equal
                intersection_count += 1;
                i += 1;
                j += 1;
            }
            union_count += 1;
        }

        // Fill remaining union slots from whichever list still has elements
        while union_count < max_size && i < self.hashes.len() {
            i += 1;
            union_count += 1;
        }
        while union_count < max_size && j < other.hashes.len() {
            j += 1;
            union_count += 1;
        }

        if union_count == 0 {
            return Ok(0.0);
        }

        Ok(intersection_count as f64 / union_count as f64)
    }

    /// Estimate containment of `self` in `other`.
    ///
    /// Containment C(A, B) = |A intersect B| / |A|.
    ///
    /// # Errors
    ///
    /// Returns an error if the sketches have different `k` values.
    pub fn containment(&self, other: &MinHash) -> Result<f64> {
        if self.k != other.k {
            return Err(CyaneaError::InvalidInput(format!(
                "incompatible k-mer sizes: {} vs {}",
                self.k, other.k
            )));
        }
        if self.hashes.is_empty() {
            return Ok(1.0);
        }

        let mut shared = 0usize;
        let mut j = 0;
        for &h in &self.hashes {
            while j < other.hashes.len() && other.hashes[j] < h {
                j += 1;
            }
            if j < other.hashes.len() && other.hashes[j] == h {
                shared += 1;
            }
        }

        Ok(shared as f64 / self.hashes.len() as f64)
    }

    /// Estimate average nucleotide identity (ANI) from Jaccard similarity.
    ///
    /// Uses the Mash formula: `ANI = 1 + (2/k) * ln(2J / (1 + J))`
    ///
    /// # Errors
    ///
    /// Returns an error if the sketches have different `k` values, or if the
    /// Jaccard similarity is zero (ANI undefined).
    pub fn ani(&self, other: &MinHash) -> Result<f64> {
        let j = self.jaccard(other)?;
        ani_from_jaccard(j, self.k)
    }

    /// Number of hash values currently in the sketch.
    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    /// Whether the sketch is empty (no k-mers hashed).
    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }

    /// The k-mer size used by this sketch.
    pub fn k(&self) -> usize {
        self.k
    }

    /// The target sketch size (bottom-k parameter).
    pub fn sketch_size(&self) -> usize {
        self.sketch_size
    }

    /// Return a reference to the sorted hash values.
    pub fn hashes(&self) -> &[u64] {
        &self.hashes
    }
}

// ---------------------------------------------------------------------------
// FracMinHash — scaled sketch
// ---------------------------------------------------------------------------

/// Scaled (FracMinHash) sketch for rapid genome comparison.
///
/// Keeps all canonical k-mer hash values below `u64::MAX / scale`. This gives
/// a predictable sampling fraction regardless of genome size, making containment
/// estimation more accurate for genomes of different sizes.
#[derive(Debug, Clone)]
pub struct FracMinHash {
    /// K-mer size.
    k: usize,
    /// Scale factor: keep hashes below `u64::MAX / scale`.
    scale: u64,
    /// Maximum hash threshold: `u64::MAX / scale`.
    max_hash: u64,
    /// Sorted hash values below threshold (ascending).
    hashes: Vec<u64>,
}

impl FracMinHash {
    /// Create a new empty FracMinHash sketch.
    ///
    /// # Errors
    ///
    /// Returns an error if `k == 0` or `scale == 0`.
    pub fn new(k: usize, scale: u64) -> Result<Self> {
        if k == 0 {
            return Err(CyaneaError::InvalidInput(
                "k-mer size must be at least 1".into(),
            ));
        }
        if scale == 0 {
            return Err(CyaneaError::InvalidInput(
                "scale must be at least 1".into(),
            ));
        }
        let max_hash = u64::MAX / scale;
        Ok(Self {
            k,
            scale,
            max_hash,
            hashes: Vec::new(),
        })
    }

    /// Build a FracMinHash sketch from a DNA sequence.
    ///
    /// # Errors
    ///
    /// Returns an error if `k == 0` or `scale == 0`.
    pub fn from_sequence(seq: &[u8], k: usize, scale: u64) -> Result<Self> {
        let mut fmh = Self::new(k, scale)?;
        fmh.add_sequence(seq);
        Ok(fmh)
    }

    /// Add k-mers from a DNA sequence to this sketch.
    ///
    /// Non-ACGT bases act as k-mer break points.
    pub fn add_sequence(&mut self, seq: &[u8]) {
        if seq.len() < self.k {
            return;
        }

        let mut rc_buf = Vec::with_capacity(self.k);

        for window in seq.windows(self.k) {
            if !window.iter().all(|&b| is_dna(b)) {
                continue;
            }

            let h = canonical_kmer_hash(window, &mut rc_buf, HASH_SEED);
            if h <= self.max_hash {
                // Insert maintaining sorted order, skip duplicates
                let pos = self.hashes.binary_search(&h).unwrap_or_else(|p| p);
                if pos >= self.hashes.len() || self.hashes[pos] != h {
                    self.hashes.insert(pos, h);
                }
            }
        }
    }

    /// Estimate Jaccard similarity between this sketch and another.
    ///
    /// Jaccard = |A intersect B| / |A union B|, computed exactly over the
    /// hashes below the shared threshold.
    ///
    /// # Errors
    ///
    /// Returns an error if the sketches have different `k` or `scale` values.
    pub fn jaccard(&self, other: &FracMinHash) -> Result<f64> {
        self.check_compatible(other)?;

        if self.hashes.is_empty() && other.hashes.is_empty() {
            return Ok(1.0);
        }
        if self.hashes.is_empty() || other.hashes.is_empty() {
            return Ok(0.0);
        }

        let (intersection, union) = self.intersection_union(other);

        if union == 0 {
            return Ok(0.0);
        }

        Ok(intersection as f64 / union as f64)
    }

    /// Estimate containment of `self` in `other`.
    ///
    /// Containment C(A, B) = |A intersect B| / |A|.
    ///
    /// # Errors
    ///
    /// Returns an error if the sketches have different `k` or `scale` values.
    pub fn containment(&self, other: &FracMinHash) -> Result<f64> {
        self.check_compatible(other)?;

        if self.hashes.is_empty() {
            return Ok(1.0);
        }

        let mut shared = 0usize;
        let mut j = 0;
        for &h in &self.hashes {
            while j < other.hashes.len() && other.hashes[j] < h {
                j += 1;
            }
            if j < other.hashes.len() && other.hashes[j] == h {
                shared += 1;
            }
        }

        Ok(shared as f64 / self.hashes.len() as f64)
    }

    /// Estimate average nucleotide identity (ANI) from Jaccard similarity.
    ///
    /// # Errors
    ///
    /// Returns an error if the sketches are incompatible or Jaccard is zero.
    pub fn ani(&self, other: &FracMinHash) -> Result<f64> {
        let j = self.jaccard(other)?;
        ani_from_jaccard(j, self.k)
    }

    /// Number of hash values in the sketch.
    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    /// Whether the sketch is empty.
    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }

    /// The k-mer size used by this sketch.
    pub fn k(&self) -> usize {
        self.k
    }

    /// The scale factor.
    pub fn scale(&self) -> u64 {
        self.scale
    }

    /// Return a reference to the sorted hash values.
    pub fn hashes(&self) -> &[u64] {
        &self.hashes
    }

    /// Check that two sketches have compatible parameters.
    fn check_compatible(&self, other: &FracMinHash) -> Result<()> {
        if self.k != other.k {
            return Err(CyaneaError::InvalidInput(format!(
                "incompatible k-mer sizes: {} vs {}",
                self.k, other.k
            )));
        }
        if self.scale != other.scale {
            return Err(CyaneaError::InvalidInput(format!(
                "incompatible scale values: {} vs {}",
                self.scale, other.scale
            )));
        }
        Ok(())
    }

    /// Count intersection and union sizes via merge of sorted hash lists.
    fn intersection_union(&self, other: &FracMinHash) -> (usize, usize) {
        let mut i = 0;
        let mut j = 0;
        let mut intersection = 0;
        let mut union = 0;

        while i < self.hashes.len() && j < other.hashes.len() {
            let a = self.hashes[i];
            let b = other.hashes[j];
            if a < b {
                i += 1;
            } else if a > b {
                j += 1;
            } else {
                intersection += 1;
                i += 1;
                j += 1;
            }
            union += 1;
        }

        // Count remaining elements
        union += (self.hashes.len() - i) + (other.hashes.len() - j);

        (intersection, union)
    }
}

// ---------------------------------------------------------------------------
// Shared ANI estimation
// ---------------------------------------------------------------------------

/// Estimate ANI from Jaccard similarity and k-mer size.
///
/// Uses the Mash formula: `ANI = 1 + (2/k) * ln(2J / (1 + J))`
///
/// Returns an error if Jaccard is zero (ANI undefined due to log(0)).
fn ani_from_jaccard(j: f64, k: usize) -> Result<f64> {
    if j <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "cannot estimate ANI when Jaccard similarity is zero".into(),
        ));
    }
    if (j - 1.0).abs() < f64::EPSILON {
        return Ok(1.0);
    }

    let mash_distance = -((2.0 * j) / (1.0 + j)).ln() / (k as f64);
    let ani = 1.0 - mash_distance;

    // Clamp to [0, 1] — numerical edge cases can push slightly outside
    Ok(ani.clamp(0.0, 1.0))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Parameter validation ---

    #[test]
    fn minhash_k_zero_error() {
        assert!(MinHash::new(0, 100).is_err());
    }

    #[test]
    fn minhash_sketch_size_zero_error() {
        assert!(MinHash::new(7, 0).is_err());
    }

    #[test]
    fn fracminhash_k_zero_error() {
        assert!(FracMinHash::new(0, 10).is_err());
    }

    #[test]
    fn fracminhash_scale_zero_error() {
        assert!(FracMinHash::new(7, 0).is_err());
    }

    // --- MinHash: identical sequences ---

    #[test]
    fn minhash_identical_jaccard_one() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let a = MinHash::from_sequence(seq, 7, 100).unwrap();
        let b = MinHash::from_sequence(seq, 7, 100).unwrap();
        let j = a.jaccard(&b).unwrap();
        assert!(
            (j - 1.0).abs() < 1e-10,
            "expected Jaccard ~1.0, got {}",
            j
        );
    }

    // --- MinHash: disjoint sequences ---

    #[test]
    fn minhash_disjoint_jaccard_near_zero() {
        // Two sequences with no shared k-mers (different base composition)
        let seq_a = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let seq_b = b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        let a = MinHash::from_sequence(seq_a, 7, 100).unwrap();
        let b = MinHash::from_sequence(seq_b, 7, 100).unwrap();
        let j = a.jaccard(&b).unwrap();
        assert!(
            j < 0.01,
            "expected Jaccard ~0.0 for disjoint sequences, got {}",
            j
        );
    }

    // --- MinHash: sketch size correct ---

    #[test]
    fn minhash_sketch_size_capped() {
        // A long enough sequence to produce many distinct k-mers
        let seq = b"ACGTACGTAAACCCGGGTTTACGTACGTAAACCCGGGTTTACGTACGT";
        let mh = MinHash::from_sequence(seq, 5, 10).unwrap();
        assert!(
            mh.len() <= 10,
            "sketch should have at most 10 hashes, got {}",
            mh.len()
        );
        // With 43 5-mers from this sequence, we should fill the sketch
        assert!(
            mh.len() > 0,
            "sketch should not be empty for a non-trivial sequence"
        );
    }

    // --- MinHash: containment of subset ---

    #[test]
    fn minhash_containment_subset() {
        // seq_a is a prefix of seq_b — all k-mers of seq_a appear in seq_b
        let seq_a = b"ACGTACGTACGTACGT";
        let seq_b = b"ACGTACGTACGTACGTTTTTCCCCAAAAGGGGG";
        let a = MinHash::from_sequence(seq_a, 5, 1000).unwrap();
        let b = MinHash::from_sequence(seq_b, 5, 1000).unwrap();
        let c = a.containment(&b).unwrap();
        assert!(
            c > 0.9,
            "expected containment ~1.0 for subset, got {}",
            c
        );
    }

    // --- MinHash: ANI of identical sequences ---

    #[test]
    fn minhash_ani_identical_is_one() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let a = MinHash::from_sequence(seq, 7, 100).unwrap();
        let b = MinHash::from_sequence(seq, 7, 100).unwrap();
        let ani = a.ani(&b).unwrap();
        assert!(
            (ani - 1.0).abs() < 1e-10,
            "expected ANI ~1.0, got {}",
            ani
        );
    }

    // --- MinHash: incompatible k error ---

    #[test]
    fn minhash_incompatible_k_error() {
        let a = MinHash::from_sequence(b"ACGTACGT", 3, 10).unwrap();
        let b = MinHash::from_sequence(b"ACGTACGT", 5, 10).unwrap();
        assert!(a.jaccard(&b).is_err());
    }

    // --- Canonical k-mers: sequence and reverse complement produce same sketch ---

    #[test]
    fn canonical_kmers_same_sketch() {
        let fwd = b"ACGTACGTACGTACGTACGTACGT";
        // Manually compute reverse complement
        let rc: Vec<u8> = fwd.iter().rev().map(|&b| dna_complement(b)).collect();

        let a = MinHash::from_sequence(fwd, 7, 100).unwrap();
        let b = MinHash::from_sequence(&rc, 7, 100).unwrap();

        assert_eq!(
            a.hashes(),
            b.hashes(),
            "forward and reverse complement should produce identical sketches"
        );
    }

    // --- FracMinHash: scaling behavior ---

    #[test]
    fn fracminhash_scales_properly() {
        let seq = b"ACGTACGTAAACCCGGGTTTACGTACGTAAACCCGGGTTTACGTACGT";
        // With scale=1, keep all hashes (max_hash = u64::MAX)
        let all = FracMinHash::from_sequence(seq, 5, 1).unwrap();
        // With scale=10, keep ~1/10 of hashes
        let tenth = FracMinHash::from_sequence(seq, 5, 10).unwrap();

        assert!(
            tenth.len() <= all.len(),
            "scale=10 should keep fewer hashes than scale=1: {} vs {}",
            tenth.len(),
            all.len()
        );
    }

    // --- FracMinHash: identical sequences ---

    #[test]
    fn fracminhash_identical_jaccard_one() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let a = FracMinHash::from_sequence(seq, 7, 2).unwrap();
        let b = FracMinHash::from_sequence(seq, 7, 2).unwrap();
        let j = a.jaccard(&b).unwrap();
        assert!(
            (j - 1.0).abs() < 1e-10,
            "expected Jaccard ~1.0, got {}",
            j
        );
    }

    // --- FracMinHash: containment ---

    #[test]
    fn fracminhash_containment_subset() {
        let seq_a = b"ACGTACGTACGTACGT";
        let seq_b = b"ACGTACGTACGTACGTTTTTCCCCAAAAGGGGG";
        let a = FracMinHash::from_sequence(seq_a, 5, 1).unwrap();
        let b = FracMinHash::from_sequence(seq_b, 5, 1).unwrap();
        let c = a.containment(&b).unwrap();
        assert!(
            c > 0.9,
            "expected containment ~1.0 for subset, got {}",
            c
        );
    }

    // --- FracMinHash: incompatible parameters error ---

    #[test]
    fn fracminhash_incompatible_scale_error() {
        let a = FracMinHash::from_sequence(b"ACGTACGT", 5, 10).unwrap();
        let b = FracMinHash::from_sequence(b"ACGTACGT", 5, 20).unwrap();
        assert!(a.jaccard(&b).is_err());
    }

    // --- Round-trip: Jaccard in reasonable range ---

    #[test]
    fn jaccard_in_valid_range() {
        let seq_a = b"ACGTACGTACGTACGTAAAA";
        let seq_b = b"ACGTACGTACGTACGTCCCC";
        let a = MinHash::from_sequence(seq_a, 5, 100).unwrap();
        let b = MinHash::from_sequence(seq_b, 5, 100).unwrap();
        let j = a.jaccard(&b).unwrap();
        assert!(
            (0.0..=1.0).contains(&j),
            "Jaccard should be in [0,1], got {}",
            j
        );
        // These sequences share a prefix, so Jaccard should be non-zero
        assert!(j > 0.0, "sequences share k-mers, Jaccard should be > 0");
        assert!(j < 1.0, "sequences differ, Jaccard should be < 1");
    }

    // --- Empty sketch behavior ---

    #[test]
    fn minhash_empty_sketch() {
        let mh = MinHash::new(7, 100).unwrap();
        assert!(mh.is_empty());
        assert_eq!(mh.len(), 0);
    }

    // --- Hash function determinism ---

    #[test]
    fn murmurhash3_deterministic() {
        let key = b"ACGTACGT";
        let h1 = murmurhash3_64(key, 42);
        let h2 = murmurhash3_64(key, 42);
        assert_eq!(h1, h2);

        // Different seed gives different hash
        let h3 = murmurhash3_64(key, 99);
        assert_ne!(h1, h3);
    }
}
