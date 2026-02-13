//! Rank/select bitvectors and wavelet matrix.
//!
//! [`RankSelectBitVec`] supports O(1) rank queries and O(log n) select queries
//! using a two-level index over u64 blocks.
//!
//! [`WaveletMatrix`] uses a stack of bitvectors to support access, rank, and
//! select over integer alphabets in O(log σ) time.

use crate::{CyaneaError, Result};

/// Superblock size in bits (must be a multiple of 64).
const SUPERBLOCK_SIZE: usize = 512;
/// Number of u64 blocks per superblock.
const BLOCKS_PER_SUPER: usize = SUPERBLOCK_SIZE / 64;

/// A bitvector with O(1) rank and O(log n) select support.
///
/// Uses u64 blocks with a superblock index every 512 bits. Rank queries
/// use `u64::count_ones()` for popcount. Select uses binary search over
/// superblocks.
#[derive(Debug, Clone)]
pub struct RankSelectBitVec {
    blocks: Vec<u64>,
    /// Cumulative popcount at the start of each superblock.
    superblocks: Vec<usize>,
    len: usize,
}

impl RankSelectBitVec {
    /// Build a bitvector from a slice of booleans.
    pub fn build(bits: &[bool]) -> Self {
        let n = bits.len();
        let num_blocks = (n + 63) / 64;
        let mut blocks = vec![0u64; num_blocks];

        for (i, &b) in bits.iter().enumerate() {
            if b {
                blocks[i / 64] |= 1u64 << (i % 64);
            }
        }

        // Build superblock index
        // superblocks[i] = cumulative popcount before the i-th superblock group.
        // An extra sentinel entry stores the total count so binary search works.
        let num_super_groups = (num_blocks + BLOCKS_PER_SUPER - 1) / BLOCKS_PER_SUPER;
        let mut superblocks = vec![0usize; num_super_groups + 1];
        let mut cumulative = 0usize;
        for (i, block) in blocks.iter().enumerate() {
            if i % BLOCKS_PER_SUPER == 0 {
                superblocks[i / BLOCKS_PER_SUPER] = cumulative;
            }
            cumulative += block.count_ones() as usize;
        }
        superblocks[num_super_groups] = cumulative;

        Self {
            blocks,
            superblocks,
            len: n,
        }
    }

    /// Get the bit at position `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= len`.
    pub fn get(&self, i: usize) -> bool {
        assert!(i < self.len, "index out of bounds");
        (self.blocks[i / 64] >> (i % 64)) & 1 == 1
    }

    /// Count the number of 1-bits in positions `[0, i)`.
    ///
    /// Returns 0 if `i == 0`. Panics if `i > len`.
    pub fn rank1(&self, i: usize) -> usize {
        assert!(i <= self.len, "rank1: index out of bounds");
        if i == 0 {
            return 0;
        }

        let block_idx = (i - 1) / 64;
        let super_idx = block_idx / BLOCKS_PER_SUPER;
        let mut count = self.superblocks[super_idx];

        // Count blocks within the superblock
        let first_block = super_idx * BLOCKS_PER_SUPER;
        for b in first_block..block_idx {
            count += self.blocks[b].count_ones() as usize;
        }

        // Count bits within the last block
        let bit_pos = i % 64;
        if bit_pos == 0 {
            count += self.blocks[block_idx].count_ones() as usize;
        } else {
            let mask = (1u64 << bit_pos) - 1;
            count += (self.blocks[block_idx] & mask).count_ones() as usize;
        }

        count
    }

    /// Count the number of 0-bits in positions `[0, i)`.
    pub fn rank0(&self, i: usize) -> usize {
        i - self.rank1(i)
    }

    /// Find the position of the `k`-th 1-bit (1-indexed).
    ///
    /// Returns `None` if there are fewer than `k` set bits.
    pub fn select1(&self, k: usize) -> Option<usize> {
        if k == 0 || k > self.count_ones() {
            return None;
        }

        // Binary search on superblocks
        let mut lo = 0;
        let mut hi = self.superblocks.len() - 1;
        while lo < hi {
            let mid = lo + (hi - lo + 1) / 2;
            if self.superblocks[mid] < k {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }

        let mut remaining = k - self.superblocks[lo];
        let first_block = lo * BLOCKS_PER_SUPER;

        // Scan blocks within the superblock
        for b in first_block..self.blocks.len() {
            let popcnt = self.blocks[b].count_ones() as usize;
            if popcnt >= remaining {
                // Find the exact bit position within this block
                let mut word = self.blocks[b];
                for _ in 1..remaining {
                    word &= word - 1; // clear lowest set bit
                }
                let bit_in_block = word.trailing_zeros() as usize;
                let pos = b * 64 + bit_in_block;
                return if pos < self.len { Some(pos) } else { None };
            }
            remaining -= popcnt;
        }

        None
    }

    /// Find the position of the `k`-th 0-bit (1-indexed).
    ///
    /// Returns `None` if there are fewer than `k` zero bits.
    pub fn select0(&self, k: usize) -> Option<usize> {
        if k == 0 || k > self.count_zeros() {
            return None;
        }

        // Linear scan (acceptable for most bioinformatics use cases)
        let mut remaining = k;
        for (b, &block) in self.blocks.iter().enumerate() {
            let zeros_in_block = if (b + 1) * 64 <= self.len {
                64 - block.count_ones() as usize
            } else {
                let valid_bits = self.len - b * 64;
                let mask = if valid_bits >= 64 {
                    u64::MAX
                } else {
                    (1u64 << valid_bits) - 1
                };
                valid_bits - (block & mask).count_ones() as usize
            };

            if zeros_in_block >= remaining {
                // Find the exact bit position
                let mut word = if (b + 1) * 64 <= self.len {
                    !block
                } else {
                    let valid_bits = self.len - b * 64;
                    let mask = (1u64 << valid_bits) - 1;
                    !block & mask
                };
                for _ in 1..remaining {
                    word &= word - 1;
                }
                let bit_in_block = word.trailing_zeros() as usize;
                let pos = b * 64 + bit_in_block;
                return if pos < self.len { Some(pos) } else { None };
            }
            remaining -= zeros_in_block;
        }

        None
    }

    /// Total number of bits in the bitvector.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the bitvector has zero length.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Total number of 1-bits.
    pub fn count_ones(&self) -> usize {
        self.rank1(self.len)
    }

    /// Total number of 0-bits.
    pub fn count_zeros(&self) -> usize {
        self.len - self.count_ones()
    }
}

// ── Wavelet Matrix ───────────────────────────────────────────────────────

/// A wavelet matrix over an integer alphabet `[0, σ)`.
///
/// Supports access, rank, and select queries in O(log σ) time using
/// ⌈log₂ σ⌉ levels of [`RankSelectBitVec`]s.
#[derive(Debug, Clone)]
pub struct WaveletMatrix {
    levels: Vec<RankSelectBitVec>,
    /// Number of zeros at each level (for navigation).
    num_zeros: Vec<usize>,
    sigma: usize,
    len: usize,
}

impl WaveletMatrix {
    /// Build a wavelet matrix from a sequence of symbols in `[0, sigma)`.
    ///
    /// # Errors
    ///
    /// Returns an error if any symbol is ≥ `sigma` or `sigma` is 0.
    pub fn build(symbols: &[usize], sigma: usize) -> Result<Self> {
        if sigma == 0 {
            return Err(CyaneaError::InvalidInput(
                "WaveletMatrix: sigma must be positive".into(),
            ));
        }
        if let Some(&s) = symbols.iter().find(|&&s| s >= sigma) {
            return Err(CyaneaError::InvalidInput(format!(
                "WaveletMatrix: symbol {} out of range [0, {})",
                s, sigma
            )));
        }

        let n = symbols.len();
        let num_levels = if sigma <= 1 { 1 } else { (sigma as f64).log2().ceil() as usize };

        let mut levels = Vec::with_capacity(num_levels);
        let mut num_zeros = Vec::with_capacity(num_levels);
        let mut current = symbols.to_vec();

        for level in (0..num_levels).rev() {
            let bit = 1 << level;
            let bits: Vec<bool> = current.iter().map(|&s| s & bit != 0).collect();
            let bv = RankSelectBitVec::build(&bits);
            let nz = bv.count_zeros();
            num_zeros.push(nz);
            levels.push(bv);

            // Stable partition: 0-bit symbols first, then 1-bit symbols
            let mut next = Vec::with_capacity(n);
            for &s in &current {
                if s & bit == 0 {
                    next.push(s);
                }
            }
            for &s in &current {
                if s & bit != 0 {
                    next.push(s);
                }
            }
            current = next;
        }

        Ok(Self {
            levels,
            num_zeros,
            sigma,
            len: n,
        })
    }

    /// Access the symbol at position `i`.
    ///
    /// Returns `None` if `i >= len`.
    pub fn access(&self, mut i: usize) -> Option<usize> {
        if i >= self.len {
            return None;
        }

        let mut symbol = 0;
        for (level_idx, bv) in self.levels.iter().enumerate() {
            let bit_val = 1 << (self.levels.len() - 1 - level_idx);
            if bv.get(i) {
                symbol |= bit_val;
                i = self.num_zeros[level_idx] + bv.rank1(i);
            } else {
                i = bv.rank0(i);
            }
        }

        Some(symbol)
    }

    /// Count occurrences of symbol `c` in positions `[0, i)`.
    pub fn rank(&self, c: usize, mut i: usize) -> usize {
        if c >= self.sigma || i == 0 {
            return 0;
        }
        if i > self.len {
            i = self.len;
        }

        let mut lo = 0;
        let mut hi = i;

        for (level_idx, bv) in self.levels.iter().enumerate() {
            let bit_val = 1 << (self.levels.len() - 1 - level_idx);
            if c & bit_val != 0 {
                lo = self.num_zeros[level_idx] + bv.rank1(lo);
                hi = self.num_zeros[level_idx] + bv.rank1(hi);
            } else {
                lo = bv.rank0(lo);
                hi = bv.rank0(hi);
            }
        }

        hi - lo
    }

    /// Find the position of the `k`-th occurrence of symbol `c` (1-indexed).
    ///
    /// Returns `None` if there are fewer than `k` occurrences.
    pub fn select(&self, c: usize, k: usize) -> Option<usize> {
        if c >= self.sigma || k == 0 {
            return None;
        }

        // Navigate down to find the range for symbol c
        let mut lo = 0usize;
        let mut hi = self.len;
        for (level_idx, bv) in self.levels.iter().enumerate() {
            let bit_val = 1 << (self.levels.len() - 1 - level_idx);
            if c & bit_val != 0 {
                lo = self.num_zeros[level_idx] + bv.rank1(lo);
                hi = self.num_zeros[level_idx] + bv.rank1(hi);
            } else {
                lo = bv.rank0(lo);
                hi = bv.rank0(hi);
            }
        }

        if k > hi - lo {
            return None;
        }

        // Navigate back up from position lo + k - 1
        let mut pos = lo + k - 1;
        for level_idx in (0..self.levels.len()).rev() {
            let bv = &self.levels[level_idx];
            let bit_val = 1 << (self.levels.len() - 1 - level_idx);
            if c & bit_val != 0 {
                // pos is in the 1-zone: pos = nz + rank1(original_pos)
                // We need to find original_pos such that nz + rank1(original_pos) == pos
                // i.e., rank1(original_pos) == pos - nz
                let target_rank = pos - self.num_zeros[level_idx] + 1;
                pos = bv.select1(target_rank)?;
            } else {
                let target_rank = pos + 1;
                pos = bv.select0(target_rank)?;
            }
        }

        Some(pos)
    }

    /// Length of the indexed sequence.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Alphabet size.
    pub fn sigma(&self) -> usize {
        self.sigma
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── RankSelectBitVec tests ───────────────────────────────────────

    #[test]
    fn rank_empty() {
        let bv = RankSelectBitVec::build(&[]);
        assert_eq!(bv.len(), 0);
        assert!(bv.is_empty());
        assert_eq!(bv.count_ones(), 0);
    }

    #[test]
    fn rank_basic() {
        // bits: 1 0 1 1 0 1 0 0
        let bits = [true, false, true, true, false, true, false, false];
        let bv = RankSelectBitVec::build(&bits);
        assert_eq!(bv.len(), 8);
        assert_eq!(bv.count_ones(), 4);
        assert_eq!(bv.count_zeros(), 4);

        assert_eq!(bv.rank1(0), 0);
        assert_eq!(bv.rank1(1), 1);
        assert_eq!(bv.rank1(2), 1);
        assert_eq!(bv.rank1(3), 2);
        assert_eq!(bv.rank1(4), 3);
        assert_eq!(bv.rank1(8), 4);
    }

    #[test]
    fn rank0_basic() {
        let bits = [true, false, true, true, false, true, false, false];
        let bv = RankSelectBitVec::build(&bits);
        assert_eq!(bv.rank0(0), 0);
        assert_eq!(bv.rank0(2), 1);
        assert_eq!(bv.rank0(8), 4);
    }

    #[test]
    fn get_bits() {
        let bits = [true, false, true, false];
        let bv = RankSelectBitVec::build(&bits);
        assert!(bv.get(0));
        assert!(!bv.get(1));
        assert!(bv.get(2));
        assert!(!bv.get(3));
    }

    #[test]
    fn select1_basic() {
        let bits = [true, false, true, true, false, true, false, false];
        let bv = RankSelectBitVec::build(&bits);
        assert_eq!(bv.select1(1), Some(0));
        assert_eq!(bv.select1(2), Some(2));
        assert_eq!(bv.select1(3), Some(3));
        assert_eq!(bv.select1(4), Some(5));
        assert_eq!(bv.select1(5), None);
        assert_eq!(bv.select1(0), None);
    }

    #[test]
    fn select0_basic() {
        let bits = [true, false, true, true, false, true, false, false];
        let bv = RankSelectBitVec::build(&bits);
        assert_eq!(bv.select0(1), Some(1));
        assert_eq!(bv.select0(2), Some(4));
        assert_eq!(bv.select0(3), Some(6));
        assert_eq!(bv.select0(4), Some(7));
        assert_eq!(bv.select0(5), None);
    }

    #[test]
    fn rank_large_bitvec() {
        // Test with > 512 bits to exercise superblock logic
        let n = 1000;
        let bits: Vec<bool> = (0..n).map(|i| i % 3 == 0).collect();
        let bv = RankSelectBitVec::build(&bits);

        // Verify rank against brute force
        for i in (0..=n).step_by(100) {
            let expected = bits[..i].iter().filter(|&&b| b).count();
            assert_eq!(bv.rank1(i), expected, "rank1({}) mismatch", i);
        }
    }

    #[test]
    fn select1_large() {
        let n = 1000;
        let bits: Vec<bool> = (0..n).map(|i| i % 3 == 0).collect();
        let bv = RankSelectBitVec::build(&bits);
        // First 1-bit is at position 0, second at position 3, etc.
        assert_eq!(bv.select1(1), Some(0));
        assert_eq!(bv.select1(2), Some(3));
        assert_eq!(bv.select1(3), Some(6));
    }

    #[test]
    fn all_ones() {
        let bits = vec![true; 200];
        let bv = RankSelectBitVec::build(&bits);
        assert_eq!(bv.count_ones(), 200);
        assert_eq!(bv.rank1(100), 100);
        assert_eq!(bv.select1(50), Some(49));
    }

    #[test]
    fn all_zeros() {
        let bits = vec![false; 200];
        let bv = RankSelectBitVec::build(&bits);
        assert_eq!(bv.count_zeros(), 200);
        assert_eq!(bv.rank0(100), 100);
        assert_eq!(bv.select0(50), Some(49));
        assert_eq!(bv.select1(1), None);
    }

    // ── WaveletMatrix tests ──────────────────────────────────────────

    #[test]
    fn wavelet_access() {
        let data = [3, 1, 4, 1, 5, 9, 2, 6];
        let wm = WaveletMatrix::build(&data, 10).unwrap();
        for (i, &expected) in data.iter().enumerate() {
            assert_eq!(wm.access(i), Some(expected), "access({}) failed", i);
        }
        assert_eq!(wm.access(8), None);
    }

    #[test]
    fn wavelet_rank() {
        let data = [3, 1, 4, 1, 5, 9, 2, 6];
        let wm = WaveletMatrix::build(&data, 10).unwrap();
        assert_eq!(wm.rank(1, 4), 2); // two 1s in [0, 4)
        assert_eq!(wm.rank(1, 2), 1); // one 1 in [0, 2)
        assert_eq!(wm.rank(4, 3), 1); // one 4 in [0, 3)
        assert_eq!(wm.rank(7, 8), 0); // no 7s
    }

    #[test]
    fn wavelet_select() {
        let data = [3, 1, 4, 1, 5, 9, 2, 6];
        let wm = WaveletMatrix::build(&data, 10).unwrap();
        assert_eq!(wm.select(1, 1), Some(1)); // first 1 at index 1
        assert_eq!(wm.select(1, 2), Some(3)); // second 1 at index 3
        assert_eq!(wm.select(1, 3), None);    // no third 1
        assert_eq!(wm.select(3, 1), Some(0)); // first 3 at index 0
    }

    #[test]
    fn wavelet_binary_alphabet() {
        let data = [0, 1, 0, 1, 1, 0];
        let wm = WaveletMatrix::build(&data, 2).unwrap();
        assert_eq!(wm.rank(0, 6), 3);
        assert_eq!(wm.rank(1, 6), 3);
        assert_eq!(wm.select(0, 1), Some(0));
        assert_eq!(wm.select(0, 2), Some(2));
    }

    #[test]
    fn wavelet_single_symbol() {
        let data = [0, 0, 0, 0];
        let wm = WaveletMatrix::build(&data, 1).unwrap();
        assert_eq!(wm.access(0), Some(0));
        assert_eq!(wm.rank(0, 4), 4);
    }

    #[test]
    fn wavelet_empty() {
        let wm = WaveletMatrix::build(&[], 4).unwrap();
        assert_eq!(wm.len(), 0);
        assert!(wm.is_empty());
        assert_eq!(wm.access(0), None);
    }

    #[test]
    fn wavelet_invalid() {
        assert!(WaveletMatrix::build(&[], 0).is_err());
        assert!(WaveletMatrix::build(&[5], 4).is_err());
    }

    #[test]
    fn wavelet_dna_encoded() {
        // DNA as 0=A, 1=C, 2=G, 3=T
        let dna = [0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let wm = WaveletMatrix::build(&dna, 4).unwrap();
        assert_eq!(wm.rank(0, 8), 2); // 2 A's
        assert_eq!(wm.rank(1, 8), 2); // 2 C's
        assert_eq!(wm.select(2, 1), Some(2)); // first G at index 2
        assert_eq!(wm.select(3, 2), Some(7)); // second T at index 7
    }
}
