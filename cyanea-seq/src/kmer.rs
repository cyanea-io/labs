//! Zero-allocation k-mer iterator.
//!
//! Wraps [`std::slice::Windows`] over a pre-validated uppercase sequence.
//! Since input is already validated and normalized, no per-kmer processing
//! is needed.

use cyanea_core::{CyaneaError, Result};

/// Iterator over k-mer windows of a byte slice.
///
/// Yields `&[u8]` slices of length `k`. Implements [`ExactSizeIterator`]
/// and [`DoubleEndedIterator`] for efficient bidirectional traversal.
pub struct KmerIter<'a> {
    inner: std::slice::Windows<'a, u8>,
    remaining: usize,
}

impl<'a> KmerIter<'a> {
    /// Create a new k-mer iterator.
    ///
    /// `k` must be in `[1, 31]` and `k <= seq.len()` (unless the sequence
    /// is empty, in which case `k >= 1` is accepted and yields nothing — wait,
    /// actually: if k > seq.len() we return an error per the plan).
    pub fn new(seq: &'a [u8], k: usize) -> Result<Self> {
        if k == 0 {
            return Err(CyaneaError::InvalidInput(
                "k-mer size must be at least 1".into(),
            ));
        }
        if k > 31 {
            return Err(CyaneaError::InvalidInput(
                "k-mer size must be at most 31".into(),
            ));
        }
        if k > seq.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "k-mer size {} exceeds sequence length {}",
                k,
                seq.len()
            )));
        }
        let remaining = seq.len() - k + 1;
        Ok(Self {
            inner: seq.windows(k),
            remaining,
        })
    }
}

impl<'a> Iterator for KmerIter<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        let item = self.inner.next()?;
        self.remaining -= 1;
        Some(item)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<'a> ExactSizeIterator for KmerIter<'a> {}

impl<'a> DoubleEndedIterator for KmerIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let item = self.inner.next_back()?;
        self.remaining -= 1;
        Some(item)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_k2() {
        let seq = b"ACGT";
        let kmers: Vec<&[u8]> = KmerIter::new(seq, 2).unwrap().collect();
        assert_eq!(kmers, vec![b"AC", b"CG", b"GT"]);
    }

    #[test]
    fn exact_size() {
        let seq = b"ACGTACGT";
        let iter = KmerIter::new(seq, 3).unwrap();
        assert_eq!(iter.len(), 6);
    }

    #[test]
    fn k_zero_error() {
        let result = KmerIter::new(b"ACGT", 0);
        assert!(result.is_err());
    }

    #[test]
    fn k_exceeds_len_error() {
        let result = KmerIter::new(b"AC", 3);
        assert!(result.is_err());
    }

    #[test]
    fn double_ended() {
        let seq = b"ACGT";
        let mut iter = KmerIter::new(seq, 2).unwrap();
        assert_eq!(iter.next_back(), Some(b"GT".as_slice()));
        assert_eq!(iter.next(), Some(b"AC".as_slice()));
        assert_eq!(iter.next(), Some(b"CG".as_slice()));
        assert_eq!(iter.next(), None);
    }
}
