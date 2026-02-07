//! Generic validated sequence type.
//!
//! [`ValidatedSeq<A>`] is a newtype over `Vec<u8>` parameterized by an
//! [`Alphabet`] marker type. Construction uppercases and validates every byte.
//! The inner data is always uppercase, so `Deref<Target=[u8]>` and
//! `as_bytes()` are zero-cost and safe to pass to downstream `&[u8]` APIs.

use std::fmt;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;
use std::ops::Deref;

use cyanea_core::{CyaneaError, ContentAddressable, Sequence, Summarizable};

use crate::alphabet::Alphabet;

/// A validated biological sequence parameterized by its alphabet.
///
/// `ValidatedSeq<DnaAlphabet>` is a DNA sequence, `ValidatedSeq<RnaAlphabet>`
/// is RNA, etc. The inner bytes are always uppercase.
#[derive(Clone)]
pub struct ValidatedSeq<A: Alphabet> {
    data: Vec<u8>,
    _alphabet: PhantomData<A>,
}

impl<A: Alphabet> ValidatedSeq<A> {
    /// Create a new validated sequence from raw bytes.
    ///
    /// Input is uppercased, then every byte is checked against the alphabet.
    /// Returns an error if any byte is not in the alphabet after uppercasing.
    pub fn new(bytes: impl AsRef<[u8]>) -> cyanea_core::Result<Self> {
        let data: Vec<u8> = bytes.as_ref().iter().map(|b| b.to_ascii_uppercase()).collect();
        for (i, &b) in data.iter().enumerate() {
            if !A::is_valid(b) {
                return Err(CyaneaError::InvalidInput(format!(
                    "invalid {} byte '{}' (0x{:02X}) at position {}",
                    A::NAME,
                    b as char,
                    b,
                    i
                )));
            }
        }
        Ok(Self {
            data,
            _alphabet: PhantomData,
        })
    }

    /// Create a sequence from pre-validated bytes, skipping validation.
    ///
    /// # Safety (logical)
    ///
    /// Caller must guarantee all bytes are valid uppercase members of `A`.
    pub(crate) fn from_validated(data: Vec<u8>) -> Self {
        Self {
            data,
            _alphabet: PhantomData,
        }
    }

    /// Consume the sequence and return the inner byte vector.
    pub fn into_bytes(self) -> Vec<u8> {
        self.data
    }
}

impl<A: Alphabet> Deref for ValidatedSeq<A> {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        &self.data
    }
}

impl<A: Alphabet> AsRef<[u8]> for ValidatedSeq<A> {
    fn as_ref(&self) -> &[u8] {
        &self.data
    }
}

impl<A: Alphabet> Sequence for ValidatedSeq<A> {
    fn as_bytes(&self) -> &[u8] {
        &self.data
    }
}

impl<A: Alphabet> ContentAddressable for ValidatedSeq<A> {
    fn content_hash(&self) -> String {
        cyanea_core::hash::sha256(&self.data)
    }
}

impl<A: Alphabet> Summarizable for ValidatedSeq<A> {
    fn summary(&self) -> String {
        let preview_len = self.data.len().min(20);
        let preview = std::str::from_utf8(&self.data[..preview_len]).unwrap_or("???");
        if self.data.len() > 20 {
            format!("{} sequence ({} bp): {}...", A::NAME, self.data.len(), preview)
        } else {
            format!("{} sequence ({} bp): {}", A::NAME, self.data.len(), preview)
        }
    }
}

impl<A: Alphabet> fmt::Debug for ValidatedSeq<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = std::str::from_utf8(&self.data).unwrap_or("???");
        write!(f, "{}(\"{}\")", A::NAME, s)
    }
}

impl<A: Alphabet> fmt::Display for ValidatedSeq<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = std::str::from_utf8(&self.data).unwrap_or("???");
        f.write_str(s)
    }
}

impl<A: Alphabet> PartialEq for ValidatedSeq<A> {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl<A: Alphabet> Eq for ValidatedSeq<A> {}

impl<A: Alphabet> Hash for ValidatedSeq<A> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.data.hash(state);
    }
}

#[cfg(feature = "serde")]
impl<A: Alphabet> serde::Serialize for ValidatedSeq<A> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error> {
        let s = std::str::from_utf8(&self.data).map_err(serde::ser::Error::custom)?;
        serializer.serialize_str(s)
    }
}

#[cfg(feature = "serde")]
impl<'de, A: Alphabet> serde::Deserialize<'de> for ValidatedSeq<A> {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> std::result::Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        Self::new(s.as_bytes()).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::DnaAlphabet;

    type DnaSeq = ValidatedSeq<DnaAlphabet>;

    #[test]
    fn stores_uppercase() {
        let seq = DnaSeq::new(b"acgt").unwrap();
        assert_eq!(seq.as_bytes(), b"ACGT");
    }

    #[test]
    fn empty_sequence_ok() {
        let seq = DnaSeq::new(b"").unwrap();
        assert!(seq.is_empty());
        assert_eq!(seq.len(), 0);
    }

    #[test]
    fn as_bytes_uppercase() {
        let seq = DnaSeq::new(b"AcGtN").unwrap();
        assert_eq!(seq.as_bytes(), b"ACGTN");
    }

    #[test]
    fn deref_to_slice() {
        let seq = DnaSeq::new(b"ACGT").unwrap();
        let slice: &[u8] = &*seq;
        assert_eq!(slice, b"ACGT");
        assert_eq!(seq[0], b'A');
    }

    #[test]
    fn content_addressable_deterministic() {
        let seq1 = DnaSeq::new(b"ACGT").unwrap();
        let seq2 = DnaSeq::new(b"acgt").unwrap();
        assert_eq!(seq1.content_hash(), seq2.content_hash());
    }

    #[test]
    fn rejects_invalid_bytes() {
        let result = DnaSeq::new(b"ACGX");
        assert!(result.is_err());
    }
}
