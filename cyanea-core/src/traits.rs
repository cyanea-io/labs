//! Core trait definitions for the Cyanea ecosystem.
//!
//! These traits define the contracts that domain types implement across crates.

/// A biological sequence (DNA, RNA, protein, etc.).
pub trait Sequence {
    /// The raw byte representation of the sequence.
    fn as_bytes(&self) -> &[u8];

    /// Length in residues/bases.
    fn len(&self) -> usize {
        self.as_bytes().len()
    }

    /// Whether the sequence is empty.
    fn is_empty(&self) -> bool {
        self.as_bytes().is_empty()
    }
}

/// A type whose identity can be derived from its content via cryptographic hash.
pub trait ContentAddressable {
    /// Return the content hash as a hex string (e.g. SHA-256).
    fn content_hash(&self) -> String;
}

/// A type that can be compressed and decompressed.
pub trait Compressible: Sized {
    /// Compress to bytes.
    fn compress(&self) -> crate::Result<Vec<u8>>;

    /// Decompress from bytes.
    fn decompress(data: &[u8]) -> crate::Result<Self>;
}

/// A type that carries a numeric score (alignment score, quality, etc.).
pub trait Scored {
    /// The score value.
    fn score(&self) -> f64;
}

/// A type that carries annotations (names, descriptions, metadata).
pub trait Annotated {
    /// A human-readable name or identifier.
    fn name(&self) -> &str;

    /// An optional description.
    fn description(&self) -> Option<&str> {
        None
    }
}

/// A type that can produce a summary of its contents.
pub trait Summarizable {
    /// A one-line summary suitable for display.
    fn summary(&self) -> String;
}
