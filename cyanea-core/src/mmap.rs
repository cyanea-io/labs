//! Memory-mapped file access for zero-copy I/O.
//!
//! Only available with the `std` feature (disabled for WASM).

use memmap2::Mmap;
use std::fs::File;
use std::path::Path;

use crate::{CyaneaError, Result};

/// A read-only memory-mapped file.
pub struct MappedFile {
    _file: File,
    mmap: Mmap,
}

impl MappedFile {
    /// Open and memory-map a file.
    ///
    /// # Safety
    ///
    /// The caller must ensure that the file is not modified by another process
    /// while the mapping is active. This is safe for immutable data files
    /// (FASTA, BAM, etc.) which are the primary use case.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: {}", path.display(), e),
            ))
        })?;
        // SAFETY: We hold the File open for the lifetime of the mapping.
        // The caller is responsible for ensuring no concurrent modification.
        let mmap = unsafe { Mmap::map(&file) }?;
        Ok(Self { _file: file, mmap })
    }

    /// The mapped bytes.
    pub fn as_bytes(&self) -> &[u8] {
        &self.mmap
    }

    /// Length in bytes.
    pub fn len(&self) -> usize {
        self.mmap.len()
    }

    /// Whether the mapped region is empty.
    pub fn is_empty(&self) -> bool {
        self.mmap.is_empty()
    }
}

impl AsRef<[u8]> for MappedFile {
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_mmap_roundtrip() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(b"hello mmap").unwrap();
        file.flush().unwrap();

        let mapped = MappedFile::open(file.path()).unwrap();
        assert_eq!(mapped.as_bytes(), b"hello mmap");
        assert_eq!(mapped.len(), 10);
        assert!(!mapped.is_empty());
    }

    #[test]
    fn test_mmap_not_found() {
        let result = MappedFile::open("/nonexistent/file.txt");
        assert!(result.is_err());
    }
}
