//! Compression utilities with algorithm auto-detection.

use std::io::{Read, Write};

use crate::{CyaneaError, Result};

/// Supported compression algorithms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algorithm {
    Zstd,
    Gzip,
}

/// Compress data using zstd at the given level (1–22).
pub fn zstd_compress(data: &[u8], level: i32) -> Result<Vec<u8>> {
    zstd::encode_all(data, level).map_err(|e| CyaneaError::Compression(e.to_string()))
}

/// Decompress zstd data.
pub fn zstd_decompress(data: &[u8]) -> Result<Vec<u8>> {
    zstd::decode_all(data).map_err(|e| CyaneaError::Compression(e.to_string()))
}

/// Compress data using gzip at the given level (0–9).
pub fn gzip_compress(data: &[u8], level: u32) -> Result<Vec<u8>> {
    use flate2::write::GzEncoder;
    use flate2::Compression;

    let mut encoder = GzEncoder::new(Vec::new(), Compression::new(level));
    encoder
        .write_all(data)
        .map_err(|e| CyaneaError::Compression(e.to_string()))?;
    encoder
        .finish()
        .map_err(|e| CyaneaError::Compression(e.to_string()))
}

/// Decompress gzip data.
pub fn gzip_decompress(data: &[u8]) -> Result<Vec<u8>> {
    use flate2::read::GzDecoder;

    let mut decoder = GzDecoder::new(data);
    let mut decompressed = Vec::new();
    decoder
        .read_to_end(&mut decompressed)
        .map_err(|e| CyaneaError::Compression(e.to_string()))?;
    Ok(decompressed)
}

/// Detect the compression algorithm from the magic bytes of `data`.
///
/// Returns `None` if the data does not match a known format.
pub fn detect_algorithm(data: &[u8]) -> Option<Algorithm> {
    if data.len() >= 4 && data[..4] == [0x28, 0xB5, 0x2F, 0xFD] {
        Some(Algorithm::Zstd)
    } else if data.len() >= 2 && data[..2] == [0x1F, 0x8B] {
        Some(Algorithm::Gzip)
    } else {
        None
    }
}

/// Decompress data by auto-detecting the algorithm from magic bytes.
///
/// Returns an error if the format is unrecognised.
pub fn decompress(data: &[u8]) -> Result<Vec<u8>> {
    match detect_algorithm(data) {
        Some(Algorithm::Zstd) => zstd_decompress(data),
        Some(Algorithm::Gzip) => gzip_decompress(data),
        None => Err(CyaneaError::Compression(
            "unknown compression format".into(),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zstd_roundtrip() {
        let original = b"Hello, world! This is some test data for compression.";
        let compressed = zstd_compress(original, 3).unwrap();
        let decompressed = zstd_decompress(&compressed).unwrap();
        assert_eq!(original.to_vec(), decompressed);
    }

    #[test]
    fn test_gzip_roundtrip() {
        let original = b"Hello, world! This is gzip test data.";
        let compressed = gzip_compress(original, 6).unwrap();
        let decompressed = gzip_decompress(&compressed).unwrap();
        assert_eq!(original.to_vec(), decompressed);
    }

    #[test]
    fn test_detect_zstd() {
        let compressed = zstd_compress(b"test", 3).unwrap();
        assert_eq!(detect_algorithm(&compressed), Some(Algorithm::Zstd));
    }

    #[test]
    fn test_detect_gzip() {
        let compressed = gzip_compress(b"test", 6).unwrap();
        assert_eq!(detect_algorithm(&compressed), Some(Algorithm::Gzip));
    }

    #[test]
    fn test_detect_unknown() {
        assert_eq!(detect_algorithm(b"not compressed"), None);
    }

    #[test]
    fn test_auto_decompress_zstd() {
        let original = b"auto-detect zstd";
        let compressed = zstd_compress(original, 3).unwrap();
        let decompressed = decompress(&compressed).unwrap();
        assert_eq!(original.to_vec(), decompressed);
    }

    #[test]
    fn test_auto_decompress_gzip() {
        let original = b"auto-detect gzip";
        let compressed = gzip_compress(original, 6).unwrap();
        let decompressed = decompress(&compressed).unwrap();
        assert_eq!(original.to_vec(), decompressed);
    }

    #[test]
    fn test_auto_decompress_unknown() {
        let result = decompress(b"not compressed data");
        assert!(result.is_err());
    }
}
