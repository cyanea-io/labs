//! Core utility wrappers: SHA-256 hashing and zstd compression.

use crate::error::wasm_ok;
#[cfg(feature = "compress")]
use crate::error::wasm_err;

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

/// SHA-256 hash of a string, returned as JSON hex string.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sha256(data: &str) -> String {
    wasm_ok(&cyanea_core::hash::sha256(data.as_bytes()))
}

/// Compress a string with zstd at the given level, return JSON byte array.
#[cfg(feature = "compress")]
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn zstd_compress(data: &str, level: i32) -> String {
    match cyanea_core::compress::zstd_compress(data.as_bytes(), level) {
        Ok(bytes) => wasm_ok(&bytes),
        Err(e) => wasm_err(e),
    }
}

/// Decompress zstd data from a JSON byte array, return JSON string.
#[cfg(feature = "compress")]
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn zstd_decompress(data_json: &str) -> String {
    let bytes: Vec<u8> = match serde_json::from_str(data_json) {
        Ok(b) => b,
        Err(e) => return wasm_err(format!("invalid JSON byte array: {e}")),
    };
    match cyanea_core::compress::zstd_decompress(&bytes) {
        Ok(decompressed) => match String::from_utf8(decompressed) {
            Ok(s) => wasm_ok(&s),
            Err(e) => wasm_err(format!("decompressed data is not UTF-8: {e}")),
        },
        Err(e) => wasm_err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_known() {
        let json = sha256("hello");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let hash = v["ok"].as_str().unwrap();
        assert_eq!(
            hash,
            "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
        );
    }

    #[cfg(feature = "compress")]
    #[test]
    fn zstd_roundtrip() {
        let compressed_json = zstd_compress("hello world", 3);
        let v: serde_json::Value = serde_json::from_str(&compressed_json).unwrap();
        let bytes = v["ok"].as_array().unwrap();
        assert!(!bytes.is_empty());

        // Round-trip through decompress
        let bytes_json = serde_json::to_string(&v["ok"]).unwrap();
        let decompressed_json = zstd_decompress(&bytes_json);
        let v2: serde_json::Value = serde_json::from_str(&decompressed_json).unwrap();
        assert_eq!(v2["ok"].as_str().unwrap(), "hello world");
    }
}
