# cyanea-core Usage Guide

## Installation

```toml
[dependencies]
cyanea-core = { version = "0.1", features = ["std"] }
```

Without the `std` feature, compression and memory-mapped file access are disabled. The core traits, error types, probability types, bitvectors, and Fenwick tree are always available.

## Error Handling

All crates in the Cyanea ecosystem use `CyaneaError` as their unified error type, with `Result<T>` as a convenience alias.

```rust
use cyanea_core::{CyaneaError, Result};

fn parse_quality(s: &str) -> Result<u8> {
    let q: u8 = s.parse().map_err(|e| {
        CyaneaError::Parse(format!("invalid quality score '{}': {}", s, e))
    })?;
    if q > 60 {
        return Err(CyaneaError::InvalidInput(
            format!("quality score {} exceeds maximum 60", q),
        ));
    }
    Ok(q)
}

// CyaneaError variants:
//   Io(std::io::Error)       - file/network I/O errors
//   Parse(String)            - format parsing failures
//   InvalidInput(String)     - bad arguments or data
//   Compression(String)      - zstd/gzip failures
//   Hash(String)             - hashing errors
//   Other(String)            - catch-all
```

## Log-Space Probability

`LogProb` stores probabilities as natural logarithms to avoid floating-point underflow when multiplying many small probabilities together. `PhredProb` wraps Phred quality scores and converts to/from `LogProb`.

```rust
use cyanea_core::{LogProb, PhredProb};

// Create from raw probabilities
let p = LogProb::from_prob(0.001).unwrap();   // ln(0.001)
let q = LogProb::from_prob(0.002).unwrap();   // ln(0.002)

// Log-space arithmetic avoids underflow
let product = p.ln_mul(q);                     // ln(0.001 * 0.002)
let sum = p.ln_add(q);                         // ln(0.001 + 0.002) via log-sum-exp

// Special values
let certain = LogProb::certain();              // ln(1) = 0.0
let impossible = LogProb::impossible();        // ln(0) = -inf

// Convert back
let raw: f64 = p.to_prob();                   // 0.001

// Phred quality scores
let phred = PhredProb::from_phred(30.0).unwrap();  // Q30
let error_prob = phred.to_prob();                    // 0.001

// Bidirectional conversion
let log_from_phred: LogProb = phred.into();
let phred_from_log: PhredProb = p.into();
```

## Rank/Select Bitvectors

`RankSelectBitVec` answers rank and select queries on a bitvector in O(1) and O(log n) time respectively, using a superblock index built over 64-bit popcount blocks.

```rust
use cyanea_core::RankSelectBitVec;

// Build from a boolean slice
let bits = vec![true, false, true, true, false, true, false, false];
let bv = RankSelectBitVec::build(&bits);

// rank1(i) counts 1-bits in [0, i)
assert_eq!(bv.rank1(0), 0);   // no 1-bits before position 0
assert_eq!(bv.rank1(4), 3);   // bits[0..4] has three 1-bits

// select1(k) finds the position of the k-th 1-bit (1-indexed)
assert_eq!(bv.select1(1), Some(0));  // first 1-bit is at position 0
assert_eq!(bv.select1(3), Some(3));  // third 1-bit is at position 3

// Wavelet matrix for alphabet rank/select
use cyanea_core::WaveletMatrix;

let symbols = vec![3, 1, 4, 1, 5, 9, 2, 6];
let wm = WaveletMatrix::build(&symbols, 10).unwrap();

assert_eq!(wm.access(2), Some(4));       // symbol at position 2
assert_eq!(wm.rank(1, 4), 2);            // two occurrences of '1' in [0, 4)
assert_eq!(wm.select(1, 1), Some(1));    // first '1' is at position 1
```

## Fenwick Tree

`FenwickTree<T>` provides O(log n) prefix sums and point updates, useful for cumulative frequency tables, running statistics, and coordination-compressed DP.

```rust
use cyanea_core::FenwickTree;

// Build from a slice
let values = vec![3, 1, 4, 1, 5, 9, 2, 6];
let mut tree = FenwickTree::from_slice(&values);

// Prefix sum: sum of elements in [0, i] (inclusive)
assert_eq!(tree.prefix_sum(3), 9);    // 3 + 1 + 4 + 1

// Range sum: sum of elements in [l, r] (inclusive)
assert_eq!(tree.range_sum(2, 5), 19); // 4 + 1 + 5 + 9

// Point update: add delta to element at index i
tree.update(2, 10);                    // values[2] is now 14
assert_eq!(tree.prefix_sum(3), 19);   // 3 + 1 + 14 + 1
```

## Compression

The `compress` module (requires `std` feature) provides zstd and gzip compression with auto-detection by magic bytes.

```rust
use cyanea_core::compress::{
    zstd_compress, zstd_decompress,
    gzip_compress, gzip_decompress,
    detect_algorithm, decompress, Algorithm,
};

let data = b"ACGTACGTACGTACGT".repeat(1000);

// Zstd compression (level 3 is the default, range 1-22)
let compressed = zstd_compress(&data, 3).unwrap();
let decompressed = zstd_decompress(&compressed).unwrap();
assert_eq!(data.as_slice(), decompressed.as_slice());

// Gzip compression
let gz = gzip_compress(&data, 6).unwrap();

// Auto-detect from magic bytes
assert_eq!(detect_algorithm(&compressed), Some(Algorithm::Zstd));
assert_eq!(detect_algorithm(&gz), Some(Algorithm::Gzip));

// Auto-detect and decompress
let auto = decompress(&compressed).unwrap();
assert_eq!(data.as_slice(), auto.as_slice());
```

## SHA-256 Hashing

Content-addressed hashing for data integrity and deduplication.

```rust
use cyanea_core::{sha256, sha256_file};

// Hash in-memory data
let hash = sha256(b"ACGTACGT");
assert_eq!(hash.len(), 64); // 256 bits = 64 hex chars

// Hash a file (streaming, memory-efficient)
// let file_hash = sha256_file("genome.fa").unwrap();
```
