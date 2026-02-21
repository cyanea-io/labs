# Usage Guide -- cyanea-gpu

Practical examples for GPU-accelerated compute in bioinformatics.

## Backend Selection and Auto-Detection

```rust
use cyanea_gpu::{auto_backend, Backend, BackendKind};

// Auto-detect the best available backend
let backend = auto_backend();
let info = backend.device_info();
println!("Backend: {} ({:?})", info.name, info.kind);
println!("Memory: {} bytes", info.total_memory);
println!("Parallelism: {}", info.max_parallelism);

// Check what we got
match info.kind {
    BackendKind::Metal => println!("Using Apple Metal GPU"),
    BackendKind::Cuda => println!("Using NVIDIA CUDA GPU"),
    BackendKind::Wgpu => println!("Using WebGPU"),
    BackendKind::Cpu => println!("Using CPU fallback"),
}
```

## GPU Buffer Creation and Data Transfer

```rust
use cyanea_gpu::{CpuBackend, Backend};

let backend = CpuBackend::new();

// Create a buffer from host data
let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
let buf = backend.buffer_from_slice(&data).unwrap();
println!("Buffer length: {}", backend.buffer_len(&buf));

// Read buffer contents back to host
let host_data = backend.read_buffer(&buf).unwrap();
assert_eq!(host_data, data);

// Create a zero-initialized buffer
let zeros = backend.buffer_zeros(100).unwrap();
let sum = backend.reduce_sum(&zeros).unwrap();
assert_eq!(sum, 0.0);

// Overwrite buffer contents
let new_data = vec![10.0, 20.0, 30.0, 40.0, 50.0];
backend.write_buffer(&buf, &new_data).unwrap();
```

## Pairwise Distance Computation on GPU

```rust
use cyanea_gpu::{auto_backend, ops, Backend, DistanceMetricGpu};

let backend = auto_backend();

// 4 points in 3D space (flattened row-major)
let points = vec![
    0.0, 0.0, 0.0,  // point 0
    1.0, 0.0, 0.0,  // point 1
    0.0, 1.0, 0.0,  // point 2
    1.0, 1.0, 0.0,  // point 3
];
let buf = backend.buffer_from_slice(&points).unwrap();

// Compute 4x4 Euclidean distance matrix
let dist_buf = backend.pairwise_distance_matrix(&buf, 4, 3, DistanceMetricGpu::Euclidean).unwrap();
let distances = backend.read_buffer(&dist_buf).unwrap();
println!("Distance matrix (4x4): {:?}", distances);

// For larger-than-memory matrices, use tiled computation
let large_data: Vec<f64> = (0..30000).map(|i| i as f64 * 0.001).collect();
let tiled = ops::tiled_pairwise_distance(
    &*backend, &large_data, 100, 300,
    DistanceMetricGpu::Euclidean, 25,
).unwrap();
println!("Tiled result: {} entries", tiled.len());
```

## GPU-Accelerated K-mer Counting

```rust
use cyanea_gpu::kmer::{gpu_kmer_count, gpu_kmer_count_cpu};

let sequences = vec![
    b"ATCGATCGATCGATCG".to_vec(),
    b"GCTAGCTAGCTAGCTA".to_vec(),
    b"AAAACCCCGGGGTTTT".to_vec(),
];

// Auto-dispatch: uses Metal/CUDA if available, falls back to CPU
let result = gpu_kmer_count(&sequences, 4).unwrap();
println!("K={}, total k-mer slots: {}", result.k, result.counts.len());

// Find most frequent k-mers
let max_count = *result.counts.iter().max().unwrap();
println!("Most frequent k-mer count: {}", max_count);

// Force CPU-only computation
let cpu_result = gpu_kmer_count_cpu(&sequences, 4).unwrap();
assert_eq!(result.counts, cpu_result.counts);
```

## GPU Batch Smith-Waterman Alignment

```rust
use cyanea_gpu::smith_waterman::{gpu_smith_waterman_batch, gpu_sw_cpu, BLOSUM62_24};

// Pairs of (query, target) protein sequences
let pairs = vec![
    (b"ACDEFGHIKLMNPQ".to_vec(), b"ACDEFGHIKLMNPQ".to_vec()),
    (b"ACDEFG".to_vec(), b"ACDXFG".to_vec()),
    (b"MKTAYIAKQRQISFVK".to_vec(), b"MKTAYIAKQRQ".to_vec()),
];

// Batch alignment with BLOSUM62, gap_open=11, gap_extend=1
let results = gpu_smith_waterman_batch(&pairs, &BLOSUM62_24, 11, 1).unwrap();

for (i, r) in results.iter().enumerate() {
    println!("Pair {}: score={}, query_end={}, target_end={}",
        i, r.score, r.query_end, r.target_end);
}
```

## GPU MinHash Sketching

```rust
use cyanea_gpu::minhash_gpu::{gpu_minhash, gpu_minhash_jaccard};

let sequences = vec![
    b"ATCGATCGATCGATCGATCG".to_vec(),
    b"ATCGATCGATCGATCGATCA".to_vec(),  // 1 base different
    b"GCTAGCTAGCTAGCTAGCTA".to_vec(),  // completely different
];

// Compute MinHash sketches with k=5, sketch_size=128
let sketches = gpu_minhash(&sequences, 5, 128).unwrap();
println!("Generated {} sketches of size {}", sketches.len(), sketches[0].hashes.len());

// Compare sketches using Jaccard similarity
let sim_01 = gpu_minhash_jaccard(&sketches[0], &sketches[1]).unwrap();
let sim_02 = gpu_minhash_jaccard(&sketches[0], &sketches[2]).unwrap();
println!("Jaccard(seq0, seq1) = {:.4} (similar)", sim_01);
println!("Jaccard(seq0, seq2) = {:.4} (different)", sim_02);
```

## CPU Fallback Backend

The CPU backend is always available and requires no feature flags:

```rust
use cyanea_gpu::{CpuBackend, Backend, DistanceMetricGpu};

let cpu = CpuBackend::new();
assert_eq!(cpu.device_info().kind, cyanea_gpu::BackendKind::Cpu);

// All operations work identically to GPU backends
let buf = cpu.buffer_from_slice(&[1.0, 2.0, 3.0]).unwrap();
let sum = cpu.reduce_sum(&buf).unwrap();
let min = cpu.reduce_min(&buf).unwrap();
let max = cpu.reduce_max(&buf).unwrap();
println!("sum={}, min={}, max={}", sum, min, max);

// Matrix multiply
let a = cpu.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();  // 2x2
let b = cpu.buffer_from_slice(&[5.0, 6.0, 7.0, 8.0]).unwrap();  // 2x2
let c = cpu.matrix_multiply(&a, &b, 2, 2, 2).unwrap();
let result = cpu.read_buffer(&c).unwrap();
println!("Matrix product: {:?}", result);
```

## Feature Flag Combinations

Configure backends via Cargo feature flags:

```toml
# CPU only (default, no features needed)
[dependencies]
cyanea-gpu = "0.1"

# Apple Metal on macOS
[dependencies]
cyanea-gpu = { version = "0.1", features = ["metal"] }

# NVIDIA CUDA
[dependencies]
cyanea-gpu = { version = "0.1", features = ["cuda"] }

# WebGPU (cross-platform: Vulkan, Metal, DX12, WebGPU)
[dependencies]
cyanea-gpu = { version = "0.1", features = ["wgpu"] }

# Multiple backends (auto_backend picks the best)
[dependencies]
cyanea-gpu = { version = "0.1", features = ["metal", "cuda", "wgpu"] }

# CPU parallelism with Rayon (no GPU)
[dependencies]
cyanea-gpu = { version = "0.1", features = ["parallel"] }

# Full: all backends + parallelism + serialization
[dependencies]
cyanea-gpu = { version = "0.1", features = ["metal", "cuda", "wgpu", "parallel", "serde"] }
```

The `auto_backend()` function selects the best available backend at runtime based on compiled features and hardware availability: Metal > CUDA > WebGPU > CPU.
