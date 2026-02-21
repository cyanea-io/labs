# cyanea-py Architecture

## Module Map

```
cyanea-py
├── lib.rs          PyO3 module root, submodule registration
├── error.rs        CyaneaError → PyErr mapping (IntoPyResult trait)
├── seq.rs          cyanea.seq — sequence types and FASTA/FASTQ
├── align.rs        cyanea.align — alignment functions
├── stats.rs        cyanea.stats — statistics and hypothesis tests
├── core_utils.rs   cyanea.core — SHA-256, compression
├── ml.rs           cyanea.ml — ML, distances, embeddings
├── chem.rs         cyanea.chem — SMILES, properties, fingerprints
├── struct_bio.rs   cyanea.struct_bio — PDB, geometry, DSSP
├── phylo.rs        cyanea.phylo — trees, distances, construction
├── io.rs           cyanea.io — file format parsing
├── omics.rs        cyanea.omics — genomic intervals, annotation
└── sc.rs           cyanea.sc — single-cell pipeline
```

## Module Registration

PyO3 uses a nested module pattern. In `lib.rs`:

1. `#[pymodule]` function `cyanea` is the root module
2. Each submodule has a `register(parent_module)` function
3. Modules are also added to `sys.modules` for direct import (`from cyanea import seq`)

## Error Mapping

The `IntoPyResult` trait in `error.rs` converts `Result<T, CyaneaError>` to `PyResult<T>`:

- `CyaneaError::Io` → `PyIOError`
- `CyaneaError::Parse` / `InvalidInput` → `PyValueError`
- `CyaneaError::Compression` / `Hash` / `Other` → `PyRuntimeError`

## Python Class Pattern

Classes use `#[pyclass(frozen)]` for thread-safe immutable objects:

```rust
#[pyclass(frozen)]
struct DnaSequence {
    inner: cyanea_seq::DnaSequence,
}

#[pymethods]
impl DnaSequence {
    #[new]
    fn new(data: &[u8]) -> PyResult<Self> { ... }
    fn gc_content(&self) -> f64 { self.inner.gc_content() }
}
```

Frozen classes are safe to share across threads and can be used as dict keys.

## NumPy Interop

Behind the `numpy` feature flag, ML functions have `_np` variants that return `numpy::PyArray`:

- Uses the `numpy` crate (0.23) for zero-copy array creation
- `PyArray::from_vec(py, data)` for 1D arrays
- `PyArray::from_vec2(py, &data)` for 2D arrays
- No copies when returning owned Vec data

## ABI Compatibility

PyO3 0.23 supports Python up to 3.13. For Python 3.14+:
- Set `PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1` environment variable
- Uses stable ABI (abi3) for forward compatibility

## GIL Management

- Short operations (property access, small computations): hold the GIL
- Long operations (alignment, PCA, clustering): could benefit from GIL release but currently hold it for simplicity
- Future optimization: `py.allow_threads(|| ...)` for expensive computations
