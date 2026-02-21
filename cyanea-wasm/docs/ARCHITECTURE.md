# cyanea-wasm Architecture

## Module Map

```
cyanea-wasm
├── lib.rs          Module declarations, re-exports, VERSION constant
├── error.rs        JSON error envelope (wasm_ok/wasm_err/wasm_result)
├── seq.rs          Wraps cyanea-seq
├── align.rs        Wraps cyanea-align
├── stats.rs        Wraps cyanea-stats
├── ml.rs           Wraps cyanea-ml
├── chem.rs         Wraps cyanea-chem
├── struct_bio.rs   Wraps cyanea-struct
├── phylo.rs        Wraps cyanea-phylo
├── io.rs           Wraps cyanea-io
├── omics.rs        Wraps cyanea-omics
└── core_utils.rs   Wraps cyanea-core (SHA-256, zstd)
```

## JSON Envelope Pattern

All WASM functions return `String` containing JSON. Three helper functions in `error.rs` provide the envelope:

- `wasm_ok(val)` — Serializes `val` via serde_json, wraps as `{"ok": ...}`
- `wasm_err(msg)` — Wraps error message as `{"error": "..."}`
- `wasm_result(r)` — Converts `Result<T, CyaneaError>` to either ok or error envelope

This pattern ensures JavaScript consumers always get parseable JSON regardless of success or failure.

## Function Binding Pattern

Each module file follows the same pattern:

1. Define `Js*` structs with `#[derive(Serialize)]` for return types
2. Write a public function that takes simple types (String, f64, etc.)
3. Parse any JSON input via `serde_json::from_str`
4. Call the corresponding domain crate function
5. Return via `wasm_ok(result)` or `wasm_result(r)`

All functions are annotated with `#[cfg_attr(feature = "wasm", wasm_bindgen)]` for conditional wasm-bindgen export.

## TypeScript Layer

The `ts/` directory contains:

- `ts/types.ts` — ~35 TypeScript interfaces matching `Js*` return structs
- `ts/index.ts` — Namespace wrappers (Seq, Align, Stats, ML, Chem, StructBio, Phylo, IO, Omics, Core) that parse JSON and provide typed APIs

## Stateless Design

All functions are pure — no global state, no initialization step. This makes them safe to call from multiple Web Workers concurrently.

## Error Mapping

`CyaneaError` variants are converted to error message strings. The JavaScript consumer sees only `{"error": "message"}` without variant information. For richer error handling, the message text encodes the error category (e.g., "parse error: ...", "invalid input: ...").

## Performance Considerations

- All data passes through JSON serialization, which adds overhead for large datasets
- For large distance matrices or PCA, consider using Web Workers to avoid blocking the main thread
- The WASM binary includes all 10 domain crates; tree-shaking is not possible at the WASM level
