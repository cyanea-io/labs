//! JSON result wrapper for WASM boundary.
//!
//! Every public function in `cyanea-wasm` returns a `String` containing JSON.
//! Success → `{"ok": <value>}`, failure → `{"error": "<message>"}`.

use std::fmt::Display;

use serde::Serialize;

/// Serialize a success value as `{"ok": val}`.
pub fn wasm_ok<T: Serialize>(val: &T) -> String {
    #[derive(Serialize)]
    struct Ok<'a, T: Serialize> {
        ok: &'a T,
    }
    serde_json::to_string(&Ok { ok: val }).unwrap_or_else(|e| wasm_err(e))
}

/// Serialize an error as `{"error": "msg"}`.
pub fn wasm_err(msg: impl Display) -> String {
    #[derive(Serialize)]
    struct Err {
        error: String,
    }
    // This serialization cannot practically fail (it's a plain string).
    serde_json::to_string(&Err {
        error: msg.to_string(),
    })
    .unwrap_or_else(|_| r#"{"error":"serialization failed"}"#.into())
}

/// Map a `cyanea_core::Result<T>` into the JSON envelope.
pub fn wasm_result<T: Serialize>(r: cyanea_core::Result<T>) -> String {
    match r {
        Ok(val) => wasm_ok(&val),
        Err(e) => wasm_err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ok_serialization() {
        let json = wasm_ok(&42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], 42);
    }

    #[test]
    fn err_serialization() {
        let json = wasm_err("something broke");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["error"], "something broke");
    }

    #[test]
    fn result_ok_variant() {
        let r: cyanea_core::Result<&str> = Ok("hello");
        let json = wasm_result(r);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], "hello");
    }

    #[test]
    fn result_err_variant() {
        let r: cyanea_core::Result<i32> =
            Err(cyanea_core::CyaneaError::InvalidInput("bad".into()));
        let json = wasm_result(r);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].as_str().unwrap().contains("bad"));
    }
}
