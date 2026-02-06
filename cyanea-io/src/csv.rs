//! CSV/TSV parsing and preview.

use std::fs::File;
use std::path::Path;

use ::csv::ReaderBuilder;
use cyanea_core::{CyaneaError, Result};
use serde_json::{json, Value};

/// Metadata about a CSV file.
#[derive(Debug, Clone)]
#[cfg_attr(
    feature = "csv",
    derive(serde::Serialize, serde::Deserialize)
)]
pub struct CsvInfo {
    pub row_count: u64,
    pub column_count: usize,
    pub columns: Vec<String>,
    pub has_headers: bool,
}

/// Parse a CSV file and return its metadata.
pub fn parse_csv_info(path: impl AsRef<Path>) -> Result<CsvInfo> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let mut reader = ReaderBuilder::new().has_headers(true).from_reader(file);

    let headers = reader
        .headers()
        .map_err(|e| CyaneaError::Parse(e.to_string()))?;
    let columns: Vec<String> = headers.iter().map(|s| s.to_string()).collect();
    let column_count = columns.len();

    let mut row_count: u64 = 0;
    for result in reader.records() {
        result.map_err(|e| CyaneaError::Parse(e.to_string()))?;
        row_count += 1;
    }

    Ok(CsvInfo {
        row_count,
        column_count,
        columns,
        has_headers: true,
    })
}

/// Parse a CSV file and return the first `limit` rows as a JSON string.
pub fn csv_preview(path: impl AsRef<Path>, limit: usize) -> Result<String> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let mut reader = ReaderBuilder::new().has_headers(true).from_reader(file);

    let headers = reader
        .headers()
        .map_err(|e| CyaneaError::Parse(e.to_string()))?
        .clone();
    let columns: Vec<&str> = headers.iter().collect();

    let mut rows: Vec<Value> = Vec::new();

    for result in reader.records().take(limit) {
        let record = result.map_err(|e| CyaneaError::Parse(e.to_string()))?;
        let mut row = serde_json::Map::new();

        for (i, field) in record.iter().enumerate() {
            if let Some(col) = columns.get(i) {
                row.insert(col.to_string(), json!(field));
            }
        }

        rows.push(Value::Object(row));
    }

    serde_json::to_string(&json!({
        "columns": columns,
        "rows": rows
    }))
    .map_err(|e| CyaneaError::Parse(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_csv_parsing() {
        let mut file = NamedTempFile::with_suffix(".csv").unwrap();
        writeln!(file, "name,age,city").unwrap();
        writeln!(file, "Alice,30,NYC").unwrap();
        writeln!(file, "Bob,25,LA").unwrap();
        file.flush().unwrap();

        let info = parse_csv_info(file.path()).unwrap();
        assert_eq!(info.row_count, 2);
        assert_eq!(info.column_count, 3);
        assert_eq!(info.columns, vec!["name", "age", "city"]);
        assert!(info.has_headers);
    }

    #[test]
    fn test_csv_preview() {
        let mut file = NamedTempFile::with_suffix(".csv").unwrap();
        writeln!(file, "name,age").unwrap();
        writeln!(file, "Alice,30").unwrap();
        writeln!(file, "Bob,25").unwrap();
        writeln!(file, "Charlie,35").unwrap();
        file.flush().unwrap();

        let json_str = csv_preview(file.path(), 2).unwrap();
        let parsed: Value = serde_json::from_str(&json_str).unwrap();
        assert_eq!(parsed["rows"].as_array().unwrap().len(), 2);
        assert_eq!(parsed["columns"].as_array().unwrap().len(), 2);
    }

    #[test]
    fn test_csv_file_not_found() {
        let result = parse_csv_info("/nonexistent/file.csv");
        assert!(result.is_err());
    }
}
