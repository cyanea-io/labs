//! Map `CyaneaError` variants to Python exceptions.

use cyanea_core::CyaneaError;
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::PyErr;

/// Extension trait for converting `Result<T, CyaneaError>` into `PyResult<T>`.
pub trait IntoPyResult<T> {
    fn into_pyresult(self) -> pyo3::PyResult<T>;
}

impl<T> IntoPyResult<T> for Result<T, CyaneaError> {
    fn into_pyresult(self) -> pyo3::PyResult<T> {
        self.map_err(to_pyerr)
    }
}

fn to_pyerr(e: CyaneaError) -> PyErr {
    match e {
        CyaneaError::Io(ref _inner) => PyIOError::new_err(e.to_string()),
        CyaneaError::Parse(_) | CyaneaError::InvalidInput(_) => {
            PyValueError::new_err(e.to_string())
        }
        CyaneaError::Compression(_) | CyaneaError::Hash(_) | CyaneaError::Other(_) => {
            PyRuntimeError::new_err(e.to_string())
        }
    }
}
