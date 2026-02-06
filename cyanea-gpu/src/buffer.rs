//! Type-safe host/device memory buffer abstraction.

use cyanea_core::Summarizable;

use crate::backend::BackendKind;

/// A buffer representing data on a compute device (or host memory for CPU).
///
/// For the CPU backend, data lives in a `Vec<f64>`. For GPU backends,
/// `host_data` may be `None` until an explicit read-back is performed.
#[derive(Debug, Clone)]
pub struct Buffer {
    /// Host-side copy of the data (always present for CPU, lazy for GPU).
    pub(crate) host_data: Option<Vec<f64>>,
    /// Number of f64 elements.
    pub(crate) len: usize,
    /// Which backend created this buffer.
    pub(crate) origin: BackendKind,
}

impl Buffer {
    /// Creates a buffer backed by host data.
    pub(crate) fn from_host(data: Vec<f64>, origin: BackendKind) -> Self {
        let len = data.len();
        Self {
            host_data: Some(data),
            len,
            origin,
        }
    }

    /// Creates a buffer marker for device-only data (no host copy yet).
    #[allow(dead_code)]
    pub(crate) fn device_only(len: usize, origin: BackendKind) -> Self {
        Self {
            host_data: None,
            len,
            origin,
        }
    }

    /// Returns the number of `f64` elements in this buffer.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if the buffer has zero elements.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns which backend created this buffer.
    pub fn origin(&self) -> BackendKind {
        self.origin
    }

    /// Returns a reference to the host-side data, if available.
    ///
    /// For CPU buffers this always returns `Some`. For GPU buffers,
    /// call `Backend::read_buffer()` first to populate host data.
    pub fn as_host_slice(&self) -> Option<&[f64]> {
        self.host_data.as_deref()
    }

    /// Consumes the buffer and returns the host-side data, if available.
    pub fn into_host_vec(self) -> Option<Vec<f64>> {
        self.host_data
    }
}

impl Summarizable for Buffer {
    fn summary(&self) -> String {
        let host_status = if self.host_data.is_some() {
            "host-resident"
        } else {
            "device-only"
        };
        format!(
            "Buffer({} f64s, {}, {})",
            self.len, self.origin, host_status
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_host_round_trip() {
        let data = vec![1.0, 2.0, 3.0];
        let buf = Buffer::from_host(data.clone(), BackendKind::Cpu);
        assert_eq!(buf.len(), 3);
        assert!(!buf.is_empty());
        assert_eq!(buf.origin(), BackendKind::Cpu);
        assert_eq!(buf.as_host_slice(), Some(data.as_slice()));
    }

    #[test]
    fn into_host_vec() {
        let data = vec![4.0, 5.0];
        let buf = Buffer::from_host(data.clone(), BackendKind::Cpu);
        assert_eq!(buf.into_host_vec(), Some(data));
    }

    #[test]
    fn device_only_returns_none() {
        let buf = Buffer::device_only(10, BackendKind::Cuda);
        assert_eq!(buf.len(), 10);
        assert_eq!(buf.origin(), BackendKind::Cuda);
        assert!(buf.as_host_slice().is_none());
        assert!(buf.into_host_vec().is_none());
    }

    #[test]
    fn empty_buffer() {
        let buf = Buffer::from_host(vec![], BackendKind::Cpu);
        assert_eq!(buf.len(), 0);
        assert!(buf.is_empty());
        assert_eq!(buf.as_host_slice(), Some([].as_slice()));
    }

    #[test]
    fn summary_host_resident() {
        let buf = Buffer::from_host(vec![1.0, 2.0], BackendKind::Cpu);
        assert_eq!(buf.summary(), "Buffer(2 f64s, CPU, host-resident)");
    }

    #[test]
    fn summary_device_only() {
        let buf = Buffer::device_only(5, BackendKind::Metal);
        assert_eq!(buf.summary(), "Buffer(5 f64s, Metal, device-only)");
    }

    #[test]
    fn clone_independence() {
        let buf = Buffer::from_host(vec![1.0, 2.0, 3.0], BackendKind::Cpu);
        let mut cloned = buf.clone();
        // Mutating the clone's internal data doesn't affect the original.
        if let Some(ref mut data) = cloned.host_data {
            data[0] = 999.0;
        }
        assert_eq!(buf.as_host_slice().unwrap()[0], 1.0);
        assert_eq!(cloned.as_host_slice().unwrap()[0], 999.0);
    }
}
