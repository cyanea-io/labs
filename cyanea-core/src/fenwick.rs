//! Generic Fenwick tree (Binary Indexed Tree) for prefix sums.
//!
//! Supports point updates and prefix/range sum queries in O(log n) time.

use core::ops::{Add, Sub};

/// A Fenwick tree (Binary Indexed Tree) supporting point updates and
/// prefix/range sum queries in O(log n) time.
///
/// Uses 1-based indexing internally. The generic type `T` must support
/// addition, subtraction, copy, and have a default (zero) value.
#[derive(Debug, Clone)]
pub struct FenwickTree<T> {
    tree: Vec<T>,
    n: usize,
}

impl<T: Copy + Default + Add<Output = T> + Sub<Output = T>> FenwickTree<T> {
    /// Create a Fenwick tree of size `n` initialized to zero.
    pub fn new(n: usize) -> Self {
        Self {
            tree: vec![T::default(); n + 1],
            n,
        }
    }

    /// Build a Fenwick tree from a slice in O(n) time.
    ///
    /// More efficient than creating an empty tree and calling `update` n times.
    pub fn from_slice(values: &[T]) -> Self {
        let n = values.len();
        let mut tree = vec![T::default(); n + 1];

        // Copy values into 1-indexed positions
        for (i, &v) in values.iter().enumerate() {
            tree[i + 1] = v;
        }

        // Build tree structure in O(n)
        for i in 1..=n {
            let parent = i + lowbit(i);
            if parent <= n {
                let child_val = tree[i];
                tree[parent] = tree[parent] + child_val;
            }
        }

        Self { tree, n }
    }

    /// Add `delta` to the element at index `i` (0-based).
    ///
    /// # Panics
    ///
    /// Panics if `i >= n`.
    pub fn update(&mut self, i: usize, delta: T) {
        assert!(i < self.n, "index out of bounds");
        let mut idx = i + 1;
        while idx <= self.n {
            self.tree[idx] = self.tree[idx] + delta;
            idx += lowbit(idx);
        }
    }

    /// Sum of elements in `[0, i]` (inclusive, 0-based).
    ///
    /// # Panics
    ///
    /// Panics if `i >= n`.
    pub fn prefix_sum(&self, i: usize) -> T {
        assert!(i < self.n, "index out of bounds");
        let mut idx = i + 1;
        let mut sum = T::default();
        while idx > 0 {
            sum = sum + self.tree[idx];
            idx -= lowbit(idx);
        }
        sum
    }

    /// Sum of elements in `[l, r]` (inclusive, 0-based).
    ///
    /// # Panics
    ///
    /// Panics if `l > r` or `r >= n`.
    pub fn range_sum(&self, l: usize, r: usize) -> T {
        assert!(l <= r, "l must be <= r");
        assert!(r < self.n, "index out of bounds");
        if l == 0 {
            self.prefix_sum(r)
        } else {
            self.prefix_sum(r) - self.prefix_sum(l - 1)
        }
    }

    /// Number of elements in the tree.
    pub fn len(&self) -> usize {
        self.n
    }

    /// Whether the tree is empty.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

/// Lowest set bit of `i` (i.e., `i & -i`).
#[inline]
fn lowbit(i: usize) -> usize {
    i & i.wrapping_neg()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_all_zeros() {
        let ft: FenwickTree<i64> = FenwickTree::new(5);
        assert_eq!(ft.prefix_sum(0), 0);
        assert_eq!(ft.prefix_sum(4), 0);
    }

    #[test]
    fn from_slice_prefix_sum() {
        let ft = FenwickTree::from_slice(&[1i64, 2, 3, 4, 5]);
        assert_eq!(ft.prefix_sum(0), 1);
        assert_eq!(ft.prefix_sum(1), 3);
        assert_eq!(ft.prefix_sum(2), 6);
        assert_eq!(ft.prefix_sum(3), 10);
        assert_eq!(ft.prefix_sum(4), 15);
    }

    #[test]
    fn update_and_query() {
        let mut ft: FenwickTree<i64> = FenwickTree::new(5);
        ft.update(0, 3);
        ft.update(1, 5);
        ft.update(2, 7);
        ft.update(3, 2);
        ft.update(4, 1);

        assert_eq!(ft.prefix_sum(0), 3);
        assert_eq!(ft.prefix_sum(2), 15);
        assert_eq!(ft.prefix_sum(4), 18);
    }

    #[test]
    fn range_sum() {
        let ft = FenwickTree::from_slice(&[1i64, 2, 3, 4, 5]);
        assert_eq!(ft.range_sum(0, 4), 15);
        assert_eq!(ft.range_sum(1, 3), 9); // 2+3+4
        assert_eq!(ft.range_sum(2, 2), 3); // just element 2
        assert_eq!(ft.range_sum(0, 0), 1);
    }

    #[test]
    fn update_then_range() {
        let mut ft = FenwickTree::from_slice(&[1i64, 2, 3, 4, 5]);
        ft.update(2, 10); // element 2 becomes 13
        assert_eq!(ft.range_sum(1, 3), 19); // 2+13+4
        assert_eq!(ft.prefix_sum(4), 25);
    }

    #[test]
    fn f64_values() {
        let ft = FenwickTree::from_slice(&[1.0_f64, 2.5, 3.0, 0.5]);
        assert!((ft.prefix_sum(3) - 7.0).abs() < 1e-10);
        assert!((ft.range_sum(1, 2) - 5.5).abs() < 1e-10);
    }

    #[test]
    fn negative_values() {
        let ft = FenwickTree::from_slice(&[-1i64, 3, -2, 5, -3]);
        assert_eq!(ft.prefix_sum(4), 2); // -1+3-2+5-3
        assert_eq!(ft.range_sum(1, 3), 6); // 3-2+5
    }

    #[test]
    fn single_element() {
        let ft = FenwickTree::from_slice(&[42i64]);
        assert_eq!(ft.prefix_sum(0), 42);
        assert_eq!(ft.range_sum(0, 0), 42);
        assert_eq!(ft.len(), 1);
    }

    #[test]
    fn len_and_is_empty() {
        let ft: FenwickTree<i64> = FenwickTree::new(0);
        assert!(ft.is_empty());
        assert_eq!(ft.len(), 0);

        let ft2 = FenwickTree::from_slice(&[1i64, 2]);
        assert!(!ft2.is_empty());
        assert_eq!(ft2.len(), 2);
    }

    #[test]
    fn large_tree() {
        let n = 1000;
        let values: Vec<i64> = (1..=n as i64).collect();
        let ft = FenwickTree::from_slice(&values);
        // Sum of 1..=1000 = 500500
        assert_eq!(ft.prefix_sum(n - 1), 500500);
        // Sum of 100..=200 = sum(1..=200) - sum(1..=99) = 20100 - 4950 = 15150
        assert_eq!(ft.range_sum(99, 199), 15150);
    }

    #[test]
    fn incremental_updates() {
        let mut ft: FenwickTree<i64> = FenwickTree::new(4);
        ft.update(0, 1);
        assert_eq!(ft.prefix_sum(3), 1);
        ft.update(1, 2);
        assert_eq!(ft.prefix_sum(3), 3);
        ft.update(2, 3);
        assert_eq!(ft.prefix_sum(3), 6);
        ft.update(3, 4);
        assert_eq!(ft.prefix_sum(3), 10);
    }
}
