//! Static augmented interval tree for fast overlap queries.
//!
//! [`IntervalTree`] stores intervals in an implicit balanced BST layout
//! (nodes in a contiguous `Vec`, children of node `i` at `2i+1`/`2i+2`).
//! Build once, query many times with O(log n + k) overlap queries.

/// A generic interval with associated data.
#[derive(Debug, Clone)]
pub struct Interval<T> {
    /// Start coordinate (inclusive).
    pub start: u64,
    /// End coordinate (exclusive).
    pub end: u64,
    /// Payload associated with this interval.
    pub data: T,
}

impl<T> Interval<T> {
    /// Create a new interval.
    pub fn new(start: u64, end: u64, data: T) -> Self {
        Self { start, end, data }
    }
}

/// Internal node in the implicit BST.
#[derive(Debug, Clone)]
struct Node<T> {
    interval: Interval<T>,
    /// Maximum end coordinate in this subtree.
    max_end: u64,
}

/// A static augmented interval tree using an implicit BST layout.
///
/// Built once from a set of intervals, then supports efficient overlap queries.
/// The tree cannot be modified after construction.
#[derive(Debug, Clone)]
pub struct IntervalTree<T> {
    nodes: Vec<Option<Node<T>>>,
}

impl<T> IntervalTree<T> {
    /// Build an interval tree from unsorted intervals. O(n log n).
    pub fn from_unsorted(mut intervals: Vec<Interval<T>>) -> Self {
        intervals.sort_by_key(|iv| iv.start);
        Self::from_sorted(intervals)
    }

    /// Build an interval tree from intervals sorted by start coordinate. O(n).
    pub fn from_sorted(intervals: Vec<Interval<T>>) -> Self {
        let n = intervals.len();
        if n == 0 {
            return Self { nodes: Vec::new() };
        }

        // Compute the required array size for the implicit BST
        let capacity = implicit_tree_size(n);
        let mut nodes: Vec<Option<Node<T>>> = (0..capacity).map(|_| None).collect();

        // Convert intervals into an indexable vec
        let mut sorted: Vec<Option<Interval<T>>> = intervals.into_iter().map(Some).collect();

        build_implicit(&mut nodes, &mut sorted, 0, 0, n);
        augment_max_end(&mut nodes, 0);

        Self { nodes }
    }

    /// Query all intervals overlapping the range `[start, end)`.
    ///
    /// Returns references to all intervals where `interval.start < end && interval.end > start`.
    pub fn query(&self, start: u64, end: u64) -> Vec<&Interval<T>> {
        let mut results = Vec::new();
        if !self.nodes.is_empty() {
            self.query_recursive(0, start, end, &mut results);
        }
        results
    }

    /// Count intervals overlapping the range `[start, end)` without allocating.
    pub fn count_overlaps(&self, start: u64, end: u64) -> usize {
        if self.nodes.is_empty() {
            return 0;
        }
        self.count_recursive(0, start, end)
    }

    /// Find the nearest interval to a point.
    ///
    /// Returns the interval whose midpoint is closest to `point`.
    /// If multiple intervals are equidistant, returns one arbitrarily.
    pub fn nearest(&self, point: u64) -> Option<&Interval<T>> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best: Option<&Interval<T>> = None;
        let mut best_dist = u64::MAX;
        self.nearest_recursive(0, point, &mut best, &mut best_dist);
        best
    }

    /// Find the nearest interval that ends at or before `point`.
    ///
    /// Returns the interval with the largest `end` that is `<= point`.
    pub fn preceding(&self, point: u64) -> Option<&Interval<T>> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best: Option<&Interval<T>> = None;
        self.preceding_recursive(0, point, &mut best);
        best
    }

    /// Find the nearest interval that starts at or after `point`.
    ///
    /// Returns the interval with the smallest `start` that is `>= point`.
    pub fn following(&self, point: u64) -> Option<&Interval<T>> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best: Option<&Interval<T>> = None;
        self.following_recursive(0, point, &mut best);
        best
    }

    /// Number of intervals in the tree.
    pub fn len(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_some()).count()
    }

    /// Whether the tree contains no intervals.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty() || self.nodes.iter().all(|n| n.is_none())
    }

    /// Iterate over all intervals in the tree (in-order traversal).
    pub fn iter(&self) -> impl Iterator<Item = &Interval<T>> {
        IntervalTreeIter {
            nodes: &self.nodes,
            stack: if self.nodes.is_empty() {
                Vec::new()
            } else {
                vec![IterState::Descend(0)]
            },
        }
    }

    fn query_recursive<'a>(
        &'a self,
        idx: usize,
        start: u64,
        end: u64,
        results: &mut Vec<&'a Interval<T>>,
    ) {
        if idx >= self.nodes.len() {
            return;
        }
        let node = match &self.nodes[idx] {
            Some(n) => n,
            None => return,
        };

        // Prune: if max_end in this subtree <= query start, no overlap possible
        if node.max_end <= start {
            return;
        }

        // Search left subtree
        let left = 2 * idx + 1;
        self.query_recursive(left, start, end, results);

        // Check current node
        if node.interval.start < end && node.interval.end > start {
            results.push(&node.interval);
        }

        // Prune right: if node.start >= end, right subtree has only larger starts
        if node.interval.start < end {
            let right = 2 * idx + 2;
            self.query_recursive(right, start, end, results);
        }
    }

    fn count_recursive(&self, idx: usize, start: u64, end: u64) -> usize {
        if idx >= self.nodes.len() {
            return 0;
        }
        let node = match &self.nodes[idx] {
            Some(n) => n,
            None => return 0,
        };

        if node.max_end <= start {
            return 0;
        }

        let mut count = 0;

        let left = 2 * idx + 1;
        count += self.count_recursive(left, start, end);

        if node.interval.start < end && node.interval.end > start {
            count += 1;
        }

        if node.interval.start < end {
            let right = 2 * idx + 2;
            count += self.count_recursive(right, start, end);
        }

        count
    }

    fn nearest_recursive<'a>(
        &'a self,
        idx: usize,
        point: u64,
        best: &mut Option<&'a Interval<T>>,
        best_dist: &mut u64,
    ) {
        if idx >= self.nodes.len() {
            return;
        }
        let node = match &self.nodes[idx] {
            Some(n) => n,
            None => return,
        };

        // Distance from point to this interval
        let dist = if point < node.interval.start {
            node.interval.start - point
        } else if point >= node.interval.end {
            point - node.interval.end + 1
        } else {
            0 // point is inside the interval
        };

        if dist < *best_dist {
            *best_dist = dist;
            *best = Some(&node.interval);
        }

        if dist == 0 {
            return; // Can't do better than overlapping
        }

        let left = 2 * idx + 1;
        let right = 2 * idx + 2;

        // Search both subtrees
        if point < node.interval.start {
            self.nearest_recursive(left, point, best, best_dist);
            if node.interval.start - point <= *best_dist {
                self.nearest_recursive(right, point, best, best_dist);
            }
        } else {
            self.nearest_recursive(right, point, best, best_dist);
            self.nearest_recursive(left, point, best, best_dist);
        }
    }

    fn preceding_recursive<'a>(
        &'a self,
        idx: usize,
        point: u64,
        best: &mut Option<&'a Interval<T>>,
    ) {
        if idx >= self.nodes.len() {
            return;
        }
        let node = match &self.nodes[idx] {
            Some(n) => n,
            None => return,
        };

        if node.interval.end <= point {
            // This interval ends before point — candidate
            let is_better = match best {
                None => true,
                Some(b) => node.interval.end > b.end
                    || (node.interval.end == b.end && node.interval.start > b.start),
            };
            if is_better {
                *best = Some(&node.interval);
            }
        }

        let left = 2 * idx + 1;
        let right = 2 * idx + 2;

        // Always check left (may have intervals ending before point)
        self.preceding_recursive(left, point, best);
        // Check right subtree too (intervals may end before point but start after current)
        self.preceding_recursive(right, point, best);
    }

    fn following_recursive<'a>(
        &'a self,
        idx: usize,
        point: u64,
        best: &mut Option<&'a Interval<T>>,
    ) {
        if idx >= self.nodes.len() {
            return;
        }
        let node = match &self.nodes[idx] {
            Some(n) => n,
            None => return,
        };

        if node.interval.start >= point {
            // This interval starts at or after point — candidate
            let is_better = match best {
                None => true,
                Some(b) => node.interval.start < b.start,
            };
            if is_better {
                *best = Some(&node.interval);
            }
        }

        let left = 2 * idx + 1;
        let right = 2 * idx + 2;

        // If current node starts after point, left subtree may have closer intervals
        if node.interval.start >= point {
            self.following_recursive(left, point, best);
        }
        // Always check right subtree
        self.following_recursive(right, point, best);
    }
}

// ---------------------------------------------------------------------------
// Implicit BST construction helpers
// ---------------------------------------------------------------------------

/// Compute the array size needed for an implicit BST with `n` elements.
fn implicit_tree_size(n: usize) -> usize {
    if n == 0 {
        return 0;
    }
    // Height of the tree
    let height = (n as f64).log2().ceil() as u32 + 1;
    (1usize << height) - 1
}

/// Recursively build the implicit BST by placing the median at each node.
fn build_implicit<T>(
    nodes: &mut [Option<Node<T>>],
    sorted: &mut [Option<Interval<T>>],
    node_idx: usize,
    lo: usize,
    hi: usize,
) {
    if lo >= hi || node_idx >= nodes.len() {
        return;
    }

    let mid = lo + (hi - lo) / 2;

    if let Some(interval) = sorted[mid].take() {
        let max_end = interval.end;
        nodes[node_idx] = Some(Node {
            interval,
            max_end,
        });

        let left = 2 * node_idx + 1;
        let right = 2 * node_idx + 2;

        build_implicit(nodes, sorted, left, lo, mid);
        build_implicit(nodes, sorted, right, mid + 1, hi);
    }
}

/// Post-order traversal to compute augmented max_end values.
fn augment_max_end<T>(nodes: &mut [Option<Node<T>>], idx: usize) -> u64 {
    if idx >= nodes.len() {
        return 0;
    }

    let node = match &nodes[idx] {
        Some(n) => n,
        None => return 0,
    };

    let own_end = node.interval.end;
    let left_max = augment_max_end(nodes, 2 * idx + 1);
    let right_max = augment_max_end(nodes, 2 * idx + 2);

    let max_end = own_end.max(left_max).max(right_max);

    if let Some(ref mut n) = nodes[idx] {
        n.max_end = max_end;
    }

    max_end
}

// ---------------------------------------------------------------------------
// Iterator
// ---------------------------------------------------------------------------

enum IterState {
    Descend(usize),
    Visit(usize),
}

struct IntervalTreeIter<'a, T> {
    nodes: &'a [Option<Node<T>>],
    stack: Vec<IterState>,
}

impl<'a, T> Iterator for IntervalTreeIter<'a, T> {
    type Item = &'a Interval<T>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let state = self.stack.pop()?;
            match state {
                IterState::Descend(idx) => {
                    if idx >= self.nodes.len() {
                        continue;
                    }
                    if self.nodes[idx].is_none() {
                        continue;
                    }
                    // Push right, then visit, then left (so left is processed first)
                    self.stack.push(IterState::Descend(2 * idx + 2));
                    self.stack.push(IterState::Visit(idx));
                    self.stack.push(IterState::Descend(2 * idx + 1));
                }
                IterState::Visit(idx) => {
                    if let Some(node) = &self.nodes[idx] {
                        return Some(&node.interval);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(start: u64, end: u64) -> Interval<()> {
        Interval::new(start, end, ())
    }

    fn iv_data(start: u64, end: u64, data: usize) -> Interval<usize> {
        Interval::new(start, end, data)
    }

    #[test]
    fn empty_tree() {
        let tree: IntervalTree<()> = IntervalTree::from_unsorted(vec![]);
        assert!(tree.is_empty());
        assert_eq!(tree.len(), 0);
        assert_eq!(tree.query(0, 100).len(), 0);
        assert_eq!(tree.count_overlaps(0, 100), 0);
        assert!(tree.nearest(50).is_none());
        assert!(tree.preceding(50).is_none());
        assert!(tree.following(50).is_none());
        assert_eq!(tree.iter().count(), 0);
    }

    #[test]
    fn single_interval() {
        let tree = IntervalTree::from_unsorted(vec![iv(10, 20)]);
        assert_eq!(tree.len(), 1);
        assert!(!tree.is_empty());

        assert_eq!(tree.query(5, 15).len(), 1);
        assert_eq!(tree.query(15, 25).len(), 1);
        assert_eq!(tree.query(10, 20).len(), 1);
        assert_eq!(tree.query(0, 10).len(), 0); // abutting
        assert_eq!(tree.query(20, 30).len(), 0); // abutting
        assert_eq!(tree.query(25, 30).len(), 0);
    }

    #[test]
    fn many_intervals() {
        let tree = IntervalTree::from_unsorted(vec![
            iv(0, 10),
            iv(5, 15),
            iv(20, 30),
            iv(25, 35),
            iv(50, 60),
        ]);
        assert_eq!(tree.len(), 5);

        // Query overlapping first two
        let hits = tree.query(8, 12);
        assert_eq!(hits.len(), 2);

        // Query overlapping middle two
        let hits = tree.query(22, 28);
        assert_eq!(hits.len(), 2);

        // Query overlapping none (gap)
        let hits = tree.query(40, 45);
        assert_eq!(hits.len(), 0);

        // Query overlapping all in first cluster
        let hits = tree.query(0, 35);
        assert_eq!(hits.len(), 4);
    }

    #[test]
    fn nested_intervals() {
        let tree = IntervalTree::from_unsorted(vec![
            iv(0, 100),
            iv(10, 90),
            iv(20, 80),
            iv(30, 70),
            iv(40, 60),
        ]);

        // Point query in center should hit all
        assert_eq!(tree.query(45, 55).len(), 5);

        // Point query at edge
        assert_eq!(tree.query(0, 1).len(), 1);
        assert_eq!(tree.query(95, 100).len(), 1);
    }

    #[test]
    fn adjacent_intervals() {
        let tree = IntervalTree::from_unsorted(vec![
            iv(0, 10),
            iv(10, 20),
            iv(20, 30),
        ]);

        // Abutting intervals don't overlap in half-open semantics
        assert_eq!(tree.query(10, 20).len(), 1);
        assert_eq!(tree.query(9, 11).len(), 2);
    }

    #[test]
    fn all_same_start() {
        let tree = IntervalTree::from_unsorted(vec![
            iv(10, 20),
            iv(10, 30),
            iv(10, 40),
            iv(10, 50),
        ]);

        assert_eq!(tree.query(10, 11).len(), 4);
        assert_eq!(tree.query(25, 26).len(), 3);
        assert_eq!(tree.query(35, 36).len(), 2);
        assert_eq!(tree.query(45, 46).len(), 1);
    }

    #[test]
    fn count_overlaps() {
        let tree = IntervalTree::from_unsorted(vec![
            iv(0, 10),
            iv(5, 15),
            iv(20, 30),
        ]);
        assert_eq!(tree.count_overlaps(8, 12), 2);
        assert_eq!(tree.count_overlaps(25, 35), 1);
        assert_eq!(tree.count_overlaps(16, 19), 0);
    }

    #[test]
    fn nearest_basic() {
        let tree = IntervalTree::from_unsorted(vec![
            iv(10, 20),
            iv(30, 40),
            iv(60, 70),
        ]);

        // Point inside an interval
        let n = tree.nearest(15).unwrap();
        assert_eq!(n.start, 10);

        // Point between intervals — closer to [30,40)
        let n = tree.nearest(28).unwrap();
        assert_eq!(n.start, 30);

        // Point before all intervals
        let n = tree.nearest(0).unwrap();
        assert_eq!(n.start, 10);

        // Point after all intervals
        let n = tree.nearest(100).unwrap();
        assert_eq!(n.start, 60);
    }

    #[test]
    fn preceding_basic() {
        let tree = IntervalTree::from_unsorted(vec![
            iv(10, 20),
            iv(30, 40),
            iv(60, 70),
        ]);

        // Before first interval
        assert!(tree.preceding(5).is_none());

        // After first interval
        let p = tree.preceding(25).unwrap();
        assert_eq!(p.start, 10);

        // After second interval
        let p = tree.preceding(50).unwrap();
        assert_eq!(p.start, 30);

        // After all intervals
        let p = tree.preceding(100).unwrap();
        assert_eq!(p.start, 60);
    }

    #[test]
    fn following_basic() {
        let tree = IntervalTree::from_unsorted(vec![
            iv(10, 20),
            iv(30, 40),
            iv(60, 70),
        ]);

        // Before first interval
        let f = tree.following(0).unwrap();
        assert_eq!(f.start, 10);

        // Between intervals
        let f = tree.following(25).unwrap();
        assert_eq!(f.start, 30);

        // At interval start
        let f = tree.following(30).unwrap();
        assert_eq!(f.start, 30);

        // After all intervals
        assert!(tree.following(75).is_none());
    }

    #[test]
    fn preceding_at_boundary() {
        let tree = IntervalTree::from_unsorted(vec![iv(10, 20)]);

        // End == point: preceding should find it (end <= point)
        let p = tree.preceding(20).unwrap();
        assert_eq!(p.start, 10);

        // End > point: not preceding
        assert!(tree.preceding(15).is_none());
    }

    #[test]
    fn following_at_boundary() {
        let tree = IntervalTree::from_unsorted(vec![iv(10, 20)]);

        // start == point
        let f = tree.following(10).unwrap();
        assert_eq!(f.start, 10);

        // start < point
        assert!(tree.following(15).is_none());
    }

    #[test]
    fn iter_in_order() {
        let tree = IntervalTree::from_unsorted(vec![
            iv(30, 40),
            iv(10, 20),
            iv(50, 60),
            iv(0, 5),
        ]);

        let starts: Vec<u64> = tree.iter().map(|i| i.start).collect();
        // In-order traversal should yield sorted by start
        assert_eq!(starts, vec![0, 10, 30, 50]);
    }

    #[test]
    fn from_sorted() {
        let sorted = vec![iv(0, 10), iv(10, 20), iv(20, 30)];
        let tree = IntervalTree::from_sorted(sorted);
        assert_eq!(tree.len(), 3);
        assert_eq!(tree.query(5, 25).len(), 3); // all three overlap [5, 25)
        assert_eq!(tree.query(5, 15).len(), 2); // [0,10) and [10,20)
    }

    #[test]
    fn data_preserved() {
        let tree = IntervalTree::from_unsorted(vec![
            iv_data(10, 20, 42),
            iv_data(30, 40, 99),
        ]);

        let hits = tree.query(15, 35);
        assert_eq!(hits.len(), 2);
        let mut data: Vec<usize> = hits.iter().map(|h| h.data).collect();
        data.sort();
        assert_eq!(data, vec![42, 99]);
    }

    #[test]
    fn large_tree() {
        let intervals: Vec<Interval<usize>> = (0..1000)
            .map(|i| iv_data(i * 10, i * 10 + 5, i as usize))
            .collect();
        let tree = IntervalTree::from_unsorted(intervals);
        assert_eq!(tree.len(), 1000);

        // Query a small range
        let hits = tree.query(500, 510);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].data, 50);

        // Query a wider range
        let hits = tree.query(0, 10000);
        assert_eq!(hits.len(), 1000);
    }

    #[test]
    fn query_matches_linear_scan() {
        // Property test: tree query results should match linear scan
        let intervals = vec![
            iv(5, 15),
            iv(10, 25),
            iv(20, 35),
            iv(30, 45),
            iv(40, 55),
            iv(0, 100),
            iv(50, 60),
            iv(70, 80),
        ];

        let tree = IntervalTree::from_unsorted(intervals.clone());

        for start in (0..100).step_by(7) {
            for end in (start + 1..110).step_by(11) {
                let tree_count = tree.count_overlaps(start, end);
                let linear_count = intervals
                    .iter()
                    .filter(|iv| iv.start < end && iv.end > start)
                    .count();
                assert_eq!(
                    tree_count, linear_count,
                    "mismatch for query [{}, {}): tree={}, linear={}",
                    start, end, tree_count, linear_count
                );
            }
        }
    }

    #[test]
    fn two_intervals() {
        let tree = IntervalTree::from_unsorted(vec![iv(0, 10), iv(20, 30)]);
        assert_eq!(tree.len(), 2);
        assert_eq!(tree.query(5, 25).len(), 2);
        assert_eq!(tree.query(5, 15).len(), 1);
        assert_eq!(tree.query(25, 35).len(), 1);
        assert_eq!(tree.query(12, 18).len(), 0);
    }
}
