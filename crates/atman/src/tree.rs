//! Average-linkage trees with bootstrap clade support, shared by
//! `axes tree` and `concordance --tree`.
//!
//! Conventions (matching scipy `linkage` and the Block N fixtures): leaves
//! are ids `0..n`, the i-th merge creates node `n + i`, a node's height is
//! the merge distance, a Newick branch length is `parent height − child
//! height` (leaves sit at height 0), and clade support is the fraction of
//! bootstrap trees containing the same leaf set, written as an integer
//! percentage label on the internal node.

use anyhow::Result;
use atman_core::modules_discover::upgma;
use std::collections::BTreeSet;
use std::path::Path;

use crate::io::{atomic_write, format_f64};

#[derive(Debug, Clone)]
pub struct Tree {
    pub labels: Vec<String>,
    /// `(left, right, height, size)` in merge order; node id = `labels.len() + index`.
    pub merges: Vec<(usize, usize, f64, usize)>,
}

impl Tree {
    /// Average-linkage tree over a symmetric dissimilarity matrix.
    pub fn from_dissimilarity(labels: Vec<String>, dissim: &[Vec<f64>]) -> Tree {
        let merges = upgma(dissim)
            .into_iter()
            .map(|(l, r, h, n)| (l.min(r), l.max(r), h, n))
            .collect();
        Tree { merges, labels }
    }

    pub fn n_leaves(&self) -> usize {
        self.labels.len()
    }

    /// Leaf indices under `node` (a leaf id or an internal id).
    pub fn leaf_set(&self, node: usize) -> BTreeSet<usize> {
        let n = self.n_leaves();
        if node < n {
            return [node].into_iter().collect();
        }
        let (l, r, _, _) = self.merges[node - n];
        let mut set = self.leaf_set(l);
        set.extend(self.leaf_set(r));
        set
    }

    fn height(&self, node: usize) -> f64 {
        let n = self.n_leaves();
        if node < n {
            0.0
        } else {
            self.merges[node - n].2
        }
    }

    /// Leaf sets of every internal node, in merge order.
    pub fn clades(&self) -> Vec<BTreeSet<usize>> {
        let n = self.n_leaves();
        (0..self.merges.len())
            .map(|i| self.leaf_set(n + i))
            .collect()
    }

    /// Newick string; `support[i]` labels internal node `n + i` as an
    /// integer percentage. Labels are sanitized (`·` and `/` → `-`,
    /// whitespace and Newick metacharacters → `_`) unless `quote` is set.
    pub fn newick(&self, support: &[Option<f64>], quote: bool) -> String {
        let n = self.n_leaves();
        if n == 0 {
            return ";".to_string();
        }
        if n == 1 {
            return format!("({});", leaf_label(&self.labels[0], quote));
        }
        let root = n + self.merges.len() - 1;
        let mut out = String::new();
        self.write_node(root, None, support, quote, &mut out);
        out.push(';');
        out
    }

    fn write_node(
        &self,
        node: usize,
        parent_height: Option<f64>,
        support: &[Option<f64>],
        quote: bool,
        out: &mut String,
    ) {
        let n = self.n_leaves();
        if node < n {
            out.push_str(&leaf_label(&self.labels[node], quote));
        } else {
            let (l, r, h, _) = self.merges[node - n];
            out.push('(');
            self.write_node(l, Some(h), support, quote, out);
            out.push(',');
            self.write_node(r, Some(h), support, quote, out);
            out.push(')');
            if let Some(Some(s)) = support.get(node - n) {
                out.push_str(&format!("{}", (s * 100.0).round() as i64));
            }
        }
        if let Some(ph) = parent_height {
            out.push(':');
            out.push_str(&branch_length(ph - self.height(node)));
        }
    }
}

/// Six decimals with trailing zeros trimmed (Newick is a figure input, not
/// a round-trip value column); tiny negatives from floating error become 0.
fn branch_length(v: f64) -> String {
    let v = if v.abs() < 5e-7 { 0.0 } else { v };
    let text = format!("{v:.6}");
    let trimmed = text.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-0" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

fn leaf_label(label: &str, quote: bool) -> String {
    if quote {
        format!("'{}'", label.replace('\'', "''"))
    } else {
        label
            .replace(['·', '/'], "-")
            .chars()
            .map(|c| {
                if c.is_whitespace() || "(),:;[]'\"".contains(c) {
                    '_'
                } else {
                    c
                }
            })
            .collect()
    }
}

/// Fraction of `replicates` that contain each internal clade of `reference`.
pub fn clade_support(reference: &Tree, replicates: &[Tree]) -> Vec<f64> {
    let ref_clades = reference.clades();
    if replicates.is_empty() {
        return vec![f64::NAN; ref_clades.len()];
    }
    let rep_clades: Vec<BTreeSet<BTreeSet<usize>>> = replicates
        .iter()
        .map(|t| t.clades().into_iter().collect())
        .collect();
    ref_clades
        .iter()
        .map(|c| {
            rep_clades.iter().filter(|set| set.contains(c)).count() as f64 / replicates.len() as f64
        })
        .collect()
}

/// scipy-style linkage table: `left, right, distance, n`.
pub fn write_linkage(path: &Path, tree: &Tree) -> Result<()> {
    let mut out = String::from("left\tright\tdistance\tn\n");
    for (l, r, h, size) in &tree.merges {
        out.push_str(&format!("{}\t{}\t{}\t{}\n", l, r, format_f64(*h), size));
    }
    atomic_write(path, out.as_bytes())
}

/// `node` (sorted leaf labels joined by ` | `), `n_leaves`, `support`, plus
/// the internal node id and height for cross-reference with the linkage.
pub fn write_support(path: &Path, tree: &Tree, support: &[f64], n_boot: usize) -> Result<()> {
    let mut out = String::from("node\tn_leaves\tsupport\tnode_id\theight\tn_boot\n");
    let n = tree.n_leaves();
    for (i, clade) in tree.clades().iter().enumerate() {
        let names: Vec<&str> = clade.iter().map(|&k| tree.labels[k].as_str()).collect();
        let mut sorted = names.clone();
        sorted.sort();
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            sorted.join(" | "),
            clade.len(),
            if support[i].is_finite() {
                format_f64(support[i])
            } else {
                String::new()
            },
            n + i,
            format_f64(tree.merges[i].2),
            n_boot
        ));
    }
    atomic_write(path, out.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tree4() -> Tree {
        // A–B close (0.1), C–D close (0.2), the two pairs far apart (1.0).
        let d = vec![
            vec![0.0, 0.1, 1.0, 1.0],
            vec![0.1, 0.0, 1.0, 1.0],
            vec![1.0, 1.0, 0.0, 0.2],
            vec![1.0, 1.0, 0.2, 0.0],
        ];
        Tree::from_dissimilarity(vec!["A".into(), "B".into(), "C".into(), "D".into()], &d)
    }

    #[test]
    fn upgma_tree_merges_clades_and_writes_newick() {
        let t = tree4();
        assert_eq!(t.merges.len(), 3);
        assert_eq!(t.merges[0].2, 0.1);
        assert_eq!(t.merges[1].2, 0.2);
        assert!((t.merges[2].2 - 1.0).abs() < 1e-12);
        let clades = t.clades();
        assert_eq!(clades[0], [0usize, 1].into_iter().collect());
        assert_eq!(clades[1], [2usize, 3].into_iter().collect());
        assert_eq!(clades[2].len(), 4);
        let support = vec![Some(0.9), Some(0.5), Some(1.0)];
        assert_eq!(
            t.newick(&support, false),
            "((A:0.1,B:0.1)90:0.9,(C:0.2,D:0.2)50:0.8)100;"
        );
    }

    #[test]
    fn clade_support_counts_replicates_containing_the_clade() {
        let t = tree4();
        // A replicate where B and C pair instead.
        let d = vec![
            vec![0.0, 1.0, 1.0, 0.1],
            vec![1.0, 0.0, 0.1, 1.0],
            vec![1.0, 0.1, 0.0, 1.0],
            vec![0.1, 1.0, 1.0, 0.0],
        ];
        let other = Tree::from_dissimilarity(t.labels.clone(), &d);
        let s = clade_support(&t, &[t.clone(), other]);
        assert_eq!(s, vec![0.5, 0.5, 1.0]);
    }

    #[test]
    fn labels_are_sanitized_or_quoted() {
        assert_eq!(
            leaf_label("CNS autoimmune / NMOSD · MS PILOT", false),
            "CNS_autoimmune_-_NMOSD_-_MS_PILOT"
        );
        assert_eq!(leaf_label("CSF leak, lateral", true), "'CSF leak, lateral'");
    }
}
