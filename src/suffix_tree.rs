use std::collections::HashMap;
use std::fs::File;
use std::io::{Result as IoResult, Write};

// ------------------------ Suffix Tree Data Structures ---------------------------
#[derive(Clone)]
pub struct StEdge {
    pub start: usize,
    pub end: usize,
    pub target_node: usize,
}

pub struct StNode {
    pub id: usize,
    pub parent: Option<usize>,
    pub suffix_link: Option<usize>,
    pub children: Vec<(char, StEdge)>,
    pub depth: usize,
    pub suffix_start: Option<usize>,
}

pub struct SuffixTree {
    pub nodes: Vec<StNode>,
    pub root: usize,
}

// ------------------------ Implementation ---------------------------
impl SuffixTree {
    pub fn new() -> Self {
        let mut nodes = Vec::new();
        nodes.push(StNode {
            id: 0,
            parent: None,
            suffix_link: None,
            children: Vec::new(),
            depth: 0,
            suffix_start: None,
        });
        Self { nodes, root: 0 }
    }

    pub fn add_node(&mut self, parent: usize, depth: usize) -> usize {
        let id = self.nodes.len();
        self.nodes.push(StNode {
            id,
            parent: Some(parent),
            suffix_link: None,
            children: Vec::new(),
            depth,
            suffix_start: None,
        });
        id
    }

    pub fn node(&self, id: usize) -> &StNode {
        &self.nodes[id]
    }

    pub fn node_mut(&mut self, id: usize) -> &mut StNode {
        &mut self.nodes[id]
    }
}

// ------------------------------ Traversal and Operations  ------------------------------------

pub fn dfs_depths(tree: &SuffixTree, alphabet_order: &HashMap<char, usize>) -> Vec<usize> {
    let mut result = Vec::new();
    let mut stack = vec![tree.root];

    while let Some(node_id) = stack.pop() {
        result.push(tree.node(node_id).depth);

        let mut children = tree.node(node_id).children.clone();
        children.sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));

        for (_, edge) in children.into_iter().rev() {
            stack.push(edge.target_node);
        }
    }

    result
}

pub fn postorder_ids(tree: &SuffixTree, alphabet_order: &HashMap<char, usize>) -> Vec<usize> {
    let mut result = Vec::new();
    let mut stack = vec![tree.root];
    let mut visited = Vec::new();

    while let Some(node_id) = stack.pop() {
        visited.push(node_id);

        let mut children = tree.node(node_id).children.clone();
        children.sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));

        for (_, edge) in children {
            stack.push(edge.target_node);
        }
    }

    // Reverse for postorder
    result.extend(visited.into_iter().rev());
    result
}

pub fn compute_bwt(tree: &SuffixTree, s: &[u8]) -> Vec<u8> {
    // Ensure we have the complete sequence with $ terminator
    let seq = if s.last() == Some(&b'$') {
        s.to_vec()
    } else {
        let mut s_copy = s.to_vec();
        s_copy.push(b'$');
        s_copy
    };

    let mut suffixes: Vec<usize> = Vec::new();

    // Traverse all nodes to find leaves with their suffix start positions
    for node in &tree.nodes {
        if node.children.is_empty() && node.suffix_start.is_some() {
            suffixes.push(node.suffix_start.unwrap());
        }
    }

    // Sort suffixes lexicographically
    suffixes.sort_by(|&pos1, &pos2| {
        let suffix1 = &seq[pos1..];
        let suffix2 = &seq[pos2..];

        // Compare the suffixes character by character
        for (c1, c2) in suffix1.iter().zip(suffix2.iter()) {
            if c1 != c2 {
                return c1.cmp(c2);
            }
        }

        // If one is a prefix of the other, shorter comes first
        suffix1.len().cmp(&suffix2.len())
    });

    let mut bwt: Vec<u8> = Vec::with_capacity(seq.len());

    for &pos in &suffixes {
        if pos == 0 {
            bwt.push(b'$');
        } else {
            bwt.push(seq[pos - 1]);
        }
    }

    bwt
}

pub fn find_longest_repeat(tree: &SuffixTree, s: &[u8]) -> (usize, Vec<usize>) {
    let mut longest_length = 0;
    let mut start_positions = Vec::new();

    // For each internal node
    for (node_id, node) in tree.nodes.iter().enumerate() {
        // Skip leaves and root
        if node.children.is_empty() || node_id == tree.root {
            continue;
        }

        // Find all leaf descendants
        let mut leaf_positions = collect_leaf_positions(tree, node_id);

        if leaf_positions.len() >= 2 {
            let repeated_length = node.depth;

            if repeated_length > longest_length {
                longest_length = repeated_length;
                start_positions = leaf_positions;
            }
        }
    }

    (longest_length, start_positions)
}

// Helper function to collect all leaf positions under a node
fn collect_leaf_positions(tree: &SuffixTree, node_id: usize) -> Vec<usize> {
    let mut positions = Vec::new();
    let mut stack = vec![node_id];

    while let Some(current) = stack.pop() {
        let node = tree.node(current);

        if node.children.is_empty() {
            if let Some(pos) = node.suffix_start {
                positions.push(pos);
            }
        } else {
            for (_, edge) in &node.children {
                stack.push(edge.target_node);
            }
        }
    }

    positions
}
// -------------------------------------------------------------------------

// ------------------------------ IO ---------------------------------------

pub fn write_tree_stats(tree: &SuffixTree, fasta_file: &str, file: &mut File) -> IoResult<()> {
    let internal_count = tree
        .nodes
        .iter()
        .filter(|node| !node.children.is_empty())
        .count();
    let leaf_count = tree.nodes.len() - internal_count;
    let total_depth: usize = tree
        .nodes
        .iter()
        .filter(|node| !node.children.is_empty())
        .map(|node| node.depth)
        .sum();
    let max_depth = tree
        .nodes
        .iter()
        .filter(|node| !node.children.is_empty())
        .map(|node| node.depth)
        .max()
        .unwrap_or(0);

    writeln!(file, "Suffix Tree Statistics for {}", fasta_file)?;
    writeln!(file, "Internal nodes: {}", internal_count)?;
    writeln!(file, "Leaves: {}", leaf_count)?;
    writeln!(file, "Total nodes: {}", tree.nodes.len())?;
    writeln!(
        file,
        "Avg internal string depth: {:.2}",
        total_depth as f64 / internal_count.max(1) as f64
    )?;
    writeln!(file, "Max string depth: {}", max_depth)?;

    Ok(())
}

// Function to write DFS string depths to a file
pub fn dfs_to_file(
    tree: &SuffixTree,
    alphabet_order: &HashMap<char, usize>,
    file: &mut File,
) -> IoResult<()> {
    let dfs_depths = dfs_depths(tree, alphabet_order);
    crate::utils::write_vector_formatted(&dfs_depths, 10, file)
}

// Function to write postorder traversal to a file
pub fn postorder_to_file(
    tree: &SuffixTree,
    alphabet_order: &HashMap<char, usize>,
    file: &mut File,
) -> IoResult<()> {
    let postorder = postorder_ids(tree, alphabet_order);
    crate::utils::write_vector_formatted(&postorder, 10, file)
}

// Enhanced function to write BWT to a file with better error handling
pub fn bwt_to_file(tree: &SuffixTree, s: &[u8], file: &mut File) -> IoResult<()> {
    let bwt = compute_bwt(tree, s);

    for (i, &b) in bwt.iter().enumerate() {
        write!(file, "{}", b as char)?;
        if (i + 1) % 60 == 0 {
            writeln!(file)?;
        }
    }

    if bwt.len() % 60 != 0 {
        writeln!(file)?;
    }

    Ok(())
}

// Function to write longest repeat information to a file
pub fn longest_repeat_to_file(tree: &SuffixTree, seq: &[u8], file: &mut File) -> IoResult<()> {
    let (repeat_length, repeat_positions) = find_longest_repeat(tree, seq);

    writeln!(file, "\nLongest exact matching repeat:")?;
    writeln!(file, "Length: {}", repeat_length)?;
    writeln!(file, "Start positions: {:?}", repeat_positions)?;

    if repeat_length > 0
        && !repeat_positions.is_empty()
        && repeat_positions[0] + repeat_length <= seq.len()
    {
        if let Ok(repeat_str) =
            std::str::from_utf8(&seq[repeat_positions[0]..repeat_positions[0] + repeat_length])
        {
            writeln!(file, "Repeat sequence: {}", repeat_str)?;
        }
    }

    Ok(())
}

// ------------------------------------------------------------------
