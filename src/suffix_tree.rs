// suffix_tree.rs

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Result as IoResult, Write};

// Represents an edge in the suffix tree
#[derive(Clone, Debug)]
pub struct StEdge {
    pub start: usize,       // Start index in concatenated string
    pub end: usize,         // Inclusive end index
    pub target_node: usize, // Child node ID
}

// Represents a node in the suffix tree
#[derive(Debug)]
pub struct StNode {
    pub id: usize,
    pub parent: Option<usize>,
    pub suffix_link: Option<usize>,     // Used in construction algorithms (McCreight)
    pub children: Vec<(char, StEdge)>,  // Outgoing edges, sorted by first char
    pub depth: usize,                   // String depth from root
    pub suffix_start: Option<usize>,    // Start index if leaf node represents suffix[i..]
    pub origin_sequences: HashSet<usize>, // For GST: Which original sequences pass through node {0, 1}
}

// The Suffix Tree itself
#[derive(Debug)]
pub struct SuffixTree {
    pub nodes: Vec<StNode>,
    pub root: usize, // ID of the root
}

impl SuffixTree {
    pub fn new() -> Self {
        let mut nodes = Vec::new();
        nodes.push(StNode { // Node 0: Root
            id: 0,
            parent: None,
            suffix_link: None,
            children: Vec::new(),
            depth: 0,
            suffix_start: None,
            origin_sequences: HashSet::new(),
        });
        Self { nodes, root: 0 }
    }

    // Add a new node to the tree.
    pub fn add_node(&mut self, parent_id: usize, depth: usize) -> usize {
        let new_id = self.nodes.len();
        self.nodes.push(StNode {
            id: new_id,
            parent: Some(parent_id),
            suffix_link: None,
            children: Vec::new(),
            depth,
            suffix_start: None,
            origin_sequences: HashSet::new(),
        });
        new_id
    }

    // Get immutable ref to node.
    pub fn node(&self, id: usize) -> &StNode {
        &self.nodes[id]
    }

    // Get mutable ref to node.
    pub fn node_mut(&mut self, id: usize) -> &mut StNode {
        &mut self.nodes[id]
    }
}

// DFS traversal, returns node depths.
pub fn dfs_depths(tree: &SuffixTree, alphabet_order: &HashMap<char, usize>) -> Vec<usize> {
    let mut result = Vec::new();
    let mut stack = vec![tree.root];
    let mut visited = HashSet::new();

    while let Some(node_id) = stack.pop() {
        if !visited.insert(node_id) { continue; }
        result.push(tree.node(node_id).depth);

        // Get children sorted according to alphabet
        let mut children_to_visit: Vec<_> = tree.node(node_id).children.iter()
            .map(|(c, edge)| (alphabet_order.get(c).copied().unwrap_or(usize::MAX), edge.target_node))
            .collect();
        children_to_visit.sort_by_key(|(order, _)| *order);

        // Push onto stack in reverse sorted order for DFS
        for (_, child_id) in children_to_visit.into_iter().rev() {
            stack.push(child_id);
        }
    }
    result
}

// Post-order traversal, returns node IDs.
pub fn postorder_ids(tree: &SuffixTree, alphabet_order: &HashMap<char, usize>) -> Vec<usize> {
    let mut result = Vec::new();
    let mut stack1 = vec![tree.root]; // Main traversal stack
    let mut stack2 = Vec::new();      // Reverse post-order result stack

    while let Some(node_id) = stack1.pop() {
        stack2.push(node_id);

        // Get children sorted by alphabet
        let mut sorted_children = tree.node(node_id).children.clone();
        sorted_children.sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));

        // Push sorted children onto stack1
        for (_, edge) in sorted_children {
            stack1.push(edge.target_node);
        }
    }
    // Pop from stack2 to get final post-order
    result = stack2.into_iter().rev().collect();
    result
}

// Computes the Burrows-Wheeler Transform. Assumes s has $.
pub fn compute_bwt(tree: &SuffixTree, s: &[u8]) -> Vec<u8> {
    // Ensure s has $ and make a local copy if needed 
    let seq = if s.last() == Some(&b'$') {
        std::borrow::Cow::Borrowed(s) // Borrow if already ends with $
    } else {
        let mut s_copy = s.to_vec();
        s_copy.push(b'$');
        std::borrow::Cow::Owned(s_copy) // Own if $ was added
    };
    let n = seq.len();

    // Get suffix start positions from all leaves
    let mut suffix_starts: Vec<usize> = tree.nodes.iter()
        .filter_map(|node| if node.children.is_empty() { node.suffix_start } else { None })
        .collect();

    if suffix_starts.len() != n {
         eprintln!("Warning: BWT expected {} leaves, found {}. May be incomplete.", n, suffix_starts.len());
    }

    suffix_starts.sort_by_key(|&start_index| &seq[start_index..]);

    // Build BWT: for each sorted suffix start `pos`, take seq[pos-1] (or $ if pos=0)
    let bwt: Vec<u8> = suffix_starts.iter().map(|&pos| {
        if pos == 0 { seq[n - 1] } 
        else { seq[pos - 1] }
    }).collect();

    bwt
}

// Finds the longest repeated substring (deepest internal node).
pub fn find_longest_repeat(tree: &SuffixTree, _s: &[u8]) -> (usize, Vec<usize>) { // s not really needed if leaves are correct
    let mut longest_len = 0;
    let mut best_node_id = tree.root;

    // Find internal node with max depth
    for node_id in 1..tree.nodes.len() { // Skip root
        let node = tree.node(node_id);
        // Internal nodes represent repeats. Find deepest one.
        if !node.children.is_empty() && node.depth > longest_len {
            longest_len = node.depth;
            best_node_id = node_id;
        }
    }

    if longest_len == 0 { return (0, Vec::new()); } // No repeats

    // Collect leaves below the deepest internal node
    let start_positions = rake_leaves(tree, best_node_id);
    (longest_len, start_positions)
}


fn rake_leaves(tree: &SuffixTree, start_node_id: usize) -> Vec<usize> { //Lol
    let mut positions = Vec::new();
    let mut stack = vec![start_node_id];
    let mut visited = HashSet::new();

    while let Some(current_id) = stack.pop() {
        if !visited.insert(current_id) { continue; }
        let node = tree.node(current_id);

        if node.children.is_empty() { // Leaf
            if let Some(pos) = node.suffix_start { positions.push(pos); }
        } else { // Internal
            for (_, edge) in &node.children { stack.push(edge.target_node); }
        }
    }
    positions.sort();
    positions
}


// Calculates origin_sequences for GST nodes using post-order traversal.
// Assumes tree built from s0 + # + s1 + $.
pub fn compute_and_store_origin_info(tree: &mut SuffixTree, len_s0: usize) {
    // Manual post-order walk (doesn't need alphabet)
    let mut post_order_nodes = Vec::with_capacity(tree.nodes.len());
    let mut visited = vec![false; tree.nodes.len()];
    let mut stack = vec![(tree.root, false)]; // (node_id, children_processed)

    while let Some((node_id, children_processed)) = stack.pop() {
        if children_processed {
            // Visited this node after its children
            post_order_nodes.push(node_id);
        } else {
            // First time visiting: mark for reprocessing, push children
            if visited[node_id] { continue; }
            visited[node_id] = true;
            stack.push((node_id, true));

            let children_targets: Vec<usize> = tree.node(node_id).children.iter()
                                                 .map(|(_, edge)| edge.target_node)
                                                 .collect();

            for &child_id in children_targets.iter().rev() {
                if !visited[child_id] {
                    stack.push((child_id, false));
                }
            }
        }
    }

    // Process nodes in post-order
    for &node_id in &post_order_nodes {
        let is_leaf;
        let suffix_start;
        let children_targets: Vec<usize>;

        // Scope to make borrow checker happy...
        {
            let node = tree.node(node_id);
            is_leaf = node.children.is_empty();
            suffix_start = node.suffix_start;
            children_targets = node.children.iter().map(|(_, edge)| edge.target_node).collect();
        }

        if is_leaf {
            if let Some(start_pos) = suffix_start {
                let origin = if start_pos < len_s0 { 0 } // Belongs to s0
                             else if start_pos > len_s0 { 1 } // Belongs to s1
                             else { continue }; // Starts exactly at #, ignore? Yes.
                tree.node_mut(node_id).origin_sequences.insert(origin);
            }
        } else { // Internal Node
            // Union origin sets from children
            let mut combined_origins = HashSet::new();
            for child_id in children_targets {
                // Borrow child origins temporarily
                let child_origins = &tree.node(child_id).origin_sequences;
                combined_origins.extend(child_origins);
            }
            tree.node_mut(node_id).origin_sequences = combined_origins;
        }
    }
}

// Write basic tree stats to file.
pub fn write_tree_stats(_tree: &SuffixTree, fasta_file: &str, file: &mut File) -> IoResult<()> {

    let tree = _tree; // Use the passed tree
    let internal_count = tree.nodes.iter().filter(|n| !n.children.is_empty() && n.id != tree.root).count();
    let leaf_count = tree.nodes.iter().filter(|n| n.children.is_empty()).count();
    let total_nodes = tree.nodes.len();

    let mut internal_depth_sum: usize = 0;
    let mut max_internal_depth: usize = 0;
    let mut internal_node_count_for_avg: usize = 0;

    for node in &tree.nodes {
        if !node.children.is_empty() && node.id != tree.root { // Proper internal node
            internal_depth_sum += node.depth;
            max_internal_depth = max_internal_depth.max(node.depth);
            internal_node_count_for_avg += 1;
        }
    }

    writeln!(file, "--- Tree Stats ({}) ---", fasta_file)?;
    writeln!(file, "Total nodes: {}", total_nodes)?;
    writeln!(file, " Internal nodes: {}", internal_count)?;
    writeln!(file, " Leaf nodes: {}", leaf_count)?;
    if internal_node_count_for_avg > 0 {
        writeln!(file, " Avg internal depth: {:.2}", internal_depth_sum as f64 / internal_node_count_for_avg as f64)?;
        writeln!(file, " Max internal depth: {}", max_internal_depth)?;
    } else {
        writeln!(file, " No internal nodes found.")?;
    }
    Ok(())
}

// Write DFS node depths to file.
pub fn dfs_to_file(tree: &SuffixTree, alphabet_order: &HashMap<char, usize>, file: &mut File) -> IoResult<()> {
    writeln!(file, "--- DFS Node Depths ---")?;
    let depths = dfs_depths(tree, alphabet_order);
    crate::utils::write_vector_formatted(&depths, 10, file)
}

// Write Postorder node IDs to file.
pub fn postorder_to_file(tree: &SuffixTree, alphabet_order: &HashMap<char, usize>, file: &mut File) -> IoResult<()> {
    writeln!(file, "--- Postorder Node IDs ---")?;
    let ids = postorder_ids(tree, alphabet_order);
    crate::utils::write_vector_formatted(&ids, 10, file)
}

// Write BWT to file.
pub fn bwt_to_file(tree: &SuffixTree, s: &[u8], file: &mut File) -> IoResult<()> {
    writeln!(file, "--- BWT ---")?;
    let bwt = compute_bwt(tree, s);
    // Write BWT, wrap lines
    for (i, &byte) in bwt.iter().enumerate() {
        write!(file, "{}", byte as char)?;
        if (i + 1) % 60 == 0 && i != bwt.len() - 1 { writeln!(file)?; }
    }
    writeln!(file)?; // Final newline
    Ok(())
}

// Write longest repeat info to file.
pub fn longest_repeat_to_file(tree: &SuffixTree, seq: &[u8], file: &mut File) -> IoResult<()> {
    writeln!(file, "\n--- Longest Repeat ---")?;
    let (repeat_length, mut repeat_positions) = find_longest_repeat(tree, seq);

    writeln!(file, "Length: {}", repeat_length)?;
    repeat_positions.sort();
    writeln!(file, "Starts ({}): {:?}", repeat_positions.len(), repeat_positions)?;

    // Print the repeat sequence itself (if reasonable length)
    if repeat_length > 0 && !repeat_positions.is_empty() {
        let first_pos = repeat_positions[0];
        if first_pos + repeat_length <= seq.len() {
            let repeat_slice = &seq[first_pos .. first_pos + repeat_length];
            match std::str::from_utf8(repeat_slice) {
                 Ok(repeat_str) => {
                     const MAX_PRINT_LEN: usize = 100;
                     if repeat_length <= MAX_PRINT_LEN {
                         writeln!(file, "Sequence: {}", repeat_str)?;
                     } else {
                          writeln!(file, "Sequence (first {} chars): {}...", MAX_PRINT_LEN, &repeat_str[..MAX_PRINT_LEN])?;
                     }
                 },
                 Err(_) => writeln!(file, "Sequence: (non-UTF8 bytes)")?,
             }
        } else {
             writeln!(file, "Sequence: (invalid coords)")?;
         }
    }
    Ok(())
}