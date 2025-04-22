// mccreight.rs

use crate::suffix_tree::{StEdge, SuffixTree};
use std::collections::HashMap;

/* McCreight's Suffix Tree Construction */
// Add suffixes one by one using suffix links and edge splits.

pub fn construct_st_mc(seq: &[u8], alphabet_order: &HashMap<char, usize>) -> SuffixTree {
    let mut s_vec = seq.to_vec();
    if s_vec.last() != Some(&b'$') {
        s_vec.push(b'$');
    } // Ensure terminator
    let s = &s_vec[..];
    let n = s.len();

    let mut tree = SuffixTree::new();
    let root = tree.root;
    tree.node_mut(root).suffix_link = Some(root); // Root SL -> root

    add_first_suffix(&mut tree, s, alphabet_order);

    let first_leaf_id = 1;
    let node1_is_leaf =
        tree.nodes.len() > first_leaf_id && tree.node(first_leaf_id).children.is_empty();

    // leaf_ids[i-1] stores node ID where Suffix(i-1) ended.
    let mut leaf_ids: Vec<usize> = vec![first_leaf_id];

    for i in 1..n {

        // Find ancestor of previous leaf's node with a suffix link (or root)
        let u = find_valid_u(&tree, leaf_ids[i - 1]);

        // Find v = SL(u) (follow or compute)
        let v = follow_or_calc_link(&mut tree, s, i, u);

        // Walk down from v matching s[i..], split edge if needed
        let parent_for_new_leaf = descend_and_split(&mut tree, &s, v, i, alphabet_order);

        // Determine where new leaf edge starts in s
        let (_, matched_len_from_v, _) = find_path(&tree, s, v, i);
        let leaf_edge_start_index = i + matched_len_from_v;

        let current_suffix_node_id;

        if leaf_edge_start_index >= n {
            let (w_end_node, _, _) = find_path(&tree, s, v, i);
            current_suffix_node_id = w_end_node; // Use the existing node
        } else {
            let new_leaf_id = add_leaf(
                &mut tree,
                s,
                parent_for_new_leaf,
                leaf_edge_start_index,
                n,
                i,
                alphabet_order,
            );
            current_suffix_node_id = new_leaf_id;
        }
        leaf_ids.push(current_suffix_node_id);
    }

    return tree;
}

fn add_first_suffix(tree: &mut SuffixTree, s: &[u8], alphabet_order: &HashMap<char, usize>) {
    let n = s.len();
    let parent_depth = tree.node(tree.root).depth;
    let edge_len = n;
    let leaf_depth = parent_depth + edge_len;

    let leaf_id = tree.add_node(tree.root, leaf_depth);
    tree.node_mut(leaf_id).suffix_start = Some(0); // Leaf for Suffix(0)

    let first_char = s[0] as char;
    let edge = StEdge {
        start: 0,
        end: n - 1,
        target_node: leaf_id,
    };

    // Add edge to root
    let children = &mut tree.node_mut(tree.root).children;
    children.push((first_char, edge));
    children.sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));
}

// Find ancestor with SL or root.
fn find_valid_u(tree: &SuffixTree, node_representing_suffix_i_1: usize) -> usize {
    let mut current_node_id = node_representing_suffix_i_1;

    if let Some(parent_id) = tree.node(current_node_id).parent {
        if tree.node(parent_id).suffix_link.is_some() {
            return parent_id;
        }
        current_node_id = parent_id;
    } else {
        return tree.root;
    }

    // Climb until SL found or we reach root
    while tree.node(current_node_id).suffix_link.is_none() && current_node_id != tree.root {
        current_node_id = tree.node(current_node_id).parent.unwrap_or(tree.root);
    }
    return current_node_id;
}

// Follow SL(u) if known, otherwise compute it.
fn follow_or_calc_link(tree: &mut SuffixTree, s: &[u8], i: usize, u: usize) -> usize {
    if u == tree.root {
        return tree.root;
    }
    if let Some(v) = tree.node(u).suffix_link {
        return v;
    }

    // SL(u)
    let u_parent_id;
    let u_depth;
    let u_parent_depth;
    let beta_start;

    {
        // Read u properties
        let u_node = tree.node(u);
        u_parent_id = u_node.parent.expect("Non-root must have parent");
        u_depth = u_node.depth;
        let parent_node = tree.node(u_parent_id);
        u_parent_depth = parent_node.depth;
        let edge_to_u = parent_node
            .children
            .iter()
            .find(|(_, edge)| edge.target_node == u)
            .expect("Edge to u not found");
        beta_start = edge_to_u.1.start;
    }

    // Recursively compute v' = SL(parent(u))
    let v_prime = follow_or_calc_link(tree, s, i - 1, u_parent_id); // Note: i-1

    // Hop down from v' using edge label beta
    let beta_len = u_depth - u_parent_depth;
    let v = node_hops(tree, v_prime, s, beta_start, beta_len);

    tree.node_mut(u).suffix_link = Some(v); // Cache SL(u) = v
    return v;
}

// Descend from v matching s[i], split if needed. Returns parent for new leaf.
fn descend_and_split(
    tree: &mut SuffixTree,
    s: &[u8],
    v: usize, // Start node
    i: usize, // Suffix index
    alphabet_order: &HashMap<char, usize>,
) -> usize {
    let (w_end_node, _, split_info) = find_path(tree, s, v, i);

    match split_info {
        // Might be able to get rid of edge start now
        Some((edge_char, edge_start, split_idx, original_target_id)) => {

            let original_edge = tree
                .node(w_end_node)
                .children
                .iter()
                .find(|(c, edge)| *c == edge_char && edge.target_node == original_target_id)
                .expect("Split info but edge not found")
                .1
                .clone();

            if original_edge.start > original_edge.end {
                eprintln!(
                    "CRITICAL ERROR: Original edge start={} > end={}",
                    original_edge.start, original_edge.end
                );
                panic!("Corrupted edge before split.");
            }

            let match_len_on_edge = (split_idx + 1).saturating_sub(original_edge.start);

            let parent_depth = tree.node(w_end_node).depth;
            let mid_node_depth = parent_depth + match_len_on_edge;
            let mid_node_id = tree.add_node(w_end_node, mid_node_depth);

            // Populate new edges
            let mid_to_target_start = split_idx + 1;
            let mid_to_target_end = original_edge.end;
            let parent_to_mid_start = original_edge.start;
            let parent_to_mid_end = split_idx;

            tree.node_mut(original_target_id).parent = Some(mid_node_id);

            if mid_to_target_start <= mid_to_target_end {
                let edge_mid_to_target = StEdge {
                    start: mid_to_target_start,
                    end: mid_to_target_end,
                    target_node: original_target_id,
                };
                let next_char = s[mid_to_target_start] as char;
                let mid_children = &mut tree.node_mut(mid_node_id).children;
                mid_children.push((next_char, edge_mid_to_target));
                mid_children
                    .sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));
            }

            let edge_parent_to_mid = StEdge {
                start: parent_to_mid_start,
                end: parent_to_mid_end,
                target_node: mid_node_id,
            };
            let parent_children = &mut tree.node_mut(w_end_node).children;
            parent_children
                .retain(|(c, edge)| !(*c == edge_char && edge.target_node == original_target_id));
            parent_children.push((edge_char, edge_parent_to_mid));
            parent_children
                .sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));

            return mid_node_id; // New leaf parent is the mid-node
        }
        None => {
            // No split needed, match ended at w_end_node.
            return w_end_node; // New leaf parent is w_end_node
        }
    }
}

fn add_leaf(
    tree: &mut SuffixTree,
    s: &[u8],
    parent_id: usize,
    edge_start: usize,
    s_len: usize,
    suffix_start_index: usize,
    alphabet_order: &HashMap<char, usize>,
) -> usize {
    let parent_depth = tree.node(parent_id).depth;
    let edge_len = s_len - edge_start;
    let leaf_depth = parent_depth + edge_len;

    let leaf_id = tree.add_node(parent_id, leaf_depth);
    tree.node_mut(leaf_id).suffix_start = Some(suffix_start_index);

    let first_char = s[edge_start] as char;
    let edge = StEdge {
        start: edge_start,
        end: s_len - 1,
        target_node: leaf_id,
    };

    if edge.start > edge.end {
        eprintln!(
            "ERROR: add_leaf invalid edge start={} > end={}",
            edge.start, edge.end
        );
        panic!("Invalid edge in add_leaf.");
    }

    // Add edge to parent, sort children
    let parent_children = &mut tree.node_mut(parent_id).children;
    parent_children.push((first_char, edge));
    parent_children.sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));

    return leaf_id;
}

// Follow path from start_node_id.
fn find_path(
    tree: &SuffixTree,
    s: &[u8],
    start_node_id: usize,
    i: usize, // Pattern start index in s
) -> (usize, usize, Option<(char, usize, usize, usize)>) {
    let mut current_node_id = start_node_id;
    let mut total_matched_from_start = 0;

    loop {
        let pattern_char_idx = i + total_matched_from_start;
        if pattern_char_idx >= s.len() {
            return (current_node_id, total_matched_from_start, None);
        } // Whole pattern matched
        let pattern_char = s[pattern_char_idx] as char;

        // Find edge starting with pattern_char
        let edge_info = tree
            .node(current_node_id)
            .children
            .iter()
            .find(|(label_char, _)| *label_char == pattern_char);

        if edge_info.is_none() {
            return (current_node_id, total_matched_from_start, None);
        }

        let (_, edge) = edge_info.unwrap();
        let edge_start = edge.start;
        let edge_end = edge.end;
        let edge_target_id = edge.target_node;
        let edge_first_char = s[edge_start] as char;

        // Compare with edge
        let mut matched_on_edge = 0;
        loop {
            let pattern_idx = pattern_char_idx + matched_on_edge;
            let edge_idx = edge_start + matched_on_edge;

            if edge_idx > edge_end {
                break;
            } 
            if pattern_idx >= s.len() {
                break;
            }

            if s[pattern_idx] != s[edge_idx] {
                let split_end_idx = edge_idx - 1;
                let split_info = (edge_first_char, edge_start, split_end_idx, edge_target_id);
                return (
                    current_node_id,
                    total_matched_from_start + matched_on_edge,
                    Some(split_info),
                );
            }
            matched_on_edge += 1;
        }

        total_matched_from_start += matched_on_edge;

        let edge_consumed = (edge_start + matched_on_edge) > edge_end;
        let pattern_consumed = (pattern_char_idx + matched_on_edge) >= s.len();

        if edge_consumed {
            current_node_id = edge_target_id; // Move to target node
        } else if pattern_consumed {
            let split_end_idx = edge_start + matched_on_edge - 1;
            let split_info = (edge_first_char, edge_start, split_end_idx, edge_target_id);
            return (current_node_id, total_matched_from_start, Some(split_info));
        } else {
            panic!("Unexpected state in find_path edge comparison.");
        }
    }
}

// Traverse using string segment for SL calculation. Splits edges.
fn node_hops(
    tree: &mut SuffixTree,
    start_node_id: usize,
    s: &[u8],
    i: usize,   // Segment start index
    len: usize, // Segment length
) -> usize {
    let mut current_node_id = start_node_id;
    let mut current_s_idx = i;
    let mut remaining_len = len;

    if remaining_len == 0 {
        return start_node_id;
    }

    while remaining_len > 0 {
        if current_s_idx >= s.len() {
            // Index check
            eprintln!("ERROR: node_hops index out of bounds.");
            return current_node_id;
        }
        let next_char = s[current_s_idx] as char;

        // Find edge, clone for split logic
        let edge_info = tree
            .node(current_node_id)
            .children
            .iter()
            .find(|(c, _)| *c == next_char)
            .map(|(_, e)| e.clone());

        let edge = match edge_info {
            Some(e) => e,
            None => {
                eprintln!("ERROR: node_hops path broken.");
                return current_node_id; // Error upstream
            }
        };

        let edge_len = (edge.end + 1).saturating_sub(edge.start);

        if remaining_len < edge_len {
            // Hop ends mid-edge -> Split
            let parent_id = current_node_id;
            let original_target_id = edge.target_node;

            // Create mid node
            let parent_depth = tree.node(parent_id).depth;
            let mid_node_depth = parent_depth + remaining_len;
            let mid_node_id = tree.add_node(parent_id, mid_node_depth);

            let mid_to_target_start = edge.start + remaining_len;
            let mid_to_target_end = edge.end;
            let parent_to_mid_start = edge.start;
            let parent_to_mid_end = edge.start + remaining_len - 1;

            tree.node_mut(original_target_id).parent = Some(mid_node_id);

            let edge_mid_to_target = StEdge {
                start: mid_to_target_start,
                end: mid_to_target_end,
                target_node: original_target_id,
            };
            let mid_target_first_char = s[mid_to_target_start] as char;
            tree.node_mut(mid_node_id)
                .children
                .push((mid_target_first_char, edge_mid_to_target));

            // Replace parent
            let edge_parent_to_mid = StEdge {
                start: parent_to_mid_start,
                end: parent_to_mid_end,
                target_node: mid_node_id,
            };
            let parent_children = &mut tree.node_mut(parent_id).children;
            parent_children
                .retain(|(c, e)| !(*c == next_char && e.target_node == original_target_id));
            parent_children.push((next_char, edge_parent_to_mid));

            return mid_node_id;
        } else {
            // remaining_len >= edge_len
            // Hop covers this entire edge
            current_node_id = edge.target_node;
            current_s_idx += edge_len;
            remaining_len -= edge_len;
        }
    }

    return current_node_id;
}
