use std::collections::HashMap;
use crate::suffix_tree::{SuffixTree, StEdge};

/* McCreight's Algorithm:
   This algorithm builds a suffix tree in linear time O(n) by using suffix links.
   
   The algorithm works as follows:
   1. Start by adding the first suffix (the entire string) to the tree
   2. For each subsequent suffix:
      a. Find a valid node u (a node with a suffix link or the root)
      b. Follow or compute the suffix link to find node v
      c. Descend from v and potentially split an edge if needed
      d. Add a leaf representing the current suffix
   
   The key insight of McCreight's algorithm is to avoid redundant work by using
   suffix links to navigate between related positions in the tree, rather than
   starting each new suffix insertion from the root.
*/

pub fn construct_st_mc(seq: &[u8], alphabet_order: &HashMap<char, usize>) -> SuffixTree {
    let mut s = seq.to_vec();
    // Add the $ terminator if its not already there
    if s.last() != Some(&b'$') {
        s.push(b'$');
    }

    let mut tree = SuffixTree::new();
    let root = tree.root;

    // Set roots suffix link to itself
    tree.node_mut(root).suffix_link = Some(root);

    // Add the first suffix (the entire string)
    add_first_suffix(&mut tree, &s, alphabet_order);

    let mut leaf_ids: Vec<usize> = vec![tree.node(root).children.first().unwrap().1.target_node];

    // Set to track added suffixes
    let mut added_suffixes = vec![0];

    // Add remaining suffixes
    for i in 1..s.len() {
        let u = find_valid_u(&tree, leaf_ids[i - 1]);
        let v = follow_or_compute_suffix_link(&mut tree, &s, i, u);
        let parent = descend_and_split(&mut tree, &s, v, i, alphabet_order);
        let leaf = add_leaf(&mut tree, &s, parent, i, s.len(), i, alphabet_order);
        leaf_ids.push(leaf);
        added_suffixes.push(i);
    }

    return tree;
}

// Adds the first suffix to an empty tree
fn add_first_suffix(tree: &mut SuffixTree, s: &[u8], alphabet_order: &HashMap<char, usize>) {
    let leaf = tree.add_node(tree.root, s.len());
    tree.node_mut(leaf).suffix_start = Some(0);

    let ch = s[0] as char;
    let edge = StEdge { start: 0, end: s.len() - 1, target_node: leaf };
    let children = &mut tree.node_mut(tree.root).children;
    children.push((ch, edge));
    children.sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));
}

// Finds a node with a valid suffix link by traversing up from the previous leaf
fn find_valid_u(tree: &SuffixTree, leaf_i_1: usize) -> usize {
    let mut u = tree.node(leaf_i_1).parent.expect("Leaf should have a parent");
    
    // Follow parent links until we find a node with a suffix link
    while tree.node(u).suffix_link.is_none() && u != tree.root {
        u = tree.node(u).parent.unwrap_or(tree.root);
    }
    
    return u;
}

// Either follows an existing suffix link or computes a new one
fn follow_or_compute_suffix_link(tree: &mut SuffixTree, s: &[u8], i: usize, u: usize) -> usize {
    let u_node = tree.node(u);
    
    if let Some(sl) = u_node.suffix_link {
        // Simply follow the suffix link if it exists
        return sl;
    } else if u == tree.root {
        // Root is its own suffix link
        return tree.root;
    } else {
        let u_prime = tree.node(u).parent.expect("u should have a parent if not root");
        let u_prime_node = tree.node(u_prime);
        let v_prime = u_prime_node.suffix_link.unwrap_or(tree.root);

        // Calculate beta length
        let beta_len = u_node.string_depth - u_prime_node.string_depth;

        // Determine the correct start index
        let start_index;
        if u_prime == tree.root {
            // If u' is the root, beta starts at the beginning of the current suffix
            start_index = i - beta_len + 1;
        } else {
            // Otherwise, beta starts after alpha in the current suffix
            start_index = i - beta_len;
        }

        let v = node_hops(tree, v_prime, s, start_index, beta_len);
        
        // Set the computed suffix link
        tree.node_mut(u).suffix_link = Some(v);
        
        return v;
    }
}

// Descends from node v and splits an edge if necessary
fn descend_and_split(
    tree: &mut SuffixTree,
    s: &[u8],
    v: usize,
    i: usize,
    alphabet_order: &HashMap<char, usize>,
) -> usize {
    let (w, _, split) = find_path(tree, s, v, i);

    match split {
        Some((ch, start, split_end, target)) => {
            let original_edge = tree.node(v)
                .children
                .iter()
                .find(|(c, _)| *c == ch)
                .unwrap()
                .1
                .clone();

            let original_end = original_edge.end;

            // If we've matched the entire edge, return the target
            if split_end + 1 > original_end {
                return target;
            }

            // Create a new internal node
            let parent_depth = tree.node(v).string_depth;
            let mid_depth = parent_depth + (split_end - start + 1);
            let mid = tree.add_node(v, mid_depth);

            // Update parent of the target
            tree.node_mut(target).parent = Some(mid);

            // Create the edge from mid to target
            let edge_mid_to_target = StEdge {
                start: split_end + 1,
                end: original_end,
                target_node: target,
            };

            tree.node_mut(mid).children.push((s[split_end + 1] as char, edge_mid_to_target));
            tree.node_mut(mid).children.sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));

            // Create the edge from v to mid
            let edge_v_to_mid = StEdge {
                start,
                end: split_end,
                target_node: mid,
            };

            // Update v's children
            let v_children = &mut tree.node_mut(v).children;
            v_children.retain(|(c, _)| *c != ch);
            v_children.push((ch, edge_v_to_mid));
            v_children.sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));

            return mid;
        }
        None => return w,
    }
}

// Adds a leaf node to the tree representing a suffix
fn add_leaf(
    tree: &mut SuffixTree,
    s: &[u8],
    parent: usize,
    start: usize,
    end: usize,
    suffix_start: usize,
    alphabet_order: &HashMap<char, usize>
) -> usize {
    // Create the new leaf node
    let leaf = tree.add_node(parent, end - start);
    tree.node_mut(leaf).suffix_start = Some(suffix_start);
    
    // Get the character for the edge label
    let ch = s[start] as char;

    // Create the edge from parent to leaf
    let edge = StEdge {
        start,
        end: end - 1,
        target_node: leaf,
    };

    // Add the edge to the parents children list
    let children = &mut tree.node_mut(parent).children;

    // Remove any duplicate edges
    let existing_edge = children.iter().position(|(c, e)| 
        *c == ch && e.start == start && e.target_node == leaf);
        
    if let Some(pos) = existing_edge {
        children.remove(pos);
    }

    // Add the new edge and sort the children
    children.push((ch, edge));
    children.sort_by_key(|(c, _)| alphabet_order.get(c).copied().unwrap_or(usize::MAX));

    return leaf;
}

// Finds the path in the tree for a given suffix and identifies where to split if needed
fn find_path(
    tree: &SuffixTree,
    s: &[u8],
    mut node_id: usize,
    mut i: usize,
) -> (usize, usize, Option<(char, usize, usize, usize)>) {
    let mut total_matched = 0;

    while i < s.len() {
        let ch = s[i] as char;

        // Find the edge for the current character
        let edge = match tree.node(node_id).children.iter().find(|(label, _)| *label == ch) {
            Some((_, e)) => e,
            None => break, // No matching edge found
        };

        // Compare the current suffix with the edge label
        let edge_slice = &s[edge.start..=edge.end];
        let mut j = 0;

        while i < s.len() && j < edge_slice.len() {
            if s[i] != edge_slice[j] {
                // Split on mismatch
                if j == 0 {
                    return (node_id, total_matched, None);
                }

                let split_start = edge.start;
                let split_end = edge.start + j - 1;
                
                let matched = total_matched + j;
                return (
                    node_id,
                    matched,
                    Some((ch, split_start, split_end, edge.target_node)),
                );
            }
            i += 1;
            j += 1;
        }

        total_matched += edge_slice.len();
        node_id = edge.target_node;
    }

    return (node_id, total_matched, None);
}

// Performs node hops (walking the tree) for a given string segment
fn node_hops(tree: &mut SuffixTree, mut node_id: usize, s: &[u8], mut i: usize, mut len: usize) -> usize {
    // Handle zero length
    if len == 0 {
        return node_id;
    }

    while len > 0 {
        let ch = s[i] as char;
        
        // Find the matching edge
        let edge = match tree.node(node_id).children.iter().find(|(c, _)| *c == ch) {
            Some((_, e)) => e.clone(),
            None => {
                return node_id;
            }
        };

        // Calculate how much of the edge we can walk
        let edge_len = edge.end - edge.start + 1;

        if len == edge_len {
            // Perfect match follow the edge completely
            node_id = edge.target_node;
            i += edge_len;
            len = 0;
        } else if len < edge_len {
            // Need to split the edge
            let mid = tree.add_node(node_id, tree.node(node_id).string_depth + len);

            // Edge from mid to the original target
            let edge_mid_to_target = StEdge {
                start: edge.start + len,
                end: edge.end,
                target_node: edge.target_node,
            };

            // Edge from the original source to mid
            let edge_parent_to_mid = StEdge {
                start: edge.start,
                end: edge.start + len - 1,
                target_node: mid,
            };

            // Update the connections
            tree.node_mut(mid).children.push((s[edge.start + len] as char, edge_mid_to_target));
            tree.node_mut(edge.target_node).parent = Some(mid);

            let children = &mut tree.node_mut(node_id).children;
            children.retain(|(c, _)| *c != ch);
            children.push((ch, edge_parent_to_mid));

            return mid;
        } else {
            // Edge is shorter than the remaining length - follow it and continue
            node_id = edge.target_node;
            i += edge_len;
            len -= edge_len;
        }
    }

    return node_id;
}