use crate::alignment::align_variant_affine;
use crate::mccreight;
use crate::suffix_tree::{self, SuffixTree};
use crate::utils;
use std::collections::{HashMap, HashSet};
use std::time::{Duration, Instant};

/// Calculates similarity for (s_i, s_j) using GST.
pub fn calculate_similarity(
    seq_i: &[u8],
    seq_j: &[u8],
    alphabet_order: &HashMap<char, usize>,
    term1: u8, // # separator
    term2: u8, // $ terminator
    match_s: i32,
    mismatch: i32,
    gap_open: i32,
    gap_extend: i32,
) -> Result<(usize, usize, Duration, Duration), String> {

    // Step 1) Build the GST
    let gst_start = Instant::now();
    let combined_seq = [seq_i, &[term1], seq_j, &[term2]].concat();
    let len_s0 = seq_i.len(); // Length of s_i
    let mut tree = mccreight::construct_st_mc(&combined_seq, alphabet_order);
    suffix_tree::compute_and_store_origin_info(&mut tree, len_s0);
    let gst_elapsed = gst_start.elapsed();

    // Step 2) Use GST to find LCS('b')
    let lcs_result_coords = find_lcs_via_gst_from_tree(&tree, len_s0)?;

    // x0 (start in s_i)
    // y0 (end in s_i)
    // x1 (start in s_j)
    // y1 (end in s_j)
    let (b_lcs, x0, y0, x1, y1) = match lcs_result_coords {
        Some(coords) => coords,
        // No common substring, b=0
        None => (0, 0, usize::MAX, 0, usize::MAX),
    };
    let b = b_lcs; // LCS length

    // Start recording time
    let alignment_start_time = Instant::now();

    // Step 3) Prefix Alignment (Reversed) to get a
    // a) Extract prefixes s_i[1...x0] and s_j[1...x1]
    let prefix_i = if x0 > 0 { &seq_i[0..x0] } else { &[] };
    let prefix_j = if x1 > 0 { &seq_j[0..x1] } else { &[] };

    // b) Reverse prefixes
    let a = if prefix_i.is_empty() && prefix_j.is_empty() {
        0 // No prefix alignment needed
    } else {
        let prefix_i_rev = utils::reverse_sequence(prefix_i);
        let prefix_j_rev = utils::reverse_sequence(prefix_j);
        // c) Compute variant global alignment
        align_variant_affine(&prefix_i_rev, &prefix_j_rev, match_s, mismatch, gap_open, gap_extend)?
    };

    // Step 4) Suffix Alignment to get c
    // a) Extract suffixes s_i[y0+1...m] and s_j[y1+1...n]
    let suffix_i_start = if y0 == usize::MAX { seq_i.len() } else { y0 + 1 };
    let suffix_j_start = if y1 == usize::MAX { seq_j.len() } else { y1 + 1 };
    let suffix_i = if suffix_i_start < seq_i.len() { &seq_i[suffix_i_start..] } else { &[] };
    let suffix_j = if suffix_j_start < seq_j.len() { &seq_j[suffix_j_start..] } else { &[] };

    // b) Compute variant global alignment
    let c = if suffix_i.is_empty() && suffix_j.is_empty() {
        0 // No suffix alignment needed
    } else {
        align_variant_affine(suffix_i, suffix_j, match_s, mismatch, gap_open, gap_extend)?
    };
    let alignment_elapsed = alignment_start_time.elapsed();

    // Step 5) Calculate Similarity = a + b + c
    let total_similarity = a + b + c;

    Ok((total_similarity, b, gst_elapsed, alignment_elapsed))
}


// Finds LCS coords via GST
fn find_lcs_via_gst_from_tree(
    tree: &SuffixTree,
    len_s0: usize,
) -> Result<Option<(usize, usize, usize, usize, usize)>, String> {

    // Find deepest internal node u covering origins {0, 1}.
    let mut best_lcs_node_id: Option<usize> = None;
    let mut max_depth = 0;
    let required_origins: HashSet<usize> = [0, 1].iter().cloned().collect();

    for node_id in 0..tree.nodes.len() {
        let node = tree.node(node_id);
        if !node.children.is_empty()
            && node.id != tree.root
            && node.origin_sequences == required_origins
            && node.depth > max_depth
        {
            max_depth = node.depth;
            best_lcs_node_id = Some(node_id);
        }
    }

    if best_lcs_node_id.is_none() { return Ok(None); } // No LCS

    let lcs_node_id = best_lcs_node_id.unwrap();
    let b = max_depth; // LCS length

    // Find representative leaves for start positions
    let mut rep_s0_start: Option<usize> = None;
    let mut rep_s1_start: Option<usize> = None;

    let mut stack = vec![lcs_node_id];
    let mut visited_for_coords = HashSet::new();

    while let Some(curr_id) = stack.pop() {
        if !visited_for_coords.insert(curr_id) { continue; }
        let node = tree.node(curr_id);

        if node.children.is_empty() { // Leaf
            if let Some(leaf_start_pos) = node.suffix_start {

                if leaf_start_pos < len_s0 { 
                    // From s_i
                    if rep_s0_start.is_none() { rep_s0_start = Some(leaf_start_pos); }
                } else if leaf_start_pos > len_s0 { 
                    // From s_j
                    if rep_s1_start.is_none() {
                        rep_s1_start = Some(leaf_start_pos - (len_s0 + 1));
                    }
                }

            }
        } else { // Internal
            for (_, edge) in node.children.iter() {
                let child_node = tree.node(edge.target_node);
                let needs_s0 = rep_s0_start.is_none() && child_node.origin_sequences.contains(&0);
                let needs_s1 = rep_s1_start.is_none() && child_node.origin_sequences.contains(&1);
                if needs_s0 || needs_s1 { stack.push(edge.target_node); }
            }
        }

        // If both are found
        if rep_s0_start.is_some() && rep_s1_start.is_some() { break; }
    }

    // Calculate end coordinates (y0, y1)
    match (rep_s0_start, rep_s1_start) {
        (Some(x0), Some(x1)) => {
            if b == 0 { return Err(format!("LCS length b=0 but found leaves from node {}", lcs_node_id)); }
            let y0 = x0 + b - 1; // End index in s_i
            let y1 = x1 + b - 1; // End index in s_j

            Ok(Some((b, x0, y0, x1, y1))) // Return (length, s_i start, s_i end, s_j start, s_j end)
        }
        _ => {
            Err(format!("Couldn't find rep leaves for both origins below LCS node {}", lcs_node_id))
        }
    }
}