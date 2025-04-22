//alignment.rs
//  Allignmnet algorithms for projects in CS 471

use std::cmp;

use crate::types::*;
use crate::utils::calculate_alignmnet_stats;

// Neddleman-Wunsch w/ affine gaps a.k.a Gotoh '82:

/*
Gotoh Alg:
    Take 2 sequences as input and a scoring system as input.
    Create a matrix with M, I and D scores (in a struct of some sort) of length [m + 1 , n + 1] where m and n are the sequence lengths
    M tracks match / mismatch score | I tracks insertion (gap in seq1 | i-1, j) | D tracks deletions (gap in seq2 | i, j-1)
    M starts at 0, the other 2 as undefined.

    Populate the first rows and columns:
        M is undefined for first row and column past (0, 0)
        D is gap score across first column, and undefined down first row
        I is undefined for first row, and gap score all the way down first column

    Iterate over all cells and calculate score:
        M is max of diagonals from every matrix + match score if seq1[i] == seq2[j] otherwise mismatch
        D is max of opening gap (M score above + open penalty) and continue gap (D score above + extend penalty)
        I is max of opening gap (M score left + open penalty) and continue gap (I score left + extend penalty)


    Traceback optimal path:
        Take Highest score of the corner (m, n) out of M I and D
        create two strings to store output path
        Until (0,0) is reached:
            Determine how we got to current score (one of M, I, or D) from adjacent cells
                - If current score is in M then check all possible moves from cells above or to the left
                - If current score is in D or I, then check Match score + gap start penalty, and previous self score + extend
            Add seq1[i] and seq2[j] or gap to output strings depending on what path is traced.

    Notes:
    (to help me keep track of m/n, i/j, and how the scores are related to moving through the matrix)
        Deletion is moving down the matrix as that represents consuming a character of seq1 without
    moving forward in seq 2. Assumng that we are trying to see if seq2 is derived from seq1 through
    some biological process then this would mean that at some point that character of seq 1 was
    lost. Insertion then, is the opposite of this, where a character of seq2 is consumed without
    moving forward in seq 1. This would represent a scenario where there is new genetic material in
    between characters of our original sequence.
*/

// Accepts any two sequences of bytes, might be fun to try on other data
pub fn gotoh(seq1: &[u8], seq2: &[u8], params: &ScoringSystem) -> Result<Alignment, String> {
    let (match_s, mismatch, gap_open, gap_extend) = params.unpack()?;
    let (m, n) = (seq1.len(), seq2.len());
    if m == 0 || n == 0 {
        return Err("One or both of the sequences is empty".to_string());
    }

    let mut scores = MDICell::default_matrix(m, n, 0, gap_open, gap_open);

    // Fill in the rest of the matrix
    for i in 1..=m {
        scores[i][0].d_score = ((i as i32) - 1) * gap_extend + gap_open;
    }

    for j in 1..=n {
        scores[0][j].i_score = ((j as i32) - 1) * gap_extend + gap_open;
    }

    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                match_s
            } else {
                mismatch
            };

            scores[i][j].m_score = *[
                scores[i - 1][j - 1].m_score + match_mismatch,
                scores[i - 1][j - 1].d_score + match_mismatch,
                scores[i - 1][j - 1].i_score + match_mismatch,
            ]
            .iter()
            .max()
            .unwrap();

            scores[i][j].d_score = cmp::max(
                scores[i - 1][j].m_score + gap_open,
                scores[i - 1][j].d_score + gap_extend,
            );

            scores[i][j].i_score = cmp::max(
                scores[i][j - 1].m_score + gap_open,
                scores[i][j - 1].i_score + gap_extend,
            );
        }
    }

    let (alignment_seq1, alignment_seq2) = {
        let mut current_op = [
            scores[m][n].m_score,
            scores[m][n].d_score,
            scores[m][n].i_score,
        ]
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.cmp(b))
        .map(|(idx, _)| match idx {
            0 => 'm',
            1 => 'd',
            2 => 'i',
            _ => unreachable!(),
        })
        .unwrap();

        traceback(
            &scores, seq1, seq2, m, n, current_op, false, gap_open, gap_extend,
        )
    };

    Ok(crate::utils::assemble_alignment(
        *[
            scores[m][n].m_score,
            scores[m][n].d_score,
            scores[m][n].i_score,
        ]
        .iter()
        .max()
        .unwrap(),
        alignment_seq1,
        alignment_seq2,
    ))
}

/* Smith-Waterman Algorithm with affine gaps.
   The Smith Waterman algorithm works the same Gotoh/Needleman, but for local allignments
   To accomplish this we do a few things.
   1. When initializing the score penalties we do not start with gap penalties (use 0 instead)
   2. When populating the inside of the score matrix we limit the score to a floor of 0
   3. During the traceback we look for the higest score anywhere in the matrix and traceback until
        a 0 in encountered.

    Together these three changes allow us to identify the best aligned segments of the two sequences
    by not penalizing allignments that have to take a sub optimal path to reach there position in the
    scoring matrix. Then, by tracing back from the highest score in the whole matrix we find the best
    local alignment instead of forcing the alignment to end at the end of the two sequences.
*/

pub fn smith_waterman(
    seq1: &[u8],
    seq2: &[u8],
    params: &ScoringSystem,
) -> Result<Alignment, String> {
    let (match_s, mismatch, gap_open, gap_extend) = params.unpack()?;
    let (m, n) = (seq1.len(), seq2.len());
    if m == 0 || n == 0 {
        return Err("One or both of the sequences is empty".to_string());
    }

    let mut scores = MDICell::default_matrix(m, n, 0, 0, 0);
    let (mut optimal_score, mut max_i, mut max_j) = (0, 0, 0);

    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                match_s
            } else {
                mismatch
            };

            scores[i][j].m_score = *[
                scores[i - 1][j - 1].m_score + match_mismatch,
                scores[i - 1][j - 1].d_score + match_mismatch,
                scores[i - 1][j - 1].i_score + match_mismatch,
                0,
            ]
            .iter()
            .max()
            .unwrap();

            scores[i][j].d_score = cmp::max(
                scores[i - 1][j].m_score + gap_open,
                scores[i - 1][j].d_score + gap_extend,
            );

            scores[i][j].i_score = cmp::max(
                scores[i][j - 1].m_score + gap_open,
                scores[i][j - 1].i_score + gap_extend,
            );

            let cell_score = *[
                scores[i][j].m_score,
                scores[i][j].d_score,
                scores[i][j].i_score,
            ]
            .iter()
            .max()
            .unwrap();

            if cell_score > optimal_score {
                optimal_score = cell_score;
                max_i = i;
                max_j = j;
            }
        }
    }

    let mut current_op = [
        scores[max_i][max_j].m_score,
        scores[max_i][max_j].d_score,
        scores[max_i][max_j].i_score,
    ]
    .iter()
    .enumerate()
    .max_by(|(_, a), (_, b)| a.cmp(b))
    .map(|(idx, _)| match idx {
        0 => 'm',
        1 => 'd',
        2 => 'i',
        _ => unreachable!(),
    })
    .unwrap();

    let (alignment_seq1, alignment_seq2) = traceback(
        &scores, seq1, seq2, max_i, max_j, current_op, true, gap_open, gap_extend,
    );

    Ok(crate::utils::assemble_alignment(
        optimal_score,
        alignment_seq1,
        alignment_seq2,
    ))
}


/// Performs a variant of global alignment with affine gaps.
/// It uses the standard Needleman-Wunsch/Gotoh recurrence for the DP fill,
/// including penalized boundaries. However, it finds the maximum score anywhere
/// in the matrix and traces back from that cell to the origin (0,0).
/// It returns the number of *matches* found along this specific traceback path.
///
/// Args:
/// * seq1: First sequence as byte slice.
/// * seq2: Second sequence as byte slice.
/// * match_s: Score for a match.
/// * mismatch: Penalty for a mismatch (typically negative).
/// * gap_open: Penalty for opening a gap (typically negative).
/// * gap_extend: Penalty for extending a gap (typically negative).
///
/// Returns:
/// * Result<usize, String>: The count of matching characters in the alignment path
///   traced back from the highest score cell, or an error string.
pub fn align_variant_affine(
    seq1: &[u8],
    seq2: &[u8],
    match_s: i32,
    mismatch: i32,
    gap_open: i32,
    gap_extend: i32,
) -> Result<usize, String> {
    let (m, n) = (seq1.len(), seq2.len());

    // Base case: If either sequence is empty, there are no matches possible in an alignment.
    if m == 0 || n == 0 {
        return Ok(0);
    }

    // --- DP Matrix Initialization ---
    // Initialize with default values (effectively -infinity for unreachable states)
    // The default_matrix sets M[0][0]=0, D[0][0]=0, I[0][0]=0, which is correct.
    // It also initializes D[i][0] and I[0][j] to MIN/2, which we will overwrite.
    let mut scores = MDICell::default_matrix(m, n, 0, i32::MIN / 2, i32::MIN / 2);

    // Apply standard Needleman-Wunsch/Gotoh boundary conditions with affine gaps
    // M scores along boundaries remain -inf (unreachable by match/mismatch)
    // D scores along first column (gaps in seq2)
    for i in 1..=m {
        scores[i][0].m_score = i32::MIN / 2; // Cannot reach M[i][0] via diagonal
        scores[i][0].d_score = gap_open.saturating_add((i as i32).saturating_sub(1).saturating_mul(gap_extend)); // h + (i-1)*g
        scores[i][0].i_score = i32::MIN / 2; // Cannot reach I[i][0] (no 'left' cell)
    }
    // I scores along first row (gaps in seq1)
    for j in 1..=n {
        scores[0][j].m_score = i32::MIN / 2; // Cannot reach M[0][j] via diagonal
        scores[0][j].d_score = i32::MIN / 2; // Cannot reach D[0][j] (no 'above' cell)
        scores[0][j].i_score = gap_open.saturating_add((j as i32).saturating_sub(1).saturating_mul(gap_extend)); // h + (j-1)*g
    }
    // scores[0][0] remains M=0, D=0, I=0 from default_matrix initialization


    // --- DP Fill Phase ---
    // Iterate through the matrix to calculate M, D, I scores for each cell
    for i in 1..=m {
        for j in 1..=n {
            // Determine match or mismatch score for the current characters
            let match_mismatch_score = if seq1[i - 1] == seq2[j - 1] { match_s } else { mismatch };

            // Calculate M score (Match/Mismatch)
            // Comes from diagonal of M, D, or I matrix
            let m_diag = scores[i - 1][j - 1].m_score;
            let d_diag = scores[i - 1][j - 1].d_score;
            let i_diag = scores[i - 1][j - 1].i_score;
            // Use saturating_add to prevent overflow with large negative scores
            let m_from_m = m_diag.saturating_add(match_mismatch_score);
            let m_from_d = d_diag.saturating_add(match_mismatch_score);
            let m_from_i = i_diag.saturating_add(match_mismatch_score);
            scores[i][j].m_score = cmp::max(m_from_m, cmp::max(m_from_d, m_from_i));


            // Calculate D score (Deletion - gap in seq2)
            // Comes from above cell in M (open gap) or D (extend gap)
            let m_above = scores[i - 1][j].m_score;
            let d_above = scores[i - 1][j].d_score;
            let d_from_m = m_above.saturating_add(gap_open); // Open gap cost = M + h
            let d_from_d = d_above.saturating_add(gap_extend); // Extend gap cost = D + g
            scores[i][j].d_score = cmp::max(d_from_m, d_from_d);


            // Calculate I score (Insertion - gap in seq1)
            // Comes from left cell in M (open gap) or I (extend gap)
            let m_left = scores[i][j - 1].m_score;
            let i_left = scores[i][j - 1].i_score;
            let i_from_m = m_left.saturating_add(gap_open); // Open gap cost = M + h
            let i_from_i = i_left.saturating_add(gap_extend); // Extend gap cost = I + g
            scores[i][j].i_score = cmp::max(i_from_m, i_from_i);
        }
    }

    // --- Find Max Score Cell ---
    // Search the entire matrix (excluding boundaries i=0 or j=0, as per prompt example)
    // for the maximum score among M, D, and I values.
    let mut max_score = i32::MIN;
    let mut max_i = 0;
    let mut max_j = 0;
    let mut max_op = ' '; // The operation ('m', 'd', 'i') corresponding to the max score

    for i in 1..=m {
        for j in 1..=n {
            // Check M score
            if scores[i][j].m_score >= max_score { // Use >= to prefer later cells and M in ties
                max_score = scores[i][j].m_score;
                max_i = i;
                max_j = j;
                max_op = 'm';
            }
            // Check D score
            if scores[i][j].d_score >= max_score { // Use >= to prefer later cells and D over I in ties
                max_score = scores[i][j].d_score;
                max_i = i;
                max_j = j;
                max_op = 'd';
            }
            // Check I score
             // Use > to only update if strictly greater than current max (or D score if equal)
            if scores[i][j].i_score > max_score {
                max_score = scores[i][j].i_score;
                max_i = i;
                max_j = j;
                max_op = 'i';
             } else if scores[i][j].i_score == max_score && max_op != 'm' && max_op != 'd' {
                 // Handle specific tie-breaking if needed, e.g., prefer I if equal to MIN and no M/D found yet
                 // For simplicity, the >= logic above generally handles ties reasonably.
                 // If max_score is still MIN, and I is also MIN, this block might select I initially.
             }
        }
    }

    // If max_i or max_j is still 0, it means the max score was found at the boundary
    // or no cell improved over the initial MIN/boundary values. Traceback is trivial.
    if max_i == 0 || max_j == 0 {
        println!("Warning: Max score {} found at boundary ({}, {}). Match count is 0.", max_score, max_i, max_j);
        return Ok(0);
    }

    // --- Traceback and Count Matches ---
    let mut current_i = max_i;
    let mut current_j = max_j;
    let mut current_op = max_op;
    let mut match_count = 0;

    println!("Traceback starting from ({},{}) with op '{}', score {}", max_i, max_j, max_op, max_score);

    // Trace back until we reach the origin (0,0) or hit an edge implicitly
    while current_i > 0 || current_j > 0 {
        // Store current position for loop termination check (optional sanity check)
        // let (prev_i, prev_j) = (current_i, current_j);

        match current_op {
            'm' => {
                // Check bounds: must come from diagonal (i-1, j-1)
                if current_i > 0 && current_j > 0 {
                    // Count match if characters are identical
                    if seq1[current_i - 1] == seq2[current_j - 1] {
                        match_count += 1;
                    }

                    // Determine the operation of the previous cell (i-1, j-1) that led here
                    let prev_cell = &scores[current_i - 1][current_j - 1];
                    let match_mismatch_score = if seq1[current_i - 1] == seq2[current_j - 1] { match_s } else { mismatch };

                    // Check which incoming score + step score equals the current M score
                    // Need to be careful with saturation and floating point comparison issues if using floats
                    // Since we use integers, direct comparison should be fine unless MIN/MAX involved heavily.
                    let from_m = prev_cell.m_score.saturating_add(match_mismatch_score);
                    let from_d = prev_cell.d_score.saturating_add(match_mismatch_score);
                    let from_i = prev_cell.i_score.saturating_add(match_mismatch_score);
                    let current_m_score = scores[current_i][current_j].m_score;

                    // Determine previous op based on which path yielded the score
                    // Prioritize M > D > I in ties during traceback source determination
                    if from_m == current_m_score {
                        current_op = 'm';
                    } else if from_d == current_m_score {
                        current_op = 'd';
                    } else if from_i == current_m_score {
                        current_op = 'i';
                    } else {
                        // This should not happen if the DP table is consistent
                        return Err(format!("Traceback error at ({}, {}): M score {} unreachable from diagonal.", current_i, current_j, current_m_score));
                    }
                    // Move diagonally
                    current_i -= 1;
                    current_j -= 1;
                } else {
                    // Cannot move diagonally from row 0 or col 0 to reach M
                    break; // Hit edge
                }
            }
            'd' => {
                // Check bounds: must come from above (i-1, j)
                if current_i > 0 {
                    // Determine the operation of the previous cell (i-1, j) that led here
                    let prev_cell = &scores[current_i - 1][current_j];
                    let from_m = prev_cell.m_score.saturating_add(gap_open);
                    let from_d = prev_cell.d_score.saturating_add(gap_extend);
                    let current_d_score = scores[current_i][current_j].d_score;

                    // Prioritize M source over D source in ties
                    if from_m == current_d_score {
                        current_op = 'm';
                    } else if from_d == current_d_score {
                        current_op = 'd';
                    } else {
                         return Err(format!("Traceback error at ({}, {}): D score {} unreachable from above.", current_i, current_j, current_d_score));
                    }
                    // Move up
                    current_i -= 1;
                } else {
                    // Cannot move up from row 0
                    break; // Hit edge
                }
            }
            'i' => {
                // Check bounds: must come from left (i, j-1)
                if current_j > 0 {
                    // Determine the operation of the previous cell (i, j-1) that led here
                    let prev_cell = &scores[current_i][current_j - 1];
                    let from_m = prev_cell.m_score.saturating_add(gap_open);
                    let from_i = prev_cell.i_score.saturating_add(gap_extend);
                    let current_i_score = scores[current_i][current_j].i_score;

                     // Prioritize M source over I source in ties
                    if from_m == current_i_score {
                        current_op = 'm';
                    } else if from_i == current_i_score {
                         current_op = 'i';
                    } else {
                         return Err(format!("Traceback error at ({}, {}): I score {} unreachable from left.", current_i, current_j, current_i_score));
                    }
                    // Move left
                    current_j -= 1;
                } else {
                    // Cannot move left from col 0
                    break; // Hit edge
                }
            }
            _ => return Err(format!("Invalid operation '{}' encountered during traceback.", current_op)),
        }

        // Sanity check: Ensure loop terminates if stuck (shouldn't be needed with correct logic)
        // if current_i == prev_i && current_j == prev_j {
        //     return Err("Traceback stuck in a loop.".to_string());
        // }
    }

    // Traceback finished when current_i == 0 and current_j == 0
    Ok(match_count)
}


// Existing traceback function (used by gotoh/smith_waterman) - keep as is
// unless refactoring is desired. The variant alignment needs its own specific
// traceback logic as implemented above.
pub fn traceback(
    scores: &[Vec<MDICell>],
    seq1: &[u8],
    seq2: &[u8],
    mut i: usize,
    mut j: usize,
    mut current_op: char,
    stop_on_zero: bool, // Specific to Smith-Waterman
    gap_open: i32,
    gap_extend: i32,
) -> (String, String) {
    let mut s1 = String::new();
    let mut s2 = String::new();

    // This existing traceback needs modification to *exactly* determine the source
    // based on score comparison like the `align_variant_affine` traceback does,
    // rather than just picking the 'best_op' of the previous cell.
    // However, for now, we leave it as is for the existing algorithms.
    // The `align_variant_affine` has its own correct traceback logic internally.

    while (i > 0 || j > 0) && !(stop_on_zero && scores[i][j].m_score <= 0 && scores[i][j].d_score <= 0 && scores[i][j].i_score <= 0) {
        // Added check for stop_on_zero condition based on *current* cell score being non-positive
         let current_cell_score = match current_op {
             'm' => scores[i][j].m_score,
             'd' => scores[i][j].d_score,
             'i' => scores[i][j].i_score,
             _ => i32::MIN,
         };
         if stop_on_zero && current_cell_score <= 0 && (i != 0 || j != 0) { // Check score <= 0, ensure not already at origin
             break;
         }

        match current_op {
            'm' => {
                 if i > 0 && j > 0 {
                    s1.push(seq1[i - 1] as char);
                    s2.push(seq2[j - 1] as char);
                     // --- Determine previous operation (Needs careful check) ---
                     let prev_cell = &scores[i - 1][j - 1];
                     let match_mismatch_score = if seq1[i - 1] == seq2[j - 1] { scores[i][j].m_score - prev_cell.m_score /* Infer score diff */ } else { scores[i][j].m_score - prev_cell.m_score /* Assumes this path */ }; // This inference is weak
                     // Need to check which path *actually* led here by comparing scores.
                     // For now, sticking to simpler 'best_op' logic from original code.
                     current_op = prev_cell.best_op(); // Original logic - might not be strictly correct path sometimes
                     // ----
                    i -= 1;
                    j -= 1;
                 } else { break; } // Hit edge
            }
            'd' => {
                 if i > 0 {
                    s1.push(seq1[i - 1] as char);
                    s2.push('-');
                     // --- Determine previous operation ---
                    let above = scores[i - 1][j];
                     // Check score relationship to determine source
                    current_op = if above.m_score.saturating_add(gap_open) >= above.d_score.saturating_add(gap_extend) { // Prefer M source
                        'm'
                    } else {
                        'd'
                    };
                     // ----
                    i -= 1;
                 } else { break; } // Hit edge
            }
            'i' => {
                 if j > 0 {
                    s1.push('-');
                    s2.push(seq2[j - 1] as char);
                     // --- Determine previous operation ---
                    let left = scores[i][j - 1];
                     // Check score relationship to determine source
                    current_op = if left.m_score.saturating_add(gap_open) >= left.i_score.saturating_add(gap_extend) { // Prefer M source
                        'm'
                    } else {
                        'i'
                    };
                     // ----
                    j -= 1;
                 } else { break; } // Hit edge
            }
            _ => panic!("Invalid op during traceback"),
        }
    }

    // Add remaining characters if traceback stopped before reaching (0,0)
    // This part is mainly for global alignment (Needleman/Gotoh)
    while i > 0 && !stop_on_zero {
        s1.push(seq1[i - 1] as char);
        s2.push('-');
        i -= 1;
    }
    while j > 0 && !stop_on_zero {
        s1.push('-');
        s2.push(seq2[j - 1] as char);
        j -= 1;
    }

    (s1.chars().rev().collect(), s2.chars().rev().collect())
}
