//alignment.rs
use std::cmp;

use crate::types::*;
use crate::utils; // For calculating stats, etc.

// Neddleman-Wunsch w/ affine gaps a.k.a Gotoh 82:

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
pub fn gotoh(seq1: &[u8], seq2: &[u8], params: &ScoringSystem) -> Result<Alignment, String> {
    let (match_s, mismatch, gap_open, gap_extend) = params.unpack()?;
    let (m, n) = (seq1.len(), seq2.len());

    // Handle empty sequences - global alignment means aligning gaps
    if m == 0 || n == 0 {
        let (s1_aligned, s2_aligned) = if m == 0 && n == 0 {
            (String::new(), String::new())
        } else if m == 0 {
            (
                "-".repeat(n),
                String::from_utf8(seq2.to_vec()).unwrap_or_default(),
            )
        } else { // n == 0
            (
                String::from_utf8(seq1.to_vec()).unwrap_or_default(),
                "-".repeat(m),
            )
        };
        // Score is just the gap penalty
        let score = if m == 0 && n == 0 { 0 }
                    else if m == 0 { gap_open + (n as i32 - 1) * gap_extend }
                    else { gap_open + (m as i32 - 1) * gap_extend };
        return Ok(utils::assemble_alignment(score, s1_aligned, s2_aligned));
    }

    // DP matrices: M, D, I
    let mut scores = MDICell::default_matrix(m, n, 0, i32::MIN / 2, i32::MIN / 2);

    // Init boundaries with affine penalties
    for i in 1..=m {
        scores[i][0].m_score = i32::MIN / 2; // Can't match
        scores[i][0].d_score = gap_open + (i as i32 - 1) * gap_extend;
        scores[i][0].i_score = i32::MIN / 2; // Can't insert
    }
    for j in 1..=n {
        scores[0][j].m_score = i32::MIN / 2;
        scores[0][j].d_score = i32::MIN / 2; // Can't delete
        scores[0][j].i_score = gap_open + (j as i32 - 1) * gap_extend;
    }

    // Fill DP table
    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] { match_s } else { mismatch };

            // M score (match/mismatch)
            let m_diag = scores[i - 1][j - 1].m_score;
            let d_diag = scores[i - 1][j - 1].d_score;
            let i_diag = scores[i - 1][j - 1].i_score;
            scores[i][j].m_score = cmp::max(m_diag, cmp::max(d_diag, i_diag)) + match_mismatch;

            // D score (deletion in seq2)
            let m_above = scores[i - 1][j].m_score;
            let d_above = scores[i - 1][j].d_score;
            scores[i][j].d_score = cmp::max(
                m_above + gap_open,  // Start gap
                d_above + gap_extend, // Extend gap
            );

            // I score (insertion in seq1)
            let m_left = scores[i][j - 1].m_score;
            let i_left = scores[i][j - 1].i_score;
            scores[i][j].i_score = cmp::max(
                m_left + gap_open,  // Start gap
                i_left + gap_extend, // Extend gap
            );
        }
    }

    // Optimal score is max of M, D, I in the bottom-right cell
    let optimal_score = *[
        scores[m][n].m_score,
        scores[m][n].d_score,
        scores[m][n].i_score,
    ]
    .iter()
    .max()
    .unwrap_or(&(i32::MIN / 2));

    // Determine starting state for traceback
    let start_op = if optimal_score == scores[m][n].m_score { 'm' }
                   else if optimal_score == scores[m][n].d_score { 'd' }
                   else { 'i' };

    // Traceback
    let (alignment_seq1, alignment_seq2) = traceback(
        &scores, seq1, seq2, m, n, start_op, false, // false: not local
         gap_open, gap_extend, match_s, mismatch
    );

    Ok(utils::assemble_alignment(
        optimal_score,
        alignment_seq1,
        alignment_seq2,
    ))
}

/* Smith-Waterman Algorithm with affine gaps. (Local)
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
        // Local alignment of empty string is empty, score 0
        return Ok(utils::assemble_alignment(0, String::new(), String::new()));
    }

    // Init DP matrices with 0s for local start
    let mut scores = MDICell::default_matrix(m, n, 0, 0, 0);
     // Boundaries stay 0

    let (mut optimal_score, mut max_i, mut max_j) = (0, 0, 0); // Track max score position

    // Fill DP table
    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] { match_s } else { mismatch };

            // M score (match/mismatch) - floor at 0 for local
            let m_diag = scores[i - 1][j - 1].m_score;
            let d_diag = scores[i - 1][j - 1].d_score;
            let i_diag = scores[i - 1][j - 1].i_score;
            scores[i][j].m_score = cmp::max(
                0, // Local alignment can start anywhere (score reset)
                cmp::max(m_diag, cmp::max(d_diag, i_diag)) + match_mismatch
            );

            // D score (deletion in seq2) - no floor here, tracks gap penalty accumulation
            let m_above = scores[i - 1][j].m_score;
            let d_above = scores[i - 1][j].d_score;
            scores[i][j].d_score = cmp::max(
                m_above + gap_open,
                d_above + gap_extend,
            );

            // I score (insertion in seq1) - no floor here
            let m_left = scores[i][j - 1].m_score;
            let i_left = scores[i][j - 1].i_score;
            scores[i][j].i_score = cmp::max(
                m_left + gap_open,
                i_left + gap_extend,
            );

            // Check if current cell's max score is the new overall max
            let current_max_cell_score = *[
                scores[i][j].m_score,
                // Don't consider D or I for *starting* the traceback? Let's use M only?
                // No, standard SW finds max of M,D,I anywhere.
                scores[i][j].d_score,
                scores[i][j].i_score,
            ]
            .iter()
            .max()
            .unwrap_or(&0); // If all negative, max is 0

            if current_max_cell_score >= optimal_score { // >= picks later cells in ties
                optimal_score = current_max_cell_score;
                max_i = i;
                max_j = j;
            }
        }
    }

     // If max score is 0, alignment is empty
     if optimal_score <= 0 {
         return Ok(utils::assemble_alignment(0, String::new(), String::new()));
     }

    // Determine starting state at the max score cell
    let start_op = if optimal_score == scores[max_i][max_j].m_score { 'm' }
                   else if optimal_score == scores[max_i][max_j].d_score { 'd' }
                   else { 'i' };

    // Traceback, stop when score hits 0
    let (alignment_seq1, alignment_seq2) = traceback(
        &scores, seq1, seq2, max_i, max_j, start_op, true, // true: stop on zero for local
         gap_open, gap_extend, match_s, mismatch
    );

    Ok(utils::assemble_alignment(
        optimal_score,
        alignment_seq1,
        alignment_seq2,
    ))
}

// Variant alignment for Project 3. NW fill, find max score anywhere, traceback to origin, count matches.
// Returns the count of matches.
pub fn align_variant_affine(
    seq1: &[u8],
    seq2: &[u8],
    match_s: i32,
    mismatch: i32,
    gap_open: i32,
    gap_extend: i32,
) -> Result<usize, String> {
    let (m, n) = (seq1.len(), seq2.len());

    if m == 0 || n == 0 {
        return Ok(0); // No matches if empty
    }

    // Init like Gotoh (NW boundaries)
    let mut scores = MDICell::default_matrix(m, n, 0, i32::MIN / 2, i32::MIN / 2);
    for i in 1..=m {
        scores[i][0].m_score = i32::MIN / 2;
        scores[i][0].d_score = gap_open + (i as i32 - 1) * gap_extend;
        scores[i][0].i_score = i32::MIN / 2;
    }
    for j in 1..=n {
        scores[0][j].m_score = i32::MIN / 2;
        scores[0][j].d_score = i32::MIN / 2;
        scores[0][j].i_score = gap_open + (j as i32 - 1) * gap_extend;
    }

    // Fill like Gotoh
    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch_score = if seq1[i - 1] == seq2[j - 1] { match_s } else { mismatch };

            // M score
            let m_diag = scores[i - 1][j - 1].m_score;
            let d_diag = scores[i - 1][j - 1].d_score;
            let i_diag = scores[i - 1][j - 1].i_score;
            scores[i][j].m_score = cmp::max(m_diag, cmp::max(d_diag, i_diag)) + match_mismatch_score;

            // D score
            let m_above = scores[i - 1][j].m_score;
            let d_above = scores[i - 1][j].d_score;
            scores[i][j].d_score = cmp::max(
                m_above + gap_open,
                d_above + gap_extend,
            );

            // I score
            let m_left = scores[i][j - 1].m_score;
            let i_left = scores[i][j - 1].i_score;
            scores[i][j].i_score = cmp::max(
                m_left + gap_open,
                i_left + gap_extend,
            );
        }
    }

    // Find Max Score Cell anywhere
    let mut max_score = i32::MIN;
    let mut max_i = 0;
    let mut max_j = 0;

    for i in 0..=m {
        for j in 0..=n {
            if i == 0 && j == 0 { continue; } // Skip origin

            let cell_max = *[
                scores[i][j].m_score,
                scores[i][j].d_score,
                scores[i][j].i_score,
            ]
            .iter()
            .max()
            .unwrap_or(&i32::MIN);

            if cell_max >= max_score { // >= picks bottom-right ties
                max_score = cell_max;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Find the operation for the max score cell
    let max_op = if max_score == scores[max_i][max_j].m_score { 'm' }
                 else if max_score == scores[max_i][max_j].d_score { 'd' }
                 else if max_score == scores[max_i][max_j].i_score { 'i' }
                 else {
                     // Handle edge cases where max score is from init, not M/D/I calculation
                     if max_i == 0 && max_j == 0 { return Ok(0); } // Max was origin
                     else if max_i == 0 { 'i' } // Max on top edge (must be insertion)
                     else if max_j == 0 { 'd' } // Max on left edge (must be deletion)
                     else {
                         return Err(format!("Max score doesn't match M/D/I scores."));
                     }
                 };

    // Traceback from max cell to origin (0,0) and count matches
    let mut current_i = max_i;
    let mut current_j = max_j;
    let mut current_op = max_op;
    let mut match_count = 0usize;


    while current_i > 0 || current_j > 0 {
        match current_op {
            'm' => {
                if current_i > 0 && current_j > 0 {
                    // Count if match
                    if seq1[current_i - 1] == seq2[current_j - 1] {
                        match_count += 1;
                    }
                    // Determine previous op
                    let prev_cell = &scores[current_i - 1][current_j - 1];
                    let current_m_score = scores[current_i][current_j].m_score;
                    let mm_score = if seq1[current_i - 1] == seq2[current_j - 1] { match_s } else { mismatch };

                    // Check which path gave current M score (prefer M > D > I)
                    if prev_cell.m_score + mm_score == current_m_score { current_op = 'm'; }
                    else if prev_cell.d_score + mm_score == current_m_score { current_op = 'd'; }
                    else if prev_cell.i_score + mm_score == current_m_score { current_op = 'i'; }
                    else { return Err(format!("Traceback error M at ({}, {})", current_i, current_j)); }

                    current_i -= 1;
                    current_j -= 1;
                } else { break; } // Hit edge
            }
            'd' => {
                if current_i > 0 {
                    let prev_cell = &scores[current_i - 1][current_j];
                    let current_d_score = scores[current_i][current_j].d_score;

                    // Check source (prefer M > D)
                    if prev_cell.m_score + gap_open == current_d_score { current_op = 'm'; }
                    else if prev_cell.d_score + gap_extend == current_d_score { current_op = 'd'; }
                    else {
                         // Handle boundary case - assume came from D above on edge
                         if current_j == 0 { current_op = 'd'; }
                         else { return Err(format!("Traceback error D at ({}, {})", current_i, current_j)); }
                    }
                    current_i -= 1;
                } else { break; } // Hit edge
            }
            'i' => {
                if current_j > 0 {
                    let prev_cell = &scores[current_i][current_j - 1];
                    let current_i_score = scores[current_i][current_j].i_score;

                    // Check source (prefer M > I)
                    if prev_cell.m_score + gap_open == current_i_score { current_op = 'm'; }
                    else if prev_cell.i_score + gap_extend == current_i_score { current_op = 'i'; }
                    else {
                         // Handle boundary case - assume came from I left on edge
                         if current_i == 0 { current_op = 'i'; }
                         else { return Err(format!("Traceback error I at ({}, {})", current_i, current_j)); }
                    }
                    current_j -= 1;
                } else { break; } // Hit edge
            }
            _ => return Err(format!("Invalid op in traceback.")),
        }
    }

    Ok(match_count)
}


// General traceback, builds aligned strings.
fn traceback(
    scores: &[Vec<MDICell>],
    seq1: &[u8],
    seq2: &[u8],
    mut i: usize,
    mut j: usize,
    mut current_op: char,
    stop_on_zero: bool, // For SW
    gap_open: i32,
    gap_extend: i32,
    match_s: i32,     // Need scores for source determination
    mismatch: i32,
) -> (String, String) {
    let mut s1_aligned = Vec::new();
    let mut s2_aligned = Vec::new();

    while i > 0 || j > 0 {
        let current_score = match current_op {
            'm' => scores[i][j].m_score,
            'd' => scores[i][j].d_score,
            'i' => scores[i][j].i_score,
            _ => i32::MIN,
        };

        // SW stops if score drops to 0
        if stop_on_zero && current_score <= 0 {
            break;
        }

        match current_op {
            'm' => {
                if i > 0 && j > 0 {
                    s1_aligned.push(seq1[i - 1]);
                    s2_aligned.push(seq2[j - 1]);
                    // Determine previous op based on which path gave current M score
                    let prev_cell = &scores[i - 1][j - 1];
                    let mm_score = if seq1[i - 1] == seq2[j - 1] { match_s } else { mismatch };
                    // Check source (M > D > I preference)
                    if prev_cell.m_score + mm_score == current_score { current_op = 'm'; }
                    else if prev_cell.d_score + mm_score == current_score { current_op = 'd'; }
                    else { current_op = 'i'; } // Must be I
                    i -= 1;
                    j -= 1;
                } else { break; } // Hit edge
            }
            'd' => { // Gap in seq2
                if i > 0 {
                    s1_aligned.push(seq1[i - 1]);
                    s2_aligned.push(b'-');
                    // Determine previous op
                    let prev_cell = &scores[i - 1][j];
                    // Check source (M > D preference)
                    if prev_cell.m_score + gap_open == current_score { current_op = 'm'; }
                    else { current_op = 'd'; } // Must be D
                    i -= 1;
                } else { break; } // Hit edge
            }
            'i' => { // Gap in seq1
                if j > 0 {
                    s1_aligned.push(b'-');
                    s2_aligned.push(seq2[j - 1]);
                    // Determine previous op
                    let prev_cell = &scores[i][j - 1];
                     // Check source (M > I preference)
                    if prev_cell.m_score + gap_open == current_score { current_op = 'm'; }
                    else { current_op = 'i'; } // Must be I
                    j -= 1;
                } else { break; } // Hit edge
            }
            _ => panic!("Invalid op during traceback"),
        }
    }

    // Reverse collected bytes and convert to String
    s1_aligned.reverse();
    s2_aligned.reverse();
    let s1_final = String::from_utf8(s1_aligned).unwrap_or_default();
    let s2_final = String::from_utf8(s2_aligned).unwrap_or_default();

    (s1_final, s2_final)
}