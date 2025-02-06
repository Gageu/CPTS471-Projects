//alignment.rs
//  Allignmnet algorithms for projects in CS 471

use std::{num::TryFromIntError, str::pattern::Pattern};

use crate::types::*;





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
fn gotoh(seq1: &str, seq2: &str, params: &ScoringSystem) -> Alignment {
    let mut alignment_seq1 = String::new();
    let mut alignment_seq2 = String::new();
    let mut optimal_score = 0;
    let (m, n) = (seq1.len(), seq2.len());

    let mut scores = vec![vec![MDICell::default(); n+1]; m+1];

    scores[0][0].m_score = 0;

    for i in 1..=m{
        scores[i][0].d_score = 
    }
    

    let stats = AllignmentStats::new(
        allignment_length,
        match_count, 
        mismatch_count, 
        gap_open_count, 
        gap_extend_count, 
        total_gaps, 
        identity_percent
    )


    let alignmnet = Alignment::new(score, sequence1, sequence2, stats);
}

fn smith_waterman(seq1: &str, seq2: &str, params: &ScoringSystem) -> Alignment {

}
