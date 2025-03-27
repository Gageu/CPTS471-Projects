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

    let mut scores = MDICell::default_matrix(
        m,
        n,
        0,
        gap_open,
        gap_open,
    );

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

        traceback(&scores, seq1, seq2, m, n, current_op, false, gap_open, gap_extend)
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

pub fn smith_waterman(seq1: &[u8], seq2: &[u8], params: &ScoringSystem) -> Result<Alignment, String> {
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
        &scores,
        seq1,
        seq2,
        max_i,
        max_j,
        current_op,
        true,
        gap_open,
        gap_extend,
    );

    Ok(crate::utils::assemble_alignment(
        optimal_score,
        alignment_seq1,
        alignment_seq2,
    ))
}


pub fn traceback(
    scores: &[Vec<MDICell>],
    seq1: &[u8],
    seq2: &[u8],
    mut i: usize,
    mut j: usize,
    mut current_op: char,
    stop_on_zero: bool,
    gap_open: i32,
    gap_extend: i32,
) -> (String, String) {
    let mut s1 = String::new();
    let mut s2 = String::new();

    while i > 0 && j > 0 {
        match current_op {
            'm' => {
                if stop_on_zero && scores[i][j].m_score == 0 { break; }
                s1.push(seq1[i - 1] as char);
                s2.push(seq2[j - 1] as char);
                current_op = scores[i - 1][j - 1].best_op();
                i -= 1;
                j -= 1;
            },
            'd' => {
                if stop_on_zero && scores[i][j].d_score == 0 { break; }
                s1.push(seq1[i - 1] as char);
                s2.push('-');
                let above = scores[i - 1][j];
                current_op = if above.m_score + gap_open >= above.d_score + gap_extend { 'm' } else { 'd' };
                i -= 1;
            },
            'i' => {
                if stop_on_zero && scores[i][j].i_score == 0 { break; }
                s1.push('-');
                s2.push(seq2[j - 1] as char);
                let left = scores[i][j - 1];
                current_op = if left.m_score + gap_open >= left.i_score + gap_extend { 'm' } else { 'i' };
                j -= 1;
            },
            _ => panic!("Invalid op")
        }
    }

    while i > 0 { s1.push(seq1[i - 1] as char); s2.push('-'); i -= 1; }
    while j > 0 { s1.push('-'); s2.push(seq2[j - 1] as char); j -= 1; }

    (s1.chars().rev().collect(), s2.chars().rev().collect())
}