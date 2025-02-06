//alignment.rs
//  Allignmnet algorithms for projects in CS 471

use std::{cmp, num::TryFromIntError, str::pattern::Pattern};

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

// Accepts any two sequences of bytes, might be fun to try on other data
fn gotoh(seq1: &[u8], seq2: &[u8], params: &ScoringSystem) -> Result<Alignment, String> {
    //TODO: switch return type to be a u8 instead of string since RUst has strings encoded in utf8
    //      which is more complexity than makes sense for fasta since it's ASCII encoded

    // Get the value of every parameter needed from the param struct
    // if is not provided return Err(String(message)) to be printed by main
    let match_s = params
        .match_score()
        .ok_or("Missing match score parameter".to_string())?;
    let mismatch = params
        .mismatch_score()
        .ok_or("Missing mismatch score parameter".to_string())?;
    let gap_open = params
        .gap_open_score()
        .ok_or("Missing gap open score parameter".to_string())?;
    let gap_extend = params
        .gap_extend_score()
        .ok_or("Missing gap extend score parameter".to_string())?;

    let mut alignment_seq1 = String::new();
    let mut alignment_seq2 = String::new();
    let mut optimal_score = 0;

    let (m, n) = (seq1.len(), seq2.len());

    let mut scores = vec![vec![MDICell::default(); n + 1]; m + 1];

    // Populate score matrix
    scores[0][0].m_score = 0;

    for i in 1..=m {
        scores[i][0].d_score = ((i as i32) - 1) * gap_extend + gap_open;
    }

    for j in 1..=n {
        scores[0][j].i_score = ((j as i32) - 1) * gap_extend + gap_open;
    }

    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch_score = if seq1[i - 1] == seq2[j - 1] {
                match_s
            } else {
                mismatch
            };

            //Somehow this is the cleanest way to get the max of threes expressions
            scores[i][j].m_score = *[
                scores[i - 1][j - 1].m_score + match_mismatch_score,
                scores[i - 1][j - 1].d_score + match_mismatch_score,
                scores[i - 1][j - 1].i_score + match_mismatch_score,
            ]
            .iter()
            .max()
            .unwrap();

            scores[i][j].d_score = cmp::max(
                scores[i - 1][j].m_score + gap_open,
                scores[i - 1][j].d_score + gap_extend,
            );

            scores[i][j].d_score = cmp::max(
                scores[i][j - 1].m_score + gap_open,
                scores[i][j - 1].d_score + gap_extend,
            );
        }
    }

    // Traceback Optimal Score

    // Find which score had the max
    // code for finding index of max elemnent from: https://stackoverflow.com/a/53908709
    // I just added the match statement to make it fit my needs
    // This works by adding all of the scores to a list, numbering them, finding the max score and
    // returning its associated number (max_by) then mapping that to the characters.
    // Kind of awful readability and should probably be split into multiple lines but this solution
    // is to fun to redo
    let mut current_op = [
        scores[m][n].m_score, //0
        scores[m][n].d_score, //1
        scores[m][n].i_score, //2
    ]
    .iter()
    .enumerate()
    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
    .map(|(index, _)| match index {
        0 => 'm',
        1 => 'd',
        2 => 'i',
        _ => panic!("You shouldn't be able to reach this"),
    })
    .unwrap();

    //don't polute the full function scope with i and j
    {
        let (mut i, mut j) = (m, n);
        while (i > 0) && (j > 0) {
            //Add letters to sequence
            //Check how we could have gotten here
            //Decrement in appropriate direction and set current_op
            match current_op {
                'm' => {
                    alignment_seq1.push(seq1[i - 1] as char);
                    alignment_seq2.push(seq2[j - 1] as char);

                    let match_mismatch_score = if seq1[i - 1] == seq2[j - 1] {
                        match_s
                    } else {
                        mismatch
                    };

                    let diag_cell = scores[i - 1][j - 1];

                    current_op = if diag_cell.m_score >= diag_cell.i_score
                        && diag_cell.m_score >= diag_cell.d_score
                    {
                        'm'
                    } else if diag_cell.d_score >= diag_cell.i_score {
                        'd'
                    } else {
                        'i'
                    };

                    i -= 1;
                    j -= 1;
                }
                'd' => {
                    alignment_seq1.push(seq1[i - 1] as char);
                    alignment_seq2.push('-');

                    let above_cell = scores[i - 1][j];

                    current_op = if above_cell.m_score >= above_cell.i_score
                        && above_cell.m_score >= above_cell.d_score
                    {
                        'm'
                    } else if above_cell.d_score >= above_cell.i_score {
                        'd'
                    } else {
                        'i'
                    };

                    i -= 1;
                }
                'i' => {}
                _ => {
                    panic!("current_op was not a deletion, insertion, or match/mismatch")
                }
            }
        }

        //Collect last row or column
        while i > 0 {
            alignment_seq1.push(seq1[i - 1] as char);
            alignment_seq2.push('-');
            i -= 1;
        }

        while j > 0 {
            alignment_seq1.push('-');
            alignment_seq2.push(seq2[j - 1] as char);
            j -= 1;
        }
    }
    alignment_seq1 = alignment_seq1.chars().rev().collect();
    alignment_seq2 = alignment_seq2.chars().rev().collect();

    let stats = AllignmentStats::new(
        allignment_length,
        match_count,
        mismatch_count,
        gap_open_count,
        gap_extend_count,
        total_gaps,
        identity_percent,
    );

    Ok(Alignment::new(
        optimal_score,
        alignment_seq1,
        alignment_seq2,
        stats,
    ))
}

fn smith_waterman(seq1: &str, seq2: &str, params: &ScoringSystem) -> Alignment {}
