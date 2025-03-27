use std::env;

use aligment::{gotoh, smith_waterman};
use types::{ProjectSelection};
use utils::{parse_sequences_from_fasta, read_score_config, print_alignment, len_without_gaps, print_sequences, print_alignment_summary};

/*
   -------- Project Structure ---------

   src/
   |-- main.rs          # Main logic for each assignment (funtions for each project that get called by main (i.e. project1(); project2(); etc... )
   |-- types.rs         # Structs and types used across projects
   |-- alignment.rs     # Alignment algorithms
   |-- suffix_trees.rs  # Suffix tree creation and algorithms that use suffix trees
   |-- utils.rs         # Helper functions (I/O stuff, file parsing, stats)
*/

// Modules ----------
mod aligment;
mod types;
mod utils;
//-------------------

// Use Statements ---
//-------------------

fn main() {
    match utils::parse_args() {
        Ok(ProjectSelection::Project1 { fasta_file, alg_select, config_file }) => {
            project1(&fasta_file, &alg_select, &config_file)
        }
        Ok(ProjectSelection::Project2 { fasta_file, alphabet_file }) => {
            //project2_Suffix_Trees(&fasta_file, &alphabet_file)
        }
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    }
}



/*
PROJECT 1:
-------- OUTLINE / TODO ---------

note:
    i-1 is above
    j-1 is left

Gotoh Alg:
    Take 2 sequences as input and a scoring system as input.
    Create a matrix with M, I and D scores (in a struct of some sort) of length [m + 1 , n + 1] where m and n are the sequence lengths
    M tracks match / mismatch score | I tracks insertion (gap in seq1 | i-1, j) | D tracks deletions (gap in seq2 | i, j-1)
    M starts at 0, the other 2 as undefined.

    Populate the first rows and columns:
        M is undefined for first row and column past (0, 0)
        I is undefined for first row, and gap score all the way down first column
        D id gap score across first row, and undefined down first column

    Iterate over all cells and calculate score:
        M is max of diagonals from every matrix + match score if seq1[i] == seq2[j] otherwise mismatch
        I is max of opening gap (M score left + open penalty) and continue gap (I score left + extend penalty)
        D is max of opening gap (M score above + open penalty) and continue gap (D score above + extend penalty)

    Traceback optimal path:
        Take Highest score of the corner (m, n) out of M I and D
        create two strings to store output path
        Until (0,0) is reached:
            Determine how we got to current score (one of M, I, or D) from adjacent cells
                - If current score is in M then check all possible moves from cells above or to the left
                - If current score is in D or I, then check Match score + gap start penalty, and previous self score + extend
            Add seq1[i] and seq2[j] or gap to output strings depending on what path is traced.

    Calculate statistics and format output:
        -sequence lengths
        -optimal path
        -optimal score
        -match, mismatch, gap open and gap extend occurences
        -identity percent
        -gap count

*/

fn project1(fasta_file: &str, alg_select: &str, config_file: &str) {


    println!("INPUT:");
    println!("******\n");

    let score_system = match read_score_config(config_file) {
        Ok(sys) => sys,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };

    let sequences = match parse_sequences_from_fasta(fasta_file) {
        Ok(seq) => seq,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };

    print_sequences(&sequences);

    println!("Alg selected {}\n", alg_select);

    let alignmnet_result = match alg_select {
        "0" => gotoh(&sequences[0], &sequences[1], &score_system),
        "1" => smith_waterman(&sequences[0], &sequences[1], &score_system),
        _ => {
            eprintln!("Invalid algorithm selection. Please provide <0: global, 1: local>");
            std::process::exit(1);
        }
    }
    .unwrap_or_else(|e| {
        eprintln!("{}", e);
        std::process::exit(1);
    });

    print_alignment_summary(&alignmnet_result, &score_system);
}



fn project2_Suffix_Trees(){

}