use std::env;

use aligment::{gotoh, smith_waterman};
use utils::{parse_sequences_from_fasta, read_score_config, print_alignment, len_without_gaps};

/*
   -------- Project Structure ---------

   src/
   |-- main.rs          # Main logic for each assignment (funtions for each project that get called by main (i.e. project1(); project2(); etc... )
   |-- types.rs         # Structs and types used across projects
   |-- alignment.rs     # Alignment algorithms
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
    project1();
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

fn project1() {

    let args: Vec<String> = env::args().collect();
    
    if args.len() < 3 || args.len() > 4 {
        eprintln!("Incorrect number of arguments provided");
    }

    // Args don't change, so use imutable borrow
    let fasta_file = &args[1];
    let alg_select = &args[2];
    let scoring_config = if args.len() == 4 {
        args[3].clone()
    } else {
        "./parameters.config".to_string()
    };

    println!("INPUT:");
    println!("******\n");
    
    let score_system = match read_score_config(&scoring_config) {
        Ok(sys) => sys,
        Err(e) => return // fix later with proper error handling
    };

    let sequences = match parse_sequences_from_fasta(&fasta_file){
        Ok(seq) => seq,
        Err(e) => return // fix later with proper error handling
    };

    for (i, sequence) in sequences.iter().enumerate() {
        println!(">s{}", i + 1);
        for &byte in sequence {
            print!("{}", byte as char);
        }
        println!("\n");
    }

    print!("Alg selected {}\n", alg_select.as_str());
    let alignmnet_result = match alg_select.as_str() {
        "0" => gotoh(&sequences[0], &sequences[1],&score_system),
        "1" => smith_waterman(&sequences[0], &sequences[1],&score_system),
        _ => {
            eprintln!("Invalid algorithm selection. Please provide <0: global, 1: local>");
            std::process::exit(1);},
    }.unwrap_or_else(|e|{
        eprintln!("{}", e);
        std::process::exit(1);
    });

    println!("\nOUTPUT:");
    println!("********\n");
    
    println!("Scores:    match = {}, mismatch = {}, h = {}, g = {}\n",
        score_system.match_score().unwrap(),
        score_system.mismatch_score().unwrap(),
        score_system.gap_open_score().unwrap(),
        score_system.gap_extend_score().unwrap()
    );

    println!("Sequence 1 = \"s1\", length = {} characters",
        len_without_gaps(alignmnet_result.sequence1()));
    println!("Sequence 2 = \"s2\", length = {} characters\n",
        len_without_gaps(alignmnet_result.sequence2()));

    print_alignment(alignmnet_result.sequence1(), alignmnet_result.sequence2());

    let stats = alignmnet_result.stats();
    println!("Report:");
    println!("Global optimal score = {}", alignmnet_result.score());
    println!("\nNumber of:  matches = {}, mismatches = {}, opening gaps = {}, gap extensions = {}",
        stats.match_count(),
        stats.mismatch_count(),
        stats.gap_open_count(),
        stats.gap_extend_count()
    );

    println!("\nIdentities = {}/{} ({}%), Gaps = {}/{} ({}%)",
    stats.match_count(),
    stats.alignment_length(),
    stats.identity_percent(),
    stats.total_gaps(),
    stats.alignment_length(),
    (stats.total_gaps() as f32 / stats.alignment_length() as f32 * 100.0).round()
    );
}
