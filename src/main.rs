use std::{fs::File};
use aligment::{gotoh, smith_waterman};
use mccreight::construct_st_mc;
use types::ProjectSelection;
use utils::{
    parse_args, parse_sequences_from_fasta, print_alignment_summary, print_sequences, read_score_config
};
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
mod suffix_tree;
mod mccreight;
mod types;
mod utils;
//-------------------

// Use Statements ---
//-------------------

fn main() {
    match parse_args() {
        Ok(ProjectSelection::Project1 { fasta_file, alg_select, config_file }) => {
            project1(&fasta_file, &alg_select, &config_file)
        }
        Ok(ProjectSelection::Project2 { fasta_file, alphabet_file }) => {
            project2_suffix_trees(&fasta_file, &alphabet_file)
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


fn project2_suffix_trees(fasta_file: &str, alphabet_file: &str) {
    // Parse the input files
    let sequences = match utils::parse_sequences_from_fasta(fasta_file) {
        Ok(seq) => seq,
        Err(e) => {
            eprintln!("Error reading FASTA file: {}", e);
            std::process::exit(1);
        }
    };

    if sequences.is_empty() {
        eprintln!("No sequence found in FASTA file.");
        std::process::exit(1);
    }

    let alphabet_order = match utils::parse_alphabet(alphabet_file) {
        Ok(map) => map,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };
    
    // Prepare output directory and file paths
    let path = std::path::Path::new(fasta_file);
    let file_stem = path.file_stem().unwrap_or_default().to_string_lossy();
    
    let output_dir = format!("output/{}_analysis", file_stem);
    
    // Create output directory structure
    if let Err(e) = std::fs::create_dir_all(&output_dir) {
        eprintln!("Error creating output directory: {}", e);
    }
    
    let report_path = format!("{}/report.txt", output_dir);
    let dfs_path = format!("{}/dfs_output.txt", output_dir);
    let postorder_path = format!("{}/postorder_output.txt", output_dir);
    let bwt_path = format!("{}/BWT.txt", output_dir);
    
    // Build the suffix tree
    println!("Building suffix tree...");
    let start_time = std::time::Instant::now();
    let tree = construct_st_mc(&sequences[0], &alphabet_order);
    let construction_time = start_time.elapsed();
    
    // Print minimal summary to console
    utils::print_minimal_summary(&tree, sequences[0].len(), construction_time);
    
    // Generate comprehensive report
    match File::create(&report_path) {
        Ok(mut file) => {
            if let Err(e) = utils::generate_comprehensive_report(
                &mut file, fasta_file, &tree, &sequences[0], construction_time
            ) {
                eprintln!("Error writing report: {}", e);
            } else {
                println!("Report written to {}", report_path);
            }
        },
        Err(e) => eprintln!("Error creating report file: {}", e)
    }
    
    // Write DFS output
    match File::create(&dfs_path) {
        Ok(mut file) => {
            if let Err(e) = suffix_tree::write_dfs_to_file(&tree, &alphabet_order, &mut file) {
                eprintln!("Error writing DFS output: {}", e);
            } else {
                println!("DFS output written to {}", dfs_path);
            }
        },
        Err(e) => eprintln!("Error creating DFS file: {}", e)
    }
    
    // Write postorder output
    match File::create(&postorder_path) {
        Ok(mut file) => {
            if let Err(e) = suffix_tree::write_postorder_to_file(&tree, &alphabet_order, &mut file) {
                eprintln!("Error writing postorder output: {}", e);
            } else {
                println!("Postorder output written to {}", postorder_path);
            }
        },
        Err(e) => eprintln!("Error creating postorder file: {}", e)
    }
    
    // Write BWT output
    match File::create(&bwt_path) {
        Ok(mut file) => {
            if let Err(e) = suffix_tree::write_bwt_to_file(&tree, &sequences[0], &mut file) {
                eprintln!("Error writing BWT: {}", e);
            } else {
                println!("BWT written to {}", bwt_path);
            }
        },
        Err(e) => eprintln!("Error creating BWT file: {}", e)
    }
    
    println!("All outputs saved to directory: {}", output_dir);
}