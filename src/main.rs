// src/main.rs

// Modules
mod alignment;
mod genome_comparison;
mod mccreight;
mod suffix_tree;
mod types;
mod utils;

// Use Statements
use alignment::{gotoh, smith_waterman};
use mccreight::construct_st_mc;
use std::fs::File;
use std::time::{Duration, Instant}; // Need Duration for project 3 timings
use types::ProjectSelection;
use utils::{
    parse_sequences_from_fasta, print_alignment_summary, print_sequences,
    read_score_config,
};


/*
    Decided to keep old project code runnable with flags (--p1, --p2, --p3).
    Defaults to project 3 if no flag given, cause that's the latest.
*/
fn main() {
    let total_start_time = Instant::now();

    match utils::parse_args() {
        // Project 1: Needs fasta (2 seqs), config file, algorithm choice (0 or 1)
        Ok(ProjectSelection::Project1 {
            fasta_file,
            alg_select,
            config_file,
        }) => project1(&fasta_file, &alg_select, &config_file),

        // Project 2: Needs fasta (1 seq), alphabet file
        Ok(ProjectSelection::Project2 {
            fasta_file,
            alphabet_file,
        }) => project2_suffix_trees(&fasta_file, &alphabet_file),

        // Project 3 (default): Needs alphabet file, 1 or more fasta files
        Ok(ProjectSelection::Project3 {
            alphabet_file,
            fasta_files,
        }) => {
            project3_similarity_matrix(&alphabet_file, &fasta_files);
        }
        Err(e) => {
            eprintln!("Arg parsing error: {}", e);
            std::process::exit(1);
        }
    }

    println!("\nTotal execution time: {:?}", total_start_time.elapsed());
}

// Project 1 alignment
fn project1(fasta_file: &str, alg_select: &str, config_file: &str) {
    println!("Running Project 1...");

    let score_system = match read_score_config(config_file) {
        Ok(sys) => sys,
        Err(e) => {
            eprintln!("Config file error: {}", e);
            std::process::exit(1);
        }
    };

    let sequences = match parse_sequences_from_fasta(fasta_file) {
        Ok(seqs) => {
            if seqs.len() < 2 {
                eprintln!(
                    "Error: Need at least 2 sequences in {} for Project 1.",
                    fasta_file
                );
                std::process::exit(1);
            }
            seqs
        }
        Err(e) => {
            eprintln!("FASTA parse error: {}", e);
            std::process::exit(1);
        }
    };

    println!("Input sequences:");
    print_sequences(&sequences[0..2]);

    // Figure out which alignment to run
    let alignment_result = match alg_select {
        "0" => {
            gotoh(&sequences[0], &sequences[1], &score_system)
        }
        "1" => {
            smith_waterman(&sequences[0], &sequences[1], &score_system)
        }
        _ => {
            eprintln!("Invalid algorithm choice {}. Use 0 for Gotoh (global) or 1 for SW (local).", alg_select);
            std::process::exit(1);
        }
    };

    match alignment_result {
        Ok(alignment) => {
            print_alignment_summary(&alignment, &score_system);
        }
        Err(e) => {
            eprintln!("\nAlignment failed: {}", e);
            std::process::exit(1);
        }
    }
}

// Project 2 suffix tree generation
fn project2_suffix_trees(fasta_file: &str, alphabet_file: &str) {
    println!("Running Project 2...");

    // Parse inputs
    let sequences = match utils::parse_sequences_from_fasta(fasta_file) {
        Ok(seqs) => {
            if seqs.is_empty() {
                eprintln!("Error: No sequences in FASTA file {}.", fasta_file);
                std::process::exit(1);
            }
            seqs
        }
        Err(e) => {
            eprintln!("FASTA read error: {}",  e);
            std::process::exit(1);
        }
    };
    let sequence = &sequences[0];

    let alphabet_order = match utils::parse_alphabet(alphabet_file) {
        Ok(map) => map,
        Err(e) => {
            eprintln!("Alphabet file error: {}", e);
            std::process::exit(1);
        }
    };

    // Setup output dir
    let path = std::path::Path::new(fasta_file);
    let file_stem = path.file_stem().unwrap_or_default().to_string_lossy();
    let output_dir = format!("output/{}_p2_analysis", file_stem);

    if let Err(e) = std::fs::create_dir_all(&output_dir) {
        eprintln!("Couldn't create output dir: {}", e);
        // Continue anyway
    }

    let report_path = format!("{}/report.txt", output_dir);
    let dfs_path = format!("{}/dfs_output.txt", output_dir);
    let postorder_path = format!("{}/postorder_output.txt", output_dir);
    let bwt_path = format!("{}/BWT.txt", output_dir);

    // Build tree
    println!("Building suffix tree...");
    let start_time = std::time::Instant::now();
    let tree = construct_st_mc(&sequence, &alphabet_order);
    let construction_time = start_time.elapsed();
    println!("Tree built in {:?}.", construction_time);

    utils::print_tree_summary(&tree, sequence.len(), construction_time);

    println!("Writing output files...");
    match File::create(&report_path) {
        Ok(mut file) => {
            if let Err(e) =
                utils::generate_report(&mut file, fasta_file, &tree, &sequence, construction_time)
            {
                eprintln!("Error writing report to {}: {}", report_path, e);
            } else {
                 println!("-> {}", report_path); // Output file message
            }
        }
        Err(e) => eprintln!("Error creating report file {}: {}", report_path, e),
    }

    match File::create(&dfs_path) {
        Ok(mut file) => {
            if let Err(e) = suffix_tree::dfs_to_file(&tree, &alphabet_order, &mut file) {
                eprintln!("Error writing DFS output to {}: {}", dfs_path, e);
            } else {
                 println!("-> {}", dfs_path); // Output file message
            }
        }
        Err(e) => eprintln!("Error creating DFS file {}: {}", dfs_path, e),
    }

    match File::create(&postorder_path) {
        Ok(mut file) => {
            if let Err(e) = suffix_tree::postorder_to_file(&tree, &alphabet_order, &mut file) {
                eprintln!("Error writing postorder output to {}: {}", postorder_path, e);
            } else {
                 println!("-> {}", postorder_path); // Output file message
            }
        }
        Err(e) => eprintln!("Error creating postorder file {}: {}", postorder_path, e),
    }

    match File::create(&bwt_path) {
        Ok(mut file) => {
            // Pass original seq, bwt func handles $
            if let Err(e) = suffix_tree::bwt_to_file(&tree, &sequence, &mut file) {
                eprintln!("Error writing BWT to {}: {}", bwt_path, e);
            } else {
                println!("-> {}", bwt_path); // Output file message
            }
        }
        Err(e) => eprintln!("Error creating BWT file {}: {}", bwt_path, e),
    }

    //println!("\nProject 2 analysis complete. Outputs in: {}", output_dir);
}

//-----------------------------------------------------------------------------------

// Project 3 Genome Similarity
fn project3_similarity_matrix(alphabet_file: &str, fasta_files: &[String]) {
    println!("Running Project 3...");
    let start_time = Instant::now();

    // Setup using # and $
    let term1 = b'#';
    let term2 = b'$';
    let mut alphabet_order = match utils::parse_alphabet(alphabet_file) {
        Ok(map) => map,
        Err(e) => {
            eprintln!("Alphabet file error: {}", e);
            std::process::exit(1);
        }
    };
    // Ensure no sepperators in the alphabet
    if alphabet_order.contains_key(&(term1 as char))
        || alphabet_order.contains_key(&(term2 as char))
    {
        eprintln!("Error: Alphabet file contains reserved characters # or $.");
        std::process::exit(1);
    }
    // Add separators to alphabet order
    let next_idx = alphabet_order.len();
    alphabet_order.insert(term1 as char, next_idx);
    alphabet_order.insert(term2 as char, next_idx + 1);

    // Load sequences
    let mut sequences = Vec::new();
    let mut seq_names = Vec::new();
    println!("Loading sequences:");
    for fasta_file in fasta_files {
        print!("  Loading {}... ", fasta_file);

        let loaded_seqs = match utils::parse_sequences_from_fasta(fasta_file) {
            Ok(seqs) => seqs,
            Err(e) => {
                println!(" failed.");
                eprintln!("\nFASTA read error {}: {}", fasta_file, e);
                std::process::exit(1);
            }
        };

        // One sequence per file
        if loaded_seqs.len() == 1 {
            let seq_data = loaded_seqs[0].clone();
            println!(" done ({} bp).", seq_data.len());
            sequences.push(seq_data);
            let path = std::path::Path::new(fasta_file);
            let name = path.file_stem().unwrap_or_default().to_string_lossy().to_string();
            seq_names.push(if name.is_empty() { fasta_file.to_string() } else { name }); // Use filename if stem fails
        } else {
            println!(" failed.");
            eprintln!("\nError: FASTA {} must have exactly one sequence.", fasta_file);
            std::process::exit(1);
        }
    }

    let k = sequences.len();
    if k < 1 {
        eprintln!("Error: No sequences loaded.");
        std::process::exit(1);
    }

    println!("\nLoaded {} sequences:", k);
    for (i, name) in seq_names.iter().enumerate() {
        println!("  {}: {} ({} bp)", i + 1, name, sequences[i].len());
    }

    // Init matrices
    let mut similarity_matrix = vec![vec![0usize; k]; k];
    let mut lcs_length_matrix = vec![vec![0usize; k]; k];

    // Default scoring
    let match_s = 1;
    let mismatch = -2;
    let gap_open = -5;
    let gap_extend = -1;

    let mut total_gst_time = Duration::ZERO;
    let mut total_alignment_time = Duration::ZERO;

    // Compare all pairs (upper triangle including diagonal)
    println!("\nComparing pairs...");
    for i in 0..k {
        // Only compute i <= j
        for j in i..k {
            print!("Comparing {} vs {} ", seq_names[i], seq_names[j]);

            let similarity_score: usize;
            let lcs_length: usize;

            if i == j {
                // For the diagonal similarity is just the length
                similarity_score = sequences[i].len();
                lcs_length = sequences[i].len();
                println!("Self, Sim = {}, LCS = {}", similarity_score, lcs_length);
            } else {
                // Off-diagonal: Calculate using GST method
                match genome_comparison::calculate_similarity(
                    &sequences[i],
                    &sequences[j],
                    &alphabet_order,
                    term1,
                    term2,
                    match_s,
                    mismatch,
                    gap_open,
                    gap_extend,
                ) {
                    Ok((sim, lcs, gst_time, align_time)) => {
                        similarity_score = sim;
                        lcs_length = lcs;
                        total_gst_time += gst_time;
                        total_alignment_time += align_time;

                        // Per-pair timing
                        println!(
                            "Sim = {}, LCS = {} (GST: {:?}, Align: {:?})",
                            similarity_score, lcs_length, gst_time, align_time
                        );
                    }
                    Err(e) => {
                        println!(" failed.");
                        eprintln!("\nError comparing pair ({}, {}): {}", i + 1, j + 1, e);
                        std::process::exit(1);
                    }
                }
            }

            // Fill both sides of the matrix since it's symmetric
            similarity_matrix[i][j] = similarity_score;
            similarity_matrix[j][i] = similarity_score;
            lcs_length_matrix[i][j] = lcs_length;
            lcs_length_matrix[j][i] = lcs_length;
        }
    }

    // Print results
    println!("\nSimilarity Matrix (D = a+b+c):");
    if let Err(e) = utils::print_similarity_matrix(&similarity_matrix, &seq_names) {
        eprintln!("Error printing similarity matrix: {}", e);
    }

    println!("\nLCS Length Matrix (b):");
    if let Err(e) = utils::print_similarity_matrix(&lcs_length_matrix, &seq_names) {
        eprintln!("Error printing LCS length matrix: {}", e);
    }

    // Timings summary
    println!("\nTimings:");
    println!("Total GST time:          {:?}", total_gst_time);
    println!("Total Alignment time:    {:?}", total_alignment_time);
}