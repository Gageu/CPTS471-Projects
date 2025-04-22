// utils.rs
use crate::suffix_tree::{StEdge, StNode, SuffixTree};
use crate::types::{Alignment, AllignmentStats, ProjectSelection, ScoringSystem};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Error as IoError, ErrorKind, Result as IoResult, Write};
use std::fs::File;
use std::path::Path;

// Parses command line args for project selection
pub fn parse_args() -> Result<ProjectSelection, String> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        return Err(
            "Not enough args. Usage: <program> [--p1|--p2|--p3] [args] \nor default (p3):\n <program> <alphabet> <fasta1> [fasta2 ...]".to_string(),
        );
    }

    match args[1].as_str() {
        "--p1" => {
            if args.len() < 4 || args.len() > 5 {
                return Err("Usage P1: cargo run -- --p1 <fasta> <alg:0|1> [config]".to_string());
            }
            let fasta_file = args[2].clone();
            let alg_select = args[3].clone();
            let config_file = args.get(4).cloned().unwrap_or("./parameters.config".to_string()); // Default config
            Ok(ProjectSelection::Project1 { fasta_file, alg_select, config_file })
        }
        "--p2" => {
            if args.len() != 4 {
                return Err("Usage P2: cargo run -- --p2 <fasta> <alphabet>".to_string());
            }
            Ok(ProjectSelection::Project2 { fasta_file: args[2].clone(), alphabet_file: args[3].clone() })
        }
        "--p3" => {
            if args.len() < 4 {
                 return Err("Usage P3: cargo run -- --p3 <alphabet> <fasta1> [fasta2 ...]".to_string());
             }
            Ok(ProjectSelection::Project3 { alphabet_file: args[2].clone(), fasta_files: args[3..].to_vec() })
        }
        _ => { // Default to Project 3
            if args.len() >= 3 {
                let alphabet_file = args[1].clone();
                let fasta_files = args[2..].to_vec();

                // Quick check if first arg looks like a flag
                if alphabet_file.starts_with('-') && !Path::new(&alphabet_file).exists() {
                    return Err(format!("Unknown option: {}.", alphabet_file));
                }

                println!("No project flag found, running Project 3.");

                Ok(ProjectSelection::Project3 { alphabet_file, fasta_files })
            } else {
                Err(format!("Invalid args for default P3."))
            }
        }
    }
}

// Reads sequences from a FASTA file. Returns Vec<Vec<u8>>.
pub fn parse_sequences_from_fasta(filename: &str) -> Result<Vec<Vec<u8>>, IoError> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut current_sequence: Option<Vec<u8>> = None;

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty() { continue; }

        if trimmed.starts_with('>') {
            // Finish previous sequence if any
            if let Some(seq) = current_sequence.take() {
                if !seq.is_empty() { sequences.push(seq); }
            }
            current_sequence = Some(Vec::new()); // Start new one
        } else if let Some(seq) = current_sequence.as_mut() {
            // Append sequence data (no whitespace)
            seq.extend(trimmed.bytes().filter(|b| !b.is_ascii_whitespace()));
        }
        // Ignore lines before first '>'
    }

    // Add the last sequence
    if let Some(seq) = current_sequence.take() {
        if !seq.is_empty() { sequences.push(seq); }
    }

    if sequences.is_empty() {
        Err(IoError::new(ErrorKind::InvalidData, format!("No sequences found in {}", filename)))
    } else {
        Ok(sequences)
    }
}

// Reverses a sequence (for P3 prefix alignment)
pub fn reverse_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().cloned().collect()
}

// Reads scores (match, mismatch, h, g) from config file.
pub fn read_score_config(filename: &str) -> Result<ScoringSystem, String> {
    let file = File::open(filename).map_err(|e| format!("Couldn't open config: {}", e))?;
    let reader = BufReader::new(file);

    let mut match_score = None;
    let mut mismatch_score = None;
    let mut gap_open_score = None;
    let mut gap_extend_score = None;

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| format!("Error reading config line {}", e))?;
        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.is_empty(){ continue; }

        if parts.len() != 2 {
            eprintln!("Warning: Skipping bad line in config");
            continue;
        }

        let param_name = parts[0];
        let value: i32 = parts[1].parse().map_err(|e| format!("Invalid score value: {}", e))?;

        match param_name {
            "match" => match_score = Some(value),
            "mismatch" => mismatch_score = Some(value),
            "h" | "gap_open" => gap_open_score = Some(value),
            "g" | "gap_extend" => gap_extend_score = Some(value),
            _ => eprintln!("Warning: Unknown param {}", param_name),
        }
    }

    // Check all required scores are present
    let m = match_score.ok_or_else(|| format!("Missing match score"))?;
    let mis = mismatch_score.ok_or_else(|| format!("Missing mismatch score"))?;
    let h = gap_open_score.ok_or_else(|| format!("Missing gap_open)score",))?;
    let g = gap_extend_score.ok_or_else(|| format!("Missing (ap_extend score"))?;


    Ok(ScoringSystem::new(Some(m), Some(mis), None, Some(h), Some(g)))
}

// Calculates stats for an alignment.
pub fn calculate_alignmnet_stats(seq1_aligned: &str, seq2_aligned: &str) -> AllignmentStats {
    let mut matches = 0;
    let mut mismatches = 0;
    let mut gap_opens = 0;
    let mut gap_extends = 0;
    let mut total_gaps = 0;
    let align_len = seq1_aligned.len();

    let mut in_gap1 = false;
    let mut in_gap2 = false;

    for (c1, c2) in seq1_aligned.chars().zip(seq2_aligned.chars()) {
        match (c1, c2) {
            (a, b) if a == b && a != '-' => { // Match
                matches += 1;
                in_gap1 = false; in_gap2 = false;
            }
            (a, b) if a != '-' && b != '-' => { // Mismatch
                mismatches += 1;
                in_gap1 = false; in_gap2 = false;
            }
            ('-', b) if b != '-' => { // Gap in Seq1
                total_gaps += 1;
                if !in_gap1 { gap_opens += 1; in_gap1 = true; }
                else { gap_extends += 1; }
                in_gap2 = false;
            }
            (a, '-') if a != '-' => { // Gap in Seq2
                total_gaps += 1;
                if !in_gap2 { gap_opens += 1; in_gap2 = true; }
                else { gap_extends += 1; }
                in_gap1 = false;
            }
            ('-', '-') => { // Parallel gap (should be rare)
                 total_gaps += 2;
                 if !in_gap1 { gap_opens += 1; in_gap1 = true; } else { gap_extends += 1; }
                 if !in_gap2 { gap_opens += 1; in_gap2 = true; } else { gap_extends += 1; }
                 eprintln!("Warning: Found parallel gap ('-','-').");
            }
            _ => { /* Should not happen */ }
        }
    }

    let identity = if align_len > 0 { (matches as f32 / align_len as f32) * 100.0 } else { 0.0 };

    AllignmentStats::new(align_len as i32, matches, mismatches, gap_opens, gap_extends, total_gaps, identity)
}

// Creates Alignment struct from parts.
pub fn assemble_alignment(score: i32, s1_aligned: String, s2_aligned: String) -> Alignment {
    let stats = calculate_alignmnet_stats(&s1_aligned, &s2_aligned);
    Alignment::new(score, s1_aligned, s2_aligned, stats)
}

// Prints alignment nicely with wrapping, indices, match bars.
pub fn print_alignment(seq1_aligned: &str, seq2_aligned: &str) {
    let len = seq1_aligned.len();
    if len == 0 { println!("(Empty Alignment)"); return; }
    let line_size = 60;
    let mut seq1_orig_idx = 0;
    let mut seq2_orig_idx = 0;

    for line_start in (0..len).step_by(line_size) {
        let line_end = std::cmp::min(line_start + line_size, len);
        let s1_slice = &seq1_aligned[line_start..line_end];
        let s2_slice = &seq2_aligned[line_start..line_end];

        // Indices for this line
        let start_idx1 = seq1_orig_idx + 1;
        let start_idx2 = seq2_orig_idx + 1;
        let end_idx1 = seq1_orig_idx + s1_slice.chars().filter(|&c| c != '-').count();
        let end_idx2 = seq2_orig_idx + s2_slice.chars().filter(|&c| c != '-').count();

        // Seq1 line
        println!("s1 {:<5} {} {}", start_idx1, s1_slice, end_idx1);

        // Match bar
        print!("         ");
        for (c1, c2) in s1_slice.chars().zip(s2_slice.chars()) {
            match (c1, c2) {
                (a, b) if a == b && a != '-' => print!("|"), // Match
                (a, b) if a != '-' && b != '-' => print!("."), // Mismatch
                _ => print!(" "), // Gap
            }
        }
        println!();

        // Seq2 line
        println!("s2 {:<5} {} {}", start_idx2, s2_slice, end_idx2);
        println!(); // Blank line

        // Update original indices for next block
        seq1_orig_idx = end_idx1;
        seq2_orig_idx = end_idx2;
    }
}

// Prints summary report for P1 alignment.
pub fn print_alignment_summary(alignment: &Alignment, scoring: &ScoringSystem) {
    println!("\nOUTPUT:");
    println!("******\n");

    let match_s = scoring.match_score().unwrap_or(0);
    let mismatch = scoring.mismatch_score().unwrap_or(0);
    let h = scoring.gap_open_score().unwrap_or(0);
    let g = scoring.gap_extend_score().unwrap_or(0);

    println!("Scores:    match = {}, mismatch = {}, h = {}, g = {}", match_s, mismatch, h, g);

    let len1_orig = alignment.sequence1().chars().filter(|&c| c != '-').count();
    let len2_orig = alignment.sequence2().chars().filter(|&c| c != '-').count();
    println!("\nSeq 1 len = {}", len1_orig);
    println!("Seq 2 len = {}\n", len2_orig);

    print_alignment(alignment.sequence1(), alignment.sequence2());

    let stats = alignment.stats();
    println!("Report:");
    println!("Score = {}", alignment.score());
    println!("Length = {}", stats.alignment_length());
    println!(
        "Matches = {}, Mismatches = {}, Gap Opens = {}, Gap Extends = {}",
        stats.match_count(), stats.mismatch_count(), stats.gap_open_count(), stats.gap_extend_count()
    );
    let align_len = stats.alignment_length();
    let gap_percent = if align_len > 0 { stats.total_gaps() as f32 / align_len as f32 * 100.0 } else { 0.0 };

    println!(
        "Identity = {}/{} ({:.1}%), Gaps = {}/{} ({:.1}%)",
        stats.match_count(), align_len, stats.identity_percent(),
        stats.total_gaps(), align_len, gap_percent
    );
}


// Parses alphabet file into map.
pub fn parse_alphabet(filename: &str) -> Result<HashMap<char, usize>, String> {
    let content = std::fs::read_to_string(filename).map_err(|e| format!("Couldn't read alphabet: {}", e))?;
    let mut order = HashMap::new();
    let mut index = 0;
    let mut seen = std::collections::HashSet::new();

    for symbol_str in content.split_whitespace() {
        if symbol_str.chars().count() == 1 {
            let symbol = symbol_str.chars().next().unwrap();
            if seen.insert(symbol) { // Only add first occurrence for ordering
                order.insert(symbol, index);
                index += 1;
            }
        } else if !symbol_str.is_empty() {
            eprintln!("Warning: Ignoring entry {} in alphabet.", symbol_str);
        }
    }

    if order.is_empty() {
        return Err(format!("No valid symbols found in alphabet."));
    }

    Ok(order)
}

// Writes a slice formatted into columns.
pub fn write_vector_formatted<T: std::fmt::Display>(
    data: &[T],
    line_width: usize, // Items per line
    mut writer: impl Write,
) -> IoResult<()> {
    if data.is_empty() { return Ok(()); }
    if line_width == 0 { return Err(IoError::new(ErrorKind::InvalidInput, "line_width cannot be zero")); }

    for (i, item) in data.iter().enumerate() {
        write!(writer, "{} ", item)?; // Simple space separation
        if (i + 1) % line_width == 0 && i != data.len() - 1 {
            writeln!(writer)?;
        }
    }
    writeln!(writer)?; // Final newline
    Ok(())
}

// Estimate of suffix tree memory usage.
pub fn calculate_tree_memory_usage(tree: &SuffixTree) -> usize {
    let mut total_size = std::mem::size_of::<SuffixTree>();
    total_size += tree.nodes.capacity() * std::mem::size_of::<StNode>();

    // Size allocated for children Vecs + HashSets within nodes
    for node in &tree.nodes {
        total_size += node.children.capacity() * std::mem::size_of::<(char, StEdge)>();
        total_size += node.origin_sequences.capacity() * std::mem::size_of::<usize>(); // Approximation
        total_size += std::mem::size_of::<std::collections::HashSet<usize>>(); // Base HashSet size
    }
    total_size
}

// Generates the P2 report file.
pub fn generate_report(
    file: &mut File,
    fasta_file: &str,
    tree: &SuffixTree,
    seq: &[u8],
    construction_time: std::time::Duration,
) -> IoResult<()> {

    writeln!(file, "SUFFIX TREE REPORT")?;
    writeln!(file, "==================")?;
    writeln!(file, "Input FASTA: {}", fasta_file)?;
    let seq_len_orig = if seq.last() == Some(&b'$') { seq.len() - 1 } else { seq.len() }; // Adjust if $ was added
    writeln!(file, "Sequence length: {} bp", seq_len_orig)?;
    writeln!(file, "Construction time: {:?}", construction_time)?;
    writeln!(file)?;

    // Tree Stats
    crate::suffix_tree::write_tree_stats(tree, fasta_file, file)?;

    // Longest Repeat
    crate::suffix_tree::longest_repeat_to_file(tree, seq, file)?;

    // Memory Estimate
    let tree_size_bytes = calculate_tree_memory_usage(tree);
    writeln!(file, "\n--- Memory Usage (Estimate) ---")?;
    writeln!(file, "Input sequence: ~{} bytes", seq_len_orig)?; // Use original length
    writeln!(file, "Suffix tree struct: ~{} bytes ({:.2} MiB)", tree_size_bytes, tree_size_bytes as f64 / (1024.0 * 1024.0))?;
    let size_per_bp = if seq_len_orig > 0 { tree_size_bytes as f64 / seq_len_orig as f64 } else { 0.0 };
    writeln!(file, "Approx size per bp: {:.2} bytes/bp", size_per_bp)?;

    writeln!(file, "\n--- End Report ---")?;
    Ok(())
}

// Prints a brief console summary of the tree.
pub fn print_tree_summary(tree: &SuffixTree, seq_len: usize, construction_time: std::time::Duration) {
     let total_nodes = tree.nodes.len();
     let internal_count = tree.nodes.iter().filter(|n| !n.children.is_empty() && n.id != tree.root).count();
     let leaf_count = tree.nodes.iter().filter(|n| n.children.is_empty()).count();

    println!(
        "Tree: {} nodes ({} internal, {} leaves). Seq Len: {} bp. Time: {:?}",
        total_nodes, internal_count, leaf_count, seq_len, construction_time
    );
}

// Prints the P3 similarity matrix.
pub fn print_similarity_matrix(matrix: &Vec<Vec<usize>>, seq_names: &[String]) -> IoResult<()> {
    let k = matrix.len();
    if k == 0 { println!("(Empty matrix)"); return Ok(()); }

    // Determine column width for formatting
    let max_name_width = seq_names.iter().map(|s| s.len()).max()
        .unwrap_or_else(|| format!("Seq{}", k).len());
    let max_score_width = matrix.iter().flat_map(|row| row.iter())
                        .map(|&score| score.to_string().len()).max().unwrap_or(1);
    let col_width = std::cmp::max(max_name_width, max_score_width) + 2; // Pad

    // Header row
    print!("{:<width$}", "", width = col_width);
    for j in 0..k {
        let name = seq_names.get(j).map_or_else(|| format!("Seq{}", j + 1), |s| s.clone());
        print!("{:<width$}", name, width = col_width);
    }
    println!();

    // Separator
    print!("{:-<width$}", "", width = col_width);
    for _ in 0..k { print!("{:-<width$}", "", width = col_width); }
    println!();

    // Matrix rows
    for i in 0..k {
        let name = seq_names.get(i).map_or_else(|| format!("Seq{}", i + 1), |s| s.clone());
        print!("{:<width$}|", name, width = col_width - 1); // Row name + separator
        for j in 0..k {
            print!("{:<width$}", matrix[i][j], width = col_width);
        }
        println!();
    }
    Ok(())
}

// Prints sequences from Vec<Vec<u8>> with simple headers.
pub fn print_sequences(sequences: &[Vec<u8>]) {
    for (i, sequence) in sequences.iter().enumerate() {
        println!(">s{}", i + 1);
        let seq_str = String::from_utf8_lossy(sequence);
        // Wrap lines roughly
        const WRAP_LEN: usize = 70;
        for (n, c) in seq_str.chars().enumerate() {
             print!("{}", c);
             if (n + 1) % WRAP_LEN == 0 && n != seq_str.len() - 1 { println!(); }
         }
         println!("\n"); // Blank line between seqs
    }
}