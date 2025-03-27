//utils.rs

use std::{cmp, fs::File, io::{BufRead, BufReader, Error}, vec};

use crate::types::{AllignmentStats, Alignment, ScoringSystem, ProjectSelection};

// Helper functions (I/O stuff, file parsing, stats)

// Read all sequences from fasta gile
pub fn parse_sequences_from_fasta(filename: &str) -> Result<Vec<Vec<u8>>, Error>{
    let file = BufReader::new(File::open(filename)?);
    let mut sequences: Vec<Vec<u8>> = vec![];
    let mut current_sequence: Vec<u8> = vec![];

    for line in file.lines() {
        let line = line?;
        let trimmed = line.trim();

        
        if trimmed.is_empty(){
            continue;
        }

        if trimmed.starts_with('>'){
            if !current_sequence.is_empty(){
                sequences.push(current_sequence);
                current_sequence = vec![];
            }
        } else {
            current_sequence.extend(
                trimmed
                .bytes()
                .filter(|byte| !byte.is_ascii_whitespace())
            );
        }
    }

    if !current_sequence.is_empty(){
        sequences.push(current_sequence);
    }

    Ok(sequences)
}

// Init scoring params
pub fn read_score_config(filename: &str) -> Result<ScoringSystem, String>{
    /* TODO:
        - Handle case where there is only a gap penalty specified and not seperate extend and open scores
     */

    let file = match File::open(filename) {
        Ok(file) => file,
        Err(e) => return Err(e.to_string())
    };
    let reader = BufReader::new(file);
    
    let mut match_score = None;
    let mut mismatch_score = None;
    let mut gap_open_score = None;
    let mut gap_extend_score = None;
    
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        
        if parts.len() < 2 {
            continue;
        }
        
        let param_name = parts[0];
        let value = match parts[1].parse::<i32>() {
            Ok(val) => val,
            Err(e) => return Err(e.to_string())
        };
        
        match param_name {
            "match" => match_score = Some(value),
            "mismatch" => mismatch_score = Some(value),
            "h" => gap_open_score = Some(value),
            "g" => gap_extend_score = Some(value),
            _ => {}
        }
    }
    
    Ok(ScoringSystem::new(
        match_score, 
        mismatch_score, 
        None, 
        gap_open_score, 
        gap_extend_score)
    )
}

// Allignment stats
pub fn calculate_alignmnet_stats(seq1: &str, seq2: &str) -> AllignmentStats {
    let mut match_count = 0;
    let mut mismatch_count = 0;
    let mut gap_open_count = 0;
    let mut gap_extend_count = 0;
    let mut total_gaps = 0;

    //Track gaps for both sequences spereately. I don't think parallel gaps could ever be part of an
    //optimal alignmnet for any of our current scoring systems, but who knows.
    let mut seq1_gap_open = false;
    let mut seq2_gap_open = false;

    for (c1, c2) in seq1.chars().zip(seq2.chars()) {
        match (c1, c2) {
            (a, b) if a == b && a != '-' => {
                match_count += 1;
                seq1_gap_open = false;
                seq2_gap_open = false;
            }

            (a, b) if a != '-' && b != '-' => {
                mismatch_count += 1;
                seq1_gap_open = false;
                seq2_gap_open = false;
            }

            ('-', _) => {
                total_gaps += 1;
                if !seq1_gap_open {
                    gap_open_count += 1;
                    seq1_gap_open = true;
                } else {
                    gap_extend_count += 1;
                }
            }

            (_, '-') => {
                total_gaps += 1;
                if !seq2_gap_open {
                    gap_open_count += 1;
                    seq2_gap_open = true;
                } else {
                    gap_extend_count += 1;
                }
            }

            _ => panic!("Alignment contained invalid characters"),
        }
    }

    let alignment_length = cmp::max(
        len_without_gaps(seq1),
        len_without_gaps(seq2));

    let identity_percent = (match_count as f32 / alignment_length as f32) * 100.0;

    AllignmentStats::new(
        alignment_length,
        match_count,
        mismatch_count,
        gap_open_count,
        gap_extend_count,
        total_gaps,
        identity_percent,
    )
}

pub fn assemble_alignment(score: i32, s1: String, s2: String) -> Alignment {
    let stats = calculate_alignmnet_stats(&s1, &s2);
    Alignment::new(score, s1, s2, stats)
}

// Print two sequences as described by the assignment description
// bars connecting matches, index at end of every line and wrpeed every 60 characters
pub fn print_alignment(seq1: &str, seq2: &str) {
    let len = seq1.len();
    let line_size = 60;
    
    for line_start in (0..len).step_by(line_size) {
        let line_end = std::cmp::min(line_start + line_size, len);
        
        println!("s1 {:<5} {} {}", 
            line_start + 1,
            &seq1[line_start..line_end],
            line_end
        );
        
        //Allign the middle bar
        print!("         ");
        for i in line_start..line_end {
            //if they match and aren;t a gap
            if seq1.chars().nth(i) == seq2.chars().nth(i) && seq1.chars().nth(i) != Some('-') {
                print!("|");
            } else {
                print!(" ");
            }
        }
        println!();
        
        println!("s2 {:<5} {} {}\n", 
            line_start + 1,
            &seq2[line_start..line_end],
            line_end
        );
    }
}

pub fn len_without_gaps(seq: &str) -> i32{
    seq.chars()
    .filter(|c| *c != '-')
    .count() as i32
}

pub fn parse_args() -> Result<ProjectSelection, String> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        return Err("Not enough arguments. Use --p1 or provide input for the latest project.".to_string());
    }

    print!("{}", args[1]);

    if args[1] == "--p1" {
        if args.len() < 4 || args.len() > 5 {
            return Err("Usage for project 1: <program> --p1 <fasta_file> <algorithm_select> [score_config_file]".to_string());
        }

        let fasta_file = args[2].clone();
        let alg_select = args[3].clone();
        let config_file = if args.len() == 5 {
            args[4].clone()
        } else {
            "./parameters.config".to_string()
        };

        Ok(ProjectSelection::Project1 {
            fasta_file,
            alg_select,
            config_file,
        })
    } else {
        if args.len() != 3 {
            return Err("Usage for project 2 (default): <program> <fasta_file> <alphabet_file>".to_string());
        }

        let fasta_file = args[1].clone();
        let alphabet_file = args[2].clone();

        Ok(ProjectSelection::Project2 {
            fasta_file,
            alphabet_file,
        })
    }
}


pub fn print_sequences(sequences: &[Vec<u8>]) {
    for (i, sequence) in sequences.iter().enumerate() {
        println!(">s{}", i + 1);
        for &byte in sequence {
            print!("{}", byte as char);
        }
        println!("\n");
    }
}


pub fn print_alignment_summary(alignment: &Alignment, scoring: &ScoringSystem) {
    println!("\nOUTPUT:");
    println!("********\n");

    println!(
        "Scores:    match = {}, mismatch = {}, h = {}, g = {}\n",
        scoring.match_score().unwrap(),
        scoring.mismatch_score().unwrap(),
        scoring.gap_open_score().unwrap(),
        scoring.gap_extend_score().unwrap()
    );

    println!(
        "Sequence 1 = \"s1\", length = {} characters",
        len_without_gaps(alignment.sequence1())
    );
    println!(
        "Sequence 2 = \"s2\", length = {} characters\n",
        len_without_gaps(alignment.sequence2())
    );

    print_alignment(alignment.sequence1(), alignment.sequence2());

    let stats = alignment.stats();
    println!("Report:");
    println!("Global optimal score = {}", alignment.score());
    println!(
        "\nNumber of:  matches = {}, mismatches = {}, opening gaps = {}, gap extensions = {}",
        stats.match_count(),
        stats.mismatch_count(),
        stats.gap_open_count(),
        stats.gap_extend_count()
    );

    println!(
        "\nIdentities = {}/{} ({}%), Gaps = {}/{} ({}%)",
        stats.match_count(),
        stats.alignment_length(),
        stats.identity_percent(),
        stats.total_gaps(),
        stats.alignment_length(),
        (stats.total_gaps() as f32 / stats.alignment_length() as f32 * 100.0).round()
    );
}

