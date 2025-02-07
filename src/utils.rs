//utils.rs

use std::{fs::File, io::{BufRead, BufReader, Error}, vec};

use crate::types::{AllignmentStats, ScoringSystem};

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
    println!("trying to open: {}", filename);

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

    let alignment_length = seq1.len() as i32;
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

// Format output
// fn format_alignment_output(result: &AlignmentResult) -> String
