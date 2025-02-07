//utils.rs

use std::{fs::File, io::{BufRead, BufReader}};

use crate::types::AllignmentStats;

// Helper functions (I/O stuff, file parsing, stats)

// Read all sequences from fasta gile
fn parse_sequences_from_fasta(filename: &str) -> Result<Vec<Vec<u8>>, Error>{
    let file = BufReader::new(File::open(filename)?);
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut current_sequence: Vec<u8> = Vec::new();

    for line in file.lines() {
        let trimmed = line?.trim();

        
        if trimmed.is_empty(){
            continue;
        }

        if trimmed.starts_with('>'){
            if current_sequence.is_empty(){
                sequences.push(current_sequence);
            }
        }
    }
}

// Init scoring params
//fn _________ (filename: Option<&str>) -> Result<ScoringParams, Error>

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
