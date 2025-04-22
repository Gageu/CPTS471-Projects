// types.rs

// Represents a computed alignment result for Project 1
pub struct Alignment {
    score: i32,
    sequence1: String, // Aligned sequence 1
    sequence2: String, // Aligned sequence 2
    stats: AllignmentStats,
}

impl Alignment {
    // Constructor requires all fields
    pub fn new(score: i32, sequence1: String, sequence2: String, stats: AllignmentStats) -> Self {
        Self { score, sequence1, sequence2, stats }
    }

    // Getters
    pub fn score(&self) -> i32 { self.score }
    pub fn sequence1(&self) -> &str { &self.sequence1 }
    pub fn sequence2(&self) -> &str { &self.sequence2 }
    pub fn stats(&self) -> &AllignmentStats { &self.stats }
}

// Statistics derived from an alignment
pub struct AllignmentStats {
    alignment_length: i32, // Length of the aligned strings (including gaps)
    match_count: i32,
    mismatch_count: i32,
    gap_open_count: i32,
    gap_extend_count: i32,
    total_gaps: i32, // Total number of - characters in both sequences
    identity_percent: f32, // (match_count / alignment_length) * 100
}

impl AllignmentStats {
    pub fn new(
        alignment_length: i32,
        match_count: i32,
        mismatch_count: i32,
        gap_open_count: i32,
        gap_extend_count: i32,
        total_gaps: i32,
        identity_percent: f32,
    ) -> Self {
        Self {
            alignment_length,
            match_count,
            mismatch_count,
            gap_open_count,
            gap_extend_count,
            total_gaps,
            identity_percent,
        }
    }

    // Getters
    pub fn alignment_length(&self) -> i32 { self.alignment_length }
    pub fn match_count(&self) -> i32 { self.match_count }
    pub fn mismatch_count(&self) -> i32 { self.mismatch_count }
    pub fn gap_open_count(&self) -> i32 { self.gap_open_count }
    pub fn gap_extend_count(&self) -> i32 { self.gap_extend_count }
    pub fn total_gaps(&self) -> i32 { self.total_gaps }
    pub fn identity_percent(&self) -> f32 { self.identity_percent }
}

// Cell structure for Gotoh/Smith-Waterman DP matrices with affine gaps
#[derive(Copy, Clone, Debug)] // Added Debug
pub(crate) struct MDICell {
    pub m_score: i32, // Score ending with Match/Mismatch
    pub d_score: i32, // Score ending with Deletion (gap in seq2)
    pub i_score: i32, // Score ending with Insertion (gap in seq1)
}

impl MDICell {
    // Default values representing negative infinity (unreachable states)
    pub fn default() -> Self {
        MDICell {
            m_score: i32::MIN / 2, // Use MIN/2 to prevent overflow during additions
            d_score: i32::MIN / 2,
            i_score: i32::MIN / 2,
        }
    }

    // Helper to determine the best operation ending at this cell (used in original traceback)
    // Note: For affine gaps, traceback needs more careful source determination based on scores.
    pub fn best_op(&self) -> char {
        if self.m_score >= self.d_score && self.m_score >= self.i_score {
            'm' // Prefer M in ties
        } else if self.d_score >= self.i_score {
            'd' // Prefer D over I in ties
        } else {
            'i'
        }
    }

    // Creates a DP matrix initialized with default (negative infinity) values.
    // Sets M[0][0] and potentially other boundary initializations based on args.
    pub fn default_matrix(
        m: usize, // rows (length of seq1)
        n: usize, // cols (length of seq2)
        init_m00: i32, // Initial M score at [0][0] (usually 0)
        init_d_boundary: i32, // Initial D score along col 0 (MIN/2 for NW, 0 for SW)
        init_i_boundary: i32, // Initial I score along row 0 (MIN/2 for NW, 0 for SW)
    ) -> Vec<Vec<MDICell>> {
        let mut matrix = vec![vec![MDICell::default(); n + 1]; m + 1];

        // Initialize top-left corner [0][0]
        matrix[0][0].m_score = init_m00;
        // D and I at [0][0] should represent starting a gap immediately.
        // For NW/Gotoh, these might be -inf or related to gap_open.
        // For SW/P3-Variant (starting from 0), they should likely be -inf.
        // Let's keep them MIN/2 by default unless specific init needed.
        // The alignment functions handle specific boundary setup after this.
        matrix[0][0].d_score = init_m00; // Or perhaps MIN/2? Check Gotoh init details. Set to MIN/2 for safety.
        matrix[0][0].i_score = init_m00; // Or perhaps MIN/2? Check Gotoh init details. Set to MIN/2 for safety.
        matrix[0][0].d_score = i32::MIN / 2;
        matrix[0][0].i_score = i32::MIN / 2;


        // Initialize boundaries with placeholder values (actual penalties applied in algos)
        // This part is mostly superseded by explicit boundary loops in alignment functions.
        for i in 1..=m {
             matrix[i][0].d_score = init_d_boundary; // Placeholder, might be overwritten
        }
        for j in 1..=n {
             matrix[0][j].i_score = init_i_boundary; // Placeholder, might be overwritten
        }

        matrix
    }
}

// Defines the scoring parameters for alignment algorithms
#[derive(Default, Debug)] // Added Debug
pub struct ScoringSystem {
    match_score: Option<i32>,
    mismatch_score: Option<i32>,
    gap_score: Option<i32>,      // For simple linear gap penalty (not used in affine)
    gap_open_score: Option<i32>, // h
    gap_extend_score: Option<i32>,// g
}

impl ScoringSystem {
    // Constructor allows setting all optional scores
    pub fn new(
        match_score: Option<i32>,
        mismatch_score: Option<i32>,
        gap_score: Option<i32>,
        gap_open_score: Option<i32>,
        gap_extend_score: Option<i32>,
    ) -> Self {
        Self { match_score, mismatch_score, gap_score, gap_open_score, gap_extend_score }
    }

    pub fn match_score(&self) -> Option<i32> { self.match_score }
    pub fn mismatch_score(&self) -> Option<i32> { self.mismatch_score }
    pub fn gap_score(&self) -> Option<i32> { self.gap_score }
    pub fn gap_open_score(&self) -> Option<i32> { self.gap_open_score }
    pub fn gap_extend_score(&self) -> Option<i32> { self.gap_extend_score }

    // Unpacks the necessary scores for the affine gap model used in P1/P3.
    // Returns Err if any required score is missing.
    pub fn unpack(&self) -> Result<(i32, i32, i32, i32), String> {
        Ok((
            self.match_score().ok_or("Missing match score ('match') in config")?,
            self.mismatch_score().ok_or("Missing mismatch score ('mismatch') in config")?,
            self.gap_open_score().ok_or("Missing gap open score ('h') in config")?,
            self.gap_extend_score().ok_or("Missing gap extend score ('g') in config")?,
        ))
    }
}

// Enum to represent the selected project and its specific arguments
pub enum ProjectSelection {
    Project1 {
        fasta_file: String,
        alg_select: String, // "0" or "1"
        config_file: String,
    },
    Project2 {
        fasta_file: String,
        alphabet_file: String,
    },
    Project3 {
        alphabet_file: String,
        fasta_files: Vec<String>, // One or more FASTA files
    },
}