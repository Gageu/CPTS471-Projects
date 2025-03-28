//types.rs

use std::{collections::BTreeMap, default};

pub enum GeneOperations {
    
}

pub struct Alignment{
    score: i32,
    sequence1: String,
    sequence2: String,
    stats: AllignmentStats,
}


impl Alignment {
    pub fn new(score: i32, sequence1: String, sequence2: String, stats: AllignmentStats) -> Self {
        Self { score, sequence1, sequence2, stats }
    }
    
    pub fn score(&self) -> i32 {
        self.score
    }

    pub fn sequence1(&self) -> &str {
        &self.sequence1
    }

    pub fn sequence2(&self) -> &str {
        &self.sequence2
    }

    pub fn stats(&self) -> &AllignmentStats {
        &self.stats
    }

}

pub struct AllignmentStats{
    alignment_length: i32,
    match_count: i32,
    mismatch_count: i32,
    gap_open_count: i32,
    gap_extend_count: i32,
    total_gaps: i32,
    identity_percent: f32,
}

impl AllignmentStats {
    pub fn new(
        alignment_length: i32, 
        match_count: i32,
        mismatch_count: i32, 
        gap_open_count: i32, 
        gap_extend_count: i32, 
        total_gaps: i32, 
        identity_percent: f32
    ) -> Self {
        Self { 
            alignment_length, 
            match_count, 
            mismatch_count, 
            gap_open_count, 
            gap_extend_count, 
            total_gaps, 
            identity_percent 
        }
    }
    
    pub fn alignment_length(&self) -> i32 {
        self.alignment_length
    }
    
    pub fn match_count(&self) -> i32 {
        self.match_count
    }
    
    pub fn mismatch_count(&self) -> i32 {
        self.mismatch_count
    }
    
    pub fn gap_open_count(&self) -> i32 {
        self.gap_open_count
    }
    
    pub fn gap_extend_count(&self) -> i32 {
        self.gap_extend_count
    }
    
    pub fn total_gaps(&self) -> i32 {
        self.total_gaps
    }
    
    pub fn identity_percent(&self) -> f32 {
        self.identity_percent
    }
}

#[derive(Copy, Clone)]
pub(crate) struct MDICell{
    pub m_score: i32,
    pub d_score: i32,
    pub i_score: i32,
}

impl MDICell{

    pub fn default() -> Self{
        MDICell{
            m_score: i32::MIN / 2,
            d_score: i32::MIN / 2,
            i_score: i32::MIN / 2,
        }
    }

    pub fn best_op(&self) -> char {
        if self.m_score >= self.d_score && self.m_score >= self.i_score {
            'm'
        } else if self.d_score >= self.i_score {
            'd'
        } else {
            'i'
        }
    }

    pub fn default_matrix(
        m: usize,
        n: usize,
        init_cell: i32,
        init_left: i32,
        init_top: i32,
    ) -> Vec<Vec<MDICell>> {
        let mut matrix = vec![vec![MDICell::default(); n + 1]; m + 1];
        matrix[0][0].m_score = init_cell;
        matrix[0][0].d_score = init_cell;
        matrix[0][0].i_score = init_cell;

        for i in 1..=m {
            matrix[i][0].d_score = init_left;
        }

        for j in 1..=n {
            matrix[0][j].i_score = init_top;
        }

        matrix
    }
}



#[derive(Default)]
pub struct ScoringSystem{
    match_score: Option<i32>,
    mismatch_score: Option<i32>,
    gap_score: Option<i32>,
    gap_open_score: Option<i32>,
    gap_extend_score: Option<i32>,
}

impl ScoringSystem {

    pub fn new(match_score: Option<i32>, mismatch_score: Option<i32>, gap_score: Option<i32>, gap_open_score: Option<i32>, gap_extend_score: Option<i32>) -> Self {
        Self { match_score, mismatch_score, gap_score, gap_open_score, gap_extend_score }
    }
    
    //--------- Getters ------------
    pub fn match_score(&self) -> Option<i32> {
        self.match_score
    }

    pub fn mismatch_score(&self) -> Option<i32> {
        self.mismatch_score
    }

    pub fn gap_score(&self) -> Option<i32> {
        self.gap_score
    }

    pub fn gap_open_score(&self) -> Option<i32> {
        self.gap_open_score
    }

    pub fn gap_extend_score(&self) -> Option<i32> {
        self.gap_extend_score
    }
    //--------------------------------


    pub fn unpack(&self) -> Result<(i32, i32, i32, i32), String> {
        Ok((
            self.match_score().ok_or("Missing match score")?,
            self.mismatch_score().ok_or("Missing mismatch score")?,
            self.gap_open_score().ok_or("Missing gap open score")?,
            self.gap_extend_score().ok_or("Missing gap extend score")?,
        ))
    }
}

pub enum ProjectSelection {
    Project1 {
        fasta_file: String,
        alg_select: String,
        config_file: String,
    },
    Project2 {
        fasta_file: String,
        alphabet_file: String,
    },
}
