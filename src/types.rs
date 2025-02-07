//types.rs

use std::default;

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
    allignment_length: i32,
    match_count: i32,
    mismatch_count: i32,
    gap_open_count: i32,
    gap_extend_count: i32,
    total_gaps: i32,
    identity_percent: f32,
}

impl AllignmentStats {
    pub fn new(
        allignment_length: i32, 
        match_count: i32,
        mismatch_count: i32, 
        gap_open_count: i32, 
        gap_extend_count: i32, 
        total_gaps: i32, 
        identity_percent: f32
    ) -> Self {
        Self { 
            allignment_length, 
            match_count, 
            mismatch_count, 
            gap_open_count, 
            gap_extend_count, 
            total_gaps, 
            identity_percent 
        }
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
}
