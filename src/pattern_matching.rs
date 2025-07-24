//! # Pattern matching algorithms BNDM and BNDMq.
//!
//! The BNDM algorithm is a backwards string matching algorithm that uses bitmasks to quickly
//! determine whether a character in the text is part of the pattern. The BNDMq algorithm is a
//! variant of BNDM that uses _q_-grams to improve performance.
//!
//! Both algorithms are based on the pseudocode found in Ďurian et al. (2009),
//! [doi:10.1137/1.9781611972894.3](https://doi.org/10.1137/1.9781611972894.3).
//!     
//! # Example
//!
//! ```
//! use pattern_matching::BNDMq;
//!
//! let pattern = b"abc";
//! let text = b"abcabcabc";
//!     
//! let bndmq = BNDMq::new(pattern, 2).unwrap();
//!
//! let matches: Vec<usize> = bndmq.find_all(text).collect();
//! assert_eq!(matches, vec![0, 3, 6]);
//! ```

use crate::pattern_preprocessing::generate_masks;
use anyhow::Result;

/// Error type for pattern matching operations
#[derive(Debug, thiserror::Error)]
pub enum PatternError {
    #[error("Invalid q-gram length: {0}. Must be between 1 and pattern length.")]
    InvalidQGramLength(usize),
    #[error("Pattern is empty.")]
    EmptyPattern,
    #[error("Pattern length {0} is too large for this architecture when using BNDM (max {1}).")]
    PatternTooLong(usize, usize),
}

/// Backwards Non-deterministic Dawg String Matching algorithm, tuned to use _q_-grams.
///
/// Implementation based on the pseudocode found in Ďurian et al. (2009), [doi:10.1137/1.9781611972894.3](https://doi.org/10.1137/1.9781611972894.3).
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BNDMq {
    m: usize,
    q: usize,
    masks: [usize; 256],
    accept: usize,
}

impl BNDMq {
    /// Creates a new BNDMq instance with the given pattern and q-gram length.
    ///
    /// # Arguments
    ///
    /// * `pattern` - The pattern to search for
    /// * `q` - The length of q-grams to use (must be between 1 and pattern length)
    ///
    /// # Errors
    ///
    /// Returns `PatternError::InvalidQGramLength` if q is 0 or greater than pattern length
    /// Returns `PatternError::EmptyPattern` if pattern is empty
    pub fn new(pattern: &[u8], q: usize) -> Result<Self, PatternError> {
        if pattern.is_empty() {
            return Err(PatternError::EmptyPattern);
        }
        if q == 0 || q > pattern.len() {
            return Err(PatternError::InvalidQGramLength(q));
        }

        let m = pattern.len();
        let (masks, accept) = generate_masks(pattern)?;

        Ok(BNDMq {
            m,
            q,
            masks,
            accept,
        })
    }

    /// Core BNDMq algorithm implementation that yields matches through a callback.
    /// Returns true if a match was found and the callback returned true.
    fn find_matches<F>(&self, text: &[u8], mut on_match: F) -> bool
    where
        F: FnMut(usize) -> bool,
    {
        // Return false if the pattern is longer than the text
        if self.m > text.len() {
            return false;
        }

        let mut i = self.m - self.q + 1;
        let n = text.len();
        let m_minus_q_plus_1 = self.m - self.q + 1;

        while i <= n - self.q + 1 {
            let mut state_mask = self.masks[text[i - 1] as usize];
            for ii in 0..self.q - 1 {
                state_mask &= self.masks[text[i + ii] as usize] << (ii + 1);
            }

            if state_mask != 0 {
                let mut j = i;
                let first = i - m_minus_q_plus_1;

                loop {
                    j -= 1;

                    if state_mask >= self.accept {
                        if j > first {
                            i = j;
                        } else if on_match(j) {
                            return true;
                        }
                    }
                    state_mask = (state_mask << 1) & self.masks[text[j - 1] as usize];

                    if state_mask == 0 {
                        break;
                    }
                }
            }
            i += m_minus_q_plus_1;
        }
        false
    }

    /// Search for the pattern and return true as soon as a match is found.
    pub fn find_match(&self, text: &[u8]) -> bool {
        self.find_matches(text, |_| true)
    }

    /// Returns an iterator over all matches of the pattern in the given text.
    pub fn find_iter<'a>(&'a self, text: &'a [u8]) -> Matches<'a> {
        Matches {
            bndmq: self,
            text,
            i: self.m - self.q + 1,
            n: text.len(),
        }
    }

    /// Returns all matches of the pattern in the given text.
    ///
    /// # Arguments
    ///
    /// * `text` - The text to search in
    ///
    /// # Returns
    ///
    /// A vector containing the starting position of each match
    pub fn find_all(&self, text: &[u8]) -> Vec<usize> {
        self.find_iter(text).collect()
    }
}

/// Iterator over the matches found by the BNDMq algorithm.
#[derive(Debug)]
pub struct Matches<'a> {
    bndmq: &'a BNDMq,
    text: &'a [u8],
    i: usize,
    n: usize,
}

impl<'a> Iterator for Matches<'a> {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // Return None if the pattern is longer than the text
        if self.bndmq.m > self.n {
            return None;
        }
        let m_minus_q_plus_1 = self.bndmq.m - self.bndmq.q + 1;

        while self.i <= self.n - self.bndmq.q + 1 {
            let mut state_mask = self.bndmq.masks[self.text[self.i - 1] as usize];
            for ii in 0..self.bndmq.q - 1 {
                state_mask &= self.bndmq.masks[self.text[self.i + ii] as usize] << (ii + 1);
            }

            if state_mask != 0 {
                let mut j = self.i;
                let first = self.i - m_minus_q_plus_1;

                loop {
                    j -= 1;

                    if state_mask >= self.bndmq.accept {
                        if j > first {
                            self.i = j;
                        } else {
                            // Advance the iterator to the next position and return the match
                            self.i += m_minus_q_plus_1;
                            return Some(j);
                        }
                    }
                    state_mask = (state_mask << 1) & self.bndmq.masks[self.text[j - 1] as usize];

                    if state_mask == 0 {
                        break;
                    }
                }
            }
            self.i += m_minus_q_plus_1;
        }
        None
    }
}

/// Tune the size of the _q_-grams for BNDMq based on the pattern length.
/// Based on the paper by Ďurian et al. (2009) and personal experience.
pub fn tune_q_value(pattern: &str) -> Result<usize> {
    let pattern_len = pattern.len();
    let q = match pattern_len {
        0..=1 => 1,
        2..=3 => 2,
        4..=8 => 3,
        9..=30 => 4,
        31..=55 => 5,
        56..=64 => 6,
        65.. => anyhow::bail!("Pattern length is too long for BNDMq."),
    };
    Ok(q)
}

/// Backwards Non-deterministic Dawg String Matching algorithm.
/// Legacy implementation that does not return matches as an iterator.
///
/// Implementation based on the pseudocode found in Ďurian et al. (2009), [doi:10.1137/1.9781611972894.3](https://doi.org/10.1137/1.9781611972894.3).
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BNDM {
    m: usize,
    masks: [usize; 256],
    accept: usize,
}

impl BNDM {
    /// Create a new BNDM instance with a given pattern.
    ///
    /// # Arguments
    ///
    /// * `pattern` - The pattern to search for
    ///
    /// # Errors
    ///
    /// Returns `PatternError::EmptyPattern` if pattern is empty
    pub fn new(pattern: &[u8]) -> Result<Self, PatternError> {
        if pattern.is_empty() {
            return Err(PatternError::EmptyPattern);
        }

        let m = pattern.len();
        let (masks, accept) = generate_masks(pattern)?;

        Ok(BNDM { m, masks, accept })
    }

    /// Find all occurrences of the pattern in the given text and return their indices.
    ///
    /// # Arguments
    ///
    /// * `text` - The text to search in
    ///
    /// # Returns
    ///
    /// A vector containing the starting position of each match
    pub fn find_all(&self, text: &[u8]) -> Vec<usize> {
        let n = text.len();
        let mut occ = Vec::new();

        // Return empty vector if the pattern is longer than the text
        if self.m > n {
            return occ;
        }

        let mut i = 0;
        while i <= n - self.m {
            let mut j = self.m;
            let mut last = self.m;
            let mut state_mask = (1 << self.m) - 1;

            while state_mask != 0 {
                state_mask &= self.masks[text[i + j - 1] as usize];
                j -= 1;

                if state_mask & self.accept != 0 {
                    if j > 0 {
                        last = j;
                    } else {
                        occ.push(i);
                        break;
                    }
                }
                state_mask <<= 1;
            }
            i += last;
        }
        occ
    }

    /// Search for the pattern and return true as soon as a match is found.
    ///
    /// # Arguments
    ///
    /// * `text` - The text to search in
    ///
    /// # Returns
    ///
    /// `true` if a match is found, `false` otherwise
    pub fn find_match(&self, text: &[u8]) -> bool {
        let n = text.len();

        // Return false if the pattern is longer than the text
        if self.m > n {
            return false;
        }

        let mut i = 0;
        while i <= n - self.m {
            let mut j = self.m;
            let mut last = self.m;
            let mut state_mask = (1 << self.m) - 1;

            while state_mask != 0 {
                state_mask &= self.masks[text[i + j - 1] as usize];
                j -= 1;

                if state_mask & self.accept != 0 {
                    if j > 0 {
                        last = j;
                    } else {
                        return true;
                    }
                }
                state_mask <<= 1;
            }
            i += last;
        }
        false
    }
}

//
// ---------------------------------- Tests ----------------------------------
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bndmq_iter() {
        let pattern = b"abc";
        let text = b"abcabcabc";
        let bndmq_iter = BNDMq::new(pattern, 2).unwrap();

        let mut occurence_iter = bndmq_iter.find_iter(text);
        assert_eq!(occurence_iter.next(), Some(0));
        assert_eq!(occurence_iter.next(), Some(3));
        assert_eq!(occurence_iter.next(), Some(6));
        assert_eq!(occurence_iter.next(), None);
    }

    #[test]
    fn test_bndmq_longer_pattern() {
        let pattern = b"1234567890";
        let text = b"123";
        let bndmq = BNDMq::new(pattern, 2).unwrap();

        let occurrences: Vec<usize> = bndmq.find_iter(text).collect();
        assert!(occurrences.is_empty());
    }

    #[test]
    fn test_bndmq_empty_text() {
        let pattern = b"abc";
        let text = b"";
        let bndmq = BNDMq::new(pattern, 2).unwrap();

        let occurrences: Vec<usize> = bndmq.find_iter(text).collect();
        assert!(occurrences.is_empty());
    }

    #[test]
    fn test_bndmq_find_all() {
        let pattern = b"abc";
        let text = b"aabcabcabc";
        let bndmq = BNDMq::new(pattern, 2).unwrap();
        let occurrences: Vec<usize> = bndmq.find_iter(text).collect();
        assert_eq!(occurrences, vec![1, 4, 7]);
    }

    #[test]
    fn test_bndm_find_all() {
        let pattern = b"abc";
        let text = b"abcabcabc";
        let bndm = BNDM::new(pattern).unwrap();
        let occurrences = bndm.find_all(text);
        assert_eq!(occurrences, vec![0, 3, 6]);
    }

    #[test]
    fn test_bndm_find_match() {
        let pattern = b"abc";
        let text = b"abcabcabc";
        let bndm = BNDM::new(pattern).unwrap();
        let has_match = bndm.find_match(text);
        assert!(has_match);
    }

    #[test]
    fn test_bndm_longer_pattern() {
        let pattern = b"1234567890";
        let text = b"123";
        let bndm = BNDM::new(pattern).unwrap();
        let occurrences = bndm.find_all(text);
        assert!(occurrences.is_empty());
    }

    #[test]
    fn test_bndm_empty_text() {
        let pattern = b"abc";
        let text = b"";
        let bndm = BNDM::new(pattern).unwrap();
        let occurrences = bndm.find_all(text);
        assert!(occurrences.is_empty());
    }

    #[test]
    fn test_bndmq_new_q_too_large() {
        let pattern = b"abc";
        assert!(matches!(
            BNDMq::new(pattern, 4),
            Err(PatternError::InvalidQGramLength(4))
        ));
    }

    #[test]
    fn test_bndmq_new_q_too_small() {
        let pattern = b"abc";
        assert!(matches!(
            BNDMq::new(pattern, 0),
            Err(PatternError::InvalidQGramLength(0))
        ));
    }

    #[test]
    fn test_bndm_empty_pattern() {
        let pattern = b"";
        assert!(matches!(
            BNDM::new(pattern),
            Err(PatternError::EmptyPattern)
        ));
    }

    #[test]
    fn test_bndmq_empty_pattern() {
        let pattern = b"";
        assert!(matches!(
            BNDMq::new(pattern, 1),
            Err(PatternError::EmptyPattern)
        ));
    }

    #[test]
    fn test_bndmq_find_match() {
        let pattern = b"abc";
        let text = b"abcabcabc";
        let bndmq = BNDMq::new(pattern, 2).unwrap();
        let has_match = bndmq.find_match(text);
        assert!(has_match);
    }

    #[test]
    fn test_bndmq_find_match_no_match() {
        let pattern = b"abc";
        let text = b"defdefdef";
        let bndmq = BNDMq::new(pattern, 2).unwrap();
        let has_match = bndmq.find_match(text);
        assert!(!has_match);
    }

    #[test]
    fn test_tune_q_value() {
        let q = tune_q_value("AAAAAAAACCCCCCCCGGGGGGGGTTTTTTT").unwrap();
        assert_eq!(q, 5);
    }
}
