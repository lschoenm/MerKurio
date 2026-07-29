//! # Exact pattern matching
//!
//! This module provides the matcher implementations used by the FASTX and
//! alignment commands:
//!
//! - [`PatternMatcher`] is the shared facade over Aho–Corasick, BNDMq, and
//!   rolling-hash matching.
//! - [`BNDMq`] searches one pattern using bit-parallel backwards matching and
//!   _q_-grams. Patterns are limited to the processor word size.
//! - [`HashMatcher`] searches many equal-length patterns in one pass using a
//!   rolling hash. It supports arbitrary bytes and verifies every hash hit, so
//!   collisions cannot produce false matches.
//! - [`BNDM`] is the legacy, non-_q_-gram implementation.
//!
//! All match positions are zero-based byte offsets. Matching is exact and
//! case-sensitive except when the Aho–Corasick variant is explicitly built
//! with ASCII case-insensitive matching. Reverse-complement and canonical
//! pattern generation happen before patterns reach this module.
//!
//! BNDM and BNDMq follow Ďurian et al. (2009),
//! [doi:10.1137/1.9781611972894.3](https://doi.org/10.1137/1.9781611972894.3).
//!
//! # Examples
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
//!
//! Rolling-hash matching accepts multiple patterns of one common length:
//!
//! ```
//! use pattern_matching::{HashMatcher, MatchHit};
//!
//! let patterns = vec!["abc".to_string(), "bca".to_string()];
//! let matcher = HashMatcher::new(&patterns).unwrap();
//! let mut matches = Vec::new();
//! matcher.for_each_match(b"abcabc", |hit| matches.push(hit));
//!
//! assert_eq!(
//!     matches,
//!     vec![
//!         MatchHit { pattern_index: 0, position: 0 },
//!         MatchHit { pattern_index: 1, position: 1 },
//!         MatchHit { pattern_index: 0, position: 3 },
//!     ]
//! );
//! ```

use crate::pattern_preprocessing::generate_masks;
use aho_corasick::AhoCorasick;
use anyhow::Result;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// Error type for pattern matching operations
#[derive(Debug, thiserror::Error)]
pub enum PatternError {
    #[error("Invalid q-gram length: {0}. Must be between 1 and pattern length.")]
    InvalidQGramLength(usize),
    #[error("Pattern is empty.")]
    EmptyPattern,
    #[error("Pattern length {0} is too large for this architecture when using BNDM (max {1}).")]
    PatternTooLong(usize, usize),
    #[error("Hash-based matching requires all patterns to have the same length.")]
    UnequalPatternLengths,
}

/// A single pattern match within one searched record.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct MatchHit {
    pub pattern_index: usize,
    pub position: usize,
}

/// Parser-independent pattern matcher used by FASTX and alignment commands.
#[derive(Debug)]
pub enum PatternMatcher {
    AhoCorasick { ac: AhoCorasick },
    Bndmq { matchers: Vec<BNDMq> },
    Hash { matcher: HashMatcher },
}

impl PatternMatcher {
    /// Build the requested matcher implementation from the preprocessed pattern list.
    pub fn new(
        pattern_list: &[String],
        use_aho_corasick: bool,
        use_hash: bool,
        case_insensitive: bool,
        q_size: Option<usize>,
    ) -> Result<Self> {
        if use_aho_corasick {
            let ac = AhoCorasick::builder()
                // Use DFA for better search performance at higher memory cost.
                .kind(Some(aho_corasick::AhoCorasickKind::DFA))
                .ascii_case_insensitive(case_insensitive)
                .build(pattern_list)?;
            Ok(Self::AhoCorasick { ac })
        } else if use_hash {
            Ok(Self::Hash {
                matcher: HashMatcher::new(pattern_list)?,
            })
        } else {
            let mut matchers = Vec::with_capacity(pattern_list.len());
            for pattern in pattern_list {
                let q = q_size.unwrap_or_else(|| tune_q_value(pattern).unwrap());
                matchers.push(BNDMq::new(pattern.as_bytes(), q)?);
            }
            Ok(Self::Bndmq { matchers })
        }
    }

    /// Search for any match and return as soon as one is found.
    pub fn find_any(&self, seq: &[u8]) -> bool {
        match self {
            Self::AhoCorasick { ac } => ac.find(seq).is_some(),
            Self::Bndmq { matchers } => matchers.iter().any(|bndmq| bndmq.find_match(seq)),
            Self::Hash { matcher } => matcher.find_match(seq),
        }
    }

    /// Stream all matches in the implementation's native order.
    pub fn for_each_match<F>(&self, seq: &[u8], mut on_match: F)
    where
        F: FnMut(MatchHit),
    {
        match self {
            Self::AhoCorasick { ac } => {
                for mat in ac.find_overlapping_iter(seq) {
                    on_match(MatchHit {
                        pattern_index: mat.pattern().as_usize(),
                        position: mat.start(),
                    });
                }
            }
            Self::Bndmq { matchers } => {
                for (pattern_index, bndmq) in matchers.iter().enumerate() {
                    for position in bndmq.find_iter(seq) {
                        on_match(MatchHit {
                            pattern_index,
                            position,
                        });
                    }
                }
            }
            Self::Hash { matcher } => matcher.for_each_match(seq, on_match),
        }
    }

    /// Mark every pattern that occurs at least once in `seq`.
    ///
    /// `matched_patterns` must have one entry per input pattern. Existing
    /// values are preserved, allowing callers to combine matches from
    /// multiple sequences without an intermediate allocation.
    pub fn mark_matching_patterns(&self, seq: &[u8], matched_patterns: &mut [bool]) {
        match self {
            Self::AhoCorasick { ac } => {
                for mat in ac.find_overlapping_iter(seq) {
                    matched_patterns[mat.pattern().as_usize()] = true;
                }
            }
            Self::Bndmq { matchers } => {
                for (pattern_index, bndmq) in matchers.iter().enumerate() {
                    if bndmq.find_match(seq) {
                        matched_patterns[pattern_index] = true;
                    }
                }
            }
            Self::Hash { matcher } => {
                matcher.for_each_match(seq, |hit| {
                    matched_patterns[hit.pattern_index] = true;
                });
            }
        }
    }

    pub fn algorithm_name(&self) -> &'static str {
        match self {
            Self::AhoCorasick { .. } => "Aho-Corasick",
            Self::Bndmq { .. } => "BNDMq",
            Self::Hash { .. } => "rolling hash",
        }
    }
}

const ROLLING_HASH_BASE: u64 = 0x9e37_79b1_85eb_ca87;

#[derive(Debug, Default)]
struct IdentityHasher(u64);

impl Hasher for IdentityHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        let mut value = 0u64;
        for &byte in bytes {
            value = value.wrapping_mul(256).wrapping_add(byte as u64);
        }
        self.0 = value;
    }

    fn write_u64(&mut self, value: u64) {
        self.0 = value;
    }
}

type IdentityBuildHasher = BuildHasherDefault<IdentityHasher>;

#[derive(Debug)]
enum HashBucket {
    One(usize),
    Many(Vec<usize>),
}

impl HashBucket {
    fn insert(&mut self, pattern_index: usize) {
        match self {
            Self::One(first) => {
                *self = Self::Many(vec![*first, pattern_index]);
            }
            Self::Many(indices) => indices.push(pattern_index),
        }
    }

    fn for_each(&self, mut visit: impl FnMut(usize) -> bool) -> bool {
        match self {
            Self::One(pattern_index) => visit(*pattern_index),
            Self::Many(pattern_indices) => {
                for &pattern_index in pattern_indices {
                    if visit(pattern_index) {
                        return true;
                    }
                }
                false
            }
        }
    }
}

/// Exact multi-pattern matcher for patterns with one common length.
///
/// A wrapping polynomial hash is updated in constant time for each text
/// window. Hash hits are always verified against the original pattern bytes,
/// so hash collisions cannot produce false matches. The algorithm treats
/// patterns and text as arbitrary bytes and is case-sensitive.
#[derive(Debug)]
pub struct HashMatcher {
    pattern_len: usize,
    highest_power: u64,
    patterns: Box<[u8]>,
    buckets: HashMap<u64, HashBucket, IdentityBuildHasher>,
}

impl HashMatcher {
    pub fn new(pattern_list: &[String]) -> Result<Self, PatternError> {
        let pattern_len = pattern_list
            .first()
            .ok_or(PatternError::EmptyPattern)?
            .len();
        if pattern_len == 0 {
            return Err(PatternError::EmptyPattern);
        }
        if pattern_list
            .iter()
            .any(|pattern| pattern.len() != pattern_len)
        {
            return Err(PatternError::UnequalPatternLengths);
        }

        let highest_power = wrapping_pow(ROLLING_HASH_BASE, pattern_len - 1);
        let mut patterns = Vec::with_capacity(pattern_list.len() * pattern_len);
        let mut buckets =
            HashMap::with_capacity_and_hasher(pattern_list.len(), IdentityBuildHasher::default());

        for (pattern_index, pattern) in pattern_list.iter().enumerate() {
            let pattern_bytes = pattern.as_bytes();
            let hash = hash_window(pattern_bytes);
            patterns.extend_from_slice(pattern_bytes);
            buckets
                .entry(hash)
                .and_modify(|bucket: &mut HashBucket| bucket.insert(pattern_index))
                .or_insert(HashBucket::One(pattern_index));
        }

        Ok(Self {
            pattern_len,
            highest_power,
            patterns: patterns.into_boxed_slice(),
            buckets,
        })
    }

    pub fn find_match(&self, text: &[u8]) -> bool {
        self.find_matches(text, |_| true)
    }

    pub fn for_each_match(&self, text: &[u8], mut on_match: impl FnMut(MatchHit)) {
        self.find_matches(text, |hit| {
            on_match(hit);
            false
        });
    }

    fn find_matches(&self, text: &[u8], mut on_match: impl FnMut(MatchHit) -> bool) -> bool {
        if text.len() < self.pattern_len {
            return false;
        }

        let last_position = text.len() - self.pattern_len;
        let mut window_hash = hash_window(&text[..self.pattern_len]);
        for position in 0..=last_position {
            if let Some(bucket) = self.buckets.get(&window_hash)
                && bucket.for_each(|pattern_index| {
                    let pattern_start = pattern_index * self.pattern_len;
                    let pattern_end = pattern_start + self.pattern_len;
                    if text[position..position + self.pattern_len]
                        == self.patterns[pattern_start..pattern_end]
                    {
                        on_match(MatchHit {
                            pattern_index,
                            position,
                        })
                    } else {
                        false
                    }
                })
            {
                return true;
            }

            if position < last_position {
                let outgoing = text[position] as u64 + 1;
                let incoming = text[position + self.pattern_len] as u64 + 1;
                window_hash = window_hash
                    .wrapping_sub(outgoing.wrapping_mul(self.highest_power))
                    .wrapping_mul(ROLLING_HASH_BASE)
                    .wrapping_add(incoming);
            }
        }
        false
    }
}

#[inline]
fn hash_window(window: &[u8]) -> u64 {
    window.iter().fold(0u64, |hash, &byte| {
        hash.wrapping_mul(ROLLING_HASH_BASE)
            .wrapping_add(byte as u64 + 1)
    })
}

#[inline]
fn wrapping_pow(mut base: u64, mut exponent: usize) -> u64 {
    let mut result = 1u64;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = result.wrapping_mul(base);
        }
        base = base.wrapping_mul(base);
        exponent >>= 1;
    }
    result
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
///
/// The ranges are a smoothed fit to the internal DNA q-value tuning benchmark
/// in `tuning/q-value`, informed by the BNDMq paper by Ďurian et al. (2009).
pub fn tune_q_value(pattern: &str) -> Result<usize> {
    let pattern_len = pattern.len();
    let q = match pattern_len {
        0..=1 => 1,
        2..=3 => 2,
        4..=5 => 3,
        6..=12 => 4,
        13..=25 => 5,
        26..=64 => 6,
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
        let cases = [
            (0, 1),
            (1, 1),
            (2, 2),
            (3, 2),
            (4, 3),
            (5, 3),
            (6, 4),
            (12, 4),
            (13, 5),
            (25, 5),
            (26, 6),
            (64, 6),
        ];

        for (pattern_len, expected_q) in cases {
            let pattern = "A".repeat(pattern_len);
            assert_eq!(
                tune_q_value(&pattern).unwrap(),
                expected_q,
                "unexpected q for pattern length {pattern_len}"
            );
        }

        assert!(tune_q_value(&"A".repeat(65)).is_err());
    }

    #[test]
    fn test_pattern_matcher_bndmq() {
        let patterns = vec!["abc".to_string()];
        let matcher = PatternMatcher::new(&patterns, false, false, false, Some(2)).unwrap();

        assert!(matcher.find_any(b"abcabc"));
        assert!(!matcher.find_any(b"defdef"));
        assert_eq!(matcher.algorithm_name(), "BNDMq");

        let mut hits = Vec::new();
        matcher.for_each_match(b"abcabc", |hit| hits.push(hit));

        assert_eq!(
            hits,
            vec![
                MatchHit {
                    pattern_index: 0,
                    position: 0,
                },
                MatchHit {
                    pattern_index: 0,
                    position: 3,
                },
            ]
        );
    }

    #[test]
    fn test_pattern_matcher_aho_corasick() {
        let patterns = vec!["abc".to_string()];
        let matcher = PatternMatcher::new(&patterns, true, false, false, None).unwrap();

        assert!(matcher.find_any(b"abcabc"));
        assert!(!matcher.find_any(b"defdef"));
        assert_eq!(matcher.algorithm_name(), "Aho-Corasick");

        let mut hits = Vec::new();
        matcher.for_each_match(b"abcabc", |hit| hits.push(hit));

        assert_eq!(
            hits,
            vec![
                MatchHit {
                    pattern_index: 0,
                    position: 0,
                },
                MatchHit {
                    pattern_index: 0,
                    position: 3,
                },
            ]
        );
    }

    #[test]
    fn test_hash_matcher_finds_overlapping_arbitrary_byte_patterns() {
        let patterns = vec!["a$a".to_string(), "$a$".to_string()];
        let matcher = HashMatcher::new(&patterns).unwrap();
        let mut hits = Vec::new();

        matcher.for_each_match(b"a$a$a", |hit| hits.push(hit));

        assert_eq!(
            hits,
            vec![
                MatchHit {
                    pattern_index: 0,
                    position: 0,
                },
                MatchHit {
                    pattern_index: 1,
                    position: 1,
                },
                MatchHit {
                    pattern_index: 0,
                    position: 2,
                },
            ]
        );
        assert!(matcher.find_match(b"xx$a$xx"));
        assert!(!matcher.find_match(b"xxAAAxx"));
    }

    #[test]
    fn test_hash_matcher_preserves_case_and_duplicate_pattern_indices() {
        let patterns = vec!["Aa".to_string(), "Aa".to_string(), "aa".to_string()];
        let matcher = HashMatcher::new(&patterns).unwrap();
        let mut hits = Vec::new();

        matcher.for_each_match(b"Aaaa", |hit| hits.push(hit));

        assert_eq!(
            hits,
            vec![
                MatchHit {
                    pattern_index: 0,
                    position: 0,
                },
                MatchHit {
                    pattern_index: 1,
                    position: 0,
                },
                MatchHit {
                    pattern_index: 2,
                    position: 1,
                },
                MatchHit {
                    pattern_index: 2,
                    position: 2,
                },
            ]
        );
    }

    #[test]
    fn test_hash_matcher_rejects_empty_and_mixed_length_patterns() {
        assert!(matches!(
            HashMatcher::new(&[]),
            Err(PatternError::EmptyPattern)
        ));
        assert!(matches!(
            HashMatcher::new(&["".to_string()]),
            Err(PatternError::EmptyPattern)
        ));
        assert!(matches!(
            HashMatcher::new(&["AA".to_string(), "AAA".to_string()]),
            Err(PatternError::UnequalPatternLengths)
        ));
    }

    #[test]
    fn test_hash_matcher_matches_naive_search_across_pattern_lengths() {
        const ALPHABET: &[u8] = b"ACGTNacgt$#";
        let mut state = 0x1234_5678_9abc_def0u64;
        let mut random_byte = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            ALPHABET[(state >> 32) as usize % ALPHABET.len()]
        };

        for pattern_len in [1, 2, 3, 8, 31, 65, 129] {
            let patterns: Vec<String> = (0..12)
                .map(|_| {
                    String::from_utf8((0..pattern_len).map(|_| random_byte()).collect()).unwrap()
                })
                .collect();
            let mut text: Vec<u8> = (0..512).map(|_| random_byte()).collect();
            for (pattern_index, position) in [0, 7, 103, 251].into_iter().enumerate() {
                text[position..position + pattern_len]
                    .copy_from_slice(patterns[pattern_index].as_bytes());
            }

            let matcher = HashMatcher::new(&patterns).unwrap();
            let mut actual = Vec::new();
            matcher.for_each_match(&text, |hit| actual.push(hit));

            let mut expected = Vec::new();
            for position in 0..=text.len() - pattern_len {
                for (pattern_index, pattern) in patterns.iter().enumerate() {
                    if text[position..position + pattern_len] == *pattern.as_bytes() {
                        expected.push(MatchHit {
                            pattern_index,
                            position,
                        });
                    }
                }
            }

            assert_eq!(actual, expected, "pattern length {pattern_len}");
        }
    }

    #[test]
    fn test_pattern_matcher_hash_marks_distinct_patterns() {
        let patterns = vec!["aba".to_string(), "bab".to_string()];
        let matcher = PatternMatcher::new(&patterns, false, true, false, None).unwrap();
        let mut matched = vec![false; patterns.len()];

        matcher.mark_matching_patterns(b"ababa", &mut matched);

        assert_eq!(matched, vec![true, true]);
        assert_eq!(matcher.algorithm_name(), "rolling hash");
    }
}
