//! # Pattern preprocessing module.
//!
//! This module contains the `generate_masks` function, which generates bitmasks for each character
//! in the pattern's alphabet and returns them along with the accept state.
//!
//! # Example
//!
//! ```
//! use pattern_preprocessing::generate_masks;
//!
//! let pattern = b"abc";
//! let (masks, accept) = generate_masks(pattern)?;
//! assert_eq!(masks[97], 4); // 'a' bitmask
//! assert_eq!(accept, 4); // accept state
//! ```

use crate::pattern_matching::PatternError;
use std::mem;

/// Generate the bitmasks for each character in the pattern's alphabet and
/// Return the bitmasks for each character in the pattern's alphabet and
/// the accept state.
/// Pattern length must be less than or equal to the processor word size.
pub fn generate_masks(pattern: &[u8]) -> Result<([usize; 256], usize), PatternError> {
    let m = pattern.len();
    let mut masks: [usize; 256] = [0; 256];

    // Check if the pattern length is not greater than the register width.
    // The BNDM algorithm only works for patterns with length less than or
    // equal to the maximum value of processor word size.
    if mem::size_of::<usize>() == 8 && pattern.len() > 64 {
        return Err(PatternError::PatternTooLong(pattern.len(), 64));
    } else if mem::size_of::<usize>() == 4 && pattern.len() > 32 {
        return Err(PatternError::PatternTooLong(pattern.len(), 32));
    }

    for j in 0..m {
        masks[pattern[j] as usize] |= 1 << (m - j - 1);
    }
    let accept = 1 << (m - 1);

    Ok((masks, accept))
}

//
// ---------------------------------- Tests ----------------------------------
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_masks_accept() {
        let pattern = b"abc";
        let (_, accept) = generate_masks(pattern).unwrap();
        assert_eq!(accept, 4); // accept state
    }

    #[test]
    fn test_generate_masks_bitmasks() {
        let pattern = b"abc";
        let (masks, _) = generate_masks(pattern).unwrap();
        assert_eq!(masks[97], 4); // 'a' bitmask
        assert_eq!(masks[98], 2); // 'b' bitmask
        assert_eq!(masks[99], 1); // 'c' bitmask
    }

    #[test]
    fn test_generate_masks_complex() {
        let pattern = b"3$$X3";
        let (masks, accept) = generate_masks(pattern).unwrap();
        assert_eq!(masks[51], 17); // '3' bitmask
        assert_eq!(masks[36], 12); // '$' bitmask
        assert_eq!(masks[88], 2); // 'X' bitmask
        assert_eq!(accept, 16); // accept state
    }

    #[test]
    fn test_masks_too_large() {
        let pattern = b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890!@#$%^&*()_+";
        let result = generate_masks(pattern);
        assert!(matches!(result, Err(PatternError::PatternTooLong(_, _))));
    }
}
