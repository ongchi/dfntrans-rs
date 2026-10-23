//! The random number generator of DFNTrans: `srand48` and `drand48` of the
//! C library. The same seed gives the same numbers as the C code, thus a
//! run can be compared with a run of DFNTrans.

const MULTIPLIER: u64 = 0x5_DEEC_E66D;
const INCREMENT: u64 = 0xB;
const MASK: u64 = (1 << 48) - 1;

/// 48-bit linear congruential generator of `drand48`
#[derive(Clone, Debug)]
pub struct Drand48 {
    state: u64,
}

impl Drand48 {
    /// `srand48(seed)`: the 32 bits of the seed are the high bits of the state.
    pub fn new(seed: u32) -> Self {
        Self {
            state: (u64::from(seed) << 16) | 0x330E,
        }
    }

    /// `drand48()`: uniform in [0, 1)
    pub fn next_f64(&mut self) -> f64 {
        self.state = (MULTIPLIER.wrapping_mul(self.state).wrapping_add(INCREMENT)) & MASK;
        // The state has 48 bits, thus the division is exact
        self.state as f64 / (1u64 << 48) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Values of `drand48` of the C library (macOS) after `srand48(seed)`,
    /// printed with 17 digits
    #[test]
    #[allow(clippy::excessive_precision)]
    fn gives_the_numbers_of_the_c_library() {
        let expected = [
            (
                12345,
                [
                    0.225_328_512_796_298_95,
                    0.919_183_068_533_556_02,
                    0.206_841_253_248_182_26,
                    0.724_779_720_275_314_79,
                ],
            ),
            (
                0,
                [
                    0.170_828_036_106_289_72,
                    0.749_901_980_484_963_81,
                    0.096_371_655_623_567_42,
                    0.870_465_227_027_075_63,
                ],
            ),
        ];
        for (seed, values) in expected {
            let mut random = Drand48::new(seed);
            for value in values {
                assert_eq!(random.next_f64(), value);
            }
        }
    }
}
