//! One reviewed deterministic PRNG primitive for all of atman.
//!
//! Before this module, ~six hand-rolled splitmix/LCG generators lived
//! across the CLI and core crates. They had *already diverged* — some
//! special-cased `seed == 0`, some XOR-folded the golden-ratio constant,
//! ratio.rs mixed an LCG step into the splitmix finaliser — and every one
//! of them drew bounded integers with a modulo reduction
//! (`next_u64() as usize % upper`), which is biased: when `upper` does not
//! divide `2^64`, the low residue values are slightly over-represented.
//!
//! This module replaces all of them with two reviewed generators and one
//! unbiased bounded-integer draw:
//!
//! - [`SplitMix64`]: a fast 64-bit splitmix stream. Used by the CLI
//!   resampling/permutation commands (`null`, `bootstrap`, `ratio`).
//! - [`Xoshiro256pp`]: xoshiro256++ seeded from a splitmix init. Used by
//!   the core decomposition family (ICA, NMF, VCA/FCLS unmixing, GSEA,
//!   archetype null, alignment bootstrap) and by Dunnett-Hsu Monte Carlo.
//!
//! Both expose [`next_u64`](SplitMix64::next_u64) and an unbiased
//! [`bounded`](SplitMix64::bounded) draw (Lemire's nearly-divisionless
//! method with rejection). `bounded` is the *only* sanctioned way to draw
//! a uniform integer in `0..upper`; there are no remaining
//! `next_u64() % upper` call sites.
//!
//! ## Seed rule (one rule, applied centrally)
//!
//! A raw seed of `0` is a degenerate state for both generators (splitmix
//! is fine at 0, but xoshiro requires a non-zero state vector, and a
//! shared all-zero convention keeps the two interchangeable). We therefore
//! fold the all-zero seed to the golden-ratio constant
//! `0x9E3779B97F4A7C15` exactly once, in [`fold_zero_seed`], and every
//! constructor routes through it. This subsumes the divergent ad-hoc
//! `if seed == 0 { .. }` checks that some of the old generators had and
//! others lacked.
//!
//! ## Determinism contract
//!
//! Byte-identical output for a fixed seed is atman's core paper claim.
//! Both generators are pure integer arithmetic with no platform-dependent
//! floating point in the integer stream, so two runs at the same seed
//! produce identical `u64` sequences and therefore identical bounded
//! draws and identical downstream output.

/// Golden-ratio odd constant; the canonical splitmix64 increment and the
/// value we fold the all-zero seed to.
pub const GOLDEN_GAMMA: u64 = 0x9E3779B97F4A7C15;

/// Apply the single, central seed rule: fold the degenerate all-zero seed
/// to [`GOLDEN_GAMMA`]; pass every other seed through unchanged.
#[inline]
pub fn fold_zero_seed(seed: u64) -> u64 {
    if seed == 0 {
        GOLDEN_GAMMA
    } else {
        seed
    }
}

/// One step of the splitmix64 finaliser on a mutable state word. Shared by
/// both [`SplitMix64`] and the seeding of [`Xoshiro256pp`].
#[inline]
fn splitmix64_step(state: &mut u64) -> u64 {
    *state = state.wrapping_add(GOLDEN_GAMMA);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// Derive a deterministic per-iteration sub-seed from `(seed, iter)`.
///
/// Used by the bootstrap / null families so the same top-level `--seed`
/// drives the same per-iteration stream everywhere. This is the historical
/// `derive_sub_seed` / `bootstrap_sub_seed` convention, centralized:
/// `splitmix64` finaliser over `seed + GOLDEN_GAMMA * (iter + 1)`.
#[inline]
pub fn derive_sub_seed(seed: u64, iter: usize) -> u64 {
    let mut z = seed.wrapping_add(GOLDEN_GAMMA.wrapping_mul(iter as u64 + 1));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// Lemire's unbiased bounded draw: a uniform integer in `0..upper` with no
/// modulo bias, given a closure that produces uniform `u64`s.
///
/// Panics if `upper == 0`. Uses the nearly-divisionless algorithm
/// (Lemire 2019) with a rejection step that is taken with probability
/// `< upper / 2^64`, so in practice it almost never rejects.
#[inline]
fn lemire_bounded(upper: u64, mut next: impl FnMut() -> u64) -> u64 {
    debug_assert!(upper > 0, "bounded(0) is undefined");
    let mut x = next();
    let mut m = (x as u128) * (upper as u128);
    let mut l = m as u64;
    if l < upper {
        // threshold = (2^64 - upper) % upper, computed without 128-bit div.
        let threshold = upper.wrapping_neg() % upper;
        while l < threshold {
            x = next();
            m = (x as u128) * (upper as u128);
            l = m as u64;
        }
    }
    (m >> 64) as u64
}

/// 64-bit splitmix generator.
///
/// Seedable, deterministic, fast. The integer stream is exactly the
/// classic splitmix64. Use [`bounded`](Self::bounded) for uniform integers
/// in a range and [`gen_bool`](Self::gen_bool) for a fair coin.
#[derive(Clone, Debug)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Construct from a seed, applying the central [`fold_zero_seed`] rule.
    #[inline]
    pub fn new(seed: u64) -> Self {
        Self {
            state: fold_zero_seed(seed),
        }
    }

    /// Next uniform `u64` from the stream.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        splitmix64_step(&mut self.state)
    }

    /// Uniform integer in `0..upper`, unbiased (Lemire). Panics if
    /// `upper == 0`.
    #[inline]
    pub fn bounded(&mut self, upper: usize) -> usize {
        assert!(upper > 0, "SplitMix64::bounded requires upper > 0");
        lemire_bounded(upper as u64, || self.next_u64()) as usize
    }

    /// Fair coin flip.
    #[inline]
    pub fn gen_bool(&mut self) -> bool {
        self.next_u64() & 1 == 1
    }
}

/// xoshiro256++ generator, seeded via a splitmix64 init.
///
/// This is the primitive shared by the core decomposition algorithms. The
/// `u64` stream is the reference xoshiro256++ sequence; floating-point
/// helpers ([`next_f64`](Self::next_f64), [`next_normal`](Self::next_normal))
/// derive from that stream and so are themselves deterministic per seed.
#[derive(Clone, Debug)]
pub struct Xoshiro256pp {
    state: [u64; 4],
}

impl Xoshiro256pp {
    /// Construct from a seed and fill the 4-word state from a splitmix64
    /// stream.
    ///
    /// The splitmix init word is pre-advanced by one `GOLDEN_GAMMA` before
    /// the first `splitmix64_step`. This reproduces, byte-for-byte, the
    /// legacy `ica.rs::Xoshiro256pp` seeding (`sm = seed + GAMMA`, then a
    /// closure that adds another `GAMMA` per word) for EVERY seed. Preserving
    /// that exact stream is load-bearing: the decomposition family (ICA, NMF,
    /// VCA, GSEA, archetype null, alignment bootstrap) consumes only the raw
    /// `next_normal()` / `next_u64()` stream, so its outputs MUST be
    /// unchanged by this consolidation. Only the *bounded-integer* draws
    /// change (modulo → unbiased), which is the intentional correctness fix.
    ///
    /// NOTE: unlike [`SplitMix64`], this path deliberately does NOT apply
    /// [`fold_zero_seed`]. Legacy `ica.rs` never folded the seed, so folding
    /// here would shift the seed-0 stream by one `GOLDEN_GAMMA` and silently
    /// change ICA/NMF/VCA/decompose point estimates at `seed == 0`. The raw
    /// `seed + GAMMA` splitmix expansion cannot produce an all-zero xoshiro
    /// state for any seed (including 0), so the zero-fold is unnecessary here.
    #[inline]
    pub fn new(seed: u64) -> Self {
        let mut sm = seed.wrapping_add(GOLDEN_GAMMA);
        let state = [
            splitmix64_step(&mut sm),
            splitmix64_step(&mut sm),
            splitmix64_step(&mut sm),
            splitmix64_step(&mut sm),
        ];
        Self { state }
    }

    /// Next uniform `u64` from the reference xoshiro256++ stream.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = self.state[0]
            .wrapping_add(self.state[3])
            .rotate_left(23)
            .wrapping_add(self.state[0]);
        let t = self.state[1] << 17;
        self.state[2] ^= self.state[0];
        self.state[3] ^= self.state[1];
        self.state[1] ^= self.state[2];
        self.state[0] ^= self.state[3];
        self.state[2] ^= t;
        self.state[3] = self.state[3].rotate_left(45);
        result
    }

    /// Uniform `[0, 1)` double.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) * (1.0_f64 / ((1u64 << 53) as f64))
    }

    /// Standard normal draw via Box-Muller (one of the paired draws).
    #[inline]
    pub fn next_normal(&mut self) -> f64 {
        let mut u1 = self.next_f64();
        while u1 <= f64::MIN_POSITIVE {
            u1 = self.next_f64();
        }
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    /// Uniform integer in `0..upper`, unbiased (Lemire). Panics if
    /// `upper == 0`.
    #[inline]
    pub fn bounded(&mut self, upper: usize) -> usize {
        assert!(upper > 0, "Xoshiro256pp::bounded requires upper > 0");
        lemire_bounded(upper as u64, || self.next_u64()) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix_is_deterministic_per_seed() {
        let mut a = SplitMix64::new(20260418);
        let mut b = SplitMix64::new(20260418);
        for _ in 0..256 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn xoshiro_is_deterministic_per_seed() {
        let mut a = Xoshiro256pp::new(20260418);
        let mut b = Xoshiro256pp::new(20260418);
        for _ in 0..256 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn zero_seed_folds_to_golden_gamma_consistently() {
        // The two paths that previously diverged: some generators folded
        // seed==0, others didn't. Now both fold identically.
        let mut a = SplitMix64::new(0);
        let mut b = SplitMix64::new(GOLDEN_GAMMA);
        for _ in 0..16 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn xoshiro_seed_zero_preserves_legacy_stream() {
        // Xoshiro256pp deliberately does NOT fold seed==0: legacy ica.rs seeded
        // from `sm = seed + GAMMA` directly, so the seed-0 state must be the
        // splitmix expansion of the RAW seed (mix(2G), mix(3G), mix(4G),
        // mix(5G)), preserving ICA/NMF/decompose point estimates at seed 0.
        let mut sm = 0u64.wrapping_add(GOLDEN_GAMMA);
        let expected = [
            splitmix64_step(&mut sm),
            splitmix64_step(&mut sm),
            splitmix64_step(&mut sm),
            splitmix64_step(&mut sm),
        ];
        assert_eq!(Xoshiro256pp::new(0).state, expected);
        // And folding would have made these equal — they must NOT be, or the
        // seed-0 stream silently shifted by one GOLDEN_GAMMA (the bug fixed here).
        let mut z = Xoshiro256pp::new(0);
        let mut g = Xoshiro256pp::new(GOLDEN_GAMMA);
        assert_ne!(z.next_u64(), g.next_u64());
    }

    #[test]
    fn bounded_stays_in_range() {
        let mut r = SplitMix64::new(7);
        for _ in 0..10_000 {
            let v = r.bounded(5);
            assert!(v < 5);
        }
        let mut x = Xoshiro256pp::new(7);
        for _ in 0..10_000 {
            let v = x.bounded(1);
            assert_eq!(v, 0);
        }
    }

    #[test]
    fn bounded_is_unbiased_within_tolerance() {
        // Pick an upper that does NOT divide 2^64 so the old modulo draw
        // would be measurably biased; check the unbiased draw is flat.
        let upper = 6usize;
        let n = 600_000usize;
        let mut counts = [0usize; 6];
        let mut r = SplitMix64::new(123456789);
        for _ in 0..n {
            counts[r.bounded(upper)] += 1;
        }
        let expected = n as f64 / upper as f64;
        for c in counts {
            let dev = (c as f64 - expected).abs() / expected;
            assert!(dev < 0.02, "bucket deviation {dev} too large: {counts:?}");
        }
    }

    #[test]
    fn derive_sub_seed_is_stable_and_distinct() {
        let s0 = derive_sub_seed(42, 0);
        let s1 = derive_sub_seed(42, 1);
        assert_eq!(s0, derive_sub_seed(42, 0));
        assert_ne!(s0, s1);
    }
}
