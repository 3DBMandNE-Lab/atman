//! One matrix-multiply entry point, with a vendor BLAS behind it where
//! one exists.
//!
//! # Why this exists
//!
//! NMF's multiplicative updates are three GEMMs per iteration and
//! nothing else that matters. Written as scalar Rust over
//! `Vec<Vec<f64>>` they run at roughly 2.3 GFLOP/s on an M4, which is
//! where unvectorised code sits. Measured on the same machine at the
//! shapes NMF actually uses (n=110, p=2000, k=10), Accelerate's
//! `cblas_dgemm` reaches 126.5 GFLOP/s for `WᵀX` and 68.7 GFLOP/s for
//! `XHᵀ` — 30 to 55 times faster, in the same f64 the rest of the
//! pipeline uses.
//!
//! # Determinism
//!
//! A blocked kernel reassociates the sum, so results move against the
//! scalar path. That is a ONE-TIME RE-BASELINE, not a loss of the
//! guarantee atman actually makes, which is: the same binary, on the
//! same input, with the same seed, produces the same bytes. The run
//! sidecar records `target_triple` and `rustc_version` precisely
//! because results are scoped to a build environment.
//!
//! Accelerate was measured for bit-stability before being adopted: 200
//! repeats of each GEMM at NMF's shapes, zero bitwise-differing runs,
//! including under eight competing CPU-bound processes. It does not
//! adapt its reduction order to load, so `--threads` invariance holds.
//!
//! What would NOT be acceptable is runtime CPU-feature dispatch that
//! varies across machines sharing one `target_triple`: the sidecar
//! would then record two differently-numbered runs as identically
//! provenanced. Accelerate is a fixed system library for a given OS
//! version and target, so it does not have that property.
//!
//! # Portability
//!
//! macOS gets Accelerate. Everything else gets [`gemm_fallback`], a
//! cache-friendly scalar kernel. The two produce DIFFERENT BITS, which
//! is why `target_triple` in the sidecar is load bearing rather than
//! decorative.

/// Whether an operand is used as-is or transposed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Use the matrix as stored.
    N,
    /// Use its transpose.
    T,
}

/// Row-major `C = op(A) · op(B)`, overwriting `C`.
///
/// `m`, `n`, `k` are the dimensions of the PRODUCT: `op(A)` is `m × k`,
/// `op(B)` is `k × n`, and `C` is `m × n`. `lda`, `ldb` and `ldc` are
/// the row strides of the STORED matrices, so for `Op::T` the stride is
/// that of the untransposed operand.
///
/// # Panics
///
/// Panics if any slice is shorter than its dimensions require. The
/// unsafe call below depends on those bounds, so they are checked here
/// rather than assumed from the caller.
#[allow(clippy::too_many_arguments)]
pub fn gemm(
    m: usize,
    n: usize,
    k: usize,
    op_a: Op,
    a: &[f64],
    lda: usize,
    op_b: Op,
    b: &[f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
) {
    // Rows of A as stored: m when N, k when T. Same for B.
    let a_rows = if op_a == Op::N { m } else { k };
    let b_rows = if op_b == Op::N { k } else { n };
    assert!(lda >= if op_a == Op::N { k } else { m }, "lda too small");
    assert!(ldb >= if op_b == Op::N { n } else { k }, "ldb too small");
    assert!(ldc >= n, "ldc too small");
    assert!(a.len() >= a_rows.saturating_mul(lda), "a too short");
    assert!(b.len() >= b_rows.saturating_mul(ldb), "b too short");
    assert!(c.len() >= m.saturating_mul(ldc), "c too short");
    if m == 0 || n == 0 || k == 0 {
        for v in c.iter_mut().take(m * ldc) {
            *v = 0.0;
        }
        return;
    }

    #[cfg(target_os = "macos")]
    {
        gemm_accelerate(m, n, k, op_a, a, lda, op_b, b, ldb, c, ldc);
    }
    #[cfg(not(target_os = "macos"))]
    {
        gemm_fallback(m, n, k, op_a, a, lda, op_b, b, ldb, c, ldc);
    }
}

#[cfg(target_os = "macos")]
mod accelerate {
    use std::os::raw::c_int;

    pub const ROW_MAJOR: c_int = 101;
    pub const NO_TRANS: c_int = 111;
    pub const TRANS: c_int = 112;
    pub const UPPER: c_int = 121;

    #[link(name = "Accelerate", kind = "framework")]
    extern "C" {
        #[allow(clippy::too_many_arguments)]
        pub fn cblas_dsyrk(
            order: c_int,
            uplo: c_int,
            trans: c_int,
            n: c_int,
            k: c_int,
            alpha: f64,
            a: *const f64,
            lda: c_int,
            beta: f64,
            c: *mut f64,
            ldc: c_int,
        );

        pub fn cblas_dgemm(
            order: c_int,
            transa: c_int,
            transb: c_int,
            m: c_int,
            n: c_int,
            k: c_int,
            alpha: f64,
            a: *const f64,
            lda: c_int,
            b: *const f64,
            ldb: c_int,
            beta: f64,
            c: *mut f64,
            ldc: c_int,
        );
    }
}

#[cfg(target_os = "macos")]
mod threading {
    use std::os::raw::{c_char, c_int, c_void};
    use std::sync::Once;

    /// `BLAS_THREADING_SINGLE_THREADED` from `vecLib/thread_api.h`.
    const SINGLE_THREADED: u32 = 1;
    /// macOS `RTLD_DEFAULT`: search every loaded image.
    const RTLD_DEFAULT: *mut c_void = -2isize as *mut c_void;

    extern "C" {
        fn dlsym(handle: *mut c_void, symbol: *const c_char) -> *mut c_void;
    }

    static INIT: Once = Once::new();

    /// Put Accelerate's BLAS in single-threaded mode, once per process.
    ///
    /// MEASURED, not assumed. NMF's GEMMs are skinny — one dimension is
    /// `k`, typically 10 — so the work per call is small and Accelerate's
    /// thread dispatch costs more than it saves. At n=110, p=10000, k=10,
    /// per-iteration cost in situ:
    ///
    /// | threading | WᵀX | XHᵀ | total/iter |
    /// |---|---|---|---|
    /// | Accelerate's default | 0.743 ms | 0.864 ms | 2.15 ms |
    /// | 2 threads | 0.466 ms | 0.548 ms | 1.38 ms |
    /// | single-threaded | 0.262 ms | 0.360 ms | 0.93 ms |
    ///
    /// 2.3x faster single-threaded. The same GEMMs measured standalone in
    /// a tight loop ran at 0.179 ms and 0.304 ms, which is why the
    /// isolated benchmark did not predict the in-loop cost: back-to-back
    /// calls keep Accelerate's threads hot, and calls separated by other
    /// work do not.
    ///
    /// It also removes a determinism question rather than raising one: a
    /// single-threaded BLAS cannot vary its reduction order with thread
    /// count or machine load.
    ///
    /// This is process-global. It is right for the skinny shapes atman
    /// decomposes and would be the wrong default for large square GEMMs,
    /// so revisit it if this module gains such a caller.
    ///
    /// `BLASSetThreading` needs macOS 15. Looked up with `dlsym` rather
    /// than linked directly so that an older system silently keeps
    /// Accelerate's default instead of failing to launch.
    pub fn set_single_threaded_once() {
        // SAFETY: three things make this sound.
        //
        // `dlsym` with RTLD_DEFAULT is safe to call with any C string and
        // returns null when the symbol is absent, which is checked before
        // the pointer is used — that null check is the whole reason for
        // resolving dynamically rather than linking directly.
        //
        // The transmuted signature matches Apple's declaration in
        // `vecLib/thread_api.h`: `int BLASSetThreading(const enum
        // BLAS_THREADING)`, where the enum's underlying type is spelled
        // `unsigned int` in that header. So `extern "C" fn(u32) -> c_int`
        // is the correct ABI, and `SINGLE_THREADED` is a value the enum
        // defines rather than an arbitrary integer.
        //
        // `Once` guarantees the call happens exactly once per process
        // even under concurrent `gemm` calls from multiple threads.
        INIT.call_once(|| unsafe {
            let sym = dlsym(RTLD_DEFAULT, c"BLASSetThreading".as_ptr());
            if !sym.is_null() {
                let f: extern "C" fn(u32) -> c_int = std::mem::transmute(sym);
                let _ = f(SINGLE_THREADED);
            }
        });
    }
}

/// Every dimension handed to Accelerate must survive the `usize` to
/// `c_int` cast.
///
/// CBLAS takes 32-bit signed dimensions. A value above `i32::MAX` would
/// truncate — silently, and possibly to a negative — and the library
/// would then read outside the slices whose lengths were checked in
/// `usize`. atman's matrices are far below this (tens of thousands of
/// assays at most), so the check never fires in practice; it exists
/// because the alternative to firing is memory unsafety rather than a
/// wrong answer.
#[cfg(target_os = "macos")]
fn assert_fits_c_int(dims: [usize; 6]) {
    for d in dims {
        assert!(
            d <= i32::MAX as usize,
            "dimension {d} exceeds the 32-bit limit of the CBLAS interface"
        );
    }
}

#[cfg(target_os = "macos")]
#[allow(clippy::too_many_arguments)]
fn gemm_accelerate(
    m: usize,
    n: usize,
    k: usize,
    op_a: Op,
    a: &[f64],
    lda: usize,
    op_b: Op,
    b: &[f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
) {
    use accelerate::*;
    threading::set_single_threaded_once();
    assert_fits_c_int([m, n, k, lda, ldb, ldc]);
    let ta = if op_a == Op::N { NO_TRANS } else { TRANS };
    let tb = if op_b == Op::N { NO_TRANS } else { TRANS };
    // SAFETY: this function is private and reachable only through
    // `gemm`, which asserts every slice length against its dimensions
    // before dispatching. `assert_fits_c_int` above then rules out the
    // one way those checks could be defeated — a `usize` dimension
    // truncating in the cast to `c_int`. The pointers come from live
    // borrows that outlive the call, and `cblas_dgemm` reads `a` and `b`
    // and writes `c` only within the bounds just validated.
    unsafe {
        cblas_dgemm(
            ROW_MAJOR,
            ta,
            tb,
            m as _,
            n as _,
            k as _,
            1.0,
            a.as_ptr(),
            lda as _,
            b.as_ptr(),
            ldb as _,
            0.0,
            c.as_mut_ptr(),
            ldc as _,
        );
    }
}

/// Row-major symmetric rank-k update: `C = op(A) · op(A)ᵀ`, `C` being
/// `n × n` and symmetric.
///
/// `Op::N` computes `A Aᵀ` from an `n × k` A; `Op::T` computes `Aᵀ A`
/// from a `k × n` A. Only about half the multiply-adds of the
/// equivalent GEMM are performed, because the result is symmetric — the
/// lower triangle is mirrored from the upper afterwards so callers can
/// treat `C` as a full dense matrix.
///
/// NMF forms two of these per iteration, `H Hᵀ` and `Wᵀ W`. `H Hᵀ` was
/// 0.19 ms of a 0.93 ms iteration at n=110, p=10000, k=10 when computed
/// as a general GEMM.
///
/// # Panics
///
/// Panics if a slice is shorter than its dimensions require.
pub fn syrk(n: usize, k: usize, op: Op, a: &[f64], lda: usize, c: &mut [f64], ldc: usize) {
    let a_rows = if op == Op::N { n } else { k };
    assert!(lda >= if op == Op::N { k } else { n }, "lda too small");
    assert!(ldc >= n, "ldc too small");
    assert!(a.len() >= a_rows.saturating_mul(lda), "a too short");
    assert!(c.len() >= n.saturating_mul(ldc), "c too short");
    if n == 0 {
        return;
    }

    #[cfg(target_os = "macos")]
    {
        use accelerate::*;
        threading::set_single_threaded_once();
        assert_fits_c_int([n, n, k, lda, ldc, ldc]);
        let tr = if op == Op::N { NO_TRANS } else { TRANS };
        // SAFETY: `syrk` asserts every slice length against its
        // dimensions above, and `assert_fits_c_int` rules out truncation
        // in the cast to `c_int`. The pointers come from live borrows
        // that outlive the call, and `cblas_dsyrk` reads `a` and writes
        // the upper triangle of `c` only within those bounds.
        unsafe {
            cblas_dsyrk(
                ROW_MAJOR,
                UPPER,
                tr,
                n as _,
                k as _,
                1.0,
                a.as_ptr(),
                lda as _,
                0.0,
                c.as_mut_ptr(),
                ldc as _,
            );
        }
    }
    #[cfg(not(target_os = "macos"))]
    {
        let (op_b, ldb) = if op == Op::N {
            (Op::T, lda)
        } else {
            (Op::N, lda)
        };
        gemm_fallback(n, n, k, op, a, lda, op_b, a, ldb, c, ldc);
        return;
    }

    // Mirror the upper triangle into the lower so callers see a full
    // dense symmetric matrix.
    #[cfg(target_os = "macos")]
    for i in 0..n {
        for j in 0..i {
            c[i * ldc + j] = c[j * ldc + i];
        }
    }
}

/// Portable scalar GEMM, used where no vendor BLAS is linked.
///
/// Loop order is i-l-j so both the `b` read and the `c` write walk
/// contiguously. This is the same order the hand-written kernels in
/// `nmf.rs` use, and it is NOT bit-identical to Accelerate.
#[allow(clippy::too_many_arguments)]
pub fn gemm_fallback(
    m: usize,
    n: usize,
    k: usize,
    op_a: Op,
    a: &[f64],
    lda: usize,
    op_b: Op,
    b: &[f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
) {
    for i in 0..m {
        for v in c.iter_mut().skip(i * ldc).take(n) {
            *v = 0.0;
        }
    }
    for i in 0..m {
        for l in 0..k {
            let a_il = match op_a {
                Op::N => a[i * lda + l],
                Op::T => a[l * lda + i],
            };
            if a_il == 0.0 {
                continue;
            }
            match op_b {
                Op::N => {
                    let brow = &b[l * ldb..l * ldb + n];
                    let crow = &mut c[i * ldc..i * ldc + n];
                    for (dst, bv) in crow.iter_mut().zip(brow.iter()) {
                        *dst += a_il * bv;
                    }
                }
                Op::T => {
                    for j in 0..n {
                        c[i * ldc + j] += a_il * b[j * ldb + l];
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive(m: usize, n: usize, k: usize, op_a: Op, a: &[f64], op_b: Op, b: &[f64]) -> Vec<f64> {
        let mut c = vec![0.0; m * n];
        for i in 0..m {
            for j in 0..n {
                let mut acc = 0.0;
                for l in 0..k {
                    let av = match op_a {
                        Op::N => a[i * k + l],
                        Op::T => a[l * m + i],
                    };
                    let bv = match op_b {
                        Op::N => b[l * n + j],
                        Op::T => b[j * k + l],
                    };
                    acc += av * bv;
                }
                c[i * n + j] = acc;
            }
        }
        c
    }

    fn seeded(len: usize, mut s: u64) -> Vec<f64> {
        (0..len)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((s >> 33) as f64 / (u32::MAX as f64)) - 0.5
            })
            .collect()
    }

    /// The dispatched GEMM must agree with a textbook triple loop to
    /// floating-point tolerance across every transpose combination and
    /// the skinny shapes NMF uses. Not bit-identity: a blocked kernel
    /// reassociates, which is the whole reason this needs a tolerance
    /// and a re-baseline rather than a golden diff.
    #[test]
    fn gemm_agrees_with_a_naive_triple_loop_in_every_orientation() {
        for &(m, n, k) in &[(4, 5, 3), (110, 10, 2000), (10, 2000, 110), (1, 1, 1)] {
            for &op_a in &[Op::N, Op::T] {
                for &op_b in &[Op::N, Op::T] {
                    let (a_rows, lda) = if op_a == Op::N { (m, k) } else { (k, m) };
                    let (b_rows, ldb) = if op_b == Op::N { (k, n) } else { (n, k) };
                    let a = seeded(a_rows * lda, 11);
                    let b = seeded(b_rows * ldb, 22);
                    let want = naive(m, n, k, op_a, &a, op_b, &b);

                    let mut got = vec![0.0; m * n];
                    gemm(m, n, k, op_a, &a, lda, op_b, &b, ldb, &mut got, n);

                    let mut fb = vec![0.0; m * n];
                    gemm_fallback(m, n, k, op_a, &a, lda, op_b, &b, ldb, &mut fb, n);

                    for idx in 0..m * n {
                        let scale = want[idx].abs().max(1.0);
                        assert!(
                            (got[idx] - want[idx]).abs() / scale < 1e-10,
                            "dispatched gemm {m}x{n}x{k} {op_a:?}{op_b:?} idx {idx}: \
                             {} vs {}",
                            got[idx],
                            want[idx]
                        );
                        assert!(
                            (fb[idx] - want[idx]).abs() / scale < 1e-10,
                            "fallback gemm {m}x{n}x{k} {op_a:?}{op_b:?} idx {idx}: \
                             {} vs {}",
                            fb[idx],
                            want[idx]
                        );
                    }
                }
            }
        }
    }

    /// `syrk` must agree with the equivalent GEMM, and must fill BOTH
    /// triangles — callers read `C` as a full dense matrix, so a syrk
    /// that left the lower half untouched would silently zero half of
    /// every Gram matrix in the solve.
    #[test]
    fn syrk_matches_the_equivalent_gemm_and_fills_both_triangles() {
        for &(n, k) in &[(5usize, 7usize), (10, 10000), (10, 110), (1, 3)] {
            for &op in &[Op::N, Op::T] {
                let (rows, lda) = if op == Op::N { (n, k) } else { (k, n) };
                let a = seeded(rows * lda, 5);
                let op_b = if op == Op::N { Op::T } else { Op::N };

                let mut want = vec![0.0; n * n];
                gemm(n, n, k, op, &a, lda, op_b, &a, lda, &mut want, n);
                let mut got = vec![0.0; n * n];
                syrk(n, k, op, &a, lda, &mut got, n);

                for i in 0..n {
                    for j in 0..n {
                        let idx = i * n + j;
                        let scale = want[idx].abs().max(1.0);
                        assert!(
                            (got[idx] - want[idx]).abs() / scale < 1e-10,
                            "syrk n={n} k={k} {op:?} at ({i},{j}): {} vs {}",
                            got[idx],
                            want[idx]
                        );
                    }
                }
                // Symmetry, explicitly.
                for i in 0..n {
                    for j in 0..n {
                        assert_eq!(got[i * n + j], got[j * n + i], "not symmetric");
                    }
                }
            }
        }
    }

    /// The dispatched GEMM must be bit-stable run to run, or every
    /// byte-identity guarantee downstream of it is void. Measured
    /// separately under CPU contention before Accelerate was adopted.
    #[test]
    fn gemm_is_bit_stable_across_repeats() {
        let (m, n, k) = (110, 10, 2000);
        let a = seeded(m * k, 7);
        let b = seeded(n * k, 8);
        let mut first = vec![0.0; m * n];
        gemm(m, n, k, Op::N, &a, k, Op::T, &b, k, &mut first, n);
        for r in 0..50 {
            let mut again = vec![0.0; m * n];
            gemm(m, n, k, Op::N, &a, k, Op::T, &b, k, &mut again, n);
            assert_eq!(first, again, "gemm differed on repeat {r}");
        }
    }
}
