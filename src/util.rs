//! Utility functions.

use std::f64::consts::PI;

use special::Gamma;

use bilby::{
    QuadratureError, QuadratureResult, adaptive_integrate, integrate_infinite,
    integrate_semi_infinite_lower, integrate_semi_infinite_upper,
};

use crate::segment::Segment;

/// Dispatch for the power function $u^r$ with a fixed exponent `r`.
///
/// `powf()` costs an order of magnitude more than `sqrt()`, which is worth avoiding
/// when the same exponent is reused, as in a quadrature integrand.
#[derive(Debug, Clone, Copy)]
pub enum PowKind {
    /// $u^0 = 1$.
    One,
    /// $u^1 = u$.
    Identity,
    /// $u^{1/2}$.
    Sqrt,
    /// $u^{-1/2}$.
    InvSqrt,
    /// $u^n$ through `powi()`.
    Powi(i32),
    /// $u^r$ through `powf()`.
    Powf(f64),
}

impl PowKind {
    /// Largest integer exponent still handled by `powi()`.
    const MAX_INT_EXPONENT: f64 = 32.0;

    /// Dispatch to use for the exponent `r`.
    pub fn of(r: f64) -> PowKind {
        if r == 0.0 {
            PowKind::One
        } else if r == 1.0 {
            PowKind::Identity
        } else if r == 0.5 {
            PowKind::Sqrt
        } else if r == -0.5 {
            PowKind::InvSqrt
        } else if r.fract() == 0.0 && r.abs() <= Self::MAX_INT_EXPONENT {
            PowKind::Powi(r as i32)
        } else {
            PowKind::Powf(r)
        }
    }

    /// $u^r$ for $u \geq 0$, agreeing with `powf()` down to $u = 0$.
    pub fn eval(&self, u: f64) -> f64 {
        match *self {
            PowKind::One => 1.0,
            PowKind::Identity => u,
            PowKind::Sqrt => u.sqrt(),
            PowKind::InvSqrt => 1.0 / u.sqrt(),
            PowKind::Powi(n) => u.powi(n),
            PowKind::Powf(r) => u.powf(r),
        }
    }
}

/// Kahan-Babuška-Neumaier summation algorithm.
pub fn kahan_babushka_neumaier_sum<I: Iterator<Item = f64>>(input: I) -> f64 {
    let (mut sum, mut c) = (0.0f64, 0.0f64);
    for x in input {
        let t = sum + x;
        c += if sum.abs() >= x.abs() {
            (sum - t) + x
        } else {
            (x - t) + sum
        };
        sum = t;
    }
    sum + c
}

/// Chebyshev coefficients of `f` over $[-1, 1]$ from its values at `n` nodes.
///
/// The nodes are those of Gauss-Chebyshev quadrature of the first kind,
/// $x_j = \cos\frac{\pi(j+1/2)}{n}$, which lie strictly inside the interval, so that
/// `f` is never sampled at an end point. The $k = 0$ coefficient is returned already
/// halved, making the expansion $\sum_k c_k T_k(x)$ as [`clenshaw_chebyshev()`]
/// evaluates it.
pub fn chebyshev_coeffs<F: Fn(f64) -> f64>(n: usize, f: F) -> Vec<f64> {
    let theta: Vec<f64> = (0..n).map(|j| PI * (j as f64 + 0.5) / n as f64).collect();
    let values: Vec<f64> = theta.iter().map(|&t| f(t.cos())).collect();
    (0..n)
        .map(|k| {
            let sum = kahan_babushka_neumaier_sum(
                values
                    .iter()
                    .zip(&theta)
                    .map(|(v, &t)| v * (k as f64 * t).cos()),
            );
            let c = 2.0 * sum / n as f64;
            if k == 0 { 0.5 * c } else { c }
        })
        .collect()
}

/// Clenshaw recurrence for the Chebyshev series $\sum_k c_k T_k(x)$.
pub fn clenshaw_chebyshev(coeffs: &[f64], x: f64) -> f64 {
    let Some((c0, rest)) = coeffs.split_first() else {
        return 0.0;
    };
    let (mut b1, mut b2) = (0.0f64, 0.0f64);
    for c in rest.iter().rev() {
        let b0 = 2.0 * x * b1 - b2 + c;
        b2 = b1;
        b1 = b0;
    }
    x * b1 - b2 + c0
}

/// Call a bilby adaptive integration function depending on the integration limits.
pub fn bilby_integrate<F: Fn(f64) -> f64>(
    f: F,
    segment: Segment,
    tol: f64,
) -> Result<QuadratureResult<f64>, QuadratureError> {
    let (a, b) = (segment.min(), segment.max());
    // Choose a bilby call depending on the integration limits
    match (a.is_infinite(), b.is_infinite()) {
        (false, false) => adaptive_integrate(f, a, b, tol),
        (false, true) => {
            assert!(b > 0.0);
            integrate_semi_infinite_upper(f, a, tol)
        }
        (true, false) => {
            assert!(a < 0.0);
            integrate_semi_infinite_lower(f, b, tol)
        }
        (true, true) => {
            assert!(a < 0.0 && b > 0.0);
            integrate_infinite(f, tol)
        }
    }
}

/// Call [`bilby_integrate`] and return zero where it refuses the request.
///
/// Refusal is about the request rather than the integrand: a tolerance that is no
/// number, an interval that is no interval. An integrand the quadrature cannot resolve
/// is not refused at all - it comes back as a value, only a less accurate one.
pub fn bilby_integrate_or_0<F: Fn(f64) -> f64>(f: F, segment: Segment, tol: f64) -> f64 {
    bilby_integrate(f, segment, tol).map_or(0.0, |r| r.value)
}

/// The whole row $\binom{n}{0}, \binom{n}{1}, \ldots, \binom{n}{n}$.
///
/// Each coefficient follows from the one before it, which is cheaper than asking for
/// them one at a time and is how they are wanted wherever a binomial expansion is
/// summed over.
pub fn binomials(n: usize) -> Vec<f64> {
    let mut c_row = vec![1.0; n + 1];
    for k in 1..=n {
        c_row[k] = c_row[k - 1] * (n - k + 1) as f64 / k as f64;
    }
    c_row
}

/// Rectangular table of numbers, indexed by row and then by column.
pub type Table = Vec<Vec<f64>>;

/// Elementwise difference of two tables.
pub fn subtract_tables(x: &Table, y: &Table) -> Table {
    x.iter()
        .zip(y)
        .map(|(rx, ry)| rx.iter().zip(ry).map(|(a, b)| a - b).collect())
        .collect()
}

/// A table with every entry turned over.
pub fn negate_table(x: &Table) -> Table {
    x.iter().map(|r| r.iter().map(|v| -v).collect()).collect()
}

/// $(-1)^n$.
pub fn alternating_sign(n: usize) -> f64 {
    if n.is_multiple_of(2) { 1.0 } else { -1.0 }
}

/// Whether `x` is one of $0, 1, 2, \ldots$
pub fn is_natural(x: f64) -> bool {
    x >= 0.0 && x.fract() == 0.0
}

/// Polygamma function $\psi^{(n)}(x)$, the $n$-th derivative of the digamma function.
///
/// Diverges at the non-positive integers, where $\Gamma$ has its poles.
///
/// The recurrence $\psi^{(n)}(x) = \psi^{(n)}(x+1) - (-1)^n n!\\,x^{-n-1}$ walks the
/// argument up to where the asymptotic series converges, which is what carries the
/// negative arguments: the reflection formula is never needed.
pub fn polygamma(n: u32, x: f64) -> f64 {
    /// $B_{2k}$ for $k = 1, 2, \ldots$
    const BERNOULLI: [f64; 8] = [
        1.0 / 6.0,
        -1.0 / 30.0,
        1.0 / 42.0,
        -1.0 / 30.0,
        5.0 / 66.0,
        -691.0 / 2730.0,
        7.0 / 6.0,
        -3617.0 / 510.0,
    ];
    // The series needs a larger argument the higher the order, the terms growing as
    // (2k+n)! before the powers of x beat them down
    let large = 10.0 + 4.0 * f64::from(n);

    let nf = f64::from(n);
    let factorial = |k: f64| Gamma::gamma(k + 1.0);
    let sign = alternating_sign(n as usize);

    // Walk up to where the series converges, collecting what the recurrence sheds:
    // ψ^(n)(x) = ψ^(n)(x+1) - (-1)^n n! x^{-n-1}, and ψ(x) = ψ(x+1) - 1/x
    let (mut x, mut shed) = (x, 0.0f64);
    while x < large {
        shed -= if n == 0 {
            x.recip()
        } else {
            sign * factorial(nf) * x.powf(-nf - 1.0)
        };
        x += 1.0;
    }

    // ψ(x) ~ ln x - 1/(2x) - Σ_k B_{2k} / (2k x^{2k}), and for n >= 1
    // ψ^(n)(x) ~ (-1)^{n-1} [ (n-1)!/x^n + n!/(2x^{n+1})
    //                         + Σ_k B_{2k} (2k+n-1)! / ((2k)! x^{2k+n}) ]
    let mut series = if n == 0 {
        x.ln() - 0.5 / x
    } else {
        factorial(nf - 1.0) * x.powf(-nf) + 0.5 * factorial(nf) * x.powf(-nf - 1.0)
    };
    for (k, b) in BERNOULLI.iter().enumerate() {
        let k2 = 2.0 * (k as f64 + 1.0);
        series += if n == 0 {
            -b / (k2 * x.powf(k2))
        } else {
            b * factorial(k2 + nf - 1.0) / (factorial(k2) * x.powf(k2 + nf))
        };
    }
    shed + if n == 0 { series } else { -sign * series }
}

/// Fermi step function.
pub fn fermi(x: f64) -> f64 {
    if x <= 0.0 {
        1.0 / (1.0 + x.exp())
    } else {
        1.0 - fermi(-x)
    }
}

#[cfg(test)]
mod tests {
    use crate::util;
    use approx::{assert_abs_diff_eq, assert_relative_eq};
    use special::Gamma;

    #[test]
    fn kahan_babushka_neumaier_sum() {
        use std::f64::consts::{E, PI};
        let v = vec![10000.0f64, PI, E];
        assert_abs_diff_eq!(
            util::kahan_babushka_neumaier_sum(v.into_iter()),
            10000.0 + PI + E,
            epsilon = 1e-12
        );

        // Neumaier's example: naive summation and the plain Kahan algorithm both
        // return 0, while the compensated sum is exact.
        let v = vec![1.0f64, 1e100, 1.0, -1e100];
        assert_eq!(util::kahan_babushka_neumaier_sum(v.into_iter()), 2.0);
    }

    #[test]
    fn pow_kind() {
        use util::PowKind;

        // Every shortcut agrees with powf() to within a few ulps. It is not bitwise
        // agreement: powi() multiplies repeatedly and InvSqrt rounds twice, while
        // powf() is correctly rounded.
        for r in [0.0f64, 1.0, 0.5, -0.5, 3.0, -7.0, 4.0, 2.5, -0.25] {
            let kind = PowKind::of(r);
            for u in [0.0f64, 1e-8, 0.25, 1.0, 7.5] {
                assert_relative_eq!(kind.eval(u), u.powf(r), max_relative = 1e-14);
            }
        }

        // u^0 is one everywhere, u^1 the identity
        assert_eq!(PowKind::of(0.0).eval(0.0), 1.0);
        assert_eq!(PowKind::of(1.0).eval(0.0), 0.0);

        // A negative exponent diverges at the origin
        assert_eq!(PowKind::of(-0.5).eval(0.0), f64::INFINITY);

        // Integer exponents beyond the powi() range fall back on powf()
        assert!(matches!(PowKind::of(32.0), PowKind::Powi(32)));
        assert!(matches!(PowKind::of(33.0), PowKind::Powf(_)));
    }

    #[test]
    fn chebyshev() {
        use util::{chebyshev_coeffs, clenshaw_chebyshev};

        // T_3(x) = 4x^3 - 3x is one coefficient and nothing else
        let t3 = |x: f64| 4.0 * x.powi(3) - 3.0 * x;
        let c = chebyshev_coeffs(8, t3);
        assert_relative_eq!(c[3], 1.0, max_relative = 1e-14);
        for (k, ck) in c.iter().enumerate() {
            if k != 3 {
                assert_abs_diff_eq!(*ck, 0.0, epsilon = 1e-14);
            }
        }
        for x in [-1.0, -0.3, 0.0, 0.5, 1.0] {
            assert_abs_diff_eq!(clenshaw_chebyshev(&c, x), t3(x), epsilon = 1e-14);
        }

        // An analytic function converges geometrically, reaching machine precision
        let c = chebyshev_coeffs(24, f64::exp);
        for x in [-1.0, -0.3, 0.0, 0.5, 1.0] {
            assert_relative_eq!(clenshaw_chebyshev(&c, x), x.exp(), max_relative = 1e-14);
        }
        assert!(c[20].abs() < 1e-14, "tail did not decay: {}", c[20]);

        // The nodes stay inside the interval, leaving the end points unsampled
        let _ = chebyshev_coeffs(16, |x| {
            assert!(x.abs() < 1.0, "sampled the end point {x}");
            x
        });

        // An empty series is the zero function
        assert_eq!(clenshaw_chebyshev(&[], 0.5), 0.0);
    }

    #[test]
    fn bilby_integrate() {
        use crate::segment::Segment;
        let pi = std::f64::consts::PI;
        let inf = f64::INFINITY;
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| x.cos() * x.cos(), Segment::new(-pi, pi), 1e-12)
                .unwrap()
                .value,
            pi,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            util::bilby_integrate(
                move |x| (-x * x / 2.0).exp(),
                Segment::new(-inf, inf),
                1e-12
            )
            .unwrap()
            .value,
            (2.0 * pi).sqrt(),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| (2.0 * x).exp(), Segment::new(-inf, 0.0), 1e-12)
                .unwrap()
                .value,
            0.5,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| (-2.0 * x).exp(), Segment::new(0.0, inf), 1e-12)
                .unwrap()
                .value,
            0.5,
            epsilon = 1e-12
        );
    }

    #[test]
    fn bilby_integrate_or_0() {
        use crate::segment::Segment;
        use std::f64::consts::E;

        assert_abs_diff_eq!(
            util::bilby_integrate_or_0(f64::exp, Segment::new(0.0, 1.0), 1e-12),
            E - 1.0,
            epsilon = 1e-12
        );

        // A segment of zero length integrates to nothing
        let point = Segment::new(1.0, 1.0);
        assert_eq!(util::bilby_integrate_or_0(f64::exp, point, 1e-12), 0.0);

        // Where the quadrature refuses the request there is no value to carry back
        let unit = Segment::new(0.0, 1.0);
        assert_eq!(util::bilby_integrate_or_0(f64::exp, unit, -1.0), 0.0);
    }

    #[test]
    fn binomials() {
        assert_eq!(util::binomials(0), vec![1.0]);
        assert_eq!(util::binomials(1), vec![1.0, 1.0]);
        assert_eq!(util::binomials(5), vec![1.0, 5.0, 10.0, 10.0, 5.0, 1.0]);

        // The running product carries a row far past where it would overflow an
        // integer of any convenient width, the rounding staying in the last few digits
        let c_row = util::binomials(60);
        assert_eq!(c_row.len(), 61);
        for (k, &c) in c_row.iter().enumerate() {
            assert_relative_eq!(c, c_row[60 - k], max_relative = 1e-14);
        }
        assert_relative_eq!(
            c_row.iter().sum::<f64>(),
            2f64.powi(60),
            max_relative = 1e-14
        );
    }

    #[test]
    fn alternating_sign() {
        assert_eq!(util::alternating_sign(0), 1.0);
        assert_eq!(util::alternating_sign(1), -1.0);
        assert_eq!(util::alternating_sign(2), 1.0);
        assert_eq!(util::alternating_sign(101), -1.0);
    }

    #[test]
    fn is_natural() {
        for x in [0.0f64, -0.0, 1.0, 2.0, 1e15] {
            assert!(util::is_natural(x), "{x} is one of 0, 1, 2, ...");
        }
        for x in [-1.0f64, -0.5, 0.5, 2.5, f64::NAN, f64::INFINITY] {
            assert!(!util::is_natural(x), "{x} is not");
        }
    }

    #[test]
    fn polygamma() {
        use std::f64::consts::PI;
        const GAMMA: f64 = 0.577_215_664_901_532_9;

        // ψ(1) = -γ, ψ(1/2) = -γ - 2ln2, and the recurrence ψ(x+1) = ψ(x) + 1/x
        assert_relative_eq!(util::polygamma(0, 1.0), -GAMMA, max_relative = 1e-13);
        assert_relative_eq!(
            util::polygamma(0, 0.5),
            -GAMMA - 2.0 * 2f64.ln(),
            max_relative = 1e-13
        );
        assert_relative_eq!(util::polygamma(0, 2.0), 1.0 - GAMMA, max_relative = 1e-13);

        // ψ'(1) = π²/6, ψ'(1/2) = π²/2, ψ''(1) = -2ζ(3), ψ'''(1) = π⁴/15
        assert_relative_eq!(
            util::polygamma(1, 1.0),
            PI.powi(2) / 6.0,
            max_relative = 1e-13
        );
        assert_relative_eq!(
            util::polygamma(1, 0.5),
            PI.powi(2) / 2.0,
            max_relative = 1e-13
        );
        assert_relative_eq!(
            util::polygamma(2, 1.0),
            -2.404_113_806_319_188_4,
            max_relative = 1e-13
        );
        assert_relative_eq!(
            util::polygamma(3, 1.0),
            PI.powi(4) / 15.0,
            max_relative = 1e-13
        );

        // Walking the argument up carries the negative ones, so no reflection formula
        // is needed: ψ(-1/2) = ψ(1/2) + 2
        assert_relative_eq!(
            util::polygamma(0, -0.5),
            2.0 - GAMMA - 2.0 * 2f64.ln(),
            max_relative = 1e-13
        );
        for n in 0..=3u32 {
            for x in [-3.25f64, -0.75, 0.3, 1.7, 12.0] {
                let shed = util::alternating_sign(n as usize)
                    * Gamma::gamma(f64::from(n) + 1.0)
                    * x.powf(-f64::from(n) - 1.0);
                assert_relative_eq!(
                    util::polygamma(n, x + 1.0) - util::polygamma(n, x),
                    shed,
                    max_relative = 1e-11,
                    epsilon = 1e-13
                );
            }
        }
    }

    #[test]
    fn fermi() {
        assert_eq!(util::fermi(0.0), 0.5);
        assert_abs_diff_eq!(
            util::fermi(-1.0),
            std::f64::consts::E / (1.0 + std::f64::consts::E),
            epsilon = 1e-14
        );
        assert_abs_diff_eq!(
            util::fermi(1.0),
            1.0 / (1.0 + std::f64::consts::E),
            epsilon = 1e-14
        );
        assert_eq!(util::fermi(f64::INFINITY), 0.0);
        assert_eq!(util::fermi(f64::NEG_INFINITY), 1.0);
    }
}
