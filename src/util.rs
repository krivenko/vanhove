//! Utility functions.

use bilby::{
    QuadratureError, QuadratureResult, adaptive_integrate, integrate_infinite,
    integrate_semi_infinite_lower, integrate_semi_infinite_upper,
};

/// Dispatch for the power function $u^r$ with a fixed exponent `r`.
///
/// `powf()` costs an order of magnitude more than `sqrt()`, which is worth avoiding
/// when the same exponent is reused, as in a quadrature integrand.
#[allow(dead_code)]
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

#[allow(dead_code)]
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

/// Call a bilby adaptive integration function depending on the integration limits.
pub fn bilby_integrate<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    tol: f64,
) -> Result<QuadratureResult<f64>, QuadratureError> {
    let a_inf = a.is_infinite();
    let b_inf = b.is_infinite();
    // Choose a bilby call depending on the integration limits
    match (a_inf, b_inf) {
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
    fn bilby_integrate() {
        let pi = std::f64::consts::PI;
        let inf = f64::INFINITY;
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| x.cos() * x.cos(), -pi, pi, 1e-12)
                .unwrap()
                .value,
            pi,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| (-x * x / 2.0).exp(), -inf, inf, 1e-12)
                .unwrap()
                .value,
            (2.0 * pi).sqrt(),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| (2.0 * x).exp(), -inf, 0.0, 1e-12)
                .unwrap()
                .value,
            0.5,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| (-2.0 * x).exp(), 0.0, inf, 1e-12)
                .unwrap()
                .value,
            0.5,
            epsilon = 1e-12
        );
    }
}
