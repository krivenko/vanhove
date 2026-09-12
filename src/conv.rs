//! Convolution of the continuous parts of two spectral functions.

use bilby::adaptive_integrate_with_breaks;

use crate::singularity::Singularity;
use crate::{ContinuousSF, SpectralFunction, non_smooth};

/// Value of a single continuous contribution, $R(\omega) + \sum_p S_p(\omega)$.
fn value_at(csf: &dyn ContinuousSF, omega: f64) -> f64 {
    let singular: f64 = csf.singularities().iter().map(|s| s.value(omega)).sum();
    csf.regular(omega) + singular
}

/// $\int C_1(\nu) C_2(\omega-\nu) d\nu$ over the overlap of the two supports.
///
/// The singular parts of $C_1$ are handled the way [`SpectralFunction::integrate()`]
/// handles them, by subtracting the value of the other factor at $\Omega_p$ and
/// restoring it through $\int S_p$ in closed form. Handing the bare product to the
/// quadrature instead leaves an inverse square root sitting on an end of the interval,
/// where Gauss-Kronrod converges too slowly to notice it is not converging.
fn pair(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF, omega: f64, tol: f64) -> f64 {
    let (a1, b1) = c1.support();
    let (a2, b2) = c2.support();
    let (lo, hi) = (a1.max(omega - b2), b1.min(omega - a2));
    if lo >= hi {
        return 0.0;
    }

    // Wherever either factor stops being smooth, the quadrature starts a fresh panel
    // instead of discovering the trouble by subdividing
    let mut breaks: Vec<f64> = non_smooth(c1).collect();
    breaks.extend(non_smooth(c2).map(|mu| omega - mu));
    breaks.retain(|&p| p > lo && p < hi);
    breaks.sort_unstable_by(f64::total_cmp);
    breaks.dedup();

    let quad = |f: &dyn Fn(f64) -> f64| {
        adaptive_integrate_with_breaks(f, lo, hi, &breaks, tol).map_or(0.0, |r| r.value)
    };
    // C_2 reflected about omega, and zero where C_2 does not reach
    let other = |nu: f64| {
        let mu = omega - nu;
        if mu < a2 || mu > b2 {
            return 0.0;
        }
        let v = value_at(c2, mu);
        // The singular points of C_2 carry no weight
        if v.is_finite() { v } else { 0.0 }
    };

    let mut total = quad(&|nu| c1.regular(nu) * other(nu));
    for sing in c1.singularities() {
        if sing.is_trivial() {
            continue;
        }
        let position = sing.position;
        if position < lo || position > hi {
            // S_p stays bounded over an overlap it does not reach into
            total += quad(&|nu| sing.value(nu) * other(nu));
            continue;
        }
        let anchor = other(position);
        total += quad(&|nu| {
            if nu == position {
                0.0
            } else {
                sing.value(nu) * (other(nu) - anchor)
            }
        });
        total += anchor * sing.integral(lo, hi);
    }
    total
}

/// Convolution of the continuous parts of two spectral functions, evaluated by
/// quadrature and carrying a singular structure settled from outside.
pub struct Convolution<'a> {
    a: &'a SpectralFunction,
    b: &'a SpectralFunction,
    support: (f64, f64),
    singularities: Box<[Singularity]>,
    breakpoints: Box<[f64]>,
    tol: f64,
}

impl<'a> Convolution<'a> {
    pub fn new(
        a: &'a SpectralFunction,
        b: &'a SpectralFunction,
        singularities: Vec<Singularity>,
        breakpoints: Vec<f64>,
        tol: f64,
    ) -> Convolution<'a> {
        let (a_lo, a_hi) = continuous_hull(a);
        let (b_lo, b_hi) = continuous_hull(b);
        let support = (a_lo + b_lo, a_hi + b_hi);
        assert!(
            support.0.is_finite() && support.1.is_finite(),
            "convolution is restricted to bounded supports"
        );
        Convolution {
            a,
            b,
            support,
            singularities: singularities.into_boxed_slice(),
            breakpoints: breakpoints.into_boxed_slice(),
            tol,
        }
    }

    /// $(C_A \ast C_B)(\omega)$, summed over every pair of contributions.
    fn value(&self, omega: f64) -> f64 {
        let mut total = 0.0;
        for (c1, w1) in &self.a.continuous {
            for (c2, w2) in &self.b.continuous {
                total += w1 * w2 * pair(c1.as_ref(), c2.as_ref(), omega, self.tol);
            }
        }
        total
    }
}

/// Smallest segment holding the support of every continuous contribution.
fn continuous_hull(sf: &SpectralFunction) -> (f64, f64) {
    sf.continuous
        .iter()
        .map(|(csf, _)| csf.support())
        .reduce(|hull, sup| (hull.0.min(sup.0), hull.1.max(sup.1)))
        .expect("a convolved spectral function must have a continuous part")
}

impl ContinuousSF for Convolution<'_> {
    fn support(&self) -> (f64, f64) {
        self.support
    }
    fn regular(&self, omega: f64) -> f64 {
        let singular: f64 = self.singularities.iter().map(|s| s.value(omega)).sum();
        self.value(omega) - singular
    }
    fn singularities(&self) -> &[Singularity] {
        &self.singularities
    }
    fn breakpoints(&self) -> &[f64] {
        &self.breakpoints
    }
    fn shifted(&self, _by: f64) -> Box<dyn ContinuousSF> {
        unimplemented!("a convolution is interpolated where it is built, never displaced")
    }
}
