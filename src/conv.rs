//! Convolution of the continuous parts of two spectral functions.

use bilby::adaptive_integrate_with_breaks;

use crate::singularity::Singularity;
use crate::{ContinuousSF, SpectralFunction, non_smooth};

/// Value of a single continuous contribution, zero where it does not reach and where
/// it has no value at all.
fn sampled(csf: &dyn ContinuousSF, omega: f64) -> f64 {
    let (omega_min, omega_max) = csf.support();
    if omega < omega_min || omega > omega_max {
        return 0.0;
    }
    let singular: f64 = csf.singularities().iter().map(|s| s.value(omega)).sum();
    let value = csf.regular(omega) + singular;
    // A singular point carries no weight under the integral
    if value.is_finite() { value } else { 0.0 }
}

/// The same, with the singularity of index `skip` left out, which is finite at the
/// point that one diverges at.
fn sampled_without(csf: &dyn ContinuousSF, omega: f64, skip: usize) -> f64 {
    let (omega_min, omega_max) = csf.support();
    if omega < omega_min || omega > omega_max {
        return 0.0;
    }
    let singular: f64 = csf
        .singularities()
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != skip)
        .map(|(_, s)| s.value(omega))
        .sum();
    let value = csf.regular(omega) + singular;
    if value.is_finite() { value } else { 0.0 }
}

/// $\int C_1(\nu) C_2(\omega-\nu) d\nu$ over the overlap of the two supports.
///
/// The overlap is divided so that each singular point of either factor falls in a
/// stretch of its own, and over that stretch the singularity is subtracted the way
/// [`SpectralFunction::integrate()`] subtracts one: the other factor is anchored at the
/// singular point, leaving a bounded integrand, and the constant is restored through
/// $\int S_p$ in closed form. Handing a bare product to the quadrature instead leaves
/// an inverse square root sitting on an end of a panel, where Gauss-Kronrod converges
/// too slowly to notice it is not converging.
fn pair(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF, omega: f64, tol: f64) -> f64 {
    let (a1, b1) = c1.support();
    let (a2, b2) = c2.support();
    let (lo, hi) = (a1.max(omega - b2), b1.min(omega - a2));
    if lo >= hi {
        return 0.0;
    }

    // Wherever either factor stops being smooth without being singular, the quadrature
    // starts a fresh panel instead of discovering the trouble by subdividing
    let mut breaks: Vec<f64> = c1.breakpoints().to_vec();
    breaks.extend(c2.breakpoints().iter().map(|mu| omega - mu));
    breaks.retain(|&p| p > lo && p < hi);
    breaks.sort_unstable_by(f64::total_cmp);
    breaks.dedup();

    // Singular points of either factor, in the integration variable
    let mut points: Vec<(f64, usize, bool)> = Vec::new();
    for (index, sing) in c1.singularities().iter().enumerate() {
        points.push((sing.position(), index, true));
    }
    for (index, sing) in c2.singularities().iter().enumerate() {
        points.push((omega - sing.position(), index, false));
    }
    points.retain(|(nu, _, _)| *nu >= lo && *nu <= hi);
    points.sort_unstable_by(|x, y| x.0.total_cmp(&y.0));

    let quad = |u: f64, v: f64, f: &dyn Fn(f64) -> f64| {
        let inner: Vec<f64> = breaks
            .iter()
            .copied()
            .filter(|p| *p > u && *p < v)
            .collect();
        adaptive_integrate_with_breaks(f, u, v, &inner, tol).map_or(0.0, |r| r.value)
    };

    if points.is_empty() {
        return quad(lo, hi, &|nu| sampled(c1, nu) * sampled(c2, omega - nu));
    }

    // One stretch per singular point, meeting halfway between neighbours
    let mut total = 0.0;
    for (i, &(nu_p, index, first)) in points.iter().enumerate() {
        let u = if i == 0 {
            lo
        } else {
            0.5 * (points[i - 1].0 + nu_p)
        };
        let v = if i + 1 == points.len() {
            hi
        } else {
            0.5 * (nu_p + points[i + 1].0)
        };
        if u >= v {
            continue;
        }

        let (sing, anchor) = if first {
            (&c1.singularities()[index], sampled(c2, omega - nu_p))
        } else {
            (&c2.singularities()[index], sampled(c1, nu_p))
        };
        // Two singular points meeting leave nothing finite to anchor against
        if !anchor.is_finite() || sing.is_trivial() {
            total += quad(u, v, &|nu| sampled(c1, nu) * sampled(c2, omega - nu));
            continue;
        }

        if first {
            let rest = |nu: f64| sampled_without(c1, nu, index);
            total += quad(u, v, &|nu| rest(nu) * sampled(c2, omega - nu));
            total += quad(u, v, &|nu| {
                if nu == nu_p {
                    0.0
                } else {
                    sing.value(nu) * (sampled(c2, omega - nu) - anchor)
                }
            });
            total += anchor * sing.integral(u, v);
        } else {
            let rest = |nu: f64| sampled_without(c2, omega - nu, index);
            total += quad(u, v, &|nu| sampled(c1, nu) * rest(nu));
            total += quad(u, v, &|nu| {
                if nu == nu_p {
                    0.0
                } else {
                    sing.value(omega - nu) * (sampled(c1, nu) - anchor)
                }
            });
            // dν = -dμ turns the stretch around
            total += anchor * sing.integral(omega - v, omega - u);
        }
    }
    total
}

/// Frequencies where the continuous parts of `a` and `b` convolve into something that
/// is not smooth.
///
/// A convolution departs from smoothness where the frequencies at which either factor
/// does meet, so the set is the pairwise sums of those of the two. The ends of a
/// support count among them: that is where a factor stops contributing at all.
pub fn breakpoints(a: &SpectralFunction, b: &SpectralFunction) -> Vec<f64> {
    let points = |sf: &SpectralFunction| {
        let mut out = Vec::new();
        for (csf, _) in &sf.continuous {
            let (omega_min, omega_max) = csf.support();
            out.push(omega_min);
            out.push(omega_max);
            out.extend(non_smooth(csf.as_ref()));
        }
        out
    };
    let (points_a, points_b) = (points(a), points(b));

    let mut sums = Vec::with_capacity(points_a.len() * points_b.len());
    for x in &points_a {
        for y in &points_b {
            sums.push(x + y);
        }
    }
    sums.sort_unstable_by(f64::total_cmp);
    sums.dedup();
    sums
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
