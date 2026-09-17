//! Convolution of the continuous parts of two spectral functions.
//!
//! $C_A = R_A + \sum_p S_p$ and $C_B = R_B + \sum_q S_q$, so the convolution is four
//! kinds of term:
//! $$
//!     C_A \ast C_B = R_A \ast R_B + \sum_p S_p \ast R_B + \sum_q R_A \ast S_q
//!                  + \sum_{p,q} S_p \ast S_q,
//! $$
//! and they are taken one at a time rather than by handing the whole product to a
//! quadrature. Only the last can put two divergences at one and the same $\nu$; the
//! first needs no care at all, and the middle two need one subtraction each against a
//! factor that is smooth by construction.

use crate::segment::Segment;
use crate::singularity::Singularity;
use crate::util::bilby_integrate_or_0;
use crate::{ContinuousSF, SpectralFunction};

/// Regular part of a contribution, zero where it does not reach.
fn regular(csf: &dyn ContinuousSF, omega: f64) -> f64 {
    if !csf.support().contains(omega) {
        return 0.0;
    }
    let value = csf.regular(omega);
    if value.is_finite() { value } else { 0.0 }
}

/// One singular part as it enters the integral over $\nu$.
///
/// The second factor of a convolution is met at $\omega - \nu$, so which variable a
/// singularity is written in has to travel with it.
struct Part<'a> {
    sing: &'a Singularity,
    /// Support of the contribution it belongs to, in that contribution's own variable.
    support: Segment,
    /// $\omega$ where the argument is $\omega - \nu$, absent where it is $\nu$ itself.
    reflected: Option<f64>,
}

impl Part<'_> {
    /// Argument the singularity is asked for when the integration variable is `nu`.
    fn argument(&self, nu: f64) -> f64 {
        match self.reflected {
            Some(omega) => omega - nu,
            None => nu,
        }
    }

    /// Where the singular point sits in $\nu$.
    fn at(&self) -> f64 {
        match self.reflected {
            Some(omega) => omega - self.sing.position,
            None => self.sing.position,
        }
    }

    /// Value at `nu`, zero outside the support and at the singular point itself, which
    /// carries no weight under the integral.
    fn value(&self, nu: f64) -> f64 {
        let argument = self.argument(nu);
        if !self.support.contains(argument) {
            return 0.0;
        }
        let value = self.sing.value(argument);
        if value.is_finite() { value } else { 0.0 }
    }

    /// $\int S$ over a stretch of the $\nu$ axis, in closed form.
    fn integral(&self, segment: Segment) -> f64 {
        let own = match self.reflected {
            Some(omega) => segment.mirrored(omega),
            None => segment,
        };
        match own.intersection(self.support) {
            Some(within) => self.sing.integral(within),
            None => 0.0,
        }
    }
}

/// $\int S(\nu) g(\nu) d\nu$ over `segment`, with the divergence of `S` subtracted.
///
/// `g` is anchored at the singular point, leaving $S(\nu)[g(\nu) - g(\Omega)]$ bounded,
/// and the constant is restored through $\int S$ in closed form. The subtraction is
/// worth nothing unless `g` is continuous at the point, which is what confines this to
/// a stretch holding one divergence and no more.
fn subtracted<G: Fn(f64) -> f64>(part: &Part, g: G, segment: Segment, tol: f64) -> f64 {
    let point = part.at();
    if !segment.contains(point) {
        return bilby_integrate_or_0(|nu| part.value(nu) * g(nu), segment, tol);
    }
    let anchor = g(point);
    if !anchor.is_finite() {
        return bilby_integrate_or_0(|nu| part.value(nu) * g(nu), segment, tol);
    }
    let bounded = bilby_integrate_or_0(
        |nu| {
            if nu == point {
                0.0
            } else {
                part.value(nu) * (g(nu) - anchor)
            }
        },
        segment,
        tol,
    );
    bounded + anchor * part.integral(segment)
}

/// $\int S_p(\nu) S_q(\omega-\nu) d\nu$ over `segment`.
///
/// Each divergence is given a stretch of its own, the two meeting halfway between them,
/// and is subtracted there against the other factor. Where the two points coincide there
/// is no such split to make and no finite anchor to take: the integrand then goes as
/// $|\nu - \Omega|^{r_1 + r_2}$ and the integral itself need not converge.
fn singular_singular(first: &Part, second: &Part, segment: Segment, tol: f64) -> f64 {
    let (p, q) = (first.at(), second.at());
    let (inside_p, inside_q) = (segment.contains(p), segment.contains(q));

    if inside_p && inside_q && p != q {
        let (lower, upper) = segment.split_at(0.5 * (p + q));
        let (near_p, near_q) = if p < q {
            (lower, upper)
        } else {
            (upper, lower)
        };
        return subtracted(first, |nu| second.value(nu), near_p, tol)
            + subtracted(second, |nu| first.value(nu), near_q, tol);
    }
    if inside_p && !inside_q {
        return subtracted(first, |nu| second.value(nu), segment, tol);
    }
    if inside_q && !inside_p {
        return subtracted(second, |nu| first.value(nu), segment, tol);
    }
    bilby_integrate_or_0(|nu| first.value(nu) * second.value(nu), segment, tol)
}

/// $\int C_1(\nu) C_2(\omega-\nu) d\nu$ over the overlap of the two supports.
fn pair(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF, omega: f64, tol: f64) -> f64 {
    let (s1, s2) = (c1.support(), c2.support());
    let Some(overlap) = s1.intersection(s2.mirrored(omega)) else {
        return 0.0;
    };
    if overlap.is_degenerate() {
        return 0.0;
    }

    let parts1: Vec<Part> = c1
        .singularities()
        .iter()
        .map(|sing| Part {
            sing,
            support: s1,
            reflected: None,
        })
        .collect();
    let parts2: Vec<Part> = c2
        .singularities()
        .iter()
        .map(|sing| Part {
            sing,
            support: s2,
            reflected: Some(omega),
        })
        .collect();

    // R_A ⊛ R_B, both smooth across the overlap with nothing to subtract
    let mut total =
        bilby_integrate_or_0(|nu| regular(c1, nu) * regular(c2, omega - nu), overlap, tol);

    // S_p ⊛ R_B and R_A ⊛ S_q, one divergence each against a factor that has none
    for part in &parts1 {
        total += subtracted(part, |nu| regular(c2, omega - nu), overlap, tol);
    }
    for part in &parts2 {
        total += subtracted(part, |nu| regular(c1, nu), overlap, tol);
    }

    // S_p ⊛ S_q, the only term where two divergences can meet
    for first in &parts1 {
        for second in &parts2 {
            total += singular_singular(first, second, overlap, tol);
        }
    }
    total
}

/// $(C_A \ast C_B)(\omega)$, summed over every pair of contributions.
pub fn value(a: &SpectralFunction, b: &SpectralFunction, omega: f64, tol: f64) -> f64 {
    let mut total = 0.0;
    for (c1, w1) in &a.continuous {
        for (c2, w2) in &b.continuous {
            total += w1 * w2 * pair(c1.as_ref(), c2.as_ref(), omega, tol);
        }
    }
    total
}

/// Support of $C_A \ast C_B$, which reaches as far as the two supports added together.
pub fn support(a: &SpectralFunction, b: &SpectralFunction) -> Segment {
    let ends = |sf: &SpectralFunction| {
        Segment::hull(sf.continuous.iter().map(|(c, _)| c.support()))
            .expect("a convolution needs a continuous contribution on either side")
    };
    let (x, y) = (ends(a), ends(b));
    assert!(
        x.is_bounded() && y.is_bounded(),
        "convolution of continuous parts is restricted to bounded supports"
    );
    Segment::new(x.min() + y.min(), x.max() + y.max())
}

/// Frequencies at which a contribution stops being smooth: its singular points, and
/// the ends of its support, where it stops contributing at all.
fn features(csf: &dyn ContinuousSF) -> Vec<f64> {
    let support = csf.support();
    let mut out: Vec<f64> = csf.singularities().iter().map(|s| s.position).collect();
    out.push(support.min());
    out.push(support.max());
    out.sort_unstable_by(f64::total_cmp);
    out.dedup();
    out
}

/// Singular structure of $C_A \ast C_B$.
///
/// $\int C_1(\nu) C_2(\omega-\nu) d\nu$ stops being smooth in $\omega$ where both
/// factors stop being smooth at one and the same $\nu$, that is where $\nu = \Omega_1$
/// and $\omega - \nu = \Omega_2$ hold together. The positions are therefore the
/// pairwise sums.
///
/// Every position carries an empty singularity for now, which says where the
/// interpolation must start a fresh panel and nothing more.
pub fn singularities(a: &SpectralFunction, b: &SpectralFunction) -> Vec<Singularity> {
    let reach = support(a, b);
    let mut positions: Vec<f64> = Vec::new();
    for (c1, _) in &a.continuous {
        for (c2, _) in &b.continuous {
            for p1 in features(c1.as_ref()) {
                for p2 in features(c2.as_ref()) {
                    positions.push(p1 + p2);
                }
            }
        }
    }
    positions.retain(|p| reach.strictly_contains(*p));
    positions.sort_unstable_by(f64::total_cmp);
    positions.dedup();
    positions
        .into_iter()
        .map(|p| Singularity::new(p, 1.0, Vec::new()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::*;
    use approx::assert_relative_eq;

    /// Two boxes convolve into a triangle, which is worth knowing exactly.
    #[test]
    fn two_boxes_make_a_triangle() {
        let (d, tol) = (1.5f64, 1e-12);
        let (a, b) = (flat(0.0, d, 0.0), flat(0.0, d, 0.0));
        let triangle = |w: f64| (2.0 * d - w.abs()).max(0.0) / (4.0 * d * d);
        for i in 0..=40 {
            let w = -2.0 * d + 4.0 * d * (i as f64) / 40.0;
            assert_relative_eq!(value(&a, &b, w, tol), triangle(w), epsilon = 1e-12);
        }
    }

    /// A chain against itself is the square lattice, singular parts and all.
    #[test]
    fn two_chains_make_a_square_lattice() {
        let t = 1.0f64;
        let (a, b) = (chain(0.0, t), chain(0.0, t));
        let reference = square(0.0, t);
        for w in [-3.9f64, -2.5, -1.0, -0.3, 0.3, 1.0, 2.5, 3.9] {
            assert_relative_eq!(
                value(&a, &b, w, 1e-12),
                reference.continuous_at(w),
                max_relative = 1e-12
            );
        }
    }

    /// The support adds, and the convolution vanishes beyond it.
    #[test]
    fn support_adds() {
        let (a, b) = (chain(0.0, 1.0), flat(0.5, 1.0, 0.0));
        let reach = support(&a, &b);
        assert_relative_eq!(reach.min(), -2.5);
        assert_relative_eq!(reach.max(), 3.5);
        assert_eq!(value(&a, &b, 4.0, 1e-10), 0.0);
        assert_eq!(value(&a, &b, -3.0, 1e-10), 0.0);
    }
}
