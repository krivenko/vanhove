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

use crate::beta::incomplete_beta_derivatives;
use crate::segment::Segment;
use crate::singularity::{LocalTerm, Singularity};
use crate::util::{bilby_integrate_or_0, binomials};
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

/// $\int_{s_1}^{s_2} s^{r_1}\ln^{m_1}\\!s\;(\Delta+s)^{r_2}\ln^{m_2}\\!(\Delta+s)\\,ds$,
/// for $\Delta > 0$ and $0 \le s_1 \le s_2$.
///
/// This is an outer stretch, the one running away from the pair of singular points
/// rather than between them. Scaling by $s = \Delta u$ and then $u = t/(1-t)$ turns it
/// into an incomplete Beta with
/// $$
///     a = r_1 + 1 > 0, \qquad b = -\rho, \qquad z = \frac{s}{\Delta+s},
/// $$
/// where $\rho = r_1 + r_2 + 1$. The logarithms come out of the same substitutions:
/// $\ln s = \ln\Delta + \ln t - \ln(1-t)$ and $\ln(\Delta+s) = \ln\Delta - \ln(1-t)$, so
/// expanding both binomially leaves nothing but entries of the derivative table.
fn outer_stretch(r1: f64, m1: usize, r2: f64, m2: usize, delta: f64, s1: f64, s2: f64) -> f64 {
    let rho = r1 + r2 + 1.0;
    let (a, b) = (r1 + 1.0, -rho);
    let z = |s: f64| s / (delta + s);
    let degree = m1 + m2;
    let upper = incomplete_beta_derivatives(a, b, z(s2), m1, degree);
    let lower = incomplete_beta_derivatives(a, b, z(s1), m1, degree);

    let ln_delta = delta.ln();
    let (rows1, rows2) = (binomials(m1), binomials(m2));
    let mut total = 0.0;
    for (i, &outer1) in rows1.iter().enumerate() {
        let inner_row = binomials(i);
        for (l, &outer2) in rows2.iter().enumerate() {
            let mut inner = 0.0;
            for (v, &weight) in inner_row.iter().enumerate() {
                let sign = if (v + l).is_multiple_of(2) { 1.0 } else { -1.0 };
                inner += weight * sign * (upper[i - v][v + l] - lower[i - v][v + l]);
            }
            total += outer1 * outer2 * ln_delta.powi((degree - i - l) as i32) * inner;
        }
    }
    delta.powf(rho) * total
}

/// $\int_{x_1}^{x_2} x^{r_1}\ln^{m_1}\\!x\;(\Delta-x)^{r_2}\ln^{m_2}\\!(\Delta-x)\\,dx$,
/// for $\Delta > 0$ and $0 \le x_1 \le x_2 \le \Delta$.
///
/// The stretch between the two singular points. Here $x = \Delta t$ is the whole
/// substitution, and the two exponents stay apart: $a = r_1+1$ and $b = r_2+1$ are both
/// positive, so this is the one stretch whose Beta never meets a pole.
fn middle_stretch(r1: f64, m1: usize, r2: f64, m2: usize, delta: f64, x1: f64, x2: f64) -> f64 {
    let rho = r1 + r2 + 1.0;
    let (a, b) = (r1 + 1.0, r2 + 1.0);
    let degree = m1 + m2;
    let upper = incomplete_beta_derivatives(a, b, x2 / delta, m1, m2);
    let lower = incomplete_beta_derivatives(a, b, x1 / delta, m1, m2);

    let ln_delta = delta.ln();
    let (rows1, rows2) = (binomials(m1), binomials(m2));
    let mut total = 0.0;
    for (i, &outer1) in rows1.iter().enumerate() {
        for (l, &outer2) in rows2.iter().enumerate() {
            let step = upper[i][l] - lower[i][l];
            total += outer1 * outer2 * ln_delta.powi((degree - i - l) as i32) * step;
        }
    }
    delta.powf(rho) * total
}

/// What $\int S_p S_q$ comes to where the two singular points sit at the same $\nu$.
///
/// There is no stretch between them and no closed form to reach for. The integrand goes
/// as $|x|^{r_1+r_2}$, so the integral diverges wherever a pair of terms brings that to
/// $-1$ or below, and the sign is the one the strongest such pair carries. Anything
/// weaker stays finite and is left to the quadrature.
///
/// The second factor is met at $-x$, so a side of it pairs with the opposite side of the
/// first: above the point it is the first from above against the second from below.
fn coincident(first: &Part, second: &Part, segment: Segment, point: f64) -> Option<f64> {
    let (reaches_above, reaches_below) = (segment.max() > point, segment.min() < point);
    let mut strongest = f64::INFINITY;
    let mut weight = 0.0;
    for t1 in &first.sing.local_form() {
        for t2 in &second.sing.local_form() {
            let order = t1.exponent + t2.exponent;
            if order > -1.0 || order > strongest {
                continue;
            }
            let mut c = 0.0;
            if reaches_above {
                c += t1.c_above * t2.c_below;
            }
            if reaches_below {
                c += t1.c_below * t2.c_above;
            }
            if c == 0.0 {
                continue;
            }
            if order < strongest {
                strongest = order;
                weight = c;
            } else {
                weight += c;
            }
        }
    }
    (weight != 0.0).then(|| f64::INFINITY * weight.signum())
}

/// $\int S_p(\nu) S_q(\omega-\nu) d\nu$ over `segment`, without a quadrature.
///
/// Measured from the first singular point, $x = \nu - \Omega_p$, the second factor is
/// met at $\Delta - x$ where $\Delta$ is the distance between the two points along the
/// integration axis. The integrand changes form where either factor turns, at $x = 0$
/// and $x = \Delta$, so the range is cut there into at most three stretches and each is
/// a Beta function.
///
/// A negative $\Delta$ is the same problem reflected: $x \mapsto -x$ carries it to a
/// positive one with both factors' sides exchanged.
///
/// [`None`] where the two points coincide, which leaves no stretch between them and an
/// integral that need not converge.
fn singular_singular_closed(first: &Part, second: &Part, segment: Segment) -> Option<f64> {
    let (p, q) = (first.at(), second.at());
    let raw = q - p;
    if !raw.is_finite() {
        return None;
    }
    if raw == 0.0 {
        return coincident(first, second, segment, p);
    }

    // Everything measured from the first singular point, then turned to face right
    let (x_lo, x_hi) = (segment.min() - p, segment.max() - p);
    let (delta, x_lo, x_hi, flipped) = if raw < 0.0 {
        (-raw, -x_hi, -x_lo, true)
    } else {
        (raw, x_lo, x_hi, false)
    };

    let sided = |t: &LocalTerm| {
        if flipped {
            (t.c_above, t.c_below)
        } else {
            (t.c_below, t.c_above)
        }
    };

    let mut total = 0.0;
    for t1 in &first.sing.local_form() {
        let (m1, r1) = (usize::from(t1.log_power), t1.exponent);
        let (below1, above1) = sided(t1);
        for t2 in &second.sing.local_form() {
            let (m2, r2) = (usize::from(t2.log_power), t2.exponent);
            let (below2, above2) = sided(t2);

            // x below zero: the first factor is met from below, the second from above
            let cut = x_hi.min(0.0);
            if x_lo < cut && below1 != 0.0 && above2 != 0.0 {
                total += below1 * above2 * outer_stretch(r1, m1, r2, m2, delta, -cut, -x_lo);
            }
            // between the two points, both factors met from above
            let (lo, hi) = (x_lo.max(0.0), x_hi.min(delta));
            if lo < hi && above1 != 0.0 && above2 != 0.0 {
                total += above1 * above2 * middle_stretch(r1, m1, r2, m2, delta, lo, hi);
            }
            // past the second point: the first from above, the second from below
            let cut = x_lo.max(delta);
            if cut < x_hi && above1 != 0.0 && below2 != 0.0 {
                total += above1
                    * below2
                    * outer_stretch(r2, m2, r1, m1, delta, cut - delta, x_hi - delta);
            }
        }
    }
    Some(total)
}

/// The same by quadrature, each divergence given a stretch of its own and subtracted
/// there against the other factor.
///
/// Kept as the reference the closed form is checked against, and as what answers where
/// the two points coincide: there is then no stretch between them to speak of, no finite
/// anchor to take, and an integrand going as $|\nu-\Omega|^{r_1+r_2}$ whose integral need
/// not converge at all.
fn singular_singular_by_quadrature(first: &Part, second: &Part, segment: Segment, tol: f64) -> f64 {
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
            total += singular_singular_closed(first, second, overlap)
                .unwrap_or_else(|| singular_singular_by_quadrature(first, second, overlap, tol));
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

    /// The closed form and the quadrature it replaces, on every pair of singular parts
    /// a chain and a square lattice have between them.
    #[test]
    fn the_closed_form_agrees_with_quadrature() {
        let (a, b) = (chain(0.0, 1.0), square(0.0, 1.0));
        let (c1, _) = &a.continuous[0];
        let (c2, _) = &b.continuous[0];
        let (s1, s2) = (c1.support(), c2.support());
        let mut checked = 0;
        for omega in [-5.5f64, -4.5, -3.3, -1.7, -0.5, 0.7, 2.2, 3.5, 4.4, 5.5] {
            let Some(overlap) = s1.intersection(s2.mirrored(omega)) else {
                continue;
            };
            if overlap.is_degenerate() {
                continue;
            }
            for sing1 in c1.singularities() {
                for sing2 in c2.singularities() {
                    let first = Part {
                        sing: sing1,
                        support: s1,
                        reflected: None,
                    };
                    let second = Part {
                        sing: sing2,
                        support: s2,
                        reflected: Some(omega),
                    };
                    let Some(closed) = singular_singular_closed(&first, &second, overlap) else {
                        continue;
                    };
                    let quadrature =
                        singular_singular_by_quadrature(&first, &second, overlap, 1e-13);
                    assert_relative_eq!(closed, quadrature, max_relative = 1e-8, epsilon = 1e-12);
                    checked += 1;
                }
            }
        }
        assert!(checked >= 20, "only {checked} pairs were reached");
    }

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

    /// Where two singular points meet the value is not a number but an infinity, and
    /// saying so is the point: two inverse square roots colliding is the van Hove peak
    /// of the square lattice.
    #[test]
    fn colliding_points_diverge() {
        let (a, b) = (chain(0.0, 1.0), chain(0.0, 1.0));
        assert_eq!(value(&a, &b, 0.0, 1e-12), f64::INFINITY);
        assert_eq!(square(0.0, 1.0).continuous_at(0.0), f64::INFINITY);

        // Either side of it the value is finite and right
        for w in [-1e-6f64, 1e-6] {
            assert_relative_eq!(
                value(&a, &b, w, 1e-12),
                square(0.0, 1.0).continuous_at(w),
                max_relative = 1e-9
            );
        }

        // Two boxes have no singular parts to collide, so their band centre is finite
        let flat_pair = flat(0.0, 1.5, 0.0);
        assert!(value(&flat_pair, &flat_pair, 0.0, 1e-12).is_finite());
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
