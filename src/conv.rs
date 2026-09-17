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

use crate::ContinuousSF;
use crate::beta::{beta_derivatives, incomplete_beta_derivatives, outer_stretch_logs};
use crate::segment::Segment;
use crate::singularity::{AsymptTerm, LocalTerm, Singularity};
use crate::util::{bilby_integrate_or_0, binomials, is_natural};

/// Exponent above which a derived term is smooth enough to leave to the interpolation.
///
/// $|\Delta|^2$ has two derivatives, which is more than a panel boundary asks for.
const MAX_DERIVED_EXPONENT: f64 = 2.0;

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
            Some(omega) => omega - self.sing.position(),
            None => self.sing.position(),
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
    ///
    /// The stretch has to be carried into the singularity's own variable, and for a
    /// reflected one its ends are measured off $\Omega_q$ rather than mirrored. The
    /// singular point sits at $\omega - \Omega_q$ along the $\nu$ axis, so mirroring that
    /// back returns $\Omega_q$ only to within rounding, and an end may land the wrong
    /// side of the very point it is there to enclose.
    fn integral(&self, segment: Segment) -> f64 {
        let own = match self.reflected {
            Some(_) => {
                let point = self.at();
                Segment::new(
                    self.sing.position() - (segment.max() - point),
                    self.sing.position() + (point - segment.min()),
                )
            }
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
pub fn pair_value(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF, omega: f64, tol: f64) -> f64 {
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

/// Support of $C_1 \ast C_2$, which reaches as far as the two supports added together.
pub fn pair_support(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF) -> Segment {
    let (x, y) = (c1.support(), c2.support());
    assert!(
        x.is_bounded() && y.is_bounded(),
        "convolution of continuous parts is restricted to bounded supports"
    );
    Segment::new(x.min() + y.min(), x.max() + y.max())
}

/// What the three stretches contribute where $r$ is a non-negative integer $n$.
///
/// $\Delta^n$ is analytic, and the outer stretches no longer converge: their integrands
/// go as $u^{n-1}$ at large $u$, so one term of the expansion of $(1+1/u)^{r_2}$
/// integrates to a logarithm rather than a power. That logarithm is the whole of what an
/// outer stretch leaves behind, the rest of it being analytic in $\Delta$ and carried by
/// the closed form.
///
/// So a convolution gains a logarithm exactly where the generic formula loses one to a
/// pole of $\Gamma(-r)$. Two locally constant factors are the case $r_1 = r_2 = 0$,
/// $n = 1$, where the binomials ask for more than they have and vanish: that is why two
/// smooth constants convolve into no logarithm at all, and why a van Hove saddle in
/// three dimensions is a square root rather than a peak.
fn degenerate_stretches(t1: LocalTerm, t2: LocalTerm, n: usize) -> [Vec<f64>; 3] {
    let (r1, r2) = (t1.exponent, t2.exponent);
    let (m1, m2) = (usize::from(t1.log_power), usize::from(t2.log_power));
    let degree = m1 + m2 + 1;
    // Δ^n carrying no logarithm is analytic, so it is no part of the singular structure.
    // It is not dropped but handed over: the closed form computes the convolution whole,
    // and what the asymptotics does not claim stays in the regular part exactly.
    let pad = |mut v: Vec<f64>| {
        v.resize(degree + 1, 0.0);
        v[0] = 0.0;
        v
    };
    [
        pad(middle_coefficients(r1, m1, r2, m2)),
        pad(outer_stretch_logs(r1, m1, r2, m2, n)),
        pad(outer_stretch_logs(r2, m2, r1, m1, n)),
    ]
}

/// What the middle stretch contributes to the coefficient of
/// $|\Delta|^r\ln^k|\Delta|$, indexed by $k$.
///
/// $\alpha$ and $\beta$ are both positive whatever $r$ comes to, so this is the one
/// stretch that never meets a pole, and it is used at a whole-number $r$ unchanged.
fn middle_coefficients(r1: f64, m1: usize, r2: f64, m2: usize) -> Vec<f64> {
    let degree = m1 + m2;
    let (c1, c2) = (binomials(m1), binomials(m2));
    let table = beta_derivatives(r1 + 1.0, r2 + 1.0, m1, m2);
    let mut mid = vec![0.0; degree + 1];
    for a in 0..=m1 {
        for b in 0..=m2 {
            mid[degree - a - b] += c1[a] * c2[b] * table[a][b];
        }
    }
    mid
}

/// What the three stretches contribute to the coefficient of $|\Delta|^r\ln^k|\Delta|$,
/// in the order middle, lower, upper, each indexed by $k$.
///
/// Substituting the length of the stretch out of the integral turns every logarithm into
/// $\ln|\Delta| + \ln(\text{something of order one})$, and expanding those binomials
/// leaves $|\Delta|^r$ times a polynomial in $\ln|\Delta|$ of degree $m_1 + m_2$. Its
/// coefficients are the integrals
/// $$
///     \int_0^1 t^{r_1}(1-t)^{r_2}\ln^a t \ln^b(1-t)\\,dt
///         = \partial_\alpha^a \partial_\beta^b B(\alpha, \beta)
/// $$
/// over the middle stretch, at $\alpha = r_1+1$, $\beta = r_2+1$, and
/// $$
///     \int_0^\infty u^{r_1}(1+u)^{r_2}\ln^a u \ln^b(1+u)\\,du
///         = (\partial_\alpha - \partial_\gamma)^a(-\partial_\gamma)^b B(\alpha, \gamma)
/// $$
/// over an outer one, at $\gamma = -r$: there $\ln(1+u)$ is $-\partial_\gamma$ and
/// $\ln u$ the difference of the two, the integrand carrying $\alpha$ in both factors.
fn stretches(t1: LocalTerm, t2: LocalTerm) -> [Vec<f64>; 3] {
    let (r1, r2) = (t1.exponent, t2.exponent);
    let (m1, m2) = (usize::from(t1.log_power), usize::from(t2.log_power));
    let r = r1 + r2 + 1.0;
    let degree = m1 + m2;
    let (c1, c2) = (binomials(m1), binomials(m2));

    let mid = middle_coefficients(r1, m1, r2, m2);

    // An outer stretch runs from its cut out to the end of the overlap, a length that
    // enters only analytically in Δ; what survives is the finite part below
    let outer = |ra: f64, rb: f64, ca: &[f64], cb: &[f64], ma: usize, mb: usize| {
        let table = beta_derivatives(ra + 1.0, -r, ma, ma + mb);
        let mut out = vec![0.0; degree + 1];
        for a in 0..=ma {
            let inner = binomials(a);
            for b in 0..=mb {
                // (∂_α - ∂_γ)^a (-∂_γ)^b = Σ_i C(a,i) (-1)^{i+b} ∂_α^{a-i} ∂_γ^{i+b}
                let value: f64 = (0..=a)
                    .map(|i| {
                        let sign = if (i + b).is_multiple_of(2) { 1.0 } else { -1.0 };
                        inner[i] * sign * table[a - i][i + b]
                    })
                    .sum();
                out[degree - a - b] += ca[a] * cb[b] * value;
            }
        }
        let _ = rb;
        out
    };

    [
        mid,
        outer(r1, r2, &c1, &c2, m1, m2),
        outer(r2, r1, &c2, &c1, m2, m1),
    ]
}

/// Terms of $T_1 \ast T_2$ at $\Omega_1 + \Omega_2$, empty where the pair is one this
/// does not reach.
///
/// The result spans $|\Delta|^r \ln^k|\Delta|$ for $k$ up to $m_1 + m_2$, each
/// coefficient the three stretches added together with the sides they draw on:
/// $$
///     C^+_k = c_1^+ c_2^+ X^{mid}_k + c_1^- c_2^+ X^{lo}_k + c_1^+ c_2^- X^{hi}_k,
/// $$
/// and $C^-_k$ the same with every $\pm$ flipped.
fn convolve_terms(t1: LocalTerm, t2: LocalTerm, weight: f64) -> Vec<AsymptTerm> {
    let r = t1.exponent + t2.exponent + 1.0;
    if r >= MAX_DERIVED_EXPONENT {
        return Vec::new();
    }
    let [mid, lo, hi] = if is_natural(r) {
        degenerate_stretches(t1, t2, r as usize)
    } else {
        stretches(t1, t2)
    };
    let (cb1, ca1, cb2, ca2) = (t1.c_below, t1.c_above, t2.c_below, t2.c_above);

    let mut terms = Vec::new();
    for k in (0..mid.len()).rev() {
        let combine = |x1: f64, x2: f64, x3: f64| {
            (
                weight * (cb1 * cb2 * x1 + ca1 * cb2 * x2 + cb1 * ca2 * x3),
                weight * (ca1 * ca2 * x1 + cb1 * ca2 * x2 + ca1 * cb2 * x3),
            )
        };
        let (c_below, c_above) = combine(mid[k], lo[k], hi[k]);
        if c_below == 0.0 && c_above == 0.0 {
            continue;
        }
        // A constant may differ between the sides here: where r is zero the term sits
        // under a logarithm that diverges, so there is no finite value at Ω_p for the
        // two sides to disagree on
        terms.push(AsymptTerm::sided_log(r, k as u8, c_below, c_above));
    }
    terms
}

/// Add `terms` to whatever has already been derived at `position`.
fn merge(derived: &mut Vec<(f64, Vec<AsymptTerm>)>, position: f64, terms: Vec<AsymptTerm>) {
    if terms.is_empty() {
        return;
    }
    match derived.iter_mut().find(|(p, _)| *p == position) {
        Some((_, existing)) => existing.extend(terms),
        None => derived.push((position, terms)),
    }
}

/// Frequencies at which a contribution stops being smooth: its singular points, and
/// the ends of its support, where it stops contributing at all.
fn features(csf: &dyn ContinuousSF) -> Vec<f64> {
    let support = csf.support();
    let mut out: Vec<f64> = csf.singularities().iter().map(|s| s.position()).collect();
    out.push(support.min());
    out.push(support.max());
    out.sort_unstable_by(f64::total_cmp);
    out.dedup();
    out
}

/// Singular structure of $C_1 \ast C_2$.
///
/// $\int C_1(\nu) C_2(\omega-\nu) d\nu$ stops being smooth in $\omega$ where both factors
/// stop being smooth at one and the same $\nu$, that is where $\nu = \Omega_1$ and
/// $\omega - \nu = \Omega_2$ hold together. The positions are therefore the pairwise sums,
/// and the ends of a support count: that is where a factor stops contributing at all.
///
/// Only singular parts carry terms. A support end where nothing diverges leaves a kink
/// going as $|\Delta|$ or milder, which is a polynomial either side of the panel
/// boundary it sits on and is fitted there exactly, so there is nothing to gain by
/// deriving it and a double count to be had by trying: two band edges of equal width
/// would each claim the whole of it. Such a position enters with no terms, saying where
/// the interpolation must start a fresh panel and nothing more.
/// Local form of `csf` about `position`: the terms of a singularity sitting there, and
/// the value everything else takes, kept apart.
///
/// Each is zeroed on the side the support does not reach. A singular part at an end of
/// its support has no side there, and reading a term as written would have a band edge
/// convolve as though the band ran on through it.
fn local_form_at(csf: &dyn ContinuousSF, position: f64) -> (Vec<LocalTerm>, LocalTerm) {
    let support = csf.support();
    let (reaches_below, reaches_above) = (position > support.min(), position < support.max());
    let sided = |exponent, log_power, c_below: f64, c_above: f64| LocalTerm {
        exponent,
        log_power,
        c_below: if reaches_below { c_below } else { 0.0 },
        c_above: if reaches_above { c_above } else { 0.0 },
    };

    let mut singular = Vec::new();
    let mut constant = csf.regular(position);
    for sing in csf.singularities() {
        if sing.position() == position {
            singular.extend(
                sing.local_form()
                    .into_iter()
                    .map(|t| sided(t.exponent, t.log_power, t.c_below, t.c_above)),
            );
        } else {
            constant += sing.value(position);
        }
    }
    (singular, sided(0.0, 0, constant, constant))
}

/// Singular structure of $C_1 \ast C_2$.
///
/// $\int C_1(\nu) C_2(\omega-\nu) d\nu$ stops being smooth in $\omega$ where both factors
/// stop being smooth at one and the same $\nu$, that is where $\nu = \Omega_1$ and
/// $\omega - \nu = \Omega_2$ hold together. The positions are therefore the pairwise sums,
/// and the ends of a support count: that is where a factor stops contributing at all.
///
/// Every pair of local terms is convolved but one: the two constants. That product is
/// $R \ast R$, which leaves a kink going as $|\Delta|$ or milder — a polynomial either
/// side of the panel boundary it sits on, and fitted there exactly. Deriving it would
/// gain nothing and cost the equal-width case, where two band edges of the same width
/// would each claim the whole of the same kink. Everything with a singular part on one
/// side or the other is derived.
pub fn pair_singularities(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF) -> Vec<Singularity> {
    let mut derived: Vec<(f64, Vec<AsymptTerm>)> = Vec::new();
    for p1 in features(c1) {
        let (singular1, constant1) = local_form_at(c1, p1);
        for p2 in features(c2) {
            let (singular2, constant2) = local_form_at(c2, p2);
            let mut terms = Vec::new();
            for &t1 in &singular1 {
                for &t2 in singular2.iter().chain(std::iter::once(&constant2)) {
                    terms.extend(convolve_terms(t1, t2, 1.0));
                }
            }
            for &t2 in &singular2 {
                terms.extend(convolve_terms(constant1, t2, 1.0));
            }
            merge(&mut derived, p1 + p2, terms);
        }
    }

    let mut positions: Vec<f64> = Vec::new();
    for p1 in features(c1) {
        for p2 in features(c2) {
            positions.push(p1 + p2);
        }
    }
    positions.sort_unstable_by(f64::total_cmp);
    positions.dedup();
    positions
        .into_iter()
        .map(|position| {
            let terms = match derived.iter_mut().find(|(p, _)| *p == position) {
                Some((_, terms)) => std::mem::take(terms),
                None => Vec::new(),
            };
            Singularity::new(position, 1.0, terms)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SpectralFunction;
    use crate::models::*;
    use approx::assert_relative_eq;

    /// The one continuous contribution a model carries.
    fn only(sf: &SpectralFunction) -> &dyn ContinuousSF {
        sf.continuous[0].0.as_ref()
    }

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

    /// The van Hove logarithm of the square lattice, which nothing in the library is
    /// told: it falls out of two inverse square roots meeting at $r = 0$.
    #[test]
    fn the_square_lattice_logarithm_is_derived() {
        let t = 1.0f64;
        let (a, b) = (chain(0.0, t), chain(0.0, t));
        let derived = pair_singularities(only(&a), only(&b));
        let centre = derived.iter().find(|s| s.position() == 0.0).unwrap();
        let slope = (centre.value(1e-8) - centre.value(1e-4)) / (1e-8f64.ln() - 1e-4f64.ln());
        assert_relative_eq!(
            slope,
            -1.0 / (2.0 * std::f64::consts::PI.powi(2) * t),
            max_relative = 1e-11
        );

        // The band edges meet at r = 0 too, where the term left over is a constant,
        // analytic and no part of the singular structure
        for position in [-4.0f64, 4.0] {
            let edge = derived.iter().find(|s| s.position() == position).unwrap();
            assert!(edge.is_trivial());
        }
    }

    /// What the derivation reaches, end to end: the interpolated result has to fit.
    #[test]
    fn the_regular_part_is_smooth_enough_to_fit() {
        use crate::interp::InterpolatedSF;
        let cases = [
            // A chain against itself pairs two genuine singularities at every point
            (
                "chain against chain",
                chain(0.0, 1.0),
                chain(0.0, 1.0),
                1e-12f64,
            ),
            // Two boxes have no singular parts at all, and their triangle is fitted
            // exactly because a kink at a panel boundary is a polynomial either side
            (
                "box against box",
                flat(0.0, 1.5, 0.0),
                flat(0.0, 1.5, 0.0),
                1e-13,
            ),
            // A band edge where one factor merely stops still pairs, its local form
            // being a constant, so the square root there is derived like any other
            (
                "chain against square",
                chain(0.0, 1.0),
                square(0.0, 1.0),
                1e-11,
            ),
            (
                "square against square",
                square(0.0, 1.0),
                square(0.0, 1.0),
                1e-11,
            ),
            (
                "triangular against itself",
                triangular(0.0, 1.0),
                triangular(0.0, 1.0),
                1e-11,
            ),
        ];
        for (name, a, b, want) in cases {
            let derived = pair_singularities(only(&a), only(&b));
            let tol = 1e-12;
            let regular = |omega: f64| {
                let singular: f64 = derived.iter().map(|s| s.value(omega)).sum();
                pair_value(only(&a), only(&b), omega, tol) - singular
            };
            let reach = pair_support(only(&a), only(&b));
            let fitted = InterpolatedSF::from_parts(reach, derived.clone(), regular, Some(tol));
            assert!(
                fitted.fit_error() < want,
                "{name} fitted to {:.2e}, wanted better than {want:.0e}",
                fitted.fit_error()
            );
        }
    }

    /// Two band edges of the same width meet at one frequency, and neither has a
    /// singular part: nothing pairs there, so neither can claim the kink twice.
    #[test]
    fn equal_band_edges_derive_nothing() {
        let d = 1.5f64;
        let (a, b) = (flat(0.0, d, 0.0), flat(0.0, d, 0.0));
        let derived = pair_singularities(only(&a), only(&b));
        assert_eq!(derived.len(), 3);
        assert!(derived.iter().all(|s| s.is_trivial()));

        // The kink is left to the regular part, which is where it is fitted exactly
        let h = 1e-3f64;
        let at = |w: f64| pair_value(only(&a), only(&b), w, 1e-13);
        let kink = (at(h) + at(-h) - 2.0 * at(0.0)) / (2.0 * h);
        assert_relative_eq!(kink, -1.0 / (4.0 * d * d), max_relative = 1e-9);
    }

    /// Two boxes convolve into a triangle, which is worth knowing exactly.
    #[test]
    fn two_boxes_make_a_triangle() {
        let (d, tol) = (1.5f64, 1e-12);
        let (a, b) = (flat(0.0, d, 0.0), flat(0.0, d, 0.0));
        let triangle = |w: f64| (2.0 * d - w.abs()).max(0.0) / (4.0 * d * d);
        for i in 0..=40 {
            let w = -2.0 * d + 4.0 * d * (i as f64) / 40.0;
            assert_relative_eq!(
                pair_value(only(&a), only(&b), w, tol),
                triangle(w),
                epsilon = 1e-12
            );
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
                pair_value(only(&a), only(&b), w, 1e-12),
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
        assert_eq!(pair_value(only(&a), only(&b), 0.0, 1e-12), f64::INFINITY);
        assert_eq!(square(0.0, 1.0).continuous_at(0.0), f64::INFINITY);

        // Either side of it the value is finite and right
        for w in [-1e-6f64, 1e-6] {
            assert_relative_eq!(
                pair_value(only(&a), only(&b), w, 1e-12),
                square(0.0, 1.0).continuous_at(w),
                max_relative = 1e-9
            );
        }

        // Two boxes have no singular parts to collide, so their band centre is finite
        let flat_pair = flat(0.0, 1.5, 0.0);
        assert!(pair_value(only(&flat_pair), only(&flat_pair), 0.0, 1e-12).is_finite());
    }

    /// The support adds, and the convolution vanishes beyond it.
    #[test]
    fn support_adds() {
        let (a, b) = (chain(0.0, 1.0), flat(0.5, 1.0, 0.0));
        let reach = pair_support(only(&a), only(&b));
        assert_relative_eq!(reach.min(), -2.5);
        assert_relative_eq!(reach.max(), 3.5);
        assert_eq!(pair_value(only(&a), only(&b), 4.0, 1e-10), 0.0);
        assert_eq!(pair_value(only(&a), only(&b), -3.0, 1e-10), 0.0);
    }
}
