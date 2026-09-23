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

use std::ops::{AddAssign, Mul};

use special::Gamma;

use crate::ContinuousSF;
use crate::beta::{
    beta_derivatives, beta_derivatives_laurent, inc_beta_derivatives, inc_beta_derivatives_split,
};
use crate::laurent::Laurent;
use crate::segment::Segment;
use crate::singularity::{AsymptTerm, Singularity, UnscaledAsymptTerm, power_log_integral};
use crate::util::{
    Table, alternating_sign, bilby_integrate_or_0, binomials, is_natural, negate_table,
    subtract_tables,
};

//
// Combinations both halves read a derivative table through
//

/// Coefficients of $\ln^m|\Delta|$ from a middle-stretch derivative table, indexed by
/// $m$.
///
/// $\ln x = \ln\Delta + \ln t$ and $\ln(\Delta-x) = \ln\Delta + \ln(1-t)$, so the two
/// exponents stay apart and each entry is used once.
fn middle_combine(table: &Table, m1: usize, m2: usize) -> Vec<f64> {
    let degree = m1 + m2;
    let (c_row1, c_row2) = (binomials(m1), binomials(m2));
    let mut out = vec![0.0; degree + 1];
    for (j, &c1) in c_row1.iter().enumerate() {
        for (k, &c2) in c_row2.iter().enumerate() {
            out[degree - j - k] += c1 * c2 * table[j][k];
        }
    }
    out
}

/// Coefficients of $\ln^m|\Delta|$ from an outer-stretch derivative table, indexed by
/// $m$.
///
/// There $\ln(1+u)$ is $-\partial_b$ and $\ln u$ the difference of the two, so
/// $(\partial_a - \partial_b)^j(-\partial_b)^k$ has to be expanded before the table is
/// read: $\sum_i \binom{j}{i} (-1)^{i+k} \partial_a^{j-i}\partial_b^{i+k}$.
///
/// The entries may be numbers or series in $\epsilon$: the generic formula reads a table
/// of the one and the whole-number case a table of the other. Reading both here is what
/// keeps the two from drifting apart, `zero` saying what an empty sum of entries is.
fn outer_combine<T>(table: &[Vec<T>], m1: usize, m2: usize, zero: &T) -> Vec<T>
where
    T: Clone + AddAssign<T>,
    for<'a> &'a T: Mul<f64, Output = T>,
{
    let degree = m1 + m2;
    let (c_row1, c_row2) = (binomials(m1), binomials(m2));
    let mut out = vec![zero.clone(); degree + 1];
    for (j, &c1) in c_row1.iter().enumerate() {
        let c_row_inner = binomials(j);
        for (k, &c2) in c_row2.iter().enumerate() {
            let mut inner = zero.clone();
            for (i, &c_inner) in c_row_inner.iter().enumerate() {
                inner += &table[j - i][i + k] * (c_inner * alternating_sign(i + k));
            }
            out[degree - j - k] += &inner * (c1 * c2);
        }
    }
    out
}

//
// Support of the result
//

/// Support of $C_1 \ast C_2$, which reaches as far as the two supports added together.
pub fn pair_support(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF) -> Segment {
    let (s1, s2) = (c1.support(), c2.support());
    assert!(
        s1.is_bounded() && s2.is_bounded(),
        "convolution of continuous parts is restricted to bounded supports"
    );
    s1 + s2
}

//
// Singular structure of the result
//

/// Singular structure of $C_1 \ast C_2$.
///
/// $\int C_1(\nu) C_2(\omega-\nu) d\nu$ stops being smooth in $\omega$ where both factors
/// stop being smooth at one and the same $\nu$, that is where $\nu = \Omega_1$ and
/// $\omega - \nu = \Omega_2$ hold together. The positions are therefore the pairwise sums,
/// and the ends of a support count: that is where a factor stops contributing at all.
///
/// Every pair of unscaled terms is convolved but one: the two constants. That product is
/// $R \ast R$, which leaves a kink going as $|\Delta|$ or milder - a polynomial either
/// side of the panel boundary it sits on, and fitted there exactly. Deriving it would
/// gain nothing and cost the equal-width case, where two band edges of the same width
/// would each claim the whole of the same kink. Everything with a singular part on one
/// side or the other is derived.
///
/// A pair of two singular parts derives one term more: the constant it leaves at
/// $\Delta = 0$, so that its share of the regular part vanishes there. That one takes
/// the geometry, the stretches it integrates over reaching only as far as the overlap
/// of the two supports.
pub fn pair_singularities(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF) -> Vec<Singularity> {
    let (s1, s2) = (c1.support(), c2.support());
    // How a contribution reads about one of its features does not depend on what that
    // feature is paired with, so each form is taken once rather than once per pairing
    let forms = |csf: &dyn ContinuousSF| -> Vec<(f64, LocalForm)> {
        features(csf)
            .into_iter()
            .map(|p| (p, unscaled_form_at(csf, p)))
            .collect()
    };
    let (forms1, forms2) = (forms(c1), forms(c2));

    let mut derived: Vec<(f64, Vec<AsymptTerm>)> = Vec::new();
    for (p1, (singular1, constant1)) in &forms1 {
        for (p2, (singular2, constant2)) in &forms2 {
            let reach = reach_between(s1, *p1, s2, *p2);
            let mut terms = Vec::new();
            for &t1 in singular1 {
                for &t2 in singular2 {
                    terms.extend(convolve_terms(t1, t2, Some(reach)));
                }
                terms.extend(convolve_terms(t1, *constant2, None));
            }
            for &t2 in singular2 {
                terms.extend(convolve_terms(*constant1, t2, None));
            }
            // Several pairs can land at one frequency. A pair deriving nothing is
            // recorded all the same: the position says where the interpolation must
            // start a fresh panel, whether or not anything is derived there.
            let position = p1 + p2;
            match derived.iter_mut().find(|(p, _)| *p == position) {
                Some((_, already)) => already.extend(terms),
                None => derived.push((position, terms)),
            }
        }
    }

    derived.sort_by(|x, y| x.0.total_cmp(&y.0));
    derived
        .into_iter()
        .map(|(position, terms)| Singularity::new(position, 1.0, terms))
        .collect()
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

/// The terms of a singularity sitting at a feature, and the value everything else takes
/// there.
type LocalForm = (Vec<UnscaledAsymptTerm>, UnscaledAsymptTerm);

/// How `csf` reads about `position`, in powers of $|\omega - \text{position}|$: the
/// terms of a singularity sitting there, and the value everything else takes, kept
/// apart.
///
/// Each is zeroed on the side the support does not reach. A singular part at an end of
/// its support has no side there, and reading a term as written would have a band edge
/// convolve as though the band ran on through it.
fn unscaled_form_at(csf: &dyn ContinuousSF, position: f64) -> LocalForm {
    let support = csf.support();
    let (reaches_below, reaches_above) = (position > support.min(), position < support.max());
    let sided = |exponent, log_power, c_below: f64, c_above: f64| UnscaledAsymptTerm {
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
                sing.unscaled_terms()
                    .into_iter()
                    .map(|t| sided(t.exponent, t.log_power, t.c_below, t.c_above)),
            );
        } else {
            constant += sing.value(position);
        }
    }
    (singular, sided(0.0, 0, constant, constant))
}

/// The range of $x = \nu - p_1$ the integral runs over where a feature `p1` of `s1`
/// meets a feature `p2` of `s2`, which is at $\omega = p_1 + p_2$.
///
/// Both factors turn at $x = 0$ there, the second being met at $\omega - \nu$ and so
/// turning at $\nu = \omega - p_2$, which is $p_1$ again. Being met that way is also why
/// its support enters mirrored about its own feature where the first enters displaced by
/// its own: going down in $\nu$ is going up in the second factor's variable.
///
/// Both operands hold the origin, a feature lying within its own support, so they always
/// meet and the range always reaches either side of the point.
fn reach_between(s1: Segment, p1: f64, s2: Segment, p2: f64) -> Segment {
    s1.shifted(-p1)
        .intersection(s2.mirrored(p2))
        .expect("a feature lies within its own support")
}

/// Exponent above which a derived term is smooth enough to leave to the interpolation.
///
/// $|\Delta|^2$ has two derivatives, which is more than a panel boundary asks for. The
/// constant a pair leaves at $\Delta = 0$ is derived whatever the exponent comes to,
/// that being the one term the interpolation cannot make good on its own.
const MAX_DERIVED_EXPONENT: f64 = 2.0;

/// Terms of $T_1 \ast T_2$ at $\Omega_1 + \Omega_2$.
///
/// The result spans $|\Delta|^r \ln^m|\Delta|$ for $m$ up to $m_1 + m_2$, or one higher
/// where $r$ is a whole number and the pole of $\Gamma(-r)$ trades a finite part for a
/// logarithm. Each coefficient is the three stretches added together with the sides they
/// draw on:
/// $$
///     C^+_m = c_1^+ c_2^+ X^{mid}_m + c_1^- c_2^+ X^{lo}_m + c_1^+ c_2^- X^{hi}_m,
/// $$
/// and $C^-_m$ the same with every $\pm$ flipped.
///
/// `reach` is how far the pair runs either side of the point, and is present
/// only for a pair of two singular parts. Such a pair has a value at $\Delta = 0$ that
/// is known in closed form, so it also derives the constant it leaves behind and what
/// goes to the regular part vanishes at the point. A singular part met with the other
/// factor's local constant has no such value - the truth there is
/// $\int S_p(\nu) R(\Omega_1+\Omega_2-\nu)d\nu$ rather than $c \int S_p$ - and derives
/// the $|\Delta|^r$ family alone.
///
/// Past `MAX_DERIVED_EXPONENT` the family is left to the interpolation and the constant
/// is all that comes back, so nothing at all comes back only where every coefficient
/// vanishes.
fn convolve_terms(
    t1: UnscaledAsymptTerm,
    t2: UnscaledAsymptTerm,
    reach: Option<Segment>,
) -> Vec<AsymptTerm> {
    let r = t1.exponent + t2.exponent + 1.0;
    // The constant is the ln^0 member of the family at r = 0 rather than a term of its
    // own, both sitting at |Δ|^0. Only a pair of singular parts can reach it: one drawn
    // against a local constant has r > 0.
    let constant = reach.map_or(0.0, |reach| coincident_constant(t1, t2, reach));
    debug_assert!(
        r != 0.0 || reach.is_some(),
        "r = 0 needs two singular parts"
    );

    let (cb1, ca1, cb2, ca2) = (t1.c_below, t1.c_above, t2.c_below, t2.c_above);
    let mut terms = Vec::new();

    if r < MAX_DERIVED_EXPONENT {
        let [mid, lo, hi] = if is_natural(r) {
            degenerate_stretches(t1, t2, r as usize)
        } else {
            generic_stretches(t1, t2)
        };
        for k in (0..mid.len()).rev() {
            let combine = |x1: f64, x2: f64, x3: f64| {
                (
                    cb1 * cb2 * x1 + ca1 * cb2 * x2 + cb1 * ca2 * x3,
                    ca1 * ca2 * x1 + cb1 * ca2 * x2 + ca1 * cb2 * x3,
                )
            };
            let (mut c_below, mut c_above) = combine(mid[k], lo[k], hi[k]);
            if r == 0.0 && k == 0 {
                c_below += constant;
                c_above += constant;
            }
            if c_below == 0.0 && c_above == 0.0 {
                continue;
            }
            terms.push(AsymptTerm::sided_log(r, k as u8, c_below, c_above));
        }
    }

    // Everything above vanishes at the point for r > 0, so what the pair leaves there
    // is a term of its own. It reaches both sides alike.
    if r != 0.0 && constant != 0.0 {
        terms.push(AsymptTerm::power(0.0, constant));
    }
    terms
}

/// What a pair of terms comes to at $\Delta = 0$ once the $|\Delta|^{\rho}$ family has
/// gone, with $\rho = r_1 + r_2 + 1$.
///
/// The middle stretch has no length there, and the outer two have their two factors
/// merge into one: $\int_0^L s^{\rho-1}\ln^M\\!s\\,ds$ for $M = m_1 + m_2$, which
/// [`power_log_integral()`] continues below $\rho = 0$. At $\rho = 0$ itself it has a
/// pole, and what survives of it against the pole of the complete beta is
/// $\ln^{M+1}(L)/(M+1)$.
///
/// The same number reaches both sides of the point. The stretch a side draws on and the
/// coefficients it draws with swap together, so what they come to does not.
fn coincident_constant(t1: UnscaledAsymptTerm, t2: UnscaledAsymptTerm, reach: Segment) -> f64 {
    // The two factors merge into one term of exponent $\rho - 1$, so the pole sits at
    // an exponent of $-1$
    let exponent = t1.exponent + t2.exponent;
    let m = t1.log_power + t2.log_power;
    let stretch = |c: f64, l: f64| {
        // A stretch of no length meets a coefficient of no weight, the side a support
        // does not reach having been zeroed already
        if c == 0.0 || l == 0.0 {
            0.0
        } else if exponent == -1.0 {
            c * l.ln().powi(i32::from(m) + 1) / f64::from(m + 1)
        } else {
            c * power_log_integral(exponent, m, l)
        }
    };
    stretch(t1.c_below * t2.c_above, -reach.min()) + stretch(t1.c_above * t2.c_below, reach.max())
}

/// What the three stretches contribute to the coefficient of $|\Delta|^r\ln^m|\Delta|$,
/// in the order middle, lower, upper, each indexed by $m$.
///
/// Substituting the length of the stretch out of the integral turns every logarithm into
/// $\ln|\Delta| + \ln(\text{something of order one})$, and expanding those binomials
/// leaves $|\Delta|^r$ times a polynomial in $\ln|\Delta|$ of degree $m_1 + m_2$. Its
/// coefficients are the integrals
/// $$
///     \int_0^1 t^{r_1}(1-t)^{r_2}\ln^j t \ln^k(1-t)\\,dt
///         = \partial_a^j \partial_b^k B(a, b)
/// $$
/// over the middle stretch, at $a = r_1+1$, $b = r_2+1$, and
/// $$
///     \int_0^\infty u^{r_1}(1+u)^{r_2}\ln^j u \ln^k(1+u)\\,du
///         = (\partial_a - \partial_b)^j(-\partial_b)^k B(a, b)
/// $$
/// over an outer one, at $b = -r$: there $\ln(1+u)$ is $-\partial_b$ and $\ln u$ the
/// difference of the two, the integrand carrying $a$ in both factors.
fn generic_stretches(t1: UnscaledAsymptTerm, t2: UnscaledAsymptTerm) -> [Vec<f64>; 3] {
    let (r1, r2) = (t1.exponent, t2.exponent);
    let (m1, m2) = (usize::from(t1.log_power), usize::from(t2.log_power));
    let r = r1 + r2 + 1.0;

    // The same tables and the same combinations the stretches themselves take, so that
    // what is declared singular here and what they set aside cannot drift apart
    let outer = |ra: f64, ma: usize, mb: usize| {
        outer_combine(&beta_derivatives(ra + 1.0, -r, ma, ma + mb), ma, mb, &0.0)
    };
    [
        middle_coefficients(r1, m1, r2, m2),
        outer(r1, m1, m2),
        outer(r2, m2, m1),
    ]
}

/// What the middle stretch contributes to the coefficient of
/// $|\Delta|^r\ln^m|\Delta|$, indexed by $m$.
///
/// $a$ and $b$ are both positive whatever $r$ comes to, so this is the one stretch that
/// never meets a pole, and it is used at a whole-number $r$ unchanged.
fn middle_coefficients(r1: f64, m1: usize, r2: f64, m2: usize) -> Vec<f64> {
    middle_combine(&beta_derivatives(r1 + 1.0, r2 + 1.0, m1, m2), m1, m2)
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
fn degenerate_stretches(t1: UnscaledAsymptTerm, t2: UnscaledAsymptTerm, n: usize) -> [Vec<f64>; 3] {
    let (r1, r2) = (t1.exponent, t2.exponent);
    let (m1, m2) = (usize::from(t1.log_power), usize::from(t2.log_power));
    let degree = m1 + m2 + 1;
    // Δ^n carrying no logarithm is analytic. For n ≥ 1 it vanishes at the point along
    // with everything else here, so its coefficient is left to the regular part: the
    // closed form computes the convolution whole and what the asymptotics does not claim
    // stays there exactly. At n = 0 there is nothing analytic about it - it is the
    // constant the pair leaves behind, and it is derived along with the logarithms.
    let pad = |mut v: Vec<f64>| {
        v.resize(degree + 1, 0.0);
        if n > 0 {
            v[0] = 0.0;
        }
        v
    };
    [
        pad(middle_coefficients(r1, m1, r2, m2)),
        pad(outer_stretch_logs(r1, m1, m2, n)),
        pad(outer_stretch_logs(r2, m2, m1, n)),
    ]
}

/// What an outer stretch contributes to the logarithms at a whole-number $r$, indexed by
/// the power of $\ln|\Delta|$.
///
/// The stretch integrates over the factor whose exponent is `r_own` and sees the other
/// displaced by $\Delta$. Only the other factor's count of logarithms is wanted, as
/// `m_other`: its exponent reaches the answer through $n$, which is
/// $r_{own} + r_{other} + 1$.
///
/// At a whole number the generic formula loses its finite parts to a pole of
/// $\Gamma(-r)$, but not its content: writing $r = n + \epsilon$ makes
/// $|\Delta|^r = |\Delta|^n e^{\epsilon\ln|\Delta|}$, so a pole of order $p$ meets
/// $\epsilon^p \ln^p|\Delta|/p!$ and leaves a logarithm behind. What survives is
/// $$
///     \sum_m \ln^m|\Delta| \sum_p \frac{(X_m)_{-p}}{p!} \ln^p|\Delta|,
/// $$
/// the $(X_m)_{-p}$ being the Laurent coefficients of the generic answer, which
/// [`beta_derivatives_laurent()`] supplies about $b = -n + \epsilon$ - read here at
/// $-\epsilon$, $b$ being $-r$.
///
/// The $\ln^0$ coefficient comes back as the $\epsilon^0$ part alone. What goes with it
/// is the stretch's own length, through the truncation the pole cancels against, and
/// that is [`coincident_constant()`]'s to add.
fn outer_stretch_logs(r_own: f64, m_own: usize, m_other: usize, n: usize) -> Vec<f64> {
    let degree = m_own + m_other;
    // Room for the pole, for the derivatives taken of it, and slack for the expansion
    let w = degree + 6;
    let table: Vec<Vec<Laurent>> = beta_derivatives_laurent(r_own + 1.0, n, degree, degree, w)
        .iter()
        .map(|row| row.iter().map(Laurent::reflected_in_epsilon).collect())
        .collect();

    // X_m, each a series in its own right: only sums and scalings of series enter
    let x = outer_combine(&table, m_own, m_other, &Laurent::zero(w));

    let mut out = vec![0.0; degree + 2];
    for (m, xm) in x.iter().enumerate() {
        for p in 0..=w {
            let c = xm.at(-(p as isize));
            if c == 0.0 {
                continue;
            }
            debug_assert!(
                m + p < out.len(),
                "a pole deeper than the logarithms it can leave behind"
            );
            out[m + p] += c / Gamma::gamma(p as f64 + 1.0);
        }
    }
    out
}

//
// Value of the result at one frequency
//

/// $\int C_1(\nu) C_2(\omega-\nu) d\nu$ over the overlap of the two supports.
///
/// The whole of it: the value the convolution takes at `omega`, not the part of it
/// [`pair_singularities()`] accounts for.
pub fn pair_value(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF, omega: f64, tol: f64) -> f64 {
    let (s1, s2) = (c1.support(), c2.support());
    let Some(overlap) = s1.intersection(s2.mirrored(omega)) else {
        return 0.0;
    };
    if overlap.is_degenerate() {
        return 0.0;
    }

    // A singular part carrying no terms is nothing to integrate, and pairing one would
    // ask `singular_singular()` for the stretch between two points that may be the same
    let parts1: Vec<Part> = c1
        .singularities()
        .iter()
        .filter(|sing| !sing.is_trivial())
        .map(|sing| Part {
            sing,
            support: s1,
            reflected: None,
        })
        .collect();
    let parts2: Vec<Part> = c2
        .singularities()
        .iter()
        .filter(|sing| !sing.is_trivial())
        .map(|sing| Part {
            sing,
            support: s2,
            reflected: Some(omega),
        })
        .collect();

    // R_A ⊛ R_B, both smooth across the overlap with nothing to subtract
    let mut total =
        bilby_integrate_or_0(|nu| c1.regular(nu) * c2.regular(omega - nu), overlap, tol);

    // S_p ⊛ R_B and R_A ⊛ S_q, one divergence each against a factor that has none
    for part in &parts1 {
        total += subtracted(part, |nu| c2.regular(omega - nu), overlap, tol);
    }
    for part in &parts2 {
        total += subtracted(part, |nu| c1.regular(nu), overlap, tol);
    }

    // S_p ⊛ S_q, the only term where two divergences can meet
    for first in &parts1 {
        for second in &parts2 {
            total += singular_singular(first, second, overlap);
        }
    }
    total
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

    /// Value at `nu`, zero outside the support and wherever $S$ has no finite value,
    /// the one frequency it diverges at carrying no weight under the integral.
    fn value(&self, nu: f64) -> f64 {
        let argument = self.argument(nu);
        if !self.support.contains(argument) {
            return 0.0;
        }
        let value = self.sing.value(argument);
        if value.is_finite() { value } else { 0.0 }
    }

    /// The singularity's terms with the side its support does not reach zeroed.
    ///
    /// The same reading [`unscaled_form_at()`] takes for the derivation. A stretch on
    /// that side has no length, so the value does not depend on it - but the singular
    /// coefficients do not know the geometry, and they have to agree with what is
    /// derived from the very same terms.
    fn sided_terms(&self) -> Vec<UnscaledAsymptTerm> {
        let position = self.sing.position();
        let reaches_below = position > self.support.min();
        let reaches_above = position < self.support.max();
        self.sing
            .unscaled_terms()
            .into_iter()
            .map(|t| UnscaledAsymptTerm {
                c_below: if reaches_below { t.c_below } else { 0.0 },
                c_above: if reaches_above { t.c_above } else { 0.0 },
                ..t
            })
            .collect()
    }

    /// $\int S$ over a stretch of the $\nu$ axis, in closed form.
    ///
    /// A reflected stretch has to be carried into the singularity's own variable, which
    /// is measuring it off the point it holds and mirroring that about $\Omega_q$.
    /// Mirroring about $\omega$ in one step would return $\Omega_q$ only to within
    /// rounding, and could put an end the wrong side of the very point it is there to
    /// enclose.
    fn integral(&self, segment: Segment) -> f64 {
        let own = match self.reflected {
            Some(_) => segment.shifted(-self.at()).mirrored(self.sing.position()),
            None => segment,
        };
        match own.intersection(self.support) {
            Some(within) => self.sing.integral(within),
            None => 0.0,
        }
    }
}

/// $\int S(\nu) g(\nu) d\nu$ over `segment`, with `g` anchored at the singular point.
///
/// Anchoring leaves $S(\nu)[g(\nu) - g(\Omega)]$ bounded where the bare product is not,
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
/// The two points must not coincide. There the middle stretch has no length and
/// $\ln|\Delta|$ no value, and nothing asks: $\omega = \Omega_p + \Omega_q$ is a singular
/// point of the result and so a panel boundary, while the regular part is sampled
/// strictly inside a panel.
fn singular_singular(first: &Part, second: &Part, segment: Segment) -> f64 {
    let (p, q) = (first.at(), second.at());
    let raw = q - p;
    assert!(
        raw != 0.0,
        "a pair has no value where its two singular points meet"
    );

    // Everything measured from the first singular point, then turned to face right
    let reach = segment.shifted(-p);
    let (delta, reach, flipped) = if raw < 0.0 {
        (-raw, reach.mirrored(0.0), true)
    } else {
        (raw, reach, false)
    };

    let sided = |t: &UnscaledAsymptTerm| {
        if flipped {
            (t.c_above, t.c_below)
        } else {
            (t.c_below, t.c_above)
        }
    };

    // The integrand changes form where either factor turns, so each stretch is the reach
    // met with one of the three pieces the cuts at x = 0 and x = Δ leave behind
    let below_zero = Segment::new(f64::NEG_INFINITY, 0.0);
    let between = Segment::new(0.0, delta);
    let past_delta = Segment::new(delta, f64::INFINITY);

    let mut total = 0.0;
    let (terms1, terms2) = (first.sided_terms(), second.sided_terms());
    for t1 in &terms1 {
        let (m1, r1) = (usize::from(t1.log_power), t1.exponent);
        let (below1, above1) = sided(t1);
        for t2 in &terms2 {
            let (m2, r2) = (usize::from(t2.log_power), t2.exponent);
            let (below2, above2) = sided(t2);
            let rho = r1 + r2 + 1.0;
            let mut take = |weight: f64, stretch: StretchContrib| {
                total += weight * (contract(&stretch.singular, rho, delta) + stretch.regular);
            };

            // x below zero: the first factor is met from below, the second from above,
            // and an outer stretch runs away from its point rather than towards it
            if let Some(s) = stretch_within(reach, below_zero)
                && below1 != 0.0
                && above2 != 0.0
            {
                take(
                    below1 * above2,
                    outer_stretch(r1, m1, r2, m2, delta, s.mirrored(0.0)),
                );
            }
            // between the two points, both factors met from above
            if let Some(x) = stretch_within(reach, between)
                && above1 != 0.0
                && above2 != 0.0
            {
                take(above1 * above2, middle_stretch(r1, m1, r2, m2, delta, x));
            }
            // past the second point: the first from above, the second from below, the
            // stretch measured off that second point
            if let Some(s) = stretch_within(reach, past_delta)
                && above1 != 0.0
                && below2 != 0.0
            {
                take(
                    above1 * below2,
                    outer_stretch(r2, m2, r1, m1, delta, s.shifted(-delta)),
                );
            }
        }
    }
    total
}

/// The part of `reach` lying within `region`, absent where the two share no length.
fn stretch_within(reach: Segment, region: Segment) -> Option<Segment> {
    reach.intersection(region).filter(|s| !s.is_degenerate())
}

/// What a stretch comes to, in two halves.
///
/// The singular half is $\sum_m c_m |\Delta|^{\rho}\ln^m|\Delta|$ and is what reaches
/// $c^{\pm}$; the regular half is what is left at this $\Delta$, analytic through
/// $\Delta = 0$. Their sum is the integral over the stretch.
struct StretchContrib {
    /// $c_m$, indexed by the power of $\ln|\Delta|$.
    singular: Vec<f64>,
    /// Everything else, evaluated at this $\Delta$.
    regular: f64,
}

impl StretchContrib {
    /// A stretch whose beta could not be divided, all of it counted as regular.
    fn undivided(degree: usize, value: f64) -> StretchContrib {
        StretchContrib {
            singular: vec![0.0; degree + 1],
            regular: value,
        }
    }
}

/// $\int_s s^{r_1}\ln^{m_1}\\!s\;(\Delta+s)^{r_2}\ln^{m_2}\\!(\Delta+s)\\,ds$, over a stretch
/// `s` of the positive axis and for $\Delta > 0$.
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
fn outer_stretch(r1: f64, m1: usize, r2: f64, m2: usize, delta: f64, s: Segment) -> StretchContrib {
    let rho = r1 + r2 + 1.0;
    let (a, b) = (r1 + 1.0, -rho);
    let z = |t: f64| t / (delta + t);
    let degree = m1 + m2;

    // Writing $B_z = B - D_z$, the stretch is $\Delta^\rho(D_{z_1} - D_{z_2})$: the
    // complete beta cancels between the ends unless the near one sits on the singular
    // point, where $z_1 = 0$ makes $D_{z_1}$ the whole of it.
    let Some((complete, far)) = inc_beta_derivatives_split(a, b, z(s.max()), m1, degree) else {
        let upper = inc_beta_derivatives(a, b, z(s.max()), m1, degree);
        let lower = inc_beta_derivatives(a, b, z(s.min()), m1, degree);
        let step = subtract_tables(&upper, &lower);
        let coefficients = outer_combine(&step, m1, m2, &0.0);
        return StretchContrib::undivided(degree, contract(&coefficients, rho, delta));
    };
    let (singular, rest) = if s.min() == 0.0 {
        (outer_combine(&complete, m1, m2, &0.0), negate_table(&far))
    } else {
        let (_, near) = inc_beta_derivatives_split(a, b, z(s.min()), m1, degree)
            .expect("the far end already split, so the near one does too");
        (vec![0.0; degree + 1], subtract_tables(&near, &far))
    };
    let regular = contract(&outer_combine(&rest, m1, m2, &0.0), rho, delta);
    StretchContrib { singular, regular }
}

/// $\int_x x^{r_1}\ln^{m_1}\\!x\;(\Delta-x)^{r_2}\ln^{m_2}\\!(\Delta-x)\\,dx$, over a stretch
/// `x` within $[0, \Delta]$ and for $\Delta > 0$.
///
/// The stretch between the two singular points. Here $x = \Delta t$ is the whole
/// substitution, and the two exponents stay apart: $a = r_1+1$ and $b = r_2+1$ are both
/// positive, so this is the one stretch whose Beta never meets a pole.
fn middle_stretch(
    r1: f64,
    m1: usize,
    r2: f64,
    m2: usize,
    delta: f64,
    x: Segment,
) -> StretchContrib {
    let rho = r1 + r2 + 1.0;
    let (a, b) = (r1 + 1.0, r2 + 1.0);
    let degree = m1 + m2;

    // Both ends of this stretch are singular points, so neither owns the complete beta:
    // written as $B_{t_2} - B_{t_1}$ it enters at the upper limit and written as
    // $D_{t_1} - D_{t_2}$ at the lower, the two being the same sum. It belongs to the
    // stretch when the stretch runs the whole way between the points and not otherwise,
    // and there the regular half is exactly zero rather than merely small.
    let Some((complete, far)) = inc_beta_derivatives_split(a, b, x.max() / delta, m1, m2) else {
        let upper = inc_beta_derivatives(a, b, x.max() / delta, m1, m2);
        let lower = inc_beta_derivatives(a, b, x.min() / delta, m1, m2);
        let step = subtract_tables(&upper, &lower);
        let coefficients = middle_combine(&step, m1, m2);
        return StretchContrib::undivided(degree, contract(&coefficients, rho, delta));
    };
    let (singular, rest) = if x.min() == 0.0 && x.max() == delta {
        (middle_combine(&complete, m1, m2), negate_table(&far))
    } else {
        let (_, near) = inc_beta_derivatives_split(a, b, x.min() / delta, m1, m2)
            .expect("the far end already split, so the near one does too");
        (vec![0.0; degree + 1], subtract_tables(&near, &far))
    };
    let regular = contract(&middle_combine(&rest, m1, m2), rho, delta);
    StretchContrib { singular, regular }
}

/// $|\Delta|^{\rho}\sum_m c_m \ln^m|\Delta|$.
fn contract(coefficients: &[f64], rho: f64, delta: f64) -> f64 {
    let ln_delta = delta.ln();
    let sum: f64 = coefficients
        .iter()
        .enumerate()
        .map(|(m, c)| c * ln_delta.powi(m as i32))
        .sum();
    delta.powf(rho) * sum
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

        // The band edges meet at r = 0 too, and what the pair leaves there is a
        // constant: the step the square lattice band edge really is, 1/(4πt). Only one
        // side of each edge lies inside the band, and the other reaches nothing.
        for position in [-4.0f64, 4.0] {
            let edge = derived.iter().find(|s| s.position() == position).unwrap();
            let (inside, outside) = if position < 0.0 {
                (position + 1e-6, position - 1e-6)
            } else {
                (position - 1e-6, position + 1e-6)
            };
            assert_relative_eq!(
                edge.value(inside),
                1.0 / (4.0 * std::f64::consts::PI * t),
                max_relative = 1e-12
            );
            assert_eq!(edge.value(outside), 0.0);
        }
    }

    /// What a pair of singular parts leaves to the regular part goes to zero at
    /// $\Omega_p + \Omega_q$ rather than to a constant.
    ///
    /// That is the one term an interpolation cannot supply on its own: a panel boundary
    /// sits on every derived point, and the fit never asks what the regular part is
    /// there. A constant left behind would be a step the two panels straddle.
    #[test]
    fn a_pair_leaves_nothing_to_the_regular_part() {
        // Each case picks one singular part out of either factor, and names the
        // exponent ρ = r₁ + r₂ + 1 the two come to
        let cases = [
            (
                "two chain edges, ρ = 0",
                chain(0.0, 1.0),
                0,
                chain(0.0, 1.0),
                1,
            ),
            (
                "two chain edges meeting, ρ = 0",
                chain(0.0, 1.0),
                1,
                chain(0.0, 1.0),
                1,
            ),
            (
                "a chain edge and a log, ρ = 1/2",
                chain(0.0, 1.0),
                1,
                square(0.0, 1.0),
                0,
            ),
            ("two logs, ρ = 1", square(0.0, 1.0), 0, square(0.0, 1.0), 0),
            (
                "a bethe edge and a chain edge, ρ = 1",
                bethe(3, 0.0, 1.0),
                1,
                chain(0.0, 1.0),
                0,
            ),
            (
                "two Dirac points, ρ = 3",
                honeycomb(0.0, 1.0),
                1,
                honeycomb(0.0, 1.0),
                1,
            ),
        ];
        for (name, a, ia, b, ib) in cases {
            let (c1, c2) = (only(&a), only(&b));
            let (s1, s2) = (c1.support(), c2.support());
            let (sing1, sing2) = (&c1.singularities()[ia], &c2.singularities()[ib]);
            let (p1, p2) = (sing1.position(), sing2.position());
            let reach = reach_between(s1, p1, s2, p2);

            // What the derivation claims for this pair and no other
            let first = Part {
                sing: sing1,
                support: s1,
                reflected: None,
            };
            let mut terms = Vec::new();
            for &t1 in &first.sided_terms() {
                let second = Part {
                    sing: sing2,
                    support: s2,
                    reflected: Some(p1 + p2),
                };
                for &t2 in &second.sided_terms() {
                    terms.extend(convolve_terms(t1, t2, Some(reach)));
                }
            }
            let claimed = Singularity::new(p1 + p2, 1.0, terms);

            // What is left over, either side of the point and twice as near the second
            // time round. It has to shrink rather than settle on a constant.
            for side in [-1.0f64, 1.0] {
                let residual = |h: f64| {
                    let omega = p1 + p2 + side * h;
                    let second = Part {
                        sing: sing2,
                        support: s2,
                        reflected: Some(omega),
                    };
                    // Past the end of the convolved support there is nothing to
                    // integrate, and the derivation must claim nothing either
                    let value = s1
                        .intersection(s2.mirrored(omega))
                        .map_or(0.0, |overlap| singular_singular(&first, &second, overlap));
                    value - claimed.value(omega)
                };
                let (far, near) = (residual(1e-4).abs(), residual(1e-6).abs());
                assert!(
                    near < 0.05 * far.max(1e-13),
                    "{name} left {near:.3e} at 1e-6 against {far:.3e} at 1e-4"
                );
            }
        }
    }

    /// $R(\omega)$ meets itself at a derived point, which is where a panel boundary
    /// sits and the fit never samples.
    ///
    /// A chain against itself brings two inverse square roots together at the band
    /// centre, $\rho = 0$, where the constant the pair leaves differs between the two
    /// sides. Left underived it would be a step in $R$ straddling the boundary.
    #[test]
    fn the_regular_part_meets_itself_at_a_derived_point() {
        let (a, b) = (chain(0.0, 1.0), chain(0.0, 1.0));
        let (c1, c2) = (only(&a), only(&b));
        let derived = pair_singularities(c1, c2);
        let regular = |omega: f64| {
            let singular: f64 = derived.iter().map(|s| s.value(omega)).sum();
            pair_value(c1, c2, omega, 1e-13) - singular
        };
        // The band centre, where the two inverse square roots collide
        let (mut below, mut above) = (0.0, 0.0);
        for h in [1e-4f64, 1e-6] {
            (below, above) = (regular(-h), regular(h));
            assert_relative_eq!(below, above, max_relative = 1e-6, epsilon = 1e-9);
        }
        // and it is a number rather than the difference of two infinities
        assert!(below.is_finite() && above.is_finite());
    }

    /// Below $\rho = 0$ the constant is an analytic continuation rather than an
    /// integral that converges, and it has to agree with what the stretch itself sets
    /// aside.
    #[test]
    fn the_constant_continues_below_rho_zero() {
        // Two edges going as |ω|^{-3/5} bring ρ to -1/5
        let (r, l) = (-0.6f64, 1.5f64);
        let term = |c_below: f64, c_above: f64| UnscaledAsymptTerm {
            exponent: r,
            log_power: 0,
            c_below,
            c_above,
        };
        let claimed = coincident_constant(term(1.0, 1.0), term(1.0, 1.0), Segment::new(-l, l));

        // What the stretches set aside as regular, taken nearer and nearer the point.
        // The two outer ones are the same call here, the exponents and the lengths
        // being equal, and the middle one runs the whole way between the points and so
        // sets aside nothing at all.
        let mut previous = f64::INFINITY;
        for delta in [1e-3f64, 1e-5, 1e-7] {
            let outer = outer_stretch(r, 0, r, 0, delta, Segment::new(0.0, l)).regular;
            let middle = middle_stretch(r, 0, r, 0, delta, Segment::new(0.0, delta));
            assert_eq!(middle.regular, 0.0);
            let error = (2.0 * outer - claimed).abs();
            assert!(
                error < previous,
                "{error:.3e} did not improve on {previous:.3e}"
            );
            previous = error;
        }
        assert!(previous < 1e-6 * claimed.abs(), "stalled at {previous:.3e}");
    }

    /// The constant a pair leaves at $\rho = 0$ differs between the two sides, which is
    /// the whole reason an [`AsymptTerm`] may carry a bare constant per side.
    ///
    /// The middle stretch is what splits them: it takes both factors from above below
    /// the point and both from below above it, and at $\rho = 0$ it no longer vanishes
    /// there. The logarithm beside it is symmetric, the two outer stretches swapping
    /// their lengths along with their coefficients.
    #[test]
    fn the_constant_a_pair_leaves_is_sided() {
        let term = |c_below: f64, c_above: f64| UnscaledAsymptTerm {
            exponent: -0.5,
            log_power: 0,
            c_below,
            c_above,
        };
        let pair = Singularity::new(
            0.0,
            1.0,
            convolve_terms(
                term(1.0, 2.0),
                term(0.7, 0.3),
                Some(Segment::new(-1.0, 3.0)),
            ),
        );

        // Two inverse square roots meet at ρ = 0, so S(±d) = C^± + L ln(d). Two points
        // a side separate the two.
        let (d1, d2) = (1e-3f64, 1e-6f64);
        let split = |side: f64| {
            let (s1, s2) = (pair.value(side * d1), pair.value(side * d2));
            let slope = (s1 - s2) / (d1.ln() - d2.ln());
            (slope, s1 - slope * d1.ln())
        };
        let (slope_above, constant_above) = split(1.0);
        let (slope_below, constant_below) = split(-1.0);

        // -(c₁⁻c₂⁺ + c₁⁺c₂⁻) is what the two outer stretches come to, alike either side
        assert_relative_eq!(slope_above, -1.7, max_relative = 1e-12);
        assert_relative_eq!(slope_below, -1.7, max_relative = 1e-12);

        // The constant does not reach the two sides alike
        assert_relative_eq!(constant_above, 5.779_713_210_193, max_relative = 1e-11);
        assert_relative_eq!(constant_below, 6.093_872_475_552, max_relative = 1e-11);
    }

    /// A singular part with no terms is no pair, so asking where its point meets
    /// another's is a question with an answer rather than a panic.
    ///
    /// Two boxes derive three positions and no terms at any of them, so convolving
    /// their triangle again sets a termless point against a real one.
    #[test]
    fn a_termless_singular_part_is_no_pair() {
        let d = 1.5f64;
        let boxes = flat(0.0, d, 0.0).conv(&flat(0.0, d, 0.0), None);
        let band = chain(0.0, 1.0);
        let (triangle, chain_band) = (only(&boxes), only(&band));
        assert!(triangle.singularities().iter().all(|s| s.is_trivial()));

        // Where a termless point of the one meets a singular point of the other
        let omega =
            triangle.singularities()[1].position() + chain_band.singularities()[0].position();
        let at = |w: f64| pair_value(triangle, chain_band, w, 1e-11);
        let value = at(omega);
        assert!(value.is_finite(), "got {value}");

        // and it is the value the convolution really takes there
        let h = 1e-6;
        assert_relative_eq!(
            value,
            0.5 * (at(omega - h) + at(omega + h)),
            max_relative = 1e-6
        );
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
                1e-14,
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
                1e-12,
            ),
            (
                "triangular against itself",
                triangular(0.0, 1.0),
                triangular(0.0, 1.0),
                1e-12,
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

    /// Two inverse square roots colliding is the van Hove peak of the square lattice,
    /// and the derived structure is what says so: the divergence is read off the terms
    /// at $\omega = \Omega_p + \Omega_q$ rather than out of an integral there.
    #[test]
    fn colliding_points_diverge() {
        let t = 1.0f64;
        let derived = chain(0.0, t).conv(&chain(0.0, t), None);
        let known = square(0.0, t);
        assert_eq!(derived.continuous_at(0.0), f64::INFINITY);
        assert_eq!(known.continuous_at(0.0), f64::INFINITY);

        // Either side of it the value is finite and right
        for w in [-1e-6f64, 1e-6] {
            assert_relative_eq!(
                derived.continuous_at(w),
                known.continuous_at(w),
                max_relative = 1e-9
            );
        }

        // Two boxes have no singular parts to collide, so their band centre is finite
        let boxes = flat(0.0, 1.5, 0.0);
        assert!(boxes.conv(&boxes, None).continuous_at(0.0).is_finite());
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
