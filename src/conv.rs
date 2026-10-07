//! Convolution of the continuous parts of two spectral functions.
//!
//! $C_A = R_A + \sum_p S_p$ and $C_B = R_B + \sum_q S_q$, so the convolution is four
//! kinds of term:
//! $$
//!     C_A \ast C_B = R_A \ast R_B + \sum_p S_p \ast R_B + \sum_q R_A \ast S_q
//!                  + \sum_{p,q} S_p \ast S_q,
//! $$
//! and they are taken one at a time.

use special::Gamma;

use crate::ContinuousSF;
use crate::beta;
use crate::beta::{
    beta_derivatives, inc_beta_derivatives, inc_beta_derivatives_split, outer_beta_complete,
    outer_beta_derivatives, outer_beta_derivatives_split, outer_beta_entire_factor,
    outer_beta_tail,
};
use crate::polynomial::Polynomial;
use crate::segment::Segment;
use crate::singularity::{AsymptTerm, Singularity, power_log_integral};
use crate::util::{
    Table, alternating_sign, bilby_integrate_or_0, binomials, integrate_by_subtraction,
    negate_table, subtract_tables,
};

//
// General functions
//

/// Given a table of coefficients $\partial_a^j \partial_b^k T(a, b)$, extract coefficients
/// in front of $\ln^m|\Delta|$ that contribute to
/// $$
///     \partial_a^j \partial_b^k \left[ \vert\Delta\vert^{a+b-1} T(a, b) \right].
/// $$
///
/// `j_max` counts the derivatives on the table's first index and `k_max` those on its
/// second. For an outer stretch that is the order the outer Beta's parameters come in,
/// $a = r_2+1$ before $b = r_1+1$, and not the order of the two factors.
fn derivatives_to_log_coeffs(table: &Table, j_max: usize, k_max: usize) -> Polynomial {
    // All stretches scale their range by $\Delta = \omega-(\Omega_p+\Omega_q)$, so
    // restoring the coefficients of $\ln^m|\Delta|$ from the derivative table of a
    // dimensionless function is the general Leibniz rule
    let degree = j_max + k_max;
    let (c_row1, c_row2) = (binomials(j_max), binomials(k_max));
    let mut out = Polynomial::zeros(degree);
    for (j, &c1) in c_row1.iter().enumerate() {
        for (k, &c2) in c_row2.iter().enumerate() {
            out[degree - j - k] += c1 * c2 * table[j][k];
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
/// and the ends of a support count too: that is where a factor stops contributing at all.
///
/// Every pair of terms is convolved but one: the two constants. That product is
/// $R \ast R$, which leaves a kink going as $|\Delta|$ or milder - a polynomial either
/// side of the panel boundary it sits on, and fitted there exactly.
///
/// A pair of two singular parts derives one term more: the constant it leaves at
/// $\Delta \to 0$, so that its share of the regular part vanishes there. That one takes
/// the geometry, the stretches it integrates over reaching only as far as the overlap
/// of the two supports.
pub fn pair_singularities(c1: &dyn ContinuousSF, c2: &dyn ContinuousSF) -> Vec<Singularity> {
    let (s1, s2) = (c1.support(), c2.support());
    let forms = |csf: &dyn ContinuousSF| -> Vec<(f64, LocalForm)> {
        features(csf)
            .into_iter()
            .map(|p| (p, local_form_at(csf, p)))
            .collect()
    };
    let (forms1, forms2) = (forms(c1), forms(c2));

    let mut derived: Vec<(f64, Vec<AsymptTerm>)> = Vec::new();
    for (p1, (singular1, constant1)) in &forms1 {
        for (p2, (singular2, constant2)) in &forms2 {
            let reach = reach_between(s1, *p1, s2, *p2);
            let mut terms = Vec::new();
            for t1 in singular1 {
                for t2 in singular2 {
                    terms.extend(convolve_terms(t1, t2, Some(reach)));
                }
                terms.extend(convolve_terms(t1, constant2, None));
            }
            for t2 in singular2 {
                terms.extend(convolve_terms(constant1, t2, None));
            }
            // Several pairs can land at one frequency. A pair deriving no singularity is
            // recorded anyway: the position says where the interpolation must start
            // a fresh panel, whether or not anything is derived there.
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
type LocalForm = (Vec<AsymptTerm>, AsymptTerm);

/// How `csf` reads about `position`, in powers of $|\omega - \text{position}|$: the
/// terms of a singularity sitting there, and the value everything else takes, kept
/// apart.
///
/// Each is zeroed on the side the support does not reach. Without that a band edge would
/// convolve as though the band ran on through it.
fn local_form_at(csf: &dyn ContinuousSF, position: f64) -> LocalForm {
    let support = csf.support();
    let (reaches_below, reaches_above) = (position > support.min(), position < support.max());
    let sided = |t: &AsymptTerm| t.sides_kept(reaches_below, reaches_above);

    let mut singular = Vec::new();
    let mut constant = csf.regular(position);
    for sing in csf.singularities() {
        if sing.position() == position {
            for t in sing.terms() {
                singular.push(sided(t));
            }
        } else {
            constant += sing.value(position);
        }
    }
    (singular, sided(&AsymptTerm::power(0.0, constant)))
}

/// The range of $x = \nu - p_1$ where a feature `p1` of `s1`
/// meets a feature `p2` of `s2`, which is at $\omega = p_1 + p_2$.
///
/// Both operands hold the origin, a feature lying within its own support, so they always
/// meet.
fn reach_between(s1: Segment, p1: f64, s2: Segment, p2: f64) -> Segment {
    s1.shifted(-p1)
        .intersection(s2.mirrored(p2))
        .expect("a feature lies within its own support")
}

/// Asymptotic terms of $T_1 \ast T_2$ at $\omega = \Omega_1 + \Omega_2$.
///
/// The pair makes one term, $|\Delta|^\rho P^{\pm}(\ln|\Delta|)$ at
/// $\Delta=\omega-(\Omega_1 + \Omega_2)$ and $\rho = r_1 + r_2 + 1$, of degree
/// $m_1 + m_2$ or one higher where $\rho$ is a whole number and the pole of
/// $\Gamma(-\rho)$ trades a finite part for a logarithm. Every pair of powers the two
/// terms hold contributes the three stretches, weighted by the sides those powers draw
/// on:
/// $$
///     P^+ \mathrel{+}= c_1^+ c_2^+ X^{mid} + c_1^- c_2^+ X^{lo} + c_1^+ c_2^- X^{hi},
/// $$
/// and $P^-$ the same with every $\pm$ flipped.
///
/// `reach` is how far the pair runs either side of the point $\Omega_1 + \Omega_2$, and
/// is absent where one of the two is a local constant. Only a pair of singular parts has
/// a value at $\Delta \to 0$ known in closed form, so only such a pair derives the constant
/// it leaves there.
///
/// Past a $\rho$ of [`AsymptTerm::MAX_EXPONENT`] the family is left to the
/// interpolation and the constant is all that comes back, so nothing at all comes back
/// only where every coefficient vanishes.
fn convolve_terms(t1: &AsymptTerm, t2: &AsymptTerm, reach: Option<Segment>) -> Vec<AsymptTerm> {
    // The exponents alone fix $\rho$, so one pair of terms is one family however many
    // powers of the logarithm the two hold
    let (r1, r2) = (t1.exponent(), t2.exponent());
    let rho = r1 + r2 + 1.0;
    // The constant is the ln^0 member of the family at rho = 0.
    // Only a pair of singular parts can reach it: one drawn against a local constant has
    // rho > 0.
    let constant = reach.map_or(0.0, |reach| coincident_constant(t1, t2, reach));
    debug_assert!(
        rho != 0.0 || reach.is_some(),
        "rho = 0 needs two singular parts"
    );

    let mut terms = Vec::new();

    if rho < AsymptTerm::MAX_EXPONENT {
        // A band around each whole $\rho$ and not an exact test. At
        // $\rho = n + \epsilon$ the family reads
        // $|\Delta|^{\rho} = |\Delta|^n(1 + \epsilon\ln|\Delta| + \ldots)$, so the
        // logarithm surviving at $\epsilon \to 0$ has coefficient $\epsilon X_0$ and a
        // finite logarithm forces $X_0$ to go as $1/\epsilon$. The generic formula duly
        // returns that pole and the regular part the opposite of it, which is the same
        // cancellation `beta` refuses to hand out, so the width is shared with it.
        let degenerate = beta::near_a_whole_exponent(rho);
        let degree = t1.log_degree() + t2.log_degree() + usize::from(degenerate);
        let (mut below, mut above) = (Polynomial::zeros(degree), Polynomial::zeros(degree));

        // Each stretch is a Beta function of the two exponents differentiated $m_1$ and
        // $m_2$ times, so the powers of the logarithm pair off one by one and the three
        // weights are the sides those two powers draw on
        for (m1, cb1, ca1) in t1.log_coeffs() {
            for (m2, cb2, ca2) in t2.log_coeffs() {
                let [mid, lo, hi] = if degenerate {
                    degenerate_stretches(r1, m1, r2, m2, rho.round() as usize)
                } else {
                    generic_stretches(r1, m1, r2, m2)
                };
                below.add_scaled(&mid, cb1 * cb2);
                below.add_scaled(&lo, ca1 * cb2);
                below.add_scaled(&hi, cb1 * ca2);
                above.add_scaled(&mid, ca1 * ca2);
                above.add_scaled(&lo, cb1 * ca2);
                above.add_scaled(&hi, ca1 * cb2);
            }
        }
        if rho == 0.0 {
            below[0] += constant;
            above[0] += constant;
        }
        // A pair whose every coefficient cancels derives nothing at all
        let family = AsymptTerm::sided_log_poly(rho, below, above);
        if !family.is_zero() {
            terms.push(family);
        }
    }

    // Everything above vanishes at the point for rho > 0, so what the pair leaves there
    // is a term of its own. It reaches both sides alike.
    if rho != 0.0 && constant != 0.0 {
        terms.push(AsymptTerm::power(0.0, constant));
    }
    terms
}

/// What a pair of terms comes to at $\Delta \to 0$ once the $|\Delta|^{\rho}$ family has
/// gone, with $\rho = r_1 + r_2 + 1$.
fn coincident_constant(t1: &AsymptTerm, t2: &AsymptTerm, reach: Segment) -> f64 {
    // The two factors merge into one term of exponent $\rho - 1$, so the pole sits at
    // an exponent of $-1$, and that exponent is the same for every pair of powers
    let exponent = t1.exponent() + t2.exponent();
    let stretch = |c: f64, m: u8, l: f64| {
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

    let mut total = 0.0;
    for (m1, cb1, ca1) in t1.log_coeffs() {
        for (m2, cb2, ca2) in t2.log_coeffs() {
            let m = (m1 + m2) as u8;
            total += stretch(cb1 * ca2, m, -reach.min()) + stretch(ca1 * cb2, m, reach.max());
        }
    }
    total
}

/// What the three stretches contribute to the coefficient of $|\Delta|^\rho\ln^m|\Delta|$,
/// in the order middle, lower, upper, each indexed by $m$.
///
fn generic_stretches(r1: f64, m1: usize, r2: f64, m2: usize) -> [Polynomial; 3] {
    // Substituting the length of the stretch out of the integral turns every logarithm
    // into $\ln|\Delta| + O(1)$, and expanding those binomials leaves $|\Delta|^\rho$
    // times a polynomial in $\ln|\Delta|$ of degree $m_1 + m_2$, whose coefficients are
    // $$
    //     \int_0^1 t^{r_1}(1-t)^{r_2}\ln^j t \ln^k(1-t)\,dt
    //         = \partial_a^j \partial_b^k B(a, b)
    // $$
    // over the middle stretch, at $a = r_1+1$, $b = r_2+1$, and
    // $$
    //     \int_0^\infty u^{r_1}(1+u)^{r_2}\ln^j u \ln^k(1+u)\,du
    //         = \partial_b^j \partial_a^k B^{\mathrm{out}}_\infty(a, b)
    // $$
    // over an outer one, at $a = r_2+1$ and $b = r_1+1$.
    let outer = |r_own: f64, m_own: usize, r_other: f64, m_other: usize| {
        let table = outer_beta_complete(r_other + 1.0, r_own + 1.0, m_other, m_own);
        derivatives_to_log_coeffs(&table, m_other, m_own)
    };
    [
        middle_coefficients(r1, m1, r2, m2),
        outer(r1, m1, r2, m2),
        outer(r2, m2, r1, m1),
    ]
}

/// What the middle stretch contributes to the coefficient of
/// $|\Delta|^\rho\ln^m|\Delta|$, indexed by $m$.
///
/// $a$ and $b$ are both positive whatever $\rho$ comes to, so this is the one stretch that
/// never meets a pole, and it is used at a whole-number $\rho$ unchanged.
fn middle_coefficients(r1: f64, m1: usize, r2: f64, m2: usize) -> Polynomial {
    derivatives_to_log_coeffs(&beta_derivatives(r1 + 1.0, r2 + 1.0, m1, m2), m1, m2)
}

/// What the three stretches contribute where $\rho$ is a non-negative integer $n$, or
/// near enough to one that [`beta::WHOLE_EXPONENT_WIDTH`] claims it.
///
fn degenerate_stretches(r1: f64, m1: usize, r2: f64, m2: usize, n: usize) -> [Polynomial; 3] {
    // $\Delta^n$ is analytic, and the outer stretches no longer converge: their
    // integrands go as $u^{n-1}$ at large $u$, so one term of the expansion of
    // $(1+1/u)^{r_2}$ integrates to a logarithm rather than a power.
    let degree = m1 + m2 + 1;
    // Δ^n with no logarithm beside it is analytic. For n ≥ 1 it vanishes at the point along
    // with everything else here, so its coefficient is left to the regular part: the
    // closed form computes the convolution whole and what the asymptotics does not claim
    // stays there exactly. At n = 0 there is nothing analytic about it - it is the
    // constant the pair leaves behind, and it is derived along with the logarithms.
    let pad = |mut v: Polynomial| {
        v.resize(degree);
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

/// What an outer stretch contributes to the logarithms at a whole-number $\rho$, indexed by
/// the power of $\ln|\Delta|$.
///
/// The stretch integrates over the factor whose exponent is `r_own` and sees the other
/// displaced by $\Delta$. Only the other factor's count of logarithms is wanted, as
/// `m_other`: its exponent reaches the answer through $n$, which is
/// $\rho = r_{own} + r_{other} + 1$.
///
/// The $\ln^0$ coefficient comes back as the $\epsilon^0$ part alone. What goes with it
/// is the stretch's own length, and that is [`coincident_constant()`]'s to add.
fn outer_stretch_logs(r_own: f64, m_own: usize, m_other: usize, n: usize) -> Polynomial {
    // At a whole number the generic formula loses its finite parts to a pole of
    // $\Gamma(-\rho)$, but not its content: writing $\rho = n + \epsilon$ makes
    // $|\Delta|^\rho = |\Delta|^n e^{\epsilon\ln|\Delta|}$, so a pole of order $p$
    // meets $\epsilon^p \ln^p|\Delta|/p!$ and leaves a logarithm behind, giving
    // $\sum_m \ln^m|\Delta| \sum_p \frac{(X_m)_{-p}}{p!} \ln^p|\Delta|$.
    //
    // `outer_beta_entire_factor()` writes the complete outer beta as $F(\delta)/\delta$
    // with $F$ analytic, which reads off the pole and the finite part at once. Since
    // $\rho = n-\delta$ puts $\delta = -\epsilon$, the entry at
    // $\partial_b^{j}\partial_\delta^{k}$ has its whole pole at order $k+1$, worth
    // $-k!\,\partial_b^{j} f_0$, and its $\epsilon^0$ part is
    // $k!\,\partial_b^{j} f_{k+1}$. That $f_0$ is a polynomial in $b$, so no $\Gamma$
    // survives into the logarithms.
    //
    // The table is indexed by $b$ and $\delta$ and not by the outer Beta's own two
    // parameters, the pole living in $\delta$ alone, so this is the one outer stretch
    // with a combination to expand before the table is read. The expansion is the change
    // of variables itself: $\partial_{r_2+1} = -\partial_\delta$ and
    // $\partial_{r_1+1} = \partial_b - \partial_\delta$.
    let degree = m_own + m_other;
    // $\partial_b^j f_l$: the pole wants $l = 0$ and the finite part $l = i+k+1$, which
    // reaches $\text{degree}+1$ at the deepest
    let f = outer_beta_entire_factor(r_own + 1.0, n, m_own, degree + 1);
    let (c_row_own, c_row_other) = (binomials(m_own), binomials(m_other));
    let factorial = |p: usize| Gamma::gamma(p as f64 + 1.0);

    let mut out = Polynomial::zeros(degree + 1);
    for (j, &c1) in c_row_own.iter().enumerate() {
        let c_row_inner = binomials(j);
        for (k, &c2) in c_row_other.iter().enumerate() {
            for (i, &c_inner) in c_row_inner.iter().enumerate() {
                // $(\partial_b-\partial_\delta)^j(-\partial_\delta)^k$ reaches
                // $\partial_b^{j-i}\partial_\delta^{i+k}$, whose pole sits at order $p$
                let p = i + k + 1;
                let scale = c1 * c2 * c_inner * factorial(i + k);
                out[degree - j - k + p] += scale * alternating_sign(p) * f[j - i][0] / factorial(p);
                out[degree - j - k] += scale * alternating_sign(i + k) * f[j - i][p];
            }
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

    // A singular part with no terms is nothing to integrate, and pairing one would
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
    /// the one frequency it diverges at having no weight under the integral.
    fn value(&self, nu: f64) -> f64 {
        let argument = self.argument(nu);
        if !self.support.contains(argument) {
            return 0.0;
        }
        let value = self.sing.value(argument);
        if value.is_finite() { value } else { 0.0 }
    }

    /// The singularity's terms with the side its support does not reach zeroed.
    fn sided_terms(&self) -> Vec<AsymptTerm> {
        let position = self.sing.position();
        let reaches_below = position > self.support.min();
        let reaches_above = position < self.support.max();
        let mut out = Vec::with_capacity(self.sing.terms().len());
        for t in self.sing.terms() {
            out.push(t.sides_kept(reaches_below, reaches_above));
        }
        out
    }

    /// $\int S$ over a stretch of the $\nu$ axis, in closed form.
    ///
    /// A reflected stretch is measured off the point it holds and mirrored about
    /// $\Omega_q$, which puts it in the singularity's own variable.
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

/// $\int S(\nu) g(\nu) d\nu$ over `segment`, with the singularity subtracted out.
///
/// Worth nothing unless `g` is continuous at the singular point, which confines this to a
/// stretch holding one divergence and no more.
fn subtracted<G: Fn(f64) -> f64>(part: &Part, g: G, segment: Segment, tol: f64) -> f64 {
    // A refused request is worth zero here, the pair being one contribution among many
    integrate_by_subtraction(
        part.at(),
        |nu| part.value(nu),
        |stretch| part.integral(stretch),
        g,
        segment,
        tol,
    )
    .map_or(0.0, |integral| integral.value)
}

/// $\int S_p(\nu) S_q(\omega-\nu) d\nu$ over `segment`, without a quadrature.
///
fn singular_singular(first: &Part, second: &Part, segment: Segment) -> f64 {
    // Measured from the first singular point, $x = \nu - \Omega_p$, the second factor is
    // met at $\Delta - x$ where $\Delta$ is the distance between the two points along
    // the integration axis. The integrand changes form where either factor turns, at
    // $x = 0$ and $x = \Delta$, so the range is cut there into at most three stretches
    // and each is a Beta function. A negative $\Delta$ is the same problem reflected:
    // $x \mapsto -x$ takes it to a positive one with both factors' sides exchanged.
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

    // A flip exchanges the two sides of every term, the integration variable having been
    // turned to face right
    let sided = |c_below: f64, c_above: f64| {
        if flipped {
            (c_above, c_below)
        } else {
            (c_below, c_above)
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
        let r1 = t1.exponent();
        for t2 in &terms2 {
            let r2 = t2.exponent();
            // The exponents alone fix $\rho$, so every pair of powers below shares it
            let rho = r1 + r2 + 1.0;
            for (m1, cb1, ca1) in t1.log_coeffs() {
                for (m2, cb2, ca2) in t2.log_coeffs() {
                    let (below1, above1) = sided(cb1, ca1);
                    let (below2, above2) = sided(cb2, ca2);
                    let mut take = |weight: f64, stretch: StretchContrib| {
                        total += weight
                            * (contract_with_ln_m(&stretch.singular, rho, delta) + stretch.regular);
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
        }
    }
    total
}

/// The part of `reach` lying within `region`, absent where the two share no length.
fn stretch_within(reach: Segment, region: Segment) -> Option<Segment> {
    reach.intersection(region).filter(|s| !s.is_degenerate())
}

/// What a stretch comes to.
///
/// The singular part is $\sum_m c_m |\Delta|^{\rho}\ln^m|\Delta|$ and is what reaches
/// $c^{\pm}$; the regular part is what is left at this $\Delta$, analytic through
/// $\Delta = 0$. Their sum is the integral over the stretch.
struct StretchContrib {
    /// $c_m$, indexed by the power of $\ln|\Delta|$.
    singular: Polynomial,
    /// Everything else, evaluated at this $\Delta$.
    regular: f64,
}

impl StretchContrib {
    /// A stretch whose beta could not be divided, all of it counted as regular.
    fn undivided(degree: usize, value: f64) -> StretchContrib {
        StretchContrib {
            singular: Polynomial::zeros(degree),
            regular: value,
        }
    }
}

/// $\int_s s^{r_1}\ln^{m_1}\\!s\;(\Delta+s)^{r_2}\ln^{m_2}\\!(\Delta+s)\\,ds$, over a stretch
/// `s` of the positive axis and for $\Delta > 0$.
///
/// The stretch running away from the pair of singular points rather than between them.
fn outer_stretch(r1: f64, m1: usize, r2: f64, m2: usize, delta: f64, s: Segment) -> StretchContrib {
    // Scaling by $s = \Delta u$ is the whole substitution, and it leaves the outer Beta
    // at $a = r_2+1$, $b = r_1+1 > 0$ and $U = s/\Delta$, where the exponent
    // $\rho = r_1 + r_2 + 1$ the pair produces is $a+b-1$ and is asked for nowhere else
    let rho = r1 + r2 + 1.0;
    let (a, b) = (r2 + 1.0, r1 + 1.0);
    let u = |t: f64| t / delta;
    let degree = m1 + m2;

    // Writing $B^{out}_U = B^{out}_\infty - T_U$, the stretch is
    // $\Delta^\rho(T_{U_1} - T_{U_2})$: the complete value cancels between the ends
    // unless the near one sits on the singular point, where $U_1 = 0$ makes $T_{U_1}$
    // the whole of it.
    let Some((complete, far)) = outer_beta_derivatives_split(a, b, u(s.max()), m2, m1) else {
        let upper = outer_beta_derivatives(a, b, u(s.max()), m2, m1);
        let lower = outer_beta_derivatives(a, b, u(s.min()), m2, m1);
        let step = subtract_tables(&upper, &lower);
        let coefficients = derivatives_to_log_coeffs(&step, m2, m1);
        return StretchContrib::undivided(degree, contract_with_ln_m(&coefficients, rho, delta));
    };
    let (singular, rest) = if s.min() == 0.0 {
        (
            derivatives_to_log_coeffs(&complete, m2, m1),
            negate_table(&far),
        )
    } else {
        let near = outer_beta_tail(a, b, u(s.min()), &complete, m2, m1);
        (Polynomial::zeros(degree), subtract_tables(&near, &far))
    };
    let regular = contract_with_ln_m(&derivatives_to_log_coeffs(&rest, m2, m1), rho, delta);
    StretchContrib { singular, regular }
}

/// $\int_x x^{r_1}\ln^{m_1}\\!x\;(\Delta-x)^{r_2}\ln^{m_2}\\!(\Delta-x)\\,dx$, over a stretch
/// `x` within $[0, \Delta]$ and for $\Delta > 0$.
///
/// The stretch between the two singular points, and the one whose Beta never meets a
/// pole.
fn middle_stretch(
    r1: f64,
    m1: usize,
    r2: f64,
    m2: usize,
    delta: f64,
    x: Segment,
) -> StretchContrib {
    // $x = \Delta t$ is the whole substitution, and the two exponents stay apart:
    // $a = r_1+1$ and $b = r_2+1$ are both positive whatever $\rho$ comes to
    let rho = r1 + r2 + 1.0;
    let (a, b) = (r1 + 1.0, r2 + 1.0);
    let degree = m1 + m2;

    // Both ends of this stretch are singular points, so neither owns the complete beta:
    // written as $B_{t_2} - B_{t_1}$ it enters at the upper limit and written as
    // $D_{t_1} - D_{t_2}$ at the lower, the two being the same sum. It belongs to the
    // stretch when the stretch runs the whole way between the points and not otherwise,
    // and there the regular part is exactly zero rather than merely small.
    let Some((complete, far)) = inc_beta_derivatives_split(a, b, x.max() / delta, m1, m2) else {
        let upper = inc_beta_derivatives(a, b, x.max() / delta, m1, m2);
        let lower = inc_beta_derivatives(a, b, x.min() / delta, m1, m2);
        let step = subtract_tables(&upper, &lower);
        let coefficients = derivatives_to_log_coeffs(&step, m1, m2);
        return StretchContrib::undivided(degree, contract_with_ln_m(&coefficients, rho, delta));
    };
    let (singular, rest) = if x.min() == 0.0 && x.max() == delta {
        (
            derivatives_to_log_coeffs(&complete, m1, m2),
            negate_table(&far),
        )
    } else {
        let (_, near) = inc_beta_derivatives_split(a, b, x.min() / delta, m1, m2)
            .expect("the far end already split, so the near one does too");
        (Polynomial::zeros(degree), subtract_tables(&near, &far))
    };
    let regular = contract_with_ln_m(&derivatives_to_log_coeffs(&rest, m1, m2), rho, delta);
    StretchContrib { singular, regular }
}

/// $|\Delta|^{\rho}\sum_m c_m \ln^m|\Delta|$.
fn contract_with_ln_m(coefficients: &Polynomial, rho: f64, delta: f64) -> f64 {
    delta.powf(rho) * coefficients.eval(delta.ln())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SpectralFunction;
    use crate::models::*;
    use approx::assert_relative_eq;

    /// The one continuous contribution a model has.
    fn only(sf: &SpectralFunction) -> &dyn ContinuousSF {
        sf.continuous[0].0.as_ref()
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

    /// The van Hove logarithm of the square lattice, which falls out of two
    /// inverse square roots meeting at $\rho = 0$.
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

        // The band edges meet at $\rho = 0$ too, and what the pair leaves there is a
        // constant: the step the square lattice band edge really is, 1/(4πt). Only one
        // side of each edge lies inside the band, and the other reaches zero.
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

    /// A pair of singular parts leaves zero to the regular part $\Omega_p + \Omega_q$.
    #[test]
    fn a_pair_leaves_nothing_to_the_regular_part() {
        // Each case picks one singular part out of either factor, and names the
        // exponent $\rho = r_1 + r_2 + 1$ the two come to
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
            for t1 in &first.sided_terms() {
                let second = Part {
                    sing: sing2,
                    support: s2,
                    reflected: Some(p1 + p2),
                };
                for t2 in &second.sided_terms() {
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
        let term = |c_below: f64, c_above: f64| AsymptTerm::sided_log(r, 0, c_below, c_above);
        let claimed = coincident_constant(&term(1.0, 1.0), &term(1.0, 1.0), Segment::new(-l, l));

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

    /// Both parts of an outer stretch against a closed form, which names them apart.
    ///
    /// Expanding $\int_{\Delta}^{\lambda} u^{a-1}(u-\Delta)^{b-1}du$ about $\Delta = 0$
    /// gives
    /// $$
    ///     \Delta^{\rho}B(b,-\rho) + \frac{\lambda^{\rho}}{\rho}
    ///         - \frac{(b-1)\lambda^{\rho-1}}{\rho-1}\Delta + O(\Delta^2),
    /// $$
    /// at $\rho = a+b-1$: the first term is the singular part and the rest the regular
    /// one. Writing $u = \Delta + s$ puts it at $r_1 = b-1$ and $r_2 = a-1$ over a
    /// stretch reaching $\lambda - \Delta$.
    #[test]
    fn an_outer_stretch_against_a_closed_form() {
        let lambda = 1.0f64;
        for &(a, b) in &[(0.7f64, 0.6f64), (1.3, 0.25), (0.4, 1.9)] {
            let rho = a + b - 1.0;
            let singular = beta_derivatives(b, -rho, 0, 0)[0][0];
            for delta in [1e-3f64, 1e-4, 1e-5] {
                let stretch = Segment::new(0.0, lambda - delta);
                let got = outer_stretch(b - 1.0, 0, a - 1.0, 0, delta, stretch);

                // The singular part is the whole of what reaches |Δ|^ρ, and knows
                // nothing of Δ or of how far the stretch runs
                assert_eq!(got.singular.degree(), 0);
                assert_relative_eq!(got.singular[0], singular, max_relative = 1e-12);

                // The regular part is the background, exact at this Δ rather than a
                // series, so the two agree to the order the series was taken
                let background = lambda.powf(rho) / rho
                    - (b - 1.0) * lambda.powf(rho - 1.0) / (rho - 1.0) * delta;
                let slack = delta * delta;
                assert!(
                    (got.regular - background).abs() < slack,
                    "a={a} b={b} Δ={delta:.0e}: {} against {background}, off by more \
                     than {slack:.1e}",
                    got.regular
                );
            }
        }
    }

    /// The constant a pair leaves at $\rho = 0$ differs between the two sides, which is
    /// the whole reason an [`AsymptTerm`] may hold a bare constant per side.
    #[test]
    fn the_constant_a_pair_leaves_is_sided() {
        let term = |c_below: f64, c_above: f64| AsymptTerm::sided_log(-0.5, 0, c_below, c_above);
        let pair = Singularity::new(
            0.0,
            1.0,
            convolve_terms(
                &term(1.0, 2.0),
                &term(0.7, 0.3),
                Some(Segment::new(-1.0, 3.0)),
            ),
        );

        // Two inverse square roots meet at ρ = 0, so $S(\pm d) = C^\pm + L ln(d)$.
        // Two points a side separate the two.
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

    #[test]
    fn the_regular_part_is_smooth_enough_to_fit() {
        use crate::interp::InterpolatedSF;
        let cases = [
            // A chain against itself pairs two real singularities at every point
            (
                "chain against chain",
                chain(0.0, 1.0),
                chain(0.0, 1.0),
                1e-12f64,
            ),
            // Two boxes have no singular parts at all, and their triangle is fitted
            // exactly because a kink at a panel boundary separates two polynomials
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
                2e-12,
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

    /// A pair whose exponent lands just off a whole number is derived as though it were
    /// on it, and the convolution stays normalised.
    #[test]
    fn an_exponent_beside_a_whole_number() {
        // $\rho = 2r+1$, so $r \to 1/2$ walks $\rho$ up to 2 from below. Both operands
        // are normalised and their second moments add, neither of which the library is
        // told. Read through the generic formula the derived coefficient goes as
        // $1/(\rho-2)$ and the regular part takes the opposite, which costs a part in
        // $10^2$ by $r = 0.499999999$.
        for r in [0.4f64, 0.49, 0.4999, 0.49999, 0.4999999, 0.499999999, 0.5] {
            let a = crate::models::pseudogap(0.0, r, 1.0);
            let moment = |sf: &crate::SpectralFunction| {
                sf.integrate(|omega: f64| omega * omega, None)
                    .unwrap()
                    .value
            };
            let one = moment(&a);
            let c = a.conv(&a, None);
            assert_relative_eq!(
                c.integrate(|_| 1.0, None).unwrap().value,
                1.0,
                epsilon = 1e-8
            );
            assert_relative_eq!(moment(&c), 2.0 * one, max_relative = 1e-8);
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

    /// Two boxes convolve into a triangle.
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

    /// A chain against itself is the square lattice.
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

    /// Two inverse square roots colliding is the van Hove peak of the square lattice.
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
}
