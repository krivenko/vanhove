//! Closed-form description of an integrable singularity of a spectral function.

use std::cmp::Ordering;

use crate::polynomial::Polynomial;
use crate::segment::Segment;
use crate::util::{PowKind, alternating_sign};

/// One term of the singular part, $|\Delta|^r P^{\pm}(\ln|\Delta|)$ in the distance
/// $\Delta = \omega - \Omega_p$ to the singular point.
///
/// One exponent makes one term, however many powers of the logarithm go with it.
#[derive(Debug, Clone)]
pub struct AsymptTerm {
    /// Exponent $r > -1$.
    exponent: f64,
    /// $P^-$, read for $\Delta < 0$.
    below: Polynomial,
    /// $P^+$, read for $\Delta \geq 0$.
    above: Polynomial,
    /// Precomputed dispatch for $|\Delta|^r$.
    pow: PowKind,
}

impl AsymptTerm {
    /// Exponent from which a term is no longer worth describing as singular.
    ///
    /// $|\Delta|^r$ has two continuous derivatives at $\Delta=0$ from $r = 2$ up.
    pub const MAX_EXPONENT: f64 = 2.0;

    /// Highest power of the logarithm a term may hold.
    ///
    /// [`Strength`] and [`power_log_integral()`] index the powers in a `u8`.
    const MAX_LOG_POWER: usize = u8::MAX as usize;

    /// $c |\Delta|^r$.
    pub fn power(exponent: f64, c: f64) -> AsymptTerm {
        AsymptTerm::sided_log(exponent, 0, c, c)
    }

    /// $-c\ln|\Delta|$, so that $c > 0$ describes a peak.
    pub fn log(c: f64) -> AsymptTerm {
        AsymptTerm::sided_log(0.0, 1, -c, -c)
    }

    /// $c^\pm |\Delta|^r \ln^m|\Delta|$, with a coefficient of its own on either side of
    /// $\Omega_p$.
    ///
    /// The sides may differ in any term, a bare constant included.
    pub fn sided_log(exponent: f64, log_power: u8, c_below: f64, c_above: f64) -> AsymptTerm {
        let mut below = Polynomial::zeros(usize::from(log_power));
        let mut above = Polynomial::zeros(usize::from(log_power));
        below[usize::from(log_power)] = c_below;
        above[usize::from(log_power)] = c_above;
        AsymptTerm::make(exponent, below, above)
    }

    /// $|\Delta|^r P^-(\ln|\Delta|)$ below $\Omega_p$ and $|\Delta|^r P^+(\ln|\Delta|)$ at
    /// and above it.
    pub(crate) fn sided_log_poly(
        exponent: f64,
        below: Polynomial,
        above: Polynomial,
    ) -> AsymptTerm {
        AsymptTerm::make(exponent, below, above)
    }

    fn make(exponent: f64, mut below: Polynomial, mut above: Polynomial) -> AsymptTerm {
        assert!(exponent > -1.0, "asymptotics exponent must satisfy r > -1");
        below.trim();
        above.trim();
        assert!(
            below.degree().max(above.degree()) <= AsymptTerm::MAX_LOG_POWER,
            "a term reaches at most the {}th power of the logarithm",
            AsymptTerm::MAX_LOG_POWER
        );
        AsymptTerm {
            // -0.0 would sort below +0.0 in `Strength`, and the two are the same exponent
            exponent: if exponent == 0.0 { 0.0 } else { exponent },
            below,
            above,
            pow: PowKind::of(exponent),
        }
    }

    /// Exponent $r$.
    pub(crate) fn exponent(&self) -> f64 {
        self.exponent
    }

    /// Whether the term weighs nothing on either side.
    pub(crate) fn is_zero(&self) -> bool {
        self.below.is_zero() && self.above.is_zero()
    }

    /// The same term with either side dropped where it is not to be kept.
    pub(crate) fn sides_kept(&self, below: bool, above: bool) -> AsymptTerm {
        let kept = |keep: bool, p: &Polynomial| {
            if keep {
                p.clone()
            } else {
                Polynomial::zeros(0)
            }
        };
        AsymptTerm::make(
            self.exponent,
            kept(below, &self.below),
            kept(above, &self.above),
        )
    }

    /// Highest power of the logarithm the term reaches on either side.
    pub(crate) fn log_degree(&self) -> usize {
        self.below.degree().max(self.above.degree())
    }

    /// Each power of the logarithm with non-zero weight, with its coefficients below
    /// and at or above $\Omega_p$.
    pub(crate) fn log_coeffs(&self) -> impl Iterator<Item = (usize, f64, f64)> {
        (0..=self.log_degree()).filter_map(|m| {
            let (c_below, c_above) = (self.below.coeff(m), self.above.coeff(m));
            (c_below != 0.0 || c_above != 0.0).then_some((m, c_below, c_above))
        })
    }

    /// Add `other`, which must share this term's exponent, into it.
    fn absorb(&mut self, other: &AsymptTerm) {
        debug_assert_eq!(self.exponent, other.exponent, "terms of one exponent merge");
        self.below += &other.below;
        self.above += &other.above;
        self.below.trim();
        self.above.trim();
    }

    /// $P^-$ where `below` and $P^+$ otherwise.
    fn side(&self, below: bool) -> &Polynomial {
        if below { &self.below } else { &self.above }
    }

    /// Value of the term at $\Delta = \omega - \Omega_p$, whose sign picks the side.
    fn value(&self, delta: f64) -> f64 {
        let below = delta < 0.0;
        let distance = delta.abs();
        if distance == 0.0 {
            return self.value_at_the_point(below);
        }
        let side = self.side(below);
        let power = self.pow.eval(distance);
        // $|\Delta|^r$ vanishing takes the term with it: no power of the logarithm
        // outgrows it
        if power == 0.0 {
            return 0.0;
        }
        power * side.eval(distance.ln())
    }

    /// Limit of the term at $\Omega_p$ itself, on the side given by `below`.
    fn value_at_the_point(&self, below: bool) -> f64 {
        // $\ln 0$ is never formed: the limit is read off the leading power instead. A
        // positive exponent takes the whole term to zero, and otherwise the leading
        // power decides, $\ln|\Delta|$ running to $-\infty$ and so flipping the sign
        // once per power of it.
        if self.exponent > 0.0 {
            return 0.0;
        }
        let side = self.side(below);
        let leading = side.coeff(side.degree());
        if leading == 0.0 {
            return 0.0;
        }
        if self.exponent == 0.0 && side.degree() == 0 {
            return leading;
        }
        f64::INFINITY * alternating_sign(side.degree()) * leading.signum()
    }

    /// Integral of the term over the side of $\Omega_p$ given by `below`, of length `l`,
    /// $\int_0^l x^r P^{\pm}(\ln x)\\, dx$.
    fn half_integral(&self, below: bool, l: f64) -> f64 {
        let side = self.side(below);
        let mut total = 0.0;
        for m in 0..=side.degree() {
            let c = side.coeff(m);
            if c != 0.0 {
                total += c * power_log_integral(self.exponent, m as u8, l);
            }
        }
        total
    }
}

/// $\int_0^l u^r \ln^m u\\, du$, for $r > -1$.
///
pub(crate) fn power_log_integral(exponent: f64, log_power: u8, l: f64) -> f64 {
    // Each power of the logarithm is integrated by parts against the one below it,
    // $I_m = (l^{r+1}\ln^m l - m I_{m-1})/(r+1)$

    // The side is empty when Ω_p sits at that end of the support
    if l == 0.0 {
        return 0.0;
    }
    let rp1 = exponent + 1.0;
    let lp = l.powf(rp1);
    let (ln_l, mut integral) = (l.ln(), lp / rp1);
    for k in 1..=i32::from(log_power) {
        integral = (lp * ln_l.powi(k) - f64::from(k) * integral) / rp1;
    }
    integral
}

/// Rate at which a term of $S_p$ grows as $\omega \to \Omega_p$.
///
/// $|\omega-\Omega_p|^r$ outgrows $|\omega-\Omega_p|^{r'}$ for $r < r'$, and at equal
/// exponents a higher power of $\ln|\omega-\Omega_p|$ outgrows a lower one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Strength {
    exponent: f64,
    log_power: u8,
}

impl Strength {
    /// Whether the term diverges rather than staying bounded.
    fn is_divergent(&self) -> bool {
        self.exponent < 0.0 || (self.exponent == 0.0 && self.log_power > 0)
    }

    /// Order the strengths, a faster growing term comparing greater.
    pub fn order(&self, other: &Strength) -> Ordering {
        other
            .exponent
            .total_cmp(&self.exponent)
            .then(self.log_power.cmp(&other.log_power))
    }
}

/// `term` rescaled from the dimensionless variable $u = |\Delta|/s$ into $|\Delta|$ itself.
fn fold_scale(term: &AsymptTerm, s: f64) -> AsymptTerm {
    // With $\ln u = \ln|\Delta| - \ln s$, a polynomial in $\ln u$ is the same polynomial
    // in $\ln|\Delta|$ read about an origin displaced by $\ln s$, and the $s^{-r}$ in
    // $u^r = s^{-r}|\Delta|^r$ scales the whole of it. At $s = 1$ both are the identity.
    if s == 1.0 {
        return term.clone();
    }
    let (ln_s, scale) = (s.ln(), s.powf(-term.exponent));
    let mut below = term.below.shifted(ln_s);
    let mut above = term.above.shifted(ln_s);
    below *= scale;
    above *= scale;
    AsymptTerm::make(term.exponent, below, above)
}

/// Isolated integrable singularity of a continuous spectral function.
///
/// The singular part is given in closed form over the whole support as
/// $$
///     S_p(\omega) = \sum_k |\Delta|^{r_k} P^{\pm}_k(\ln|\Delta|), \qquad
///     \Delta = \omega - \Omega_p.
/// $$
#[derive(Debug, Clone)]
pub struct Singularity {
    /// Position of the singular point, $\Omega_p$.
    position: f64,
    /// Terms of $S_p(\omega)$.
    terms: Box<[AsymptTerm]>,
}

impl Singularity {
    /// `scale` is the positive length used as the energy unit, in which `terms` are written.
    pub fn new(position: f64, scale: f64, terms: Vec<AsymptTerm>) -> Singularity {
        assert!(
            position.is_finite(),
            "a singularity must sit at a finite frequency"
        );
        assert!(
            scale > 0.0 && scale.is_finite(),
            "singularity scale must be positive and finite"
        );
        // The scale is notation for whoever writes the terms down, and nothing past this
        // point asks for it. Terms of one exponent are then one term, their polynomials
        // added.
        let mut merged: Vec<AsymptTerm> = Vec::with_capacity(terms.len());
        for t in &terms {
            let folded = fold_scale(t, scale);
            match merged.iter_mut().find(|m| m.exponent == folded.exponent) {
                Some(m) => m.absorb(&folded),
                None => merged.push(folded),
            }
        }
        Singularity {
            position,
            terms: merged.into_boxed_slice(),
        }
    }

    /// Terms of $S_p(\omega)$, one element per exponent.
    pub(crate) fn terms(&self) -> &[AsymptTerm] {
        &self.terms
    }

    /// Position of the singular point, $\Omega_p$.
    pub fn position(&self) -> f64 {
        self.position
    }

    /// Whether $S_p$ vanishes identically and can be skipped altogether.
    pub fn is_trivial(&self) -> bool {
        self.terms.is_empty()
    }

    /// The same singularity displaced in frequency by `by`.
    pub fn shifted(&self, by: f64) -> Singularity {
        Singularity {
            position: self.position + by,
            ..self.clone()
        }
    }

    /// $S_p(\omega)$.
    ///
    /// Diverges at $\Omega_p$ unless every term stays bounded there.
    pub fn value(&self, omega: f64) -> f64 {
        let d = omega - self.position;
        self.terms.iter().map(|t| t.value(d)).sum()
    }

    /// $\int_{\omega_{min}}^{\omega_{max}} S_p(\omega)d\omega$ over `support`,
    /// in closed form.
    pub fn integral(&self, support: Segment) -> f64 {
        if self.is_trivial() {
            return 0.0;
        }
        assert!(
            support.is_bounded(),
            "a spectral function with singularities must have a bounded support"
        );
        debug_assert!(
            support.contains(self.position),
            "Ω_p lies outside the support"
        );
        let (lower, upper) = support.split_at(self.position);
        let (below, above) = (lower.length(), upper.length());
        self.terms
            .iter()
            .map(|t| t.half_integral(true, below) + t.half_integral(false, above))
            .sum()
    }

    /// Limit of $S_p(\omega)$ at $\Omega_p$ with every divergent term dropped, taken
    /// from above if the sides differ.
    ///
    /// A term of positive $r$ vanishes there, one of negative $r$ never stays finite,
    /// and of the rest only the bare constant survives.
    pub fn finite_limit(&self) -> f64 {
        let mut total = 0.0;
        for t in &self.terms {
            if t.exponent == 0.0 {
                total += t.side(false).coeff(0);
            }
        }
        total
    }

    /// Divergent terms at $\Omega_p$, each with the coefficient of the $\pm\infty$
    /// it tends to.
    pub fn divergences(&self) -> impl Iterator<Item = (Strength, f64)> {
        // One term holds a whole polynomial and each of its powers diverges at a rate of
        // its own, so a caller weighing two singularities against each other needs them
        // apart: the leading ones can cancel and leave a lower power standing.
        let mut found = Vec::new();
        for t in &self.terms {
            // Only the side at and above $\Omega_p$ is read, so only its own powers are
            // worth walking
            let above = t.side(false);
            for m in (0..=above.degree()).rev() {
                let strength = Strength {
                    exponent: t.exponent,
                    log_power: m as u8,
                };
                let c_above = above.coeff(m);
                // A power whose above side weighs nothing never reaches the limit
                if !strength.is_divergent() || c_above == 0.0 {
                    continue;
                }
                // ln|ω-Ω_p| is negative on the way in, so m of them flip the sign m times
                found.push((strength, alternating_sign(m) * c_above));
            }
        }
        found.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::{AsymptTerm, Singularity, Strength};
    use crate::segment::Segment;
    use crate::util::bilby_integrate;
    use approx::assert_relative_eq;
    use std::cmp::Ordering;

    #[test]
    fn value() {
        // S(ω) = 2|ω-1|^{-1/2} - 3 ln|ω-1| + 0.5 + 4|ω-1|^{3/2}
        let sing = Singularity::new(
            1.0,
            1.0,
            vec![
                AsymptTerm::power(-0.5, 2.0),
                AsymptTerm::log(3.0),
                AsymptTerm::power(0.0, 0.5),
                AsymptTerm::power(1.5, 4.0),
            ],
        );
        let reference = |omega: f64| {
            let u = (omega - 1.0f64).abs();
            2.0 / u.sqrt() - 3.0 * u.ln() + 0.5 + 4.0 * u.powf(1.5)
        };
        for omega in [-2.0, 0.5, 0.999, 1.001, 1.5, 4.0] {
            assert_relative_eq!(sing.value(omega), reference(omega), max_relative = 1e-14);
        }

        // The divergence is approached from either side
        assert_eq!(sing.value(1.0), f64::INFINITY);

        // A logarithm alone diverges downwards for a negative coefficient
        let peak = Singularity::new(0.0, 1.0, vec![AsymptTerm::log(-1.0)]);
        assert_eq!(peak.value(0.0), f64::NEG_INFINITY);
    }

    #[test]
    fn integral() {
        // ∫_0^l u^r du = l^{r+1}/(r+1) on a support with Ω_p at the lower end
        for r in [-0.5f64, 0.0, 0.5, 2.0] {
            let sing = Singularity::new(1.0, 1.0, vec![AsymptTerm::power(r, 3.0)]);
            let l = 2.0f64;
            assert_relative_eq!(
                sing.integral(Segment::new(1.0, 1.0 + l)),
                3.0 * l.powf(r + 1.0) / (r + 1.0),
                max_relative = 1e-14
            );
        }

        // ∫_0^l u^r ln u du = l^{r+1}[ln(l)/(r+1) - 1/(r+1)^2], written here for the
        // -c ln u convention of `AsymptTerm::log()`
        let sing = Singularity::new(0.0, 1.0, vec![AsymptTerm::log(1.0)]);
        let l = 4.0f64;
        assert_relative_eq!(
            sing.integral(Segment::new(0.0, l)),
            -l * (l.ln() - 1.0),
            max_relative = 1e-14
        );

        // An interior Ω_p contributes both sides
        assert_relative_eq!(
            sing.integral(Segment::new(-l, l)),
            -2.0 * l * (l.ln() - 1.0),
            max_relative = 1e-14
        );

        // An empty term list integrates to zero over any support, unbounded included
        let trivial = Singularity::new(0.0, 1.0, vec![]);
        assert!(trivial.is_trivial());
        let unbounded = Segment::new(f64::NEG_INFINITY, f64::INFINITY);
        assert_eq!(trivial.integral(unbounded), 0.0);
    }

    #[test]
    fn integral_matches_quadrature() {
        // The closed form agrees with adaptive quadrature of S_p itself
        let sing = Singularity::new(
            0.5,
            1.0,
            vec![
                AsymptTerm::power(-0.5, 1.5),
                AsymptTerm::log(0.75),
                AsymptTerm::power(0.0, -0.25),
                AsymptTerm::power(0.5, 2.0),
            ],
        );
        let (omega_min, omega_max) = (-1.0, 3.0);
        // S_p has no value at Ω_p, which the quadrature is free to sample
        let quad = bilby_integrate(
            |omega| {
                if omega == sing.position {
                    0.0
                } else {
                    sing.value(omega)
                }
            },
            Segment::new(omega_min, omega_max),
            1e-8,
        )
        .unwrap()
        .value;
        assert_relative_eq!(
            sing.integral(Segment::new(omega_min, omega_max)),
            quad,
            max_relative = 1e-6
        );
    }

    #[test]
    fn sided() {
        // 3u below Ω_p and 5u above it, with Ω_p at the origin and u = |ω|
        let sing = Singularity::new(0.0, 1.0, vec![AsymptTerm::sided_log(1.0, 0, 3.0, 5.0)]);
        for omega in [-2.0f64, -0.5, 0.5, 2.0] {
            let c = if omega < 0.0 { 3.0 } else { 5.0 };
            assert_relative_eq!(sing.value(omega), c * omega.abs(), max_relative = 1e-14);
        }

        // Each side of the integral takes its own coefficient:
        // ∫_{-2}^{0} 3|ω| dω + ∫_0^2 5ω dω = 6 + 10
        assert_relative_eq!(
            sing.integral(Segment::new(-2.0, 2.0)),
            16.0,
            max_relative = 1e-14
        );

        // The scale divides the distance on both sides alike
        let scaled = Singularity::new(0.0, 2.0, vec![AsymptTerm::sided_log(1.0, 0, 3.0, 5.0)]);
        assert_relative_eq!(scaled.value(-2.0), 3.0, max_relative = 1e-14);
        assert_relative_eq!(scaled.value(2.0), 5.0, max_relative = 1e-14);

        // A term vanishing at Ω_p leaves nothing behind and diverges nowhere
        assert_eq!(sing.finite_limit(), 0.0);
        assert_eq!(sing.divergences().count(), 0);

        // A symmetric term is the same read from either side
        let plain = Singularity::new(0.0, 1.0, vec![AsymptTerm::power(1.0, 3.0)]);
        assert_eq!(plain.value(-2.0), plain.value(2.0));
    }

    #[test]
    fn sided_constant() {
        // 3 below Ω_p and 5 at and above it. The step is legitimate: S_p is a piece of
        // the splitting, and what A(ω) does at Ω_p is decided by every piece together.
        let sing = Singularity::new(1.0, 1.0, vec![AsymptTerm::sided_log(0.0, 0, 3.0, 5.0)]);
        assert_eq!(sing.value(0.5), 3.0);
        assert_eq!(sing.value(1.5), 5.0);
        assert_eq!(sing.value(1.0), 5.0);

        // Each side of the integral takes its own coefficient
        assert_relative_eq!(
            sing.integral(Segment::new(-1.0, 4.0)),
            3.0 * 2.0 + 5.0 * 3.0,
            max_relative = 1e-14
        );
        assert_eq!(sing.finite_limit(), 5.0);
        assert_eq!(sing.divergences().count(), 0);
    }

    #[test]
    fn value_at_the_singular_point() {
        // c u^r ln^m u is zero at Ω_p for r > 0, the power vanishing faster than any
        // power of the logarithm diverges
        for m in 0u8..=3 {
            let sing = Singularity::new(2.0, 1.0, vec![AsymptTerm::sided_log(0.5, m, 3.0, 5.0)]);
            assert_eq!(sing.value(2.0), 0.0, "m = {m}");
            assert!(sing.value(2.0 + 1e-300).abs() < 1e-140);
        }

        // A term of no weight is nothing at Ω_p however fast u^r diverges
        for r in [-0.5f64, 0.0, 0.5] {
            for m in 0u8..=2 {
                let sing = Singularity::new(0.0, 1.0, vec![AsymptTerm::sided_log(r, m, 0.0, 0.0)]);
                assert_eq!(sing.value(0.0), 0.0, "r = {r}, m = {m}");
            }
        }

        // Value at Ω_p: a constant as it stands, and an infinity for anything divergent
        let constant = Singularity::new(0.0, 1.0, vec![AsymptTerm::power(0.0, 1.5)]);
        assert_eq!(constant.value(0.0), 1.5);
        let divergent = Singularity::new(0.0, 1.0, vec![AsymptTerm::power(-0.5, 2.0)]);
        assert_eq!(divergent.value(0.0), f64::INFINITY);
        let peak = Singularity::new(0.0, 1.0, vec![AsymptTerm::log(1.0)]);
        assert_eq!(peak.value(0.0), f64::INFINITY);
    }

    #[test]
    fn higher_logarithm_power() {
        // A term may have any power of the logarithm, and its integral follows the
        // recurrence instead of a closed form written out per power
        let sing = Singularity::new(0.0, 1.0, vec![AsymptTerm::sided_log(-0.5, 2, 1.0, 1.0)]);
        let u = 0.25f64;
        assert_relative_eq!(
            sing.value(u),
            u.powf(-0.5) * u.ln().powi(2),
            max_relative = 1e-14
        );

        // ∫_0^l u^r ln^2 u du against the quadrature it stands for
        let l = 0.75f64;
        let reference = bilby_integrate(
            |u: f64| {
                if u <= 0.0 {
                    0.0
                } else {
                    u.powf(-0.5) * u.ln().powi(2)
                }
            },
            Segment::new(0.0, l),
            1e-13,
        )
        .unwrap()
        .value;
        assert_relative_eq!(
            sing.integral(Segment::new(0.0, l)),
            reference,
            max_relative = 1e-9
        );
    }

    #[test]
    fn finite_limit_with_several_logarithms() {
        // c ln^m u = c (ln|Δ| - ln s)^m, so what survives Δ → 0 is c (-ln s)^m, and
        // only the m = 1 case looks like -c ln s
        let (s, c) = (3.0f64, 2.0f64);
        for m in 0u8..=3 {
            let sing = Singularity::new(0.0, s, vec![AsymptTerm::sided_log(0.0, m, c, c)]);
            assert_relative_eq!(
                sing.finite_limit(),
                c * (-s.ln()).powi(i32::from(m)),
                max_relative = 1e-14
            );
        }

        // A scale of one leaves nothing behind but the constant term
        let sing = Singularity::new(
            0.0,
            1.0,
            vec![
                AsymptTerm::sided_log(0.0, 0, c, c),
                AsymptTerm::sided_log(0.0, 2, 5.0, 5.0),
            ],
        );
        assert_relative_eq!(sing.finite_limit(), c, max_relative = 1e-14);
    }

    #[test]
    fn divergence_signs() {
        // c ln^m u tends to (-1)^m sign(c) ∞ for every m > 0
        for (r, m) in [
            (0.0f64, 1u8),
            (0.0, 2),
            (0.0, 3),
            (0.0, 4),
            (-0.5, 0),
            (-0.5, 1),
            (-0.5, 2),
            (-0.5, 3),
        ] {
            let sing = Singularity::new(0.0, 1.0, vec![AsymptTerm::sided_log(r, m, 2.0, 2.0)]);
            let reported = sing.divergences().next().unwrap().1;
            assert_eq!(
                reported.signum(),
                sing.value(1e-8).signum(),
                "r = {r}, m = {m} reported {reported}"
            );
        }
    }

    #[test]
    fn divergences_reach_below_the_leading_term() {
        // Folding the scale into c u^r ln u leaves two terms behind, c s^{-r}|Δ|^r ln|Δ|
        // and -c s^{-r} ln(s) |Δ|^r, and for r < 0 both of them diverge
        let (s, c, r) = (3.0f64, 2.0f64, -0.5f64);
        let sing = Singularity::new(0.0, s, vec![AsymptTerm::sided_log(r, 1, c, c)]);
        let found: Vec<(u8, f64)> = sing
            .divergences()
            .map(|(strength, x)| (strength.log_power, x))
            .collect();

        assert_eq!(found.len(), 2, "only the leading divergence was reported");
        let leading = c * s.powf(-r);
        assert_eq!(found[0].0, 1);
        assert_relative_eq!(found[0].1, -leading, max_relative = 1e-14);
        assert_eq!(found[1].0, 0);
        assert_relative_eq!(found[1].1, -leading * s.ln(), max_relative = 1e-14);

        // The leading sign is the one the term really takes
        assert_eq!(found[0].1.signum(), sing.value(1e-8).signum());

        // At a scale of one there is nothing below the leading term
        let unit = Singularity::new(0.0, 1.0, vec![AsymptTerm::sided_log(r, 1, c, c)]);
        assert_eq!(unit.divergences().count(), 1);
    }

    #[test]
    fn finite_limit() {
        // Only the constant terms survive at Ω_p
        let sing = Singularity::new(
            0.0,
            1.0,
            vec![
                AsymptTerm::power(-0.5, 2.0),
                AsymptTerm::power(0.0, 0.25),
                AsymptTerm::power(0.0, 0.5),
                AsymptTerm::power(1.5, 4.0),
                AsymptTerm::log(3.0),
            ],
        );
        assert_relative_eq!(sing.finite_limit(), 0.75, max_relative = 1e-14);

        // A square-root band edge leaves nothing behind
        let edge = Singularity::new(0.0, 1.0, vec![AsymptTerm::power(0.5, 1.0)]);
        assert_eq!(edge.finite_limit(), 0.0);
    }

    #[test]
    fn scale() {
        // S(ω) = 2u^{-1/2} - 3 ln u + 4u^{3/2} for u = |ω-1|/s
        let s = 8.0f64;
        let terms = || {
            vec![
                AsymptTerm::power(-0.5, 2.0),
                AsymptTerm::log(3.0),
                AsymptTerm::power(1.5, 4.0),
            ]
        };
        let sing = Singularity::new(1.0, s, terms());
        let reference = |omega: f64| {
            let u = (omega - 1.0f64).abs() / s;
            2.0 / u.sqrt() - 3.0 * u.ln() + 4.0 * u.powf(1.5)
        };
        for omega in [-2.0, 0.5, 1.001, 4.0, 20.0] {
            assert_relative_eq!(sing.value(omega), reference(omega), max_relative = 1e-14);
        }

        // The closed-form integral follows the substitution dω = s du
        let (omega_min, omega_max) = (-3.0, 9.0);
        let unit = Singularity::new(1.0, 1.0, terms());
        assert_relative_eq!(
            sing.integral(Segment::new(omega_min, omega_max)),
            s * unit.integral(Segment::new(
                1.0 + (omega_min - 1.0) / s,
                1.0 + (omega_max - 1.0) / s
            )),
            max_relative = 1e-14
        );

        // Rescaling a logarithm leaves c ln(s) behind, which is what survives at Ω_p
        // once the divergence is cancelled against another singularity
        assert_relative_eq!(sing.finite_limit(), 3.0 * s.ln(), max_relative = 1e-14);
        assert_eq!(unit.finite_limit(), 0.0);

        // A power law produces no such remainder whatever the scale
        let edge = Singularity::new(0.0, 5.0, vec![AsymptTerm::power(0.5, 1.0)]);
        assert_eq!(edge.finite_limit(), 0.0);
    }

    #[test]
    #[should_panic(expected = "singularity scale must be positive and finite")]
    fn non_positive_scale() {
        let _ = Singularity::new(0.0, 0.0, vec![AsymptTerm::log(1.0)]);
    }

    #[test]
    fn divergences() {
        let sing = Singularity::new(
            0.0,
            1.0,
            vec![
                AsymptTerm::power(0.0, 0.25),
                AsymptTerm::power(0.5, 4.0),
                AsymptTerm::log(3.0),
                AsymptTerm::power(-0.5, 2.0),
            ],
        );
        let divergences: Vec<_> = sing.divergences().collect();

        // The constant and the growing power are bounded, the other two are not
        assert_eq!(divergences.len(), 2);

        // -c ln u tends to sign(c) ∞, c u^{-a} to sign(c) ∞
        assert_eq!(divergences[0].1, 3.0);
        assert_eq!(divergences[1].1, 2.0);

        // A power law outgrows a logarithm
        assert_eq!(divergences[1].0.order(&divergences[0].0), Ordering::Greater);
    }

    #[test]
    fn strength_order() {
        let strength = |exponent, log_power| Strength {
            exponent,
            log_power,
        };
        let (log, constant) = (strength(0.0, 1), strength(0.0, 0));
        let (inv_sqrt, inv) = (strength(-0.5, 0), strength(-0.9, 0));

        // A more negative exponent outgrows a less negative one, and both outgrow a
        // logarithm; a bounded term is weaker than every divergence
        assert_eq!(inv.order(&inv_sqrt), Ordering::Greater);
        assert_eq!(inv_sqrt.order(&log), Ordering::Greater);
        assert_eq!(log.order(&constant), Ordering::Greater);
        assert_eq!(log.order(&log), Ordering::Equal);

        // A tie between equal exponents is broken by the power of the logarithm, and no
        // number of logarithms makes up for an exponent that is less negative. Listed
        // from the weakest, every pair has to order the way the list does.
        let ascending = [
            strength(0.5, 0),
            strength(0.5, 2),
            strength(0.0, 0),
            strength(0.0, 1),
            strength(0.0, 3),
            strength(-0.5, 0),
            strength(-0.5, 2),
            strength(-0.5, 7),
            strength(-0.9, 0),
            strength(-0.9, 4),
        ];
        for (i, weaker) in ascending.iter().enumerate() {
            for (j, stronger) in ascending.iter().enumerate() {
                assert_eq!(
                    weaker.order(stronger),
                    i.cmp(&j),
                    "{weaker:?} vs {stronger:?}"
                );
            }
        }

        // `order()` compares exponents with `total_cmp`, which sorts -0.0 below +0.0.
        // A term built from either has to reach the same strength all the same, which is
        // what the exponent is normalized for on the way in.
        let from = |exponent| Strength {
            exponent: Singularity::new(0.0, 1.0, vec![AsymptTerm::power(exponent, 1.0)]).terms()[0]
                .exponent,
            log_power: 0,
        };
        assert_eq!(from(-0.0).order(&from(0.0)), Ordering::Equal);
        assert_eq!(from(-0.0), constant);
    }

    #[test]
    #[should_panic(expected = "asymptotics exponent must satisfy r > -1")]
    fn non_integrable_term() {
        let _ = AsymptTerm::power(-1.0, 1.0);
    }

    #[test]
    #[should_panic(expected = "must have a bounded support")]
    fn unbounded_support() {
        let sing = Singularity::new(0.0, 1.0, vec![AsymptTerm::power(0.5, 1.0)]);
        let _ = sing.integral(Segment::new(f64::NEG_INFINITY, 1.0));
    }
}
