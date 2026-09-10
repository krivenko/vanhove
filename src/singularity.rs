//! Closed-form description of an integrable singularity of a spectral function.

use std::cmp::Ordering;

use crate::util::PowKind;

/// One term of the singular part,
/// $c|\omega - \Omega_p|^r \ln^m|\omega - \Omega_p|$.
#[derive(Debug, Clone, Copy)]
pub struct AsymptTerm {
    /// Exponent $r > -1$.
    exponent: f64,
    /// Power $m \in \\{0, 1\\}$ of the logarithmic factor.
    log_power: u8,
    /// Coefficient $c$.
    c: f64,
    /// Precomputed dispatch for $|\omega - \Omega_p|^r$.
    pow: PowKind,
}

impl AsymptTerm {
    /// $c|\omega - \Omega_p|^r$.
    pub fn power(exponent: f64, c: f64) -> AsymptTerm {
        AsymptTerm::make(exponent, 0, c)
    }

    /// $-c\ln|\omega - \Omega_p|$, so that $c > 0$ describes a peak.
    pub fn log(c: f64) -> AsymptTerm {
        AsymptTerm::make(0.0, 1, -c)
    }

    /// Constant $c$.
    pub fn constant(c: f64) -> AsymptTerm {
        AsymptTerm::make(0.0, 0, c)
    }

    fn make(exponent: f64, log_power: u8, c: f64) -> AsymptTerm {
        assert!(exponent > -1.0, "asymptotics exponent must satisfy r > -1");
        assert!(
            log_power <= 1,
            "at most one logarithmic factor is supported"
        );
        AsymptTerm {
            // -0.0 would sort below +0.0 in `Strength`, and the two are the same exponent
            exponent: if exponent == 0.0 { 0.0 } else { exponent },
            log_power,
            c,
            pow: PowKind::of(exponent),
        }
    }

    /// Value of the term at $u = |\omega - \Omega_p|$, given $\ln u$.
    fn value(&self, u: f64, ln_u: f64) -> f64 {
        let v = self.c * self.pow.eval(u);
        if self.log_power == 0 { v } else { v * ln_u }
    }

    /// Integral of the term over one side of $\Omega_p$ of length `l`,
    /// $\int_0^l c u^r \ln^m u\\, du$.
    fn half_integral(&self, l: f64) -> f64 {
        // The side is empty when Ω_p sits at that end of the support
        if l == 0.0 {
            return 0.0;
        }
        let rp1 = self.exponent + 1.0;
        let lp = l.powf(rp1);
        if self.log_power == 0 {
            self.c * lp / rp1
        } else {
            self.c * lp * (l.ln() / rp1 - 1.0 / (rp1 * rp1))
        }
    }

    /// Strength of the term as $\omega \to \Omega_p$.
    fn strength(&self) -> Strength {
        Strength {
            exponent: self.exponent,
            log_power: self.log_power,
        }
    }
}

/// Rate at which a term of $S_p$ grows as $\omega \to \Omega_p$.
///
/// $|\omega-\Omega_p|^r$ outgrows $|\omega-\Omega_p|^{r'}$ for $r < r'$, and a
/// logarithmic factor breaks the tie between equal exponents.
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

/// Isolated integrable singularity of a continuous spectral function.
///
/// The singular part is given in closed form over the whole support as
/// $$
///     S_p(\omega) = \sum_k c_k |\omega-\Omega_p|^{r_k}
///         \ln^{m_k}|\omega-\Omega_p|.
/// $$
#[derive(Debug, Clone)]
pub struct Singularity {
    /// Position of the singular point, $\Omega_p$.
    pub position: f64,
    /// Terms of $S_p(\omega)$.
    terms: Box<[AsymptTerm]>,
    /// Whether any term carries a logarithmic factor.
    has_log: bool,
}

impl Singularity {
    pub fn new(position: f64, terms: Vec<AsymptTerm>) -> Singularity {
        Singularity {
            position,
            has_log: terms.iter().any(|t| t.log_power > 0),
            terms: terms.into_boxed_slice(),
        }
    }

    /// Whether $S_p$ vanishes identically and can be skipped altogether.
    pub fn is_trivial(&self) -> bool {
        self.terms.is_empty()
    }

    /// $S_p(\omega)$.
    ///
    /// Diverges at $\Omega_p$ unless every term stays bounded there.
    pub fn value(&self, omega: f64) -> f64 {
        let u = (omega - self.position).abs();
        let ln_u = if self.has_log { u.ln() } else { 0.0 };
        self.terms.iter().map(|t| t.value(u, ln_u)).sum()
    }

    /// $\int_{\omega_{min}}^{\omega_{max}} S_p(\omega)d\omega$ in closed form.
    pub fn integral(&self, omega_min: f64, omega_max: f64) -> f64 {
        if self.is_trivial() {
            return 0.0;
        }
        assert!(
            omega_min.is_finite() && omega_max.is_finite(),
            "a spectral function with singularities must have a bounded support"
        );
        let (below, above) = (self.position - omega_min, omega_max - self.position);
        debug_assert!(below >= 0.0 && above >= 0.0, "Ω_p lies outside the support");
        self.terms
            .iter()
            .map(|t| t.half_integral(below) + t.half_integral(above))
            .sum()
    }

    /// Limit of $S_p(\omega)$ at $\Omega_p$ with every divergent term dropped.
    ///
    /// Only the constant terms survive, the ones with $r > 0$ vanishing.
    pub fn finite_limit(&self) -> f64 {
        self.terms
            .iter()
            .filter(|t| t.exponent == 0.0 && t.log_power == 0)
            .map(|t| t.c)
            .sum()
    }

    /// Divergent terms at $\Omega_p$, each with the coefficient of the $\pm\infty$
    /// it tends to.
    pub fn divergences(&self) -> impl Iterator<Item = (Strength, f64)> + '_ {
        self.terms.iter().filter_map(|t| {
            let strength = t.strength();
            // c u^r tends to sign(c) ∞ for r < 0, while c ln u tends to -sign(c) ∞
            strength
                .is_divergent()
                .then(|| (strength, if t.log_power == 0 { t.c } else { -t.c }))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::{AsymptTerm, Singularity};
    use crate::util::bilby_integrate;
    use approx::assert_relative_eq;
    use std::cmp::Ordering;

    #[test]
    fn value() {
        // S(ω) = 2|ω-1|^{-1/2} - 3 ln|ω-1| + 0.5 + 4|ω-1|^{3/2}
        let sing = Singularity::new(
            1.0,
            vec![
                AsymptTerm::power(-0.5, 2.0),
                AsymptTerm::log(3.0),
                AsymptTerm::constant(0.5),
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
        let peak = Singularity::new(0.0, vec![AsymptTerm::log(-1.0)]);
        assert_eq!(peak.value(0.0), f64::NEG_INFINITY);
    }

    #[test]
    fn integral() {
        // ∫_0^l u^r du = l^{r+1}/(r+1) on a support with Ω_p at the lower end
        for r in [-0.5f64, 0.0, 0.5, 2.0] {
            let sing = Singularity::new(1.0, vec![AsymptTerm::power(r, 3.0)]);
            let l = 2.0f64;
            assert_relative_eq!(
                sing.integral(1.0, 1.0 + l),
                3.0 * l.powf(r + 1.0) / (r + 1.0),
                max_relative = 1e-14
            );
        }

        // ∫_0^l u^r ln u du = l^{r+1}[ln(l)/(r+1) - 1/(r+1)^2], written here for the
        // -c ln u convention of `AsymptTerm::log()`
        let sing = Singularity::new(0.0, vec![AsymptTerm::log(1.0)]);
        let l = 4.0f64;
        assert_relative_eq!(
            sing.integral(0.0, l),
            -l * (l.ln() - 1.0),
            max_relative = 1e-14
        );

        // An interior Ω_p contributes both sides
        assert_relative_eq!(
            sing.integral(-l, l),
            -2.0 * l * (l.ln() - 1.0),
            max_relative = 1e-14
        );

        // An empty term list integrates to zero over any support, unbounded included
        let trivial = Singularity::new(0.0, vec![]);
        assert!(trivial.is_trivial());
        assert_eq!(trivial.integral(f64::NEG_INFINITY, f64::INFINITY), 0.0);
    }

    #[test]
    fn integral_matches_quadrature() {
        // The closed form agrees with adaptive quadrature of S_p itself
        let sing = Singularity::new(
            0.5,
            vec![
                AsymptTerm::power(-0.5, 1.5),
                AsymptTerm::log(0.75),
                AsymptTerm::constant(-0.25),
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
            omega_min,
            omega_max,
            1e-8,
        )
        .unwrap()
        .value;
        assert_relative_eq!(
            sing.integral(omega_min, omega_max),
            quad,
            max_relative = 1e-6
        );
    }

    #[test]
    fn finite_limit() {
        // Only the constant terms survive at Ω_p
        let sing = Singularity::new(
            0.0,
            vec![
                AsymptTerm::power(-0.5, 2.0),
                AsymptTerm::constant(0.25),
                AsymptTerm::constant(0.5),
                AsymptTerm::power(1.5, 4.0),
                AsymptTerm::log(3.0),
            ],
        );
        assert_relative_eq!(sing.finite_limit(), 0.75, max_relative = 1e-14);

        // A square-root band edge leaves nothing behind
        let edge = Singularity::new(0.0, vec![AsymptTerm::power(0.5, 1.0)]);
        assert_eq!(edge.finite_limit(), 0.0);
    }

    #[test]
    fn divergences() {
        let sing = Singularity::new(
            0.0,
            vec![
                AsymptTerm::constant(0.25),
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
        let strength = |term: AsymptTerm| term.strength();
        let log = strength(AsymptTerm::log(1.0));
        let constant = strength(AsymptTerm::constant(1.0));
        let inv_sqrt = strength(AsymptTerm::power(-0.5, 1.0));
        let inv = strength(AsymptTerm::power(-0.9, 1.0));

        // A more negative exponent outgrows a less negative one, and both outgrow a
        // logarithm; a bounded term is weaker than every divergence
        assert_eq!(inv.order(&inv_sqrt), Ordering::Greater);
        assert_eq!(inv_sqrt.order(&log), Ordering::Greater);
        assert_eq!(log.order(&constant), Ordering::Greater);
        assert_eq!(log.order(&log), Ordering::Equal);

        // A constant built from -0.0 is the same strength as one built from +0.0
        assert_eq!(strength(AsymptTerm::power(-0.0, 1.0)), constant);
    }

    #[test]
    #[should_panic(expected = "asymptotics exponent must satisfy r > -1")]
    fn non_integrable_term() {
        let _ = AsymptTerm::power(-1.0, 1.0);
    }

    #[test]
    #[should_panic(expected = "must have a bounded support")]
    fn unbounded_support() {
        let sing = Singularity::new(0.0, vec![AsymptTerm::power(0.5, 1.0)]);
        let _ = sing.integral(f64::NEG_INFINITY, 1.0);
    }
}
