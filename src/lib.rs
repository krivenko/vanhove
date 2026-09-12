#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod conv;
pub mod discrete;
pub mod interp;
pub mod models;
pub mod singularity;
mod util;

use std::f64::consts::PI;
use std::ops::{Add, Mul, Neg, Sub};
use std::sync::Arc;

use bilby::QuadratureError;
use num_complex::Complex64;

use crate::conv::Convolution;
use crate::discrete::DiscreteSF;
use crate::interp::Interpolated;
use crate::singularity::{Singularity, Strength};

//
// ContinuousSF
//

/// Continuous spectral function possibly containing integrable singularities.
///
/// It has the form $A(\omega) = R(\omega) + \sum_p S_p(\omega)$ for
/// $\omega \in [\omega_{min}, \omega_{max}]$ and zero otherwise.
/// $R(\omega)$ is a smooth function and each $S_p(\omega)$ has one isolated
/// integrable singularity at $\Omega_p$.
///
/// $S_p(\omega)$ is described in closed form by the corresponding [`singularity::Singularity`],
/// which fixes it over the whole support and not merely near $\Omega_p$. A support
/// carrying singularities must be bounded.
///
/// Implementations are shared behind an [`Arc`], so they must be [`Send`] and
/// [`Sync`]. A spectral function is read-only once built, which every model
/// satisfies by holding nothing but plain data.
/// # Example
///
/// A band rising as $\frac{3}{2}\sqrt{\omega-\epsilon}$ out of its lower edge and
/// cut off a unit of frequency above it.
///
/// ```
/// use vanhove::singularity::{AsymptTerm, Singularity};
/// use vanhove::{ContinuousSF, SpectralFunction};
///
/// struct SqrtBand {
///     eps: f64,
///     edge: [Singularity; 1],
/// }
///
/// impl SqrtBand {
///     fn new(eps: f64) -> SqrtBand {
///         let terms = vec![AsymptTerm::power(0.5, 1.5)];
///         SqrtBand { eps, edge: [Singularity::new(eps, 1.0, terms)] }
///     }
/// }
///
/// impl ContinuousSF for SqrtBand {
///     fn support(&self) -> (f64, f64) {
///         (self.eps, self.eps + 1.0)
///     }
///     // The whole of A(ω) is the singular part, so nothing is left over
///     fn regular(&self, _omega: f64) -> f64 {
///         0.0
///     }
///     fn singularities(&self) -> &[Singularity] {
///         &self.edge
///     }
///     fn shifted(&self, by: f64) -> Box<dyn ContinuousSF> {
///         Box::new(SqrtBand::new(self.eps + by))
///     }
/// }
///
/// let dos = SpectralFunction::from_continuous(SqrtBand::new(0.0));
///
/// // The band is normalized, and its first moment is 3/5
/// assert!((dos.integrate(|_| 1.0, None).unwrap() - 1.0).abs() < 1e-12);
/// assert!((dos.integrate(|omega| omega, None).unwrap() - 0.6).abs() < 1e-12);
/// ```
pub trait ContinuousSF: Send + Sync {
    /// Support of the spectral function specified as a segment
    /// $[\omega_{min}, \omega_{max}]$.
    fn support(&self) -> (f64, f64);
    /// Regular part, $R(\omega)$.
    fn regular(&self, omega: f64) -> f64;
    /// Singular points $\Omega_p$ along with the closed form of $S_p$ at each.
    fn singularities(&self) -> &[Singularity] {
        &[]
    }
    /// Frequencies where $R(\omega)$ is not smooth, beyond the singular points.
    fn breakpoints(&self) -> &[f64] {
        &[]
    }
    /// The same spectral function displaced in frequency by `by`.
    fn shifted(&self, by: f64) -> Box<dyn ContinuousSF>;
}

/// Spectral function as a weighted sum of discrete resonances and continuous
/// contributions.
///
/// Weights may be negative, so that $A(\omega)$ is not assumed to be sign-definite.
#[derive(Clone)]
pub struct SpectralFunction {
    // Contributions of discrete resonances.
    discrete: DiscreteSF,
    // Continuous contributions with their weights.
    // Invariant: all weights are non-zero and no two entries share a contribution.
    continuous: Vec<(Arc<dyn ContinuousSF>, f64)>,
}

/// Multiply spectral function by a real number from the right.
impl Mul<f64> for SpectralFunction {
    type Output = Self;

    fn mul(self, a: f64) -> Self {
        // A weight can vanish upon scaling, either exactly for a == 0 or by underflow.
        // Such contributions are dropped by `from_discrete_continuous()`.
        Self::from_discrete_continuous(
            self.discrete * a,
            self.continuous
                .into_iter()
                .map(|(cd, w)| (cd, w * a))
                .collect(),
        )
    }
}
/// Multiply spectral function by a real number from the left.
impl Mul<SpectralFunction> for f64 {
    type Output = SpectralFunction;
    fn mul(self, sf: SpectralFunction) -> SpectralFunction {
        sf * self
    }
}
/// Negation of a spectral function.
impl Neg for SpectralFunction {
    type Output = Self;
    fn neg(self) -> SpectralFunction {
        self * (-1.0)
    }
}
/// Addition of two spectral functions.
impl Add for SpectralFunction {
    type Output = Self;
    fn add(self, rhs: SpectralFunction) -> SpectralFunction {
        // Concatenation can repeat a contribution the operands have in common, which
        // `from_discrete_continuous()` merges back into a single weight.
        let mut continuous = self.continuous;
        continuous.extend(rhs.continuous);
        SpectralFunction::from_discrete_continuous(self.discrete + rhs.discrete, continuous)
    }
}
/// Subtraction of two spectral functions.
impl Sub for SpectralFunction {
    type Output = Self;
    fn sub(self, rhs: SpectralFunction) -> SpectralFunction {
        self + (-rhs)
    }
}
/// Frequencies where a continuous spectral function is not smooth.
fn non_smooth(csf: &dyn ContinuousSF) -> impl Iterator<Item = f64> + '_ {
    csf.singularities()
        .iter()
        .map(|s| s.position())
        .chain(csf.breakpoints().iter().copied())
}

/// Continuous contributions of the convolution of `discrete` with the continuous part
/// of `sf`.
fn conv_discrete_continuous(
    discrete: &DiscreteSF,
    sf: &SpectralFunction,
) -> Vec<(Arc<dyn ContinuousSF>, f64)> {
    let mut contributions = Vec::with_capacity(discrete.len() * sf.continuous.len());
    // A resonance displaces every band to its position and scales it by its weight
    for (csf, w) in &sf.continuous {
        for res in discrete.iter() {
            contributions.push((Arc::from(csf.shifted(res.eps)), w * res.weight));
        }
    }
    contributions
}

impl SpectralFunction {
    /// Build a `SpectralFunction` from a discrete spectral function and a list of
    /// continuous contributions with their weights.
    ///
    /// Repeated contributions are merged by summing their weights, keeping the order
    /// of first appearance. A contribution is dropped once its total weight is small
    /// enough against the total magnitude of the weights to be a cancellation
    /// artefact, by the same relative tolerance that governs the discrete weights.
    ///
    /// Two contributions count as one only when they share an allocation, so that
    /// separately built models of identical parameters stay apart.
    fn from_discrete_continuous(
        dsf: DiscreteSF,
        csf: Vec<(Arc<dyn ContinuousSF>, f64)>,
    ) -> SpectralFunction {
        // Contributions grouped by identity: (contribution, ∑ w, ∑ |w|)
        let mut merged: Vec<(Arc<dyn ContinuousSF>, f64, f64)> = Vec::with_capacity(csf.len());
        for (cd, w) in csf {
            match merged
                .iter_mut()
                .find(|(other, _, _)| Arc::ptr_eq(other, &cd))
            {
                Some((_, sum, magnitude)) => {
                    *sum += w;
                    *magnitude += w.abs();
                }
                None => merged.push((cd, w, w.abs())),
            }
        }
        SpectralFunction {
            discrete: dsf,
            continuous: merged
                .into_iter()
                .filter(|(_, sum, magnitude)| sum.abs() > DiscreteSF::WEIGHT_TOL * magnitude)
                .map(|(cd, sum, _)| (cd, sum))
                .collect(),
        }
    }

    /// Build a `SpectralFunction` out of a single continuous contribution of unit weight.
    pub fn from_continuous<C: ContinuousSF + 'static>(csf: C) -> SpectralFunction {
        SpectralFunction::from_discrete_continuous(DiscreteSF::new(), vec![(Arc::new(csf), 1.0)])
    }

    /// Discrete part of the spectral function.
    pub fn discrete(&self) -> &DiscreteSF {
        &self.discrete
    }

    /// Support of the spectral function.
    ///
    /// It is the smallest segment $[\omega_{min}, \omega_{max}]$ containing the supports
    /// of all contributions. The segment is not necessarily tight: the spectral function
    /// may vanish within the gaps between disjoint contributions. Returns [`None`] for an
    /// empty spectral function.
    pub fn support(&self) -> Option<(f64, f64)> {
        self.continuous
            .iter()
            .map(|(cd, _)| cd.support())
            .chain(self.discrete.support())
            .reduce(|hull, sup| (hull.0.min(sup.0), hull.1.max(sup.1)))
    }

    /// Replace the regular part of every continuous contribution with a Chebyshev
    /// interpolation of it.
    ///
    /// The singular parts are carried over exactly, so the result departs from the
    /// original in $R(\omega)$ alone, and only within the relative tolerance `tol` of
    /// the fit. It pays for itself when one spectral function is integrated many times
    /// over and its regular part is expensive, as for the lattice models built on an
    /// elliptic integral.
    ///
    /// Contributions of unbounded support are left as they are, there being no
    /// interval to expand them over.
    pub fn precomputed(&self, tol: Option<f64>) -> SpectralFunction {
        let continuous = self
            .continuous
            .iter()
            .map(|(csf, w)| {
                let support = csf.support();
                if support.0.is_finite() && support.1.is_finite() {
                    let interpolated = Interpolated::new(csf.as_ref(), tol);
                    (Arc::new(interpolated) as Arc<dyn ContinuousSF>, *w)
                } else {
                    (Arc::clone(csf), *w)
                }
            })
            .collect();
        SpectralFunction::from_discrete_continuous(self.discrete.clone(), continuous)
    }

    /// Convolution with another spectral function,
    /// $\int A(\nu) B(\omega - \nu) d\nu$.
    ///
    /// Two continuous parts meet to within some $10^{-9}$ of the peak of $A(\omega)$:
    /// the frequencies where the result stops being smooth are derived, the asymptotics
    /// there are not, and the regular part is interpolated across them.
    pub fn conv(&self, other: &SpectralFunction) -> SpectralFunction {
        let derived = crate::conv::singularities(self, other);
        self.conv_with(other, derived, None)
    }

    /// Convolution with another spectral function, given the singular structure
    /// $C_A \ast C_B$ is known to have.
    ///
    /// `singularities` describes the convolution of the two continuous parts alone,
    /// the other three terms carrying their own. The frequencies where the result
    /// stops being smooth are derived either way.
    pub fn conv_with(
        &self,
        other: &SpectralFunction,
        singularities: Vec<Singularity>,
        tol: Option<f64>,
    ) -> SpectralFunction {
        // Every resonance of one operand against the continuous part of the other
        let mut continuous = conv_discrete_continuous(&self.discrete, other);
        continuous.extend(conv_discrete_continuous(&other.discrete, self));

        // The two continuous parts against each other, evaluated by quadrature and
        // interpolated over the singular structure given
        if !self.continuous.is_empty() && !other.continuous.is_empty() {
            let tol = tol.unwrap_or(Interpolated::DEFAULT_TOL);
            let breakpoints = crate::conv::breakpoints(self, other);
            let convolution = Convolution::new(self, other, singularities, breakpoints, tol);
            let interpolated = Interpolated::new(&convolution, Some(tol));
            continuous.push((Arc::new(interpolated) as Arc<dyn ContinuousSF>, 1.0));
        }

        SpectralFunction::from_discrete_continuous(self.discrete.conv(&other.discrete), continuous)
    }

    /// Total spectral weight.
    pub fn total_weight(&self) -> f64 {
        self.discrete.total_weight() + self.continuous.iter().map(|(_, w)| w).sum::<f64>()
    }

    /// Value of the continuous part of the spectral function at a frequency `omega`.
    ///
    /// Returns $\pm\infty$ where the spectral function diverges, and zero outside of
    /// the support of every continuous contribution.
    ///
    /// The discrete resonances are left out, a $\delta$-function having no value at a
    /// point. Use [`SpectralFunction::discrete()`] to inspect them.
    pub fn continuous_at(&self, omega: f64) -> f64 {
        // Divergent terms encountered at omega, grouped by strength:
        // (strength, ∑ w c, ∑ |w c|).
        let mut divergent: Vec<(Strength, f64, f64)> = Vec::new();
        // Everything that stays finite at omega
        let mut finite = 0.0f64;

        for (csf, w) in &self.continuous {
            let (omega_min, omega_max) = csf.support();
            // R(ω) is not defined outside of the support
            if omega < omega_min || omega > omega_max {
                continue;
            }

            let mut value = csf.regular(omega);
            for sing in csf.singularities() {
                // Away from Ω_p the asymptotics is finite and needs no analysis
                if sing.position() != omega {
                    value += sing.value(omega);
                    continue;
                }
                value += sing.finite_limit();
                for (strength, c) in sing.divergences() {
                    let coeff = w * c;
                    match divergent.iter_mut().find(|(s, _, _)| *s == strength) {
                        Some((_, sum, magnitude)) => {
                            *sum += coeff;
                            *magnitude += coeff.abs();
                        }
                        None => divergent.push((strength, coeff, coeff.abs())),
                    }
                }
            }
            finite += w * value;
        }

        // The strongest divergence whose coefficients do not cancel fixes the value. The
        // same relative tolerance decides a cancellation here as for the discrete weights.
        match divergent
            .iter()
            .filter(|(_, sum, magnitude)| sum.abs() > DiscreteSF::WEIGHT_TOL * magnitude)
            .max_by(|x, y| x.0.order(&y.0))
        {
            Some((_, sum, _)) => sum.signum() * f64::INFINITY,
            // Every divergence has cancelled, leaving the finite limits behind
            None => finite,
        }
    }

    //
    // Spectral function integration
    //

    /// Integrate the spectral function $A(\omega)$ against a real-valued function
    /// $f(\omega)$.
    ///
    /// The integral $\int A(\omega)f(\omega)d\omega$ is computed to the absolute
    /// tolerance `tol`, which defaults to $10^{-10}$.
    pub fn integrate<F: Fn(f64) -> f64>(
        &self,
        f: F,
        tol: Option<f64>,
    ) -> Result<f64, QuadratureError> {
        let mut result = 0.0f64;

        // Discrete spectral contributions
        for r in self.discrete.iter() {
            result += r.weight * f(r.eps);
        }

        // Continuous spectral contributions
        // ∫A(ω)f(ω)dω = ∫R(ω)f(ω)dω + ∑_p ∫S_p(ω)[f(ω) - f(Ω_p)]dω + ∑_p f(Ω_p) ∫S_p(ω)dω
        let tol = tol.unwrap_or(1e-10);
        for (csf, w) in &self.continuous {
            let mut res_contrib = 0.0f64;
            let (omega_min, omega_max) = csf.support();

            // Integrate the regular part, ∫R(ω)f(ω)dω
            res_contrib += util::bilby_integrate(
                |omega| csf.regular(omega) * f(omega),
                omega_min,
                omega_max,
                tol,
            )?
            .value;

            // Add integrals of the asymptotics
            for sing in csf.singularities() {
                if sing.is_trivial() {
                    continue;
                }
                // ∫S_p(ω)[f(ω) - f(Ω_p)]dω
                let omega_p = sing.position();
                let f_p = f(omega_p);
                res_contrib += util::bilby_integrate(
                    |omega| {
                        if omega == omega_p {
                            0.0
                        } else {
                            sing.value(omega) * (f(omega) - f_p)
                        }
                    },
                    omega_min,
                    omega_max,
                    tol,
                )?
                .value;
                // ∫S_p(ω)dω f(Ω_p)
                res_contrib += sing.integral(omega_min, omega_max) * f_p;
            }

            result += w * res_contrib;
        }
        Ok(result)
    }

    /// Integrate the spectral function $A(\omega)$ against a complex-valued function
    /// $f(\omega)$.
    ///
    /// The real and imaginary parts of $\int A(\omega)f(\omega)d\omega$ are computed
    /// separately, each to the absolute tolerance `tol`.
    pub fn integrate_complex<F: Fn(f64) -> Complex64>(
        &self,
        f: F,
        tol: Option<f64>,
    ) -> Result<Complex64, QuadratureError> {
        Ok(self.integrate(|omega| f(omega).re, tol)?
            + Complex64::I * self.integrate(|omega| f(omega).im, tol)?)
    }

    /// Value of the broadened spectral function at a frequency `omega`.
    ///
    /// The broadening is a convolution with the Lorentzian line shape function of the
    /// half-width at half-maximum `delta`,
    /// $$
    ///     \int \frac{1}{\pi} \frac{\delta}{\delta^2 + (\omega-\omega')^2}
    ///         A(\omega') d\omega'.
    /// $$
    /// The discrete part contributes a Lorentzian per resonance, so that the
    /// broadened spectral function is smooth in `omega` for any $\delta > 0$.
    ///
    /// The line shape function is sharply peaked for small `delta`, and its integral
    /// is of the order $1/(\pi\delta)$. It may, therefore, be necessary to relax the
    /// absolute quadrature tolerance `tol` from its default value of $10^{-10}$.
    pub fn broadened(
        &self,
        omega: f64,
        delta: f64,
        tol: Option<f64>,
    ) -> Result<f64, QuadratureError> {
        assert!(delta > 0.0, "broadening must be positive");
        let (weight, delta_sq) = (delta / PI, delta.powi(2));
        let f = |omega_prime: f64| -> f64 { weight / (delta_sq + (omega_prime - omega).powi(2)) };
        self.integrate(f, tol)
    }
}

#[cfg(test)]
mod tests {
    use crate::SpectralFunction;
    use crate::models::*;
    use approx::{assert_abs_diff_eq, assert_relative_eq};
    use std::f64::consts::PI;

    #[test]
    fn send_sync() {
        // Spectral functions cross thread boundaries, so that frequency scans can be
        // run in parallel. `ContinuousSF` requires as much of every model.
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<SpectralFunction>();
    }

    #[test]
    fn total_weight() {
        let dos = 2.0 * discrete(&[-0.7, 1.2], &[0.25, 0.6]) + 5.0 * gaussian(1.4, 0.5);
        assert_relative_eq!(dos.total_weight(), 6.7, epsilon = 1e-12);
    }

    #[test]
    fn discrete_part() {
        let dos = 2.0 * discrete(&[-0.7, 1.2], &[0.25, 0.6]) + 5.0 * gaussian(1.4, 0.5);
        let d = dos.discrete();
        assert_eq!(d.len(), 2);
        assert_relative_eq!(d.total_weight(), 1.7, epsilon = 1e-12);
        assert_relative_eq!(d.find(1.2).unwrap().weight, 1.2, epsilon = 1e-12);

        // Continuous contributions leave no trace in the discrete part
        assert!(d.find(1.4).is_none());
        assert!(gaussian(1.4, 0.5).discrete().is_empty());
    }

    #[test]
    fn support() {
        // An empty spectral function has no support
        assert_eq!(discrete(&[], &[]).support(), None);
        assert_eq!((gaussian(1.4, 0.5) * 0.0).support(), None);

        // Discrete contributions alone: positions of the outermost resonances
        assert_eq!(
            discrete(&[-0.7, 1.2, 0.3], &[0.25, 0.6, 0.15]).support(),
            Some((-0.7, 1.2))
        );

        // A single continuous contribution: its band edges
        assert_eq!(chain(0.5, 1.0).support(), Some((-1.5, 2.5)));

        // Mixed contributions: the hull of all supports. Neither the resonance at 0.0
        // nor the narrower square lattice band widens the chain band, the resonance
        // at 5.0 does.
        let dos = discrete(&[0.0, 5.0], &[0.3, 0.3]) + chain(0.5, 1.0) + square(0.0, 0.25);
        assert_eq!(dos.support(), Some((-1.5, 5.0)));

        // An unbounded contribution makes the whole support unbounded
        let dos = discrete(&[-0.7], &[0.5]) + gaussian(1.4, 0.5);
        assert_eq!(dos.support(), Some((f64::NEG_INFINITY, f64::INFINITY)));
    }

    #[test]
    fn continuous_at() {
        // Away from the singularities the value is that of A(ω) itself
        let (eps, t) = (0.5f64, 1.0f64);
        let dos = chain(eps, t);
        for omega in [-1.0, 0.0, 0.8, 2.0] {
            let ref_value = 1.0 / (PI * ((2.0 * t).powi(2) - (omega - eps).powi(2)).sqrt());
            assert_relative_eq!(dos.continuous_at(omega), ref_value, max_relative = 1e-12);
        }

        // Outside of the support of every contribution
        assert_eq!(dos.continuous_at(10.0), 0.0);

        // Divergent singular points
        assert_eq!(dos.continuous_at(eps - 2.0 * t), f64::INFINITY);
        assert_eq!(dos.continuous_at(eps + 2.0 * t), f64::INFINITY);
        assert_eq!(square(eps, t).continuous_at(eps), f64::INFINITY);
        assert_eq!(
            (-1.0 * square(eps, t)).continuous_at(eps),
            f64::NEG_INFINITY
        );

        // A square-root band edge is a singular point where A(ω) stays finite
        assert_eq!(semicircle(eps, 2.0).continuous_at(eps - 2.0), 0.0);

        // Discrete resonances contribute nothing
        assert_eq!(discrete(&[-0.7, 1.2], &[0.25, 0.6]).continuous_at(1.2), 0.0);
    }

    #[test]
    fn continuous_at_cancelling_divergences() {
        // Two logarithmic singularities meet at ω = 0, and the weights are chosen to
        // cancel their coefficients c = 1/(2π^2 t) exactly. What is left of A(ω) there is
        // the difference of the scales the two logarithms are written with,
        // ln(16t_1)/(2π^2 t_1) - 2ln(16t_2)/(2π^2 t_2) = -ln(2)/(2π^2).
        let dos = square(0.0, 1.0) + (-2.0) * square(0.0, 2.0);
        assert_relative_eq!(
            dos.continuous_at(0.0),
            -2f64.ln() / (2.0 * PI.powi(2)),
            max_relative = 1e-12
        );

        // A power law outranks a logarithm: the chain band edge wins over the
        // logarithmic peak of the square lattice, both sitting at ω = 0.
        let dos = square(0.0, 1.0) + (-1.0) * chain(2.0, 1.0);
        assert_eq!(dos.continuous_at(0.0), f64::NEG_INFINITY);
    }

    #[test]
    fn mul() {
        let dos = 2.0 * discrete(&[-0.7, 1.2], &[0.25, 0.6]) + 5.0 * gaussian(1.4, 0.5);
        assert_eq!(dos.discrete.len(), 2);
        assert_eq!(dos.continuous.len(), 1);

        // Weights underflowing to zero are dropped
        let scaled = f64::MIN_POSITIVE * (f64::MIN_POSITIVE * dos.clone());
        assert!(scaled.discrete.is_empty());
        assert!(scaled.continuous.is_empty());
        assert_eq!(scaled.total_weight(), 0.0);

        // Scaling by zero empties the spectral function
        let scaled = dos * 0.0;
        assert!(scaled.discrete.is_empty());
        assert!(scaled.continuous.is_empty());
        assert_eq!(scaled.total_weight(), 0.0);
    }

    #[test]
    fn add_merges_repeated_contributions() {
        let dos = chain(0.5, 1.0);
        let omega = 0.8;

        // A contribution repeated by addition is listed once, with the weights summed
        let tripled = dos.clone() + dos.clone() + dos.clone();
        assert_eq!(tripled.continuous.len(), 1);
        assert_relative_eq!(tripled.total_weight(), 3.0, epsilon = 1e-12);
        assert_relative_eq!(
            tripled.continuous_at(omega),
            3.0 * dos.continuous_at(omega),
            max_relative = 1e-12
        );

        // Weights cancelling each other out drop the contribution altogether
        let zero = 2.0 * dos.clone() - 2.0 * dos.clone();
        assert!(zero.continuous.is_empty());
        assert_eq!(zero.continuous_at(omega), 0.0);

        // Contributions are told apart by their allocation rather than by their
        // parameters, so separately built models of the same band stay apart
        let pair = chain(0.5, 1.0) + chain(0.5, 1.0);
        assert_eq!(pair.continuous.len(), 2);
        assert_relative_eq!(
            pair.continuous_at(omega),
            2.0 * dos.continuous_at(omega),
            max_relative = 1e-12
        );

        // Merging keeps the order of first appearance
        let mixed = dos.clone() + square(0.0, 1.0) + dos.clone();
        assert_eq!(mixed.continuous.len(), 2);
        assert_eq!(mixed.continuous[0].1, 2.0);
        assert_eq!(mixed.continuous[1].1, 1.0);
    }

    #[test]
    fn neg() {
        let dos = 2.0 * discrete(&[-0.7, 1.2], &[0.25, 0.6]) + 5.0 * gaussian(1.4, 0.5);
        let neg_dos = -dos.clone();

        // Negation flips the sign of every weight while preserving the structure
        assert_eq!(neg_dos.discrete.len(), 2);
        assert_eq!(neg_dos.continuous.len(), 1);
        assert_eq!(neg_dos.support(), dos.support());
        assert_relative_eq!(neg_dos.total_weight(), -dos.total_weight(), epsilon = 1e-12);
        assert_relative_eq!(
            neg_dos.discrete().find(1.2).unwrap().weight,
            -1.2,
            epsilon = 1e-12
        );
        for omega in [-1.0, 0.5, 1.4, 3.0] {
            assert_relative_eq!(
                neg_dos.continuous_at(omega),
                -dos.continuous_at(omega),
                max_relative = 1e-12
            );
        }

        // A divergence changes sign along with the spectral function
        let (eps, t) = (0.5f64, 1.0f64);
        assert_eq!((-square(eps, t)).continuous_at(eps), f64::NEG_INFINITY);

        // Negating twice is the identity
        let dos2 = -(-dos.clone());
        assert_relative_eq!(dos2.total_weight(), dos.total_weight(), epsilon = 1e-12);
        assert_relative_eq!(
            dos2.continuous_at(0.5),
            dos.continuous_at(0.5),
            epsilon = 1e-12
        );

        // A spectral function and its negation cancel each other out
        let zero = dos.clone() + (-dos);
        assert!(zero.discrete.is_empty());
        assert!(zero.continuous.is_empty());
        assert_relative_eq!(zero.total_weight(), 0.0, epsilon = 1e-12);
        assert_eq!(zero.continuous_at(1.4), 0.0);
    }

    #[test]
    fn sub() {
        let dos = 2.0 * discrete(&[-0.7, 1.2], &[0.25, 0.6]) + 5.0 * gaussian(1.4, 0.5);
        let rhs = discrete(&[1.2], &[1.0]) + 2.0 * gaussian(1.4, 0.5);
        let diff = dos.clone() - rhs.clone();

        // Subtraction acts on the discrete and the continuous parts alike
        assert_eq!(diff.discrete.len(), 2);
        assert_relative_eq!(
            diff.discrete().find(1.2).unwrap().weight,
            0.2,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            diff.total_weight(),
            dos.total_weight() - rhs.total_weight(),
            epsilon = 1e-12
        );
        for omega in [-1.0, 0.5, 1.4, 3.0] {
            assert_relative_eq!(
                diff.continuous_at(omega),
                dos.continuous_at(omega) - rhs.continuous_at(omega),
                max_relative = 1e-12
            );
        }

        // Subtraction agrees with adding the negation
        let diff2 = dos.clone() + (-rhs);
        assert_relative_eq!(diff2.total_weight(), diff.total_weight(), epsilon = 1e-12);
        assert_relative_eq!(
            diff2.continuous_at(0.5),
            diff.continuous_at(0.5),
            max_relative = 1e-12
        );

        // A spectral function subtracted from itself cancels out
        let zero = dos.clone() - dos;
        assert!(zero.discrete.is_empty());
        assert_relative_eq!(zero.total_weight(), 0.0, epsilon = 1e-12);
        assert_eq!(zero.continuous_at(1.4), 0.0);
    }

    #[test]
    fn precomputed() {
        let dos = square(0.0, 1.0);
        let fast = dos.precomputed(None);

        assert_eq!(fast.support(), dos.support());
        assert_relative_eq!(fast.total_weight(), dos.total_weight(), epsilon = 1e-14);

        // The singular part is carried over as it stands, divergence included
        assert_eq!(fast.continuous_at(0.0), f64::INFINITY);

        // Away from the singular point the expansion stands in for A(ω)
        for omega in [-3.5, -1.2, 0.3, 2.0, 3.9] {
            assert_relative_eq!(
                fast.continuous_at(omega),
                dos.continuous_at(omega),
                max_relative = 1e-10
            );
        }

        // Moments of even order survive the substitution, the odd ones vanishing
        for order in [0, 2, 4] {
            assert_relative_eq!(
                fast.integrate(|omega| omega.powi(order), None).unwrap(),
                dos.integrate(|omega| omega.powi(order), None).unwrap(),
                max_relative = 1e-10
            );
        }

        // An unbounded contribution has no interval to expand over and is left alone
        let g = gaussian(1.4, 0.5);
        let fast_g = g.precomputed(None);
        for omega in [0.0, 1.4, 3.0] {
            assert_eq!(fast_g.continuous_at(omega), g.continuous_at(omega));
        }

        // A band centre that is no singular point still breaks a panel, without which
        // the kink there would cap the fit five orders short
        for dos in [honeycomb(0.0, 1.0), lieb(0.0, 1.0)] {
            let fast = dos.precomputed(None);
            let (lo, hi) = dos.support().unwrap();
            let (mut worst, mut peak) = (0.0f64, 0.0f64);
            for i in 1..500 {
                let omega = lo + (hi - lo) * (i as f64) / 500.0;
                let (a, b) = (dos.continuous_at(omega), fast.continuous_at(omega));
                if a.is_finite() && b.is_finite() {
                    worst = worst.max((b - a).abs());
                    peak = peak.max(a.abs());
                }
            }
            assert!(worst < 1e-7 * peak, "fit stalled at {:.2e}", worst / peak);
        }

        // The discrete part passes through untouched
        let mixed = 0.5 * square(0.0, 1.0) + 0.5 * discrete(&[3.0], &[1.0]);
        let fast_mixed = mixed.precomputed(None);
        assert_eq!(fast_mixed.discrete().len(), 1);
        assert_relative_eq!(fast_mixed.total_weight(), 1.0, epsilon = 1e-14);
    }

    #[test]
    fn conv_discrete_continuous() {
        // A unit resonance displaces a band to its position
        let shifted = discrete(&[1.5], &[1.0]).conv(&semicircle(0.0, 2.0));
        let reference = semicircle(1.5, 2.0);
        assert_eq!(shifted.support(), reference.support());
        for omega in [-0.5, 0.4, 1.5, 2.6, 3.5] {
            assert_relative_eq!(
                shifted.continuous_at(omega),
                reference.continuous_at(omega),
                max_relative = 1e-14
            );
        }

        // Singular points travel with the band
        let shifted = discrete(&[2.0], &[1.0]).conv(&square(0.0, 1.0));
        assert_eq!(shifted.continuous_at(2.0), f64::INFINITY);
        assert_eq!(shifted.support(), Some((-2.0, 6.0)));

        // Weight is multiplicative, and the operands need not be normalized
        let a = 0.5 * discrete(&[-1.0, 2.0], &[0.25, 0.75]);
        let b = 3.0 * chain(0.0, 1.0);
        let c = a.conv(&b);
        assert_relative_eq!(
            c.total_weight(),
            a.total_weight() * b.total_weight(),
            max_relative = 1e-14
        );

        // One band per resonance, each displaced to its own position
        assert_eq!(c.continuous.len(), 2);
        assert_eq!(c.support(), Some((-3.0, 4.0)));

        // Convolution commutes
        for omega in [-2.5, 0.0, 1.3, 3.5] {
            assert_relative_eq!(
                c.continuous_at(omega),
                b.conv(&a).continuous_at(omega),
                max_relative = 1e-14
            );
        }

        // Moments obey M_n = Σ_k C(n,k) M_k^A M_{n-k}^B
        let moment = |sf: &SpectralFunction, n: i32| sf.integrate(|w| w.powi(n), None).unwrap();
        let binomial = |n: usize, k: usize| -> f64 {
            (1..=k).map(|i| (n - k + i) as f64 / i as f64).product()
        };
        for n in 0..=4usize {
            let reference: f64 = (0..=n)
                .map(|k| binomial(n, k) * moment(&a, k as i32) * moment(&b, (n - k) as i32))
                .sum();
            assert_relative_eq!(moment(&c, n as i32), reference, max_relative = 1e-9);
        }
    }

    #[test]
    fn conv_shifts_every_model() {
        // Displacing a contribution must move all of it: support, regular part,
        // singular points and breakpoints alike
        let by = 1.25f64;
        let unit = discrete(&[by], &[1.0]);
        let models = [
            chain(0.5, 1.0),
            semicircle(0.5, 2.0),
            bethe(4, 0.5, 1.0),
            powerlaw(0.5, -0.5, 2.0),
            powerlaw(0.5, 2.5, 2.0),
            pseudogap(0.5, 0.5, 2.0),
            pseudogap(0.5, 2.5, 2.0),
            square(0.5, 1.0),
            triangular(0.5, 1.0),
            honeycomb(0.5, 1.0),
            lieb(0.5, 1.0),
            kagome(0.5, 1.0),
            gaussian(0.5, 1.0),
            flat(0.5, 2.0, 0.0),
            flat(0.5, 2.0, 0.3),
            square(0.5, 1.0).precomputed(None),
        ];
        for dos in models {
            let moved = unit.conv(&dos);
            let (lo, hi) = dos.support().unwrap();
            let (lo, hi) = (lo.max(-20.0), hi.min(20.0));
            for i in 0..=40 {
                let omega = lo + (hi - lo) * (i as f64) / 40.0;
                let (before, after) = (dos.continuous_at(omega), moved.continuous_at(omega + by));
                assert_relative_eq!(after, before, max_relative = 1e-12, epsilon = 1e-300);
            }
            assert_relative_eq!(
                moved.total_weight(),
                dos.total_weight(),
                max_relative = 1e-14
            );
        }
    }

    #[test]
    fn conv_continuous_continuous() {
        // Two sharp-edged boxes of half-width d convolve into a triangle,
        // (2d - |ω|)/(4d^2), with a kink at the centre and none elsewhere
        let d = 1.5f64;
        let box_dos = flat(0.0, d, 0.0);
        let triangle = box_dos.conv_with(&box_dos, vec![], None);
        assert_eq!(triangle.support(), Some((-2.0 * d, 2.0 * d)));

        let reference = |omega: f64| (2.0 * d - omega.abs()) / (4.0 * d * d);
        for i in 0..=40 {
            let omega = -2.0 * d + 4.0 * d * (i as f64) / 40.0;
            assert_abs_diff_eq!(
                triangle.continuous_at(omega),
                reference(omega),
                epsilon = 1e-12
            );
        }

        // Weight is multiplicative here as everywhere
        assert_relative_eq!(triangle.total_weight(), 1.0, max_relative = 1e-12);

        // The kink is derived, not supplied: it is the sum of the two boxes' edges
        assert_eq!(
            crate::conv::breakpoints(&box_dos, &box_dos),
            vec![-2.0 * d, 0.0, 2.0 * d]
        );
    }

    #[test]
    fn conv_derives_its_breakpoints() {
        // A convolution stops being smooth where the frequencies at which either factor
        // does meet. For a semicircle on [-1.5, 2.5] against a box on [-2, 0] that is
        // the four sums of their edges.
        let (a, b) = (semicircle(0.5, 2.0), flat(-1.0, 1.0, 0.0));
        assert_eq!(crate::conv::breakpoints(&a, &b), vec![-3.5, -1.5, 0.5, 2.5]);

        // The chain edges against the square lattice edges and its logarithmic peak
        // give the four van Hove points of the simple cubic band
        let (ch, sq) = (chain(0.0, 1.0), square(0.0, 1.0));
        assert_eq!(
            crate::conv::breakpoints(&ch, &sq),
            vec![-6.0, -2.0, 2.0, 6.0]
        );
    }

    #[test]
    fn conv_two_continuous() {
        // `conv()` derives the structure that `simple_cubic()` is handed
        let t = 1.0f64;
        let derived = chain(0.0, t).conv(&square(0.0, t));
        let known = simple_cubic(0.0, t);
        assert_eq!(derived.support(), known.support());
        assert_relative_eq!(derived.total_weight(), 1.0, max_relative = 1e-14);

        // Deriving the frequencies without the asymptotics that go with them costs
        // about two orders of magnitude against the closed-form cusps
        for (order, reference) in [(0i32, 1.0f64), (2, 6.0 * t * t), (4, 90.0 * t.powi(4))] {
            let moment = derived.integrate(|omega| omega.powi(order), None).unwrap();
            assert_relative_eq!(moment, reference, max_relative = 1e-8);
        }

        // The two agree pointwise to the accuracy the interpolation affords
        let peak = known.continuous_at(0.0);
        for i in 1..200 {
            let omega = -6.0 * t + 12.0 * t * (i as f64) / 200.0;
            assert_abs_diff_eq!(
                derived.continuous_at(omega),
                known.continuous_at(omega),
                epsilon = 1e-6 * peak
            );
        }
    }

    #[test]
    fn conv_derives_its_band_edges() {
        use special::Gamma;

        // Two one-sided power laws convolve into a third, exactly:
        // c_1 x^{r_1} ⊛ c_2 x^{r_2} = c_1 c_2 B(r_1+1, r_2+1) x^{r_1+r_2+1}
        let w = 2.0f64;
        for (r1, r2) in [(-0.5f64, 0.5f64), (-0.5, -0.5), (0.0, 0.0), (-0.9, -0.9)] {
            let convolved = powerlaw(0.0, r1, w).conv(&powerlaw(0.0, r2, w));
            let prefactor = (r1 + 1.0) * (r2 + 1.0) / w.powf(r1 + r2 + 2.0);
            let beta =
                Gamma::gamma(r1 + 1.0) * Gamma::gamma(r2 + 1.0) / Gamma::gamma(r1 + r2 + 2.0);

            let singularities =
                crate::conv::singularities(&powerlaw(0.0, r1, w), &powerlaw(0.0, r2, w));
            let edge = singularities
                .iter()
                .find(|s| s.position() == 0.0)
                .expect("the lower band edge is derived");

            // The derived asymptotics carries the Beta coefficient, and with it the
            // whole of A(ω) near the edge, the regular part vanishing there
            for x in [1e-4f64, 1e-3, 1e-2] {
                let exact = prefactor * beta * x.powf(r1 + r2 + 1.0);
                assert_relative_eq!(edge.value(x), exact, max_relative = 1e-12);
                assert_relative_eq!(convolved.continuous_at(x), exact, max_relative = 1e-9);
            }
        }

        // The chain against the square lattice gives the simple cubic band edges, whose
        // coefficient is known: an inverse square root meeting the value the square
        // lattice takes at its own edge
        let t = 1.0f64;
        let singularities = crate::conv::singularities(&chain(0.0, t), &square(0.0, t));
        assert_eq!(singularities.len(), 2);
        assert_eq!(singularities[0].position(), -6.0 * t);
        assert_eq!(singularities[1].position(), 6.0 * t);
        let expected = 1.0 / (4.0 * PI.powi(2) * t.powf(1.5));
        for h in [1e-5f64, 1e-7] {
            let coefficient = singularities[0].value(-6.0 * t + h) / h.sqrt();
            assert_relative_eq!(coefficient, expected, max_relative = 1e-4);
        }
    }

    #[test]
    fn conv_continuous_moments() {
        // A band against a box. The result kinks wherever an edge of one meets an
        // edge of the other, at -1.5 and 0.5 inside the support [-3.5, 2.5].
        let a = semicircle(0.5, 2.0);
        let b = flat(-1.0, 1.0, 0.0);
        let c = a.conv_with(&b, vec![], None);
        assert_eq!(c.support(), Some((-3.5, 2.5)));

        assert_relative_eq!(
            c.total_weight(),
            a.total_weight() * b.total_weight(),
            max_relative = 1e-10
        );

        let moment = |sf: &SpectralFunction, n: i32| sf.integrate(|w| w.powi(n), None).unwrap();
        let binomial = |n: usize, k: usize| -> f64 {
            (1..=k).map(|i| (n - k + i) as f64 / i as f64).product()
        };
        for n in 0..=4usize {
            let reference: f64 = (0..=n)
                .map(|k| binomial(n, k) * moment(&a, k as i32) * moment(&b, (n - k) as i32))
                .sum();
            assert_relative_eq!(moment(&c, n as i32), reference, max_relative = 1e-8);
        }
    }

    #[test]
    fn broadened_discrete() {
        // A discrete spectral function is broadened into a sum of Lorentzians,
        // which serves as an analytic reference value.
        let levels = [-0.7f64, 1.2];
        let weights = [0.25f64, 0.6];
        let dos = 2.0 * discrete(&levels, &weights);

        for delta in [1e-1, 1e-2, 1e-4] {
            for omega in [-0.7, 0.0, 0.5, 1.2, 50.0] {
                let ref_value: f64 = 2.0
                    * levels
                        .iter()
                        .zip(&weights)
                        .map(|(eps, w)| w * (delta / PI) / (delta.powi(2) + (omega - eps).powi(2)))
                        .sum::<f64>();
                assert_relative_eq!(
                    dos.broadened(omega, delta, None).unwrap(),
                    ref_value,
                    max_relative = 1e-10
                );
            }
        }
    }

    #[test]
    fn broadened_narrow_limit() {
        // As delta -> 0 the broadened spectral function approaches A(omega) itself.
        // The error of the Lorentzian broadening is linear in delta.
        let (eps, t, omega) = (0.5f64, 1.0f64, 0.8f64);
        let dos = chain(eps, t);
        let ref_value = 1.0 / (PI * ((2.0 * t).powi(2) - (omega - eps).powi(2)).sqrt());
        for delta in [1e-4, 1e-6] {
            assert_relative_eq!(
                dos.broadened(omega, delta, None).unwrap(),
                ref_value,
                max_relative = 10.0 * delta
            );
        }
    }

    #[test]
    fn broadened_log_singularity() {
        // Broadening forces the quadrature to sample the spectral function arbitrarily
        // close to the logarithmic van Hove singularity of the square lattice. At the
        // band center the broadened spectral function diverges as -ln(delta)/(2π^2 t)
        // for delta -> 0.
        let t = 2.0f64;
        let dos = square(0.0, t);
        let mut prev = dos.broadened(0.0, 1e-2, None).unwrap();
        for delta in [1e-3, 1e-4, 1e-5, 1e-6] {
            let value = dos.broadened(0.0, delta, None).unwrap();
            assert_relative_eq!(
                value - prev,
                10f64.ln() / (2.0 * PI.powi(2) * t),
                max_relative = 1e-3
            );
            prev = value;
        }
    }

    #[test]
    fn broadened_mixed() {
        // Reference value computed independently with mpmath (40 decimal digits).
        let dos = 2.0 * discrete(&[-0.7, 1.2], &[0.25, 0.6]) + 5.0 * gaussian(1.4, 0.5);
        assert_relative_eq!(
            dos.broadened(0.5, 1e-2, None).unwrap(),
            0.8138402146,
            epsilon = 1e-9
        );
    }

    #[test]
    #[should_panic(expected = "broadening must be positive")]
    fn broadened_zero_delta() {
        let _ = gaussian(1.4, 0.5).broadened(0.5, 0.0, None);
    }
}
