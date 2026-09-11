#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod discrete;
pub mod models;
mod singularity;
mod util;

use std::f64::consts::PI;
use std::ops::{Add, Mul, Neg, Sub};
use std::sync::Arc;

use bilby::QuadratureError;
use num_complex::Complex64;

use crate::discrete::DiscreteSF;
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
/// $S_p(\omega)$ is described in closed form by the corresponding [`Singularity`],
/// which fixes it over the whole support and not merely near $\Omega_p$. A support
/// carrying singularities must be bounded.
///
/// Implementations are shared behind an [`Arc`], so they must be [`Send`] and
/// [`Sync`]. A spectral function is read-only once built, which every model
/// satisfies by holding nothing but plain data.
trait ContinuousSF: Send + Sync {
    /// Support of the spectral function specified as a segment
    /// $[\omega_{min}, \omega_{max}]$.
    fn support(&self) -> (f64, f64);
    /// Regular part, $R(\omega)$.
    fn regular(&self, omega: f64) -> f64;
    /// Singular points $\Omega_p$ along with the closed form of $S_p$ at each.
    fn singularities(&self) -> &[Singularity] {
        &[]
    }
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
    fn from_continuous<C: ContinuousSF + 'static>(csf: C) -> SpectralFunction {
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
                if sing.position != omega {
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
                let omega_p = sing.position;
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
    use crate::models::{chain, discrete, gaussian, semicircle, square};
    use approx::assert_relative_eq;
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
