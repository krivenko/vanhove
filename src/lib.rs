pub mod discrete;
pub mod models;
mod util;

use std::f64::consts::PI;
use std::ops::{Add, Mul};
use std::rc::Rc;

use bilby::QuadratureError;
use num_complex::Complex64;

use crate::discrete::DiscreteDOS;

//
// ContinuousDOS
//

/// Continuous density of states possibly containing integrable singularities.
/// The DOS has the form $A(\omega) = R(\omega) + \sum_p S_p(\omega)$ for
/// $\omega \in [\omega_{min}, \omega_{max}]$ and zero otherwise.
/// $R(\omega)$ is a smooth function and each $S_p(\omega)$ has one isolated
/// integrable singularity at $\Omega_p$. The total spectral weight is expected to be 1.
///
/// The singular terms are keyed by their index p, which runs over the valid
/// indices of the slice returned by `singularities()`.
trait ContinuousDOS {
    /// Support of the DOS function specified as a segment $[\omega_{min}, \omega_{max}]$
    fn support(&self) -> (f64, f64);
    /// Regular part of the DOS, $R(\omega)$
    fn regular(&self, omega: f64) -> f64;
    /// Positions of the singular points, Ω_p
    fn singularities(&self) -> &[f64] {
        &[]
    }
    /// Asymptotic form of the DOS near the p-th singular point, $S_p(\omega)$
    fn asymptotics(&self, _p: usize, _omega: f64) -> f64 {
        unreachable!("this DOS has no singularities")
    }
    /// Analytically derived value of $\int S_p(\omega)d\omega$ over
    /// $[\omega_{min}, \omega_{max}]$
    fn asympt_int(&self, _p: usize) -> f64 {
        unreachable!("this DOS has no singularities")
    }
}

/// Density of states as a weighted sum of discrete resonances and continuous contributions.
#[derive(Clone)]
pub struct DensityOfStates {
    // Contributions of discrete resonances
    discrete: DiscreteDOS,
    // Continuous contributions with their weights
    continuous: Vec<(Rc<dyn ContinuousDOS>, f64)>,
}

/// Multiply DOS by a real number from the right.
impl Mul<f64> for DensityOfStates {
    type Output = Self;

    fn mul(self, a: f64) -> Self {
        if a == 0.0 {
            // DOS with no contributions
            Self {
                discrete: DiscreteDOS::new(),
                continuous: vec![],
            }
        } else {
            Self {
                discrete: self.discrete * a,
                continuous: self
                    .continuous
                    .into_iter()
                    .map(|(cd, w)| (cd, w * a))
                    .collect(),
            }
        }
    }
}
/// Multiply DOS by a real number from the left.
impl Mul<DensityOfStates> for f64 {
    type Output = DensityOfStates;
    fn mul(self, dos: DensityOfStates) -> DensityOfStates {
        dos * self
    }
}
/// Addition of two densities of states.
impl Add for DensityOfStates {
    type Output = Self;
    fn add(self, rhs: DensityOfStates) -> DensityOfStates {
        let mut continuous = self.continuous;
        continuous.extend(rhs.continuous);
        DensityOfStates {
            discrete: self.discrete + rhs.discrete,
            continuous,
        }
    }
}
impl DensityOfStates {
    /// Total spectral weight of the density of states.
    pub fn norm(&self) -> f64 {
        self.discrete.norm() + self.continuous.iter().map(|(_, w)| w).sum::<f64>()
    }

    //
    // DOS integration
    //

    /// Given a density of states $A(\omega)$ and a real-valued function $f(\omega)$,
    /// compute the integral $\int A(\omega)f(\omega)d\omega$.
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
        for (cdos, w) in &self.continuous {
            let mut res_contrib = 0.0f64;
            let (omega_min, omega_max) = cdos.support();

            // Integrate the regular part, ∫R(ω)f(ω)dω
            res_contrib += util::bilby_integrate(
                |omega| cdos.regular(omega) * f(omega),
                omega_min,
                omega_max,
                tol,
            )?
            .value;

            // Add integrals of the asymptotics
            for (p, &omega_p) in cdos.singularities().iter().enumerate() {
                // ∫S_p(ω)[f(ω) - f(Ω_p)]dω
                let f_p = f(omega_p);
                res_contrib += util::bilby_integrate(
                    |omega| {
                        if omega == omega_p {
                            0.0
                        } else {
                            cdos.asymptotics(p, omega) * (f(omega) - f_p)
                        }
                    },
                    omega_min,
                    omega_max,
                    tol,
                )?
                .value;
                // ∫S_p(ω)dω f(Ω_p)
                res_contrib += cdos.asympt_int(p) * f_p;
            }

            result += w * res_contrib;
        }
        Ok(result)
    }

    /// Given a density of states $A(\omega)$ and a complex-valued function $f(\omega)$,
    /// compute the integral $\int A(\omega)f(\omega)d\omega$.
    pub fn integrate_complex<F: Fn(f64) -> Complex64>(
        &self,
        f: F,
        tol: Option<f64>,
    ) -> Result<Complex64, QuadratureError> {
        Ok(self.integrate(|omega| f(omega).re, tol)?
            + Complex64::I * self.integrate(|omega| f(omega).im, tol)?)
    }

    /// Evaluate value of the broadened density of states at a frequency `omega` by
    /// computing its convolution with the Lorentzian line shape function of the
    /// half-width at half-maximum `delta`,
    /// $$
    ///     \int \frac{1}{\pi} \frac{\delta}{\delta^2 + (\omega-\omega')^2}
    ///         A(\omega') d\omega'.
    /// $$
    /// The discrete part of the DOS contributes a Lorentzian per resonance, so that
    /// the broadened DOS is a smooth function of `omega` for any $\delta > 0$.
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
    use crate::models::{chain, discrete, gaussian, square};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn norm() {
        let dos = 2.0 * discrete(&[-0.7, 1.2], &[0.25, 0.6]) + 5.0 * gaussian(1.4, 0.5);
        assert_relative_eq!(dos.norm(), 6.7, epsilon = 1e-12);
    }

    #[test]
    fn broadened_discrete() {
        // A discrete DOS is broadened into a sum of Lorentzians, which serves as
        // an analytic reference value.
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
        // As delta -> 0 the broadened DOS approaches A(omega) itself. The error of
        // the Lorentzian broadening is linear in delta.
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
        // Broadening forces the quadrature to sample the DOS arbitrarily close to
        // the logarithmic van Hove singularity of the square lattice. At the band
        // center the broadened DOS diverges as -ln(delta)/(2π^2 t) for delta -> 0.
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
