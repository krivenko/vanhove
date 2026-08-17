pub mod discrete;
pub mod models;
mod util;

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
            let mut discrete: DiscreteDOS = self.discrete.clone();
            discrete.scale(a);
            Self {
                discrete,
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
        // Discrete part
        let mut sum = self;
        sum.discrete.extend(rhs.discrete.into_resonances());
        // Continuous part
        sum.continuous.extend(rhs.continuous);
        sum
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
    ) -> Result<Complex64, QuadratureError> {
        Ok(self.integrate(|omega| f(omega).re, None)?
            + Complex64::I * self.integrate(|omega| f(omega).im, None)?)
    }
}

#[cfg(test)]
mod tests {
    use crate::models::{discrete, gaussian};
    use approx::assert_relative_eq;

    #[test]
    fn norm() {
        let dos = 2.0 * discrete(&[-0.7, 1.2], &[0.25, 0.6]) + 5.0 * gaussian(1.4, 0.5);
        assert_relative_eq!(dos.norm(), 6.7, epsilon = 1e-12);
    }
}
