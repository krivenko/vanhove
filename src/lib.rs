mod util;

use std::ops::{Add, Mul};
use std::rc::Rc;

use bilby::QuadratureError;
use num_complex::Complex64;

//
// DiscreteDOS
//

/// Single discrete resonance w δ(ω - ε)
#[derive(Debug, Clone, Copy)]
pub struct Resonance {
    // Position of the resonance, ε
    eps: f64,
    // Weight of the resonance, w
    weight: f64,
}

/// Discrete density of states, A(ω) = ∑_p w_p δ(ω - ε_p).
#[derive(Debug, Default, Clone)]
pub struct DiscreteDOS {
    /// List of resonances, (ε_p, w_p)
    resonances: Vec<Resonance>,
    /// resonances[..clean] is sorted w.r.t. ε_p, deduplicated, zero-weight-free
    clean: usize,
    /// Relative tolerance level for negligible weights
    tol: f64,
}

impl DiscreteDOS {
    pub fn default() -> Self {
        DiscreteDOS {
            resonances: vec![],
            clean: 0,
            tol: f64::EPSILON,
        }
    }

    /// Make an empty discrete DOS (zero resonances)
    pub fn new() -> Self {
        Self::default()
    }

    /// Make a discrete DOS with one resonance
    pub fn one_resonance(eps: f64, weight: f64) -> DiscreteDOS {
        DiscreteDOS {
            resonances: vec![Resonance { eps, weight }],
            clean: 1,
            tol: f64::EPSILON,
        }
    }

    /// Access the resonance list
    pub fn resonances(&mut self) -> &[Resonance] {
        if self.clean != self.resonances.len() {
            self.prune();
        }
        &self.resonances
    }

    /// Total spectral weight of the DOS
    pub fn norm(&self) -> f64 {
        self.resonances[..self.clean].iter().map(|r| r.weight).sum()
    }

    /// Add a resonance to the DOS
    pub fn add(&mut self, res: &Resonance) {
        debug_assert!(!res.eps.is_nan(), "res.eps must not be NaN");
        if res.weight != 0.0 {
            // canonicalize -0.0 so it groups with 0.0 under total_cmp
            self.resonances.push(Resonance {
                eps: if res.eps == 0.0 { 0.0 } else { res.eps },
                weight: res.weight,
            });
        }
        if self.resonances.len() > 2 * self.clean + 16 {
            self.prune();
        }
    }

    /// Sort, deduplicate and remove small-weight resonances
    fn prune(&mut self) {
        // Sort the resonances w.r.t. ε
        let res = &mut self.resonances;
        res.sort_by(|r1, r2| r1.eps.total_cmp(&r2.eps));

        let (mut w, mut i) = (0, 0);
        while i < res.len() {
            // Find a slice of self.resonances sharing the same ε
            let start = i;
            let eps = res[i].eps;
            while i < res.len() && res[i].eps == eps {
                i += 1;
            }

            // Sum all weights within the slice, with Neumaier compensation
            let group = &mut res[start..i];
            group.sort_unstable_by(|r1, r2| r1.weight.abs().total_cmp(&r2.weight.abs()));
            let total = util::kahan_babushka_neumaier_sum(group.iter().map(|r| r.weight));

            // Add a resonance to the result if the total weight is not too small compared
            // to the total magnitude of the weights
            let mag: f64 = group.iter().map(|r| r.weight.abs()).sum();
            if total.abs() > self.tol * mag {
                res[w] = Resonance { eps, weight: total };
                w += 1;
            }
        }
        // Remove unused elements of self.resonances
        res.truncate(w);
        self.clean = w;
    }
}

//
// ContinuousDOS
//

/// Continuous density of states possibly containing integrable singularities.
/// The DOS has the form A(ω) = R(ω) + ∑_p S_p(ω) for ω ∈ [ω_{min}, ω_{max}]
/// and zero otherwise. R(ω) is a smooth function and each S_p(ω) has one isolated
/// integrable singularity at Ω_p. The total spectral weight is expected to be 1.
pub trait ContinuousDOS {
    /// Support of the DOS function specified as a segment [ω_{min}, ω_{max}]
    fn support(&self) -> (f64, f64);
    /// Regular part of the DOS, R(ω)
    fn regular(&self, omega: f64) -> f64;
    /// Positions of singularities
    fn sing_pos(&self) -> &[f64];
    /// Singular contributions S_p(ω)
    fn asymptotics(&self) -> &[Rc<dyn Fn(f64) -> f64>];
    /// Analytically derived values of ∫S_p(ω)dω over [ω_{min}, ω_{max}].
    fn asympt_int(&self) -> &[f64];
}

// Density of states as a weighted sum of discrete resonances and continuous contributions
#[derive(Clone)]
pub struct DensityOfStates {
    // Contributions of discrete resonances
    discrete: DiscreteDOS,
    // Continuous contributions with their weights
    continuous: Vec<(Rc<dyn ContinuousDOS>, f64)>,
}

/// Multiply DOS by a real number from the right
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
            for res in &mut discrete.resonances {
                res.weight *= a;
            }
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

/// Multiply DOS by a real number from the left
impl Mul<DensityOfStates> for f64 {
    type Output = DensityOfStates;
    fn mul(self, dos: DensityOfStates) -> DensityOfStates {
        dos * self
    }
}

/// Addition of two densities of states
impl Add for DensityOfStates {
    type Output = Self;
    fn add(self, mut rhs: DensityOfStates) -> DensityOfStates {
        // Discrete part
        let mut sum = self.clone();
        for r in rhs.discrete.resonances() {
            sum.discrete.add(r);
        }
        for r in rhs.continuous {
            sum.continuous.push(r);
        }
        sum
    }
}

//
// DOS integration
//

impl DensityOfStates {
    /// Given a density of states A(ω) and a real-valued function f(ω),
    /// compute the integral ∫A(ω)f(ω)dω.
    pub fn integrate<F: Fn(f64) -> f64>(
        &self,
        f: F,
        tol: Option<f64>,
    ) -> Result<f64, QuadratureError> {
        let mut result = 0.0f64;

        // Discrete spectral contributions
        for r in &self.discrete.resonances {
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
            for (p, omega_p) in cdos.sing_pos().into_iter().enumerate() {
                // ∫S_p(ω)[f(ω) - f(Ω_p)]dω
                let s_p = &cdos.asymptotics()[p];
                let f_p = f(*omega_p);
                res_contrib += util::bilby_integrate(
                    |omega| s_p(omega) * (f(omega) - f_p),
                    omega_min,
                    omega_max,
                    tol,
                )?
                .value;
                // f(Ω_p) ∫S_p(ω)dω
                res_contrib += cdos.asympt_int()[p] * f_p;
            }

            result += w * res_contrib;
        }
        Ok(result)
    }

    /// Given a density of states A(ω) and a complex-valued function f(ω),
    /// compute the integral ∫A(ω)f(ω)dω
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
    use super::*;

    #[test]
    fn discrete_dos() {
        let mut dos = DiscreteDOS::new();
        assert_eq!(dos.resonances().len(), 0);
        assert_eq!(dos.norm(), 0.0);

        let mut dos = DiscreteDOS::one_resonance(2.0, 1.0);
        assert_eq!(dos.resonances().len(), 1);
        assert_eq!(dos.resonances()[0].eps, 2.0);
        assert_eq!(dos.resonances()[0].weight, 1.0);
        assert_eq!(dos.norm(), 1.0);

        dos.add(&Resonance {
            eps: -1.5,
            weight: 0.25,
        });
        assert_eq!(dos.resonances().len(), 2);
        assert_eq!(dos.norm(), 1.25);
        dos.add(&Resonance {
            eps: 2.0,
            weight: -1.0,
        });
        assert_eq!(dos.resonances().len(), 1);
        assert_eq!(dos.norm(), 0.25);
        dos.add(&Resonance {
            eps: 5.0,
            weight: 0.0,
        });
        assert_eq!(dos.resonances().len(), 1);
        assert_eq!(dos.norm(), 0.25);
    }
}
