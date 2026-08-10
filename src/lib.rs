mod util;

use std::ops::{Add, Mul};
use std::rc::Rc;

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
    pub fn one_resonance(eps: f64) -> DiscreteDOS {
        DiscreteDOS {
            resonances: vec![Resonance {
                eps,
                weight: 1.0f64,
            }],
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

/// Descriptor of an integrable singularity in DOS.
#[derive(Clone)]
pub struct Singularity {
    /// Position of the singularity, Ω_p
    position: f64,
    /// Asymptotic behavior of the DOS near the singularity, S_p(ω).
    /// To be chosen such that \lim_{ω → Ω_p} [A(ω) - S_p(ω)] = 0.
    asymptotics: Rc<dyn Fn(f64) -> f64>,
}

/// Continuous density of states possibly containing integrable singularities.
/// The DOS has the form A(ω) = R(ω) + ∑_p S_p(ω) for ω ∈ [ω_{min}, ω_{max}]
/// and zero otherwise. R(ω) is a smooth function and each S_p(ω) has one isolated
/// integrable singularity at Ω_p. The total spectral weight is expected to be 1.
pub trait ContinuousDOS {
    /// Support of the DOS function specified as a segment [ω_{min}, ω_{max}]
    fn support(&self) -> (f64, f64);
    /// Regular part of the DOS, R(ω)
    fn regular(&self, omega: f64) -> f64;
    /// Singular contributions S_p(ω), if any
    fn singularities(&self) -> &Singularity;
    /// Analytically derived values of ∫S_p(ω)dω over [ω_{min}, ω_{max}].
    fn sing_int(&self) -> [f64];
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn discrete_dos() {
        let mut dos = DiscreteDOS::new();
        assert_eq!(dos.resonances().len(), 0);
        assert_eq!(dos.norm(), 0.0);

        let mut dos = DiscreteDOS::one_resonance(2.0);
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
