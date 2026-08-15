use crate::{ContinuousDOS, DensityOfStates, DiscreteDOS, Resonance};
/// Model densities of states commonly used in physics
use std::rc::Rc;

//
// Discrete DOS
//

/// Returns density of states comprised by a finite number of
/// discrete energy levels $\varepsilon_p$ with given weights $w_p$,
/// $$
///     A(\omega) = \sum_p w_p \delta(\omega - \varepsilon_p).
/// $$
pub fn discrete(levels: &[f64], weights: &[f64]) -> DensityOfStates {
    let mut discrete = DiscreteDOS::new();
    for (eps, weight) in levels.iter().zip(weights) {
        discrete.add(&Resonance {
            eps: *eps,
            weight: *weight,
        });
    }
    DensityOfStates {
        discrete,
        continuous: vec![],
    }
}

//
// Gaussian DOS
//

/// Gaussian density of states
struct GaussianDOS {
    eps: f64,
    denom: f64,
    prefactor: f64,
}
impl GaussianDOS {
    fn new(eps: f64, sigma: f64) -> GaussianDOS {
        GaussianDOS {
            eps,
            denom: 2.0 * sigma.powi(2),
            prefactor: 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt()),
        }
    }
}
impl ContinuousDOS for GaussianDOS {
    fn support(&self) -> (f64, f64) {
        (-f64::INFINITY, f64::INFINITY)
    }
    fn regular(&self, omega: f64) -> f64 {
        self.prefactor * (-(omega - self.eps).powi(2) / self.denom).exp()
    }
}

/// Returns the normalized Gaussian density of states centered at `eps` with width `sigma`,
/// $$
///     A(\omega) = \frac{1}{\sqrt{2\pi\sigma^2}}
///         \exp\left(-\frac{(\omega - \eps)^2}{2\sigma^2}\right).
/// $$
pub fn gaussian(eps: f64, sigma: f64) -> DensityOfStates {
    assert!(sigma > 0.0, "width must be positive");
    DensityOfStates {
        discrete: DiscreteDOS::new(),
        continuous: vec![(Rc::new(GaussianDOS::new(eps, sigma)), 1.0)],
    }
}

//
// Semicircle DOS
//

/// Semicircle (Wigner) density of states
struct SemicircleDOS {
    eps: f64,
    radius: f64,
    /// Band edges. For this model they double as the positions of the two
    /// square-root singularities, which sit exactly at the edges of the support.
    edges: [f64; 2],
    prefactor: f64,
}
impl SemicircleDOS {
    fn new(eps: f64, radius: f64) -> SemicircleDOS {
        assert!(radius > 0.0, "radius must be positive");
        SemicircleDOS {
            eps,
            radius,
            edges: [eps - radius, eps + radius],
            prefactor: 2.0 / (std::f64::consts::PI * radius),
        }
    }
}
impl ContinuousDOS for SemicircleDOS {
    fn support(&self) -> (f64, f64) {
        self.edges.into()
    }
    fn regular(&self, omega: f64) -> f64 {
        if omega == self.edges[0] || omega == self.edges[1] {
            -2.0 * self.prefactor
        } else {
            let x = (omega - self.eps) / self.radius;
            self.prefactor
                * ((1.0 - x * x).sqrt() - (2.0 * (1.0 - x)).sqrt() - (2.0 * (1.0 + x)).sqrt())
        }
    }
    fn singularities(&self) -> &[f64] {
        &self.edges
    }
    fn asymptotics(&self, p: usize, omega: f64) -> f64 {
        let x = (omega - self.eps) / self.radius;
        // p == 0 is the lower edge, p == 1 the upper one
        let sign = if p == 0 { 1.0 } else { -1.0 };
        self.prefactor * (2.0 * (1.0 + sign * x)).sqrt()
    }
    fn asympt_int(&self, _p: usize) -> f64 {
        16.0 / (3.0 * std::f64::consts::PI)
    }
}

/// Returns the normalized semicircle (Wigner) density of states centered at `eps` with
/// radius `r`,
/// $$
///     A(\omega) = \frac{2}{\pi r^2} \sqrt{r^2 - \omega^2}\theta(r^2 - \omega^2).
/// $$
pub fn semicircle(eps: f64, r: f64) -> DensityOfStates {
    assert!(r > 0.0, "radius must be positive");
    DensityOfStates {
        discrete: DiscreteDOS::new(),
        continuous: vec![(Rc::new(SemicircleDOS::new(eps, r)), 1.0)],
    }
}

#[cfg(test)]
mod tests {
    use crate::{DensityOfStates, models};
    use approx::assert_relative_eq;

    fn compute_moment(dos: &DensityOfStates, order: i32) -> f64 {
        dos.integrate(|omega| omega.powi(order), None).unwrap()
    }

    #[test]
    fn discrete() {
        let levels = vec![-1.0, 0.0, 2.0];
        let weights = vec![0.5, 0.5, 0.25];
        let dos = models::discrete(&levels, &weights);
        for order in 0..=6 {
            let moment = compute_moment(&dos, order as i32);
            let moment_ref = levels
                .iter()
                .zip(&weights)
                .map(|(level, weight)| level.powi(order) * weight)
                .sum();
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn gaussian() {
        let eps = 1.0f64;
        let sigma = 3.0f64;
        let dos = models::gaussian(eps, sigma);
        let moments_ref = vec![
            1.0,
            eps,
            eps.powi(2) + sigma.powi(2),
            eps.powi(3) + 3.0 * eps * sigma.powi(2),
            eps.powi(4) + 6.0 * eps.powi(2) * sigma.powi(2) + 3.0 * sigma.powi(4),
            eps.powi(5) + 10.0 * eps.powi(3) * sigma.powi(2) + 15.0 * eps * sigma.powi(4),
            eps.powi(6)
                + 15.0 * eps.powi(4) * sigma.powi(2)
                + 45.0 * eps.powi(2) * sigma.powi(4)
                + 15.0 * sigma.powi(6),
        ];
        for order in 0..=6 {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moments_ref[order], max_relative = 1e-10);
        }
    }

    #[test]
    fn semicircle() {
        let eps = 0.5f64;
        let t = 2.0f64;
        let dos = models::semicircle(eps, 2.0 * t);
        let moments_ref = vec![
            1.0,
            eps,
            eps.powi(2) + t.powi(2),
            eps.powi(3) + 3.0 * eps * t.powi(2),
            eps.powi(4) + 6.0 * eps.powi(2) * t.powi(2) + 2.0 * t.powi(4),
            eps.powi(5) + 10.0 * eps.powi(3) * t.powi(2) + 10.0 * eps * t.powi(4),
            eps.powi(6)
                + 15.0 * eps.powi(4) * t.powi(2)
                + 30.0 * eps.powi(2) * t.powi(4)
                + 5.0 * t.powi(6),
        ];
        for order in 0..=6 {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moments_ref[order], max_relative = 1e-10);
        }
    }
}
