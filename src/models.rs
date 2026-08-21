//! Model densities of states commonly used in physics

use crate::discrete::Resonance;
use crate::util::fermi;
use crate::{ContinuousSF, SingularLaw, Singularity, SpectralFunction};

use special::Elliptic;
use std::f64::consts::PI;

//
// Discrete DOS
//

/// Returns density of states comprised by a finite number of
/// discrete energy levels $\varepsilon_p$ with given weights $w_p$,
/// $$
///     A(\omega) = \sum_p w_p \delta(\omega - \varepsilon_p).
/// $$
pub fn discrete(levels: &[f64], weights: &[f64]) -> SpectralFunction {
    SpectralFunction::from_discrete_continuous(
        levels
            .iter()
            .zip(weights)
            .map(|(&eps, &weight)| Resonance { eps, weight })
            .collect(),
        vec![],
    )
}

//
// (Soft-edged) flat DOS
//

struct FlatDOS {
    eps: f64,
    d: f64,
    delta: f64,
    prefactor: f64,
}
impl FlatDOS {
    fn new(eps: f64, d: f64, delta: f64) -> FlatDOS {
        assert!(d > 0.0, "bandwidth must be positive");
        assert!(delta >= 0.0, "edge width must be non-negative");
        FlatDOS {
            eps,
            d,
            delta,
            prefactor: if delta == 0.0 {
                1.0 / (2.0 * d)
            } else {
                let x = (d / delta).tanh();
                x / (d * (1.0 + x))
            },
        }
    }
}
impl ContinuousSF for FlatDOS {
    fn support(&self) -> (f64, f64) {
        if self.delta == 0.0 {
            (self.eps - self.d, self.eps + self.d)
        } else {
            (f64::NEG_INFINITY, f64::INFINITY)
        }
    }
    fn regular(&self, omega: f64) -> f64 {
        if self.delta == 0.0 {
            self.prefactor
        } else {
            self.prefactor
                * fermi((omega - self.eps - self.d) / self.delta)
                * fermi(-(omega - self.eps + self.d) / self.delta)
        }
    }
}

/// Returns the normalized flat density of states centered at `eps` with half-bandwidth `d`,
/// and smooth, Fermi-like band edges of width `delta`,
/// $$
///     A(\omega) = \frac{1}{d(1 + \coth(d / \delta))}
///         \frac{1}{(e^{(\omega-\epsilon-d)/\delta}+1)(e^{-(\omega-\epsilon+d)/\delta}+1)}.
/// $$
/// $\delta=0$ is a valid value that turns the edges into sharp discontinuous jumps,
/// $$
///     A(\omega) = \frac{1}{2d}\theta(d - |\omega-\epsilon|),
/// $$
/// the jump points $\omega = \epsilon \pm d$ evaluating to $1/(2d)$.
pub fn flat(eps: f64, d: f64, delta: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(FlatDOS::new(eps, d, delta))
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
        assert!(sigma > 0.0, "width must be positive");
        GaussianDOS {
            eps,
            denom: 2.0 * sigma.powi(2),
            prefactor: 1.0 / (sigma * (2.0 * PI).sqrt()),
        }
    }
}
impl ContinuousSF for GaussianDOS {
    fn support(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }
    fn regular(&self, omega: f64) -> f64 {
        self.prefactor * (-(omega - self.eps).powi(2) / self.denom).exp()
    }
}

/// Returns the normalized Gaussian density of states centered at `eps` with width `sigma`,
/// $$
///     A(\omega) = \frac{1}{\sqrt{2\pi\sigma^2}}
///         \exp\left(-\frac{(\omega - \epsilon)^2}{2\sigma^2}\right).
/// $$
pub fn gaussian(eps: f64, sigma: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(GaussianDOS::new(eps, sigma))
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
    /// $S_p$ vanishes at either edge, so neither law is divergent.
    edges: [Singularity; 2],
    prefactor: f64,
}
impl SemicircleDOS {
    fn new(eps: f64, radius: f64) -> SemicircleDOS {
        assert!(radius > 0.0, "radius must be positive");
        SemicircleDOS {
            eps,
            radius,
            edges: [
                Singularity {
                    position: eps - radius,
                    law: SingularLaw::Finite,
                },
                Singularity {
                    position: eps + radius,
                    law: SingularLaw::Finite,
                },
            ],
            prefactor: 2.0 / (PI * radius),
        }
    }
}
impl ContinuousSF for SemicircleDOS {
    fn support(&self) -> (f64, f64) {
        (self.edges[0].position, self.edges[1].position)
    }
    fn regular(&self, omega: f64) -> f64 {
        if omega == self.edges[0].position || omega == self.edges[1].position {
            -2.0 * self.prefactor
        } else {
            let x = (omega - self.eps) / self.radius;
            self.prefactor
                * ((1.0 - x * x).sqrt() - (2.0 * (1.0 - x)).sqrt() - (2.0 * (1.0 + x)).sqrt())
        }
    }
    fn singularities(&self) -> &[Singularity] {
        &self.edges
    }
    fn asymptotics(&self, p: usize, omega: f64) -> f64 {
        debug_assert!(p <= 1);
        let x = (omega - self.eps) / self.radius;
        // p == 0 is the lower edge, p == 1 the upper one
        let sign = if p == 0 { 1.0 } else { -1.0 };
        self.prefactor * (2.0 * (1.0 + sign * x)).sqrt()
    }
    fn asympt_int(&self, _p: usize) -> f64 {
        16.0 / (3.0 * PI)
    }
}

/// Returns the normalized semicircle (Wigner) density of states centered at `eps` with
/// radius `r`,
/// $$
///     A(\omega) = \frac{2}{\pi r^2} \sqrt{r^2 - \omega^2}\theta(r^2 - \omega^2).
/// $$
pub fn semicircle(eps: f64, r: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(SemicircleDOS::new(eps, r))
}

//
// Chain DOS
//

/// Density of states of a linear chain
struct ChainDOS {
    eps: f64,
    t: f64,
    /// Band edges. For this model they double as the positions of the two
    /// square-root singularities, which sit exactly at the edges of the support.
    edges: [Singularity; 2],
    prefactor: f64,
}
impl ChainDOS {
    fn new(eps: f64, t: f64) -> ChainDOS {
        assert!(t > 0.0, "hopping constant must be positive");
        let law = SingularLaw::Power {
            a: 0.5,
            c: 1.0 / (2.0 * PI * t.sqrt()),
            l: 0.0,
        };
        ChainDOS {
            eps,
            t,
            edges: [
                Singularity {
                    position: eps - 2.0 * t,
                    law,
                },
                Singularity {
                    position: eps + 2.0 * t,
                    law,
                },
            ],
            prefactor: 1.0 / (2.0 * PI * t),
        }
    }
}
impl ContinuousSF for ChainDOS {
    fn support(&self) -> (f64, f64) {
        (self.edges[0].position, self.edges[1].position)
    }
    fn regular(&self, omega: f64) -> f64 {
        if omega == self.edges[0].position || omega == self.edges[1].position {
            -0.75 * self.prefactor
        } else {
            let x = (omega - self.eps) / (2.0 * self.t);

            // 'sm' regularizes the derivative near ω = -2t to ease integration
            let rm = (2.0 * (1.0 + x)).sqrt();
            let sm = self.prefactor * (1.0 / rm + rm / 8.0);

            // 'sp' regularizes the derivative near ω = 2t to ease integration
            let rp = (2.0 * (1.0 - x)).sqrt();
            let sp = self.prefactor * (1.0 / rp + rp / 8.0);

            self.prefactor / (1.0 - x * x).sqrt() - sp - sm
        }
    }
    fn singularities(&self) -> &[Singularity] {
        &self.edges
    }
    fn asymptotics(&self, p: usize, omega: f64) -> f64 {
        debug_assert!(p <= 1);
        let x = (omega - self.eps) / (2.0 * self.t);
        // p == 0 is the lower edge, p == 1 the upper one
        let sign = if p == 0 { 1.0 } else { -1.0 };
        let s = (2.0 * (1.0 + sign * x)).sqrt();
        // The second term regularizes the derivative
        self.prefactor * (1.0 / s + s / 8.0)
    }
    fn asympt_int(&self, _p: usize) -> f64 {
        7.0 / (3.0 * PI)
    }
}

/// Returns the normalized density of states of a linear chain with the hopping constant `t`
/// and the local energy level `eps`,
/// $$
///     A(\omega) = \frac{1}{\pi} \frac{1}{\sqrt{(2t)^2 - (\omega-\epsilon)^2}}
///         \theta((2t)^2 - (\omega-\epsilon)^2).
/// $$
pub fn chain(eps: f64, t: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(ChainDOS::new(eps, t))
}

//
// Square lattice DOS
//

/// Value of $|\omega-\epsilon|/(4t)$ below which the regular part of the square lattice DOS
/// is evaluated using a series expansion. At the threshold both expressions agree to
/// within a relative error of $10^{-5}$.
const SQUARE_DOS_SERIES_THRESHOLD: f64 = 1e-3;

/// Density of states of a square lattice
struct SquareDOS {
    eps: f64,
    t: f64,
    /// Band edges
    edges: [f64; 2],
    /// Position of the logarithmic van Hove singularity, which sits at the band center.
    singularity: [Singularity; 1],
    prefactor: f64,
}
impl SquareDOS {
    fn new(eps: f64, t: f64) -> SquareDOS {
        assert!(t > 0.0, "hopping constant must be positive");
        let prefactor = 1.0 / (2.0 * PI.powi(2) * t);
        SquareDOS {
            eps,
            t,
            edges: [eps - 4.0 * t, eps + 4.0 * t],
            singularity: [Singularity {
                position: eps,
                law: SingularLaw::Log {
                    c: prefactor,
                    l: prefactor * (16.0 * t).ln(),
                },
            }],
            prefactor,
        }
    }
}
impl ContinuousSF for SquareDOS {
    fn support(&self) -> (f64, f64) {
        self.edges.into()
    }
    fn regular(&self, omega: f64) -> f64 {
        let ax = ((omega - self.eps) / (4.0 * self.t)).abs();
        if ax == 0.0 {
            0.0
        } else if ax < SQUARE_DOS_SERIES_THRESHOLD {
            // Close to the singularity, K(1 - ax^2) ≈ ln(4/ax) nearly cancels the
            // logarithm, and evaluating the difference directly is catastrophically
            // inaccurate. Worse, 1 - ax^2 rounds to 1 for ax below ~1e-8, which makes
            // elliptic_k() overflow. Use the expansion
            // K(1 - a^2) = ln(4/a) + (a^2/4)[ln(4/a) - 1] + O(a^4 ln a) instead.
            self.prefactor * 0.25 * ax.powi(2) * ((4.0 / ax).ln() - 1.0)
        } else {
            self.prefactor * ((1.0 - ax * ax).elliptic_k() + (0.25 * ax).ln())
        }
    }
    fn singularities(&self) -> &[Singularity] {
        &self.singularity
    }
    fn asymptotics(&self, p: usize, omega: f64) -> f64 {
        debug_assert!(p == 0);
        let x = (omega - self.eps) / (16.0 * self.t);
        -self.prefactor * x.abs().ln()
    }
    fn asympt_int(&self, _p: usize) -> f64 {
        4.0 / PI.powi(2) * (1.0 + 2.0 * std::f64::consts::LN_2)
    }
}

/// Returns the normalized density of states of a square lattice with the hopping constant
/// `t` and the local energy level `eps`,
/// $$
///     A(\omega) = \frac{1}{2\pi^2 t} K\left(1 - \left(\frac{\omega-\epsilon}{4t}\right)^2\right)
///         \theta((4t)^2 - (\omega-\epsilon)^2),
/// $$
/// where $K(m)$ is the complete elliptic integral of the first kind.
pub fn square(eps: f64, t: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(SquareDOS::new(eps, t))
}

#[cfg(test)]
mod tests {
    use crate::{SpectralFunction, models};
    use approx::assert_relative_eq;

    fn compute_moment(dos: &SpectralFunction, order: i32) -> f64 {
        dos.integrate(|omega| omega.powi(order), None).unwrap()
    }

    fn central_to_moments(eps: f64, mu2: f64, mu4: f64, mu6: f64) -> Vec<f64> {
        vec![
            1.0,
            eps,
            eps.powi(2) + mu2,
            eps.powi(3) + 3.0 * eps * mu2,
            eps.powi(4) + 6.0 * eps.powi(2) * mu2 + mu4,
            eps.powi(5) + 5.0 * eps * mu4 + 10.0 * eps.powi(3) * mu2,
            eps.powi(6) + 15.0 * eps.powi(4) * mu2 + 15.0 * eps.powi(2) * mu4 + mu6,
        ]
    }

    #[test]
    fn discrete() {
        let levels = vec![-1.0, 0.0, 2.0];
        let weights = vec![0.5, 0.5, 0.25];
        let dos = models::discrete(&levels, &weights);
        for order in 0..=6 {
            let moment = compute_moment(&dos, order);
            let moment_ref = levels
                .iter()
                .zip(&weights)
                .map(|(level, weight)| level.powi(order) * weight)
                .sum();
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn flat() {
        let eps = 0.5f64;
        let d = 2.0f64;
        let delta = 0.2f64;

        use std::f64::consts::PI;
        let a = d / delta;
        let mu2 = d.powi(2) * (PI.powi(2) + a.powi(2)) / (3.0 * a.powi(2));
        let mu4 = mu2 * (3.0 * a.powi(2) + 7.0 * PI.powi(2)) * delta.powi(2) / 5.0;
        let mu6 = (3.0 * a.powi(6)
            + 21.0 * a.powi(4) * PI.powi(2)
            + 49.0 * a.powi(2) * PI.powi(4)
            + 31.0 * PI.powi(6))
            * delta.powi(6)
            / 21.0;
        let moments_ref = central_to_moments(eps, mu2, mu4, mu6);

        let dos = models::flat(eps, d, delta);
        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn flat_sharp() {
        let eps = 0.5f64;
        let d = 2.0f64;

        // δ -> 0 limits of the central moments used in flat()
        let moments_ref =
            central_to_moments(eps, d.powi(2) / 3.0, d.powi(4) / 5.0, d.powi(6) / 7.0);

        let dos = models::flat(eps, d, 0.0);
        assert_eq!(dos.support(), Some((eps - d, eps + d)));
        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn gaussian() {
        let eps = 1.0f64;
        let sigma = 3.0f64;
        let moments_ref = central_to_moments(
            eps,
            sigma.powi(2),
            3.0 * sigma.powi(4),
            15.0 * sigma.powi(6),
        );

        let dos = models::gaussian(eps, sigma);
        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn semicircle() {
        let eps = 0.5f64;
        let t = 2.0f64;
        let r = 2.0 * t;
        let moments_ref = central_to_moments(
            eps,
            r.powi(2) / 4.0,
            r.powi(4) / 8.0,
            5.0 * r.powi(6) / 64.0,
        );

        let dos = models::semicircle(eps, r);
        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn chain() {
        let eps = 0.5f64;
        let t = 2.0f64;
        let moments_ref =
            central_to_moments(eps, 2.0 * t.powi(2), 6.0 * t.powi(4), 20.0 * t.powi(6));

        let dos = models::chain(eps, t);
        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn square_regular_near_singularity() {
        use crate::ContinuousSF;
        use crate::models::SquareDOS;

        let (eps, t) = (0.5f64, 2.0f64);
        let dos = SquareDOS::new(eps, t);
        let prefactor = 1.0 / (2.0 * std::f64::consts::PI.powi(2) * t);

        assert_eq!(dos.regular(eps), 0.0);
        for e in 2..=16 {
            for sign in [-1.0f64, 1.0] {
                let omega = eps + sign * 10f64.powi(-e);
                // Leading order of the expansion of R(ω) around the singularity
                let ax = 10f64.powi(-e) / (4.0 * t);
                let ref_value = prefactor * 0.25 * ax.powi(2) * ((4.0 / ax).ln() - 1.0);
                assert_relative_eq!(dos.regular(omega), ref_value, max_relative = 1e-4);
            }
        }
    }

    #[test]
    fn square() {
        let eps = 0.5f64;
        let t = 2.0f64;
        let moments_ref =
            central_to_moments(eps, 4.0 * t.powi(2), 36.0 * t.powi(4), 400.0 * t.powi(6));

        let dos = models::square(eps, t);
        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }
}
