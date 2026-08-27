//! Model densities of states commonly used in physics.

use crate::discrete::{DiscreteSF, Resonance};
use crate::util::fermi;
use crate::{ContinuousSF, SingularLaw, Singularity, SpectralFunction};

use special::Elliptic;
use std::f64::consts::{PI, SQRT_2};
use std::rc::Rc;

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

/// Gaussian density of states.
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

/// Semicircle (Wigner) density of states.
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
// Power-law band edge
//

/// Power-law threshold density of states.
struct PowerLawDOS {
    eps: f64,
    r: f64,
    /// Band edges.
    edges: [f64; 2],
    singularity: [Singularity; 1],
    prefactor: f64,
}
impl PowerLawDOS {
    fn new(eps: f64, r: f64, w: f64) -> PowerLawDOS {
        assert!(r > -1.0, "asymptotics exponent must satisfy r > -1");
        assert!(w > 0.0, "bandwidth must be positive");
        let prefactor = (r + 1.0) / (w.powf(r + 1.0));
        let sing_law = if r < 0.0 {
            SingularLaw::Power {
                a: -r,
                c: prefactor,
                l: 0.0,
            }
        } else {
            SingularLaw::Finite
        };
        PowerLawDOS {
            eps,
            r,
            edges: [eps, eps + w],
            singularity: [Singularity {
                position: eps,
                law: sing_law,
            }],
            prefactor,
        }
    }
    fn value(&self, omega: f64) -> f64 {
        self.prefactor * (omega - self.eps).powf(self.r)
    }
}
impl ContinuousSF for PowerLawDOS {
    fn support(&self) -> (f64, f64) {
        self.edges.into()
    }
    fn regular(&self, omega: f64) -> f64 {
        if self.r < 1.0 { 0.0 } else { self.value(omega) }
    }
    fn singularities(&self) -> &[Singularity] {
        &self.singularity
    }
    fn asymptotics(&self, p: usize, omega: f64) -> f64 {
        debug_assert!(p == 0);
        if self.r < 1.0 { self.value(omega) } else { 0.0 }
    }
    fn asympt_int(&self, _p: usize) -> f64 {
        if self.r < 1.0 { 1.0 } else { 0.0 }
    }
}

/// Returns the normalized density of states with a power-law band edge
/// located at `eps`. The power-law asymptotics is characterized by
/// the exponent `r > -1`,
/// $$
///     A(\omega) = \frac{r+1}{w^{r+1}} (\omega - \epsilon)^r
///         \theta(\omega-\epsilon) \theta(w - (\omega-\epsilon)),
/// $$
/// and `w` is the bandwidth.
pub fn powerlaw(eps: f64, r: f64, w: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(PowerLawDOS::new(eps, r, w))
}

//
// Pseudogap DOS
//

/// Density of states with a pseudogap.
struct PseudogapDOS {
    eps: f64,
    r: f64,
    /// Band edges.
    edges: [f64; 2],
    singularity: [Singularity; 1],
    prefactor: f64,
}
impl PseudogapDOS {
    fn new(eps: f64, r: f64, d: f64) -> PseudogapDOS {
        assert!(r > 0.0, "asymptotics exponent must be positive");
        assert!(d > 0.0, "bandwidth must be positive");
        PseudogapDOS {
            eps,
            r,
            edges: [eps - d, eps + d],
            singularity: [Singularity {
                position: eps,
                law: SingularLaw::Finite,
            }],
            prefactor: (r + 1.0) / (2.0 * d.powf(r + 1.0)),
        }
    }
    fn value(&self, omega: f64) -> f64 {
        self.prefactor * (omega - self.eps).abs().powf(self.r)
    }
}
impl ContinuousSF for PseudogapDOS {
    fn support(&self) -> (f64, f64) {
        self.edges.into()
    }
    fn regular(&self, omega: f64) -> f64 {
        if self.r < 1.0 { 0.0 } else { self.value(omega) }
    }
    fn singularities(&self) -> &[Singularity] {
        &self.singularity
    }
    fn asymptotics(&self, p: usize, omega: f64) -> f64 {
        debug_assert!(p == 0);
        if self.r < 1.0 { self.value(omega) } else { 0.0 }
    }
    fn asympt_int(&self, _p: usize) -> f64 {
        if self.r < 1.0 { 1.0 } else { 0.0 }
    }
}

/// Returns the normalized density of states with a pseudogap located at `eps` with
/// power-law asymptotics characterized by the positive exponent `r`,
/// $$
///     A(\omega) = \frac{1+r}{2 d^{1+r}} |\omega - \epsilon|^r
///         \theta(d^2 - (\omega-\epsilon)^2).
/// $$
/// and the half-bandwidth `d`.
pub fn pseudogap(eps: f64, r: f64, d: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(PseudogapDOS::new(eps, r, d))
}

//
// Chain DOS
//

/// Density of states of a linear chain.
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
/// and the local energy level `eps` (derived from the dispersion law
/// $\varepsilon(k) = \epsilon - 2t\cos(k)$),
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

/// Density of states of a square lattice.
struct SquareDOS {
    eps: f64,
    t: f64,
    /// Band edges.
    edges: [f64; 2],
    /// Position of the logarithmic van Hove singularity, which sits at the band center.
    singularity: [Singularity; 1],
    prefactor: f64,
}
impl SquareDOS {
    /// Value of $|\omega-\epsilon|/(4t)$ below which the regular part is evaluated using a
    /// series expansion. At the threshold both expressions agree to within a relative
    /// error of $10^{-5}$.
    const SERIES_THRESHOLD: f64 = 1e-3;

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
        } else if ax < Self::SERIES_THRESHOLD {
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
/// `t` and the local energy level `eps` (derived from the dispersion law
/// $\varepsilon(k) = \epsilon - 2t[\cos(k_x) + \cos(k_y)]$),
/// $$
///     A(\omega) = \frac{1}{2\pi^2 t} K\left(1 - \left(\frac{\omega-\epsilon}{4t}\right)^2\right)
///         \theta((4t)^2 - (\omega-\epsilon)^2),
/// $$
/// where $K(m)$ is the complete elliptic integral of the first kind.
pub fn square(eps: f64, t: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(SquareDOS::new(eps, t))
}

//
// Triangular lattice DOS
//

/// Density of states of a triangular lattice.
struct TriangularDOS {
    eps: f64,
    t: f64,
    /// Band edges.
    edges: [f64; 2],
    /// Position of the logarithmic van Hove singularity.
    singularity: [Singularity; 1],
    prefactor: f64,
}
impl TriangularDOS {
    /// Value of $|\omega-\epsilon-2t|/(8|t|)$ below which the regular part is evaluated
    /// using a series expansion. At the threshold both expressions agree to within a
    /// relative error of $10^{-5}$.
    const SERIES_THRESHOLD: f64 = 1e-3;

    fn new(eps: f64, t: f64) -> TriangularDOS {
        assert!(t != 0.0, "hopping constant must be non-zero");
        let omega_min = eps + (-6.0 * t).min(3.0 * t);
        let omega_max = eps + (-6.0 * t).max(3.0 * t);
        let prefactor = 1.0 / (PI.powi(2) * t.abs());
        TriangularDOS {
            eps,
            t,
            edges: [omega_min, omega_max],
            singularity: [Singularity {
                position: eps + 2.0 * t,
                law: SingularLaw::Log {
                    c: prefactor * 0.75,
                    l: prefactor * 0.75 * (8.0 * t.abs()).ln(),
                },
            }],
            prefactor,
        }
    }
}
impl ContinuousSF for TriangularDOS {
    fn support(&self) -> (f64, f64) {
        self.edges.into()
    }
    fn regular(&self, omega: f64) -> f64 {
        let ax = ((omega - self.singularity[0].position) / (8.0 * self.t)).abs();
        if ax == 0.0 {
            0.0
        } else if ax < Self::SERIES_THRESHOLD {
            // Close to the singularity, K(z_1/z_0)/\sqrt{z_0} ≈ -0.75 ln(ax) nearly cancels
            // the logarithm, and evaluating the difference directly is catastrophically
            // inaccurate. Worse, 1 - z_1/z_0 = O(d^3) drowns in the round-off of z_1/z_0 as
            // d -> 0. Use the expansion in d = (ω-ε)/t - 2 instead,
            // d(3Λ/16 - 9/32) + d^2(15Λ/128 - 99/512) + d^3(21Λ/256 - 75/512)
            // + O(d^4 ln d), where Λ = -ln(ax).
            let d = (omega - self.singularity[0].position) / self.t;
            let l = -ax.ln();
            self.prefactor
                * d
                * ((3.0 / 16.0) * l - 9.0 / 32.0
                    + d * ((15.0 / 128.0) * l - 99.0 / 512.0
                        + d * ((21.0 / 256.0) * l - 75.0 / 512.0)))
        } else {
            let x = (omega - self.eps) / self.t;
            let s = (3.0 - x).sqrt();
            let a = 3.0 + 2.0 * s - (x / 2.0).powi(2);
            let b = 4.0 * s;
            let (z0, z1) = if x >= 2.0 { (a, b) } else { (b, a) };
            self.prefactor * (((z1 / z0).elliptic_k() / z0.sqrt()) + 0.75 * ax.ln())
        }
    }
    fn singularities(&self) -> &[Singularity] {
        &self.singularity
    }
    fn asymptotics(&self, p: usize, omega: f64) -> f64 {
        debug_assert!(p == 0);
        let x = (omega - self.singularity[0].position) / (8.0 * self.t);
        -self.prefactor * 0.75 * x.abs().ln()
    }
    fn asympt_int(&self, _p: usize) -> f64 {
        9.0 / (4.0 * PI.powi(2)) * (3.0 + std::f64::consts::LN_2)
    }
}

/// Returns the normalized density of states of a triangular lattice with the hopping constant
/// `t` and the local energy level `eps` (derived from the dispersion law
/// $\varepsilon(k) = \epsilon - 2t[\cos(k_x) + 2\cos(k_x/2)\cos(\sqrt{3} k_y/2)]$),
/// $$
///     A(\omega) = \frac{1}{\pi^2|t|\sqrt{z_0(\omega-\epsilon)}}
///         K\left(\frac{z_1(\omega-\epsilon)}{z_0(\omega-\epsilon)}\right)
///         \theta(\omega - \omega_\mathrm{min}) \theta(\omega_\mathrm{max} - \omega),
/// $$
/// where $K(m)$ is the complete elliptic integral of the first kind and
/// $$
///     z_0(\nu) = \begin{cases}
///         3 + 2\sqrt{3 - \nu/t} - \frac{\nu^2}{4t^2}, & \nu/t \geq 2, \\\\
///         4\sqrt{3 - \nu/t}, & \nu/t < 2,
///     \end{cases} \\\\
///     z_1(\nu) = \begin{cases}
///         4\sqrt{3 - \nu/t}, & \nu/t \geq 2, \\\\
///         3 + 2\sqrt{3 - \nu/t} - \frac{\nu^2}{4t^2}, & \nu/t < 2.
///     \end{cases}
/// $$
/// Positions of the band edges are $\omega_\mathrm{min} = \epsilon+\min(-6t, 3t)$,
/// $\omega_\mathrm{max} = \epsilon+\max(-6t, 3t)$.
pub fn triangular(eps: f64, t: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(TriangularDOS::new(eps, t))
}

//
// Honeycomb lattice DOS
//

/// Density of states of a honeycomb lattice.
///
/// The two bands $\pm t|f(k)|$ share the single-band dispersion of the triangular lattice
/// via $|f(k)|^2 = 3 - x(k)$, $x$ being the triangular lattice energy measured from
/// $\epsilon$ in units of $t$. The DOS is, therefore, that of the triangular lattice at
/// $x = 3 - ((\omega-\epsilon)/t)^2$, reweighted by the Jacobian $|\omega-\epsilon|/t$ of
/// that substitution.
struct HoneycombDOS {
    eps: f64,
    t: f64,
    /// Band edges.
    edges: [f64; 2],
    /// Positions of the logarithmic van Hove singularities.
    singularities: [Singularity; 2],
    prefactor: f64,
}
impl HoneycombDOS {
    /// Value of $|t^2-(\omega-\epsilon)^2|/(8t^2)$ below which the regular part is
    /// evaluated using a series expansion. At the threshold both expressions agree to
    /// within a relative error of $10^{-8}$.
    const SERIES_THRESHOLD: f64 = 1e-3;

    fn new(eps: f64, t: f64) -> HoneycombDOS {
        assert!(t > 0.0, "hopping constant must be positive");
        let prefactor = 1.0 / (PI.powi(2) * t);
        HoneycombDOS {
            eps,
            t,
            edges: [eps - 3.0 * t, eps + 3.0 * t],
            singularities: [
                Singularity {
                    position: eps - t,
                    law: SingularLaw::Log {
                        c: prefactor * 0.75,
                        l: prefactor * 0.75 * (4.0 * t).ln(),
                    },
                },
                Singularity {
                    position: eps + t,
                    law: SingularLaw::Log {
                        c: prefactor * 0.75,
                        l: prefactor * 0.75 * (4.0 * t).ln(),
                    },
                },
            ],
            prefactor,
        }
    }
}
impl ContinuousSF for HoneycombDOS {
    fn support(&self) -> (f64, f64) {
        self.edges.into()
    }
    fn regular(&self, omega: f64) -> f64 {
        let ax = ((omega - self.eps) / self.t).abs();
        // Distance to the singular points in the expansion variable of the triangular
        // lattice DOS, d = x - 2. The factored form is free of cancellation for ax ≈ 1.
        let d = (1.0 - ax) * (1.0 + ax);
        // Argument of the logarithm subtracted by the triangular lattice DOS
        let y = 0.125 * d.abs();
        if y == 0.0 {
            -self.prefactor * 0.75 * std::f64::consts::LN_2
        } else if y < Self::SERIES_THRESHOLD {
            // Close to a singularity, K(z_1/z_0)/\sqrt{z_0} ≈ -0.75 ln(y) nearly cancels the
            // logarithm, and evaluating the difference directly is catastrophically
            // inaccurate. Worse, 1 - z_1/z_0 = O(d^3) drowns in the round-off of z_1/z_0 as
            // d -> 0. Use the expansion of the triangular lattice regular part,
            // d(3Λ/16 - 9/32) + d^2(15Λ/128 - 99/512) + d^3(21Λ/256 - 75/512)
            // + O(d^4 ln d), where Λ = -ln(y), and pick up the leftover logarithm through
            // ax - 1.
            let l = -y.ln();
            let s = d
                * ((3.0 / 16.0) * l - 9.0 / 32.0
                    + d * ((15.0 / 128.0) * l - 99.0 / 512.0
                        + d * ((21.0 / 256.0) * l - 75.0 / 512.0)));
            self.prefactor * (ax * s - 0.75 * (l * d / (1.0 + ax) + std::f64::consts::LN_2))
        } else {
            let a = 3.0 + 2.0 * ax - 0.25 * (3.0 - ax.powi(2)).powi(2);
            let b = 4.0 * ax;
            let (z0, z1) = if ax <= 1.0 { (a, b) } else { (b, a) };
            self.prefactor * (ax * ((z1 / z0).elliptic_k() / z0.sqrt()) + 0.75 * (0.5 * y).ln())
        }
    }
    fn singularities(&self) -> &[Singularity] {
        &self.singularities
    }
    fn asymptotics(&self, p: usize, omega: f64) -> f64 {
        debug_assert!(p <= 1);
        let x = (omega - self.singularities[p].position) / (4.0 * self.t);
        -self.prefactor * 0.75 * x.abs().ln()
    }
    fn asympt_int(&self, _p: usize) -> f64 {
        (9.0 + 3.0 * std::f64::consts::LN_2) / (2.0 * PI.powi(2))
    }
}

/// Returns the normalized density of states of a honeycomb lattice with the hopping constant
/// `t` and the local energy level `eps` (derived from the two-band dispersion law
/// $\varepsilon_\pm(k) = \epsilon \pm t\sqrt{3 + 2[\cos(k_x) + 2\cos(k_x/2)\cos(\sqrt{3}k_y/2)]}$),
/// $$
///     A(\omega) = \frac{|\omega-\epsilon|}{\pi^2 t^2\sqrt{z_0(\omega-\epsilon)}}
///         K\left(\frac{z_1(\omega-\epsilon)}{z_0(\omega-\epsilon)}\right)
///         \theta((3t)^2 - (\omega-\epsilon)^2),
/// $$
/// where $K(m)$ is the complete elliptic integral of the first kind and
/// $$
///     z_0(\nu) = \begin{cases}
///         3 + 2\frac{|\nu|}{t} - \frac{1}{4}\left(3 - \frac{\nu^2}{t^2}\right)^2, & |\nu|/t \leq 1, \\\\
///         4\frac{|\nu|}{t}, & |\nu|/t > 1,
///     \end{cases} \\\\
///     z_1(\nu) = \begin{cases}
///         4\frac{|\nu|}{t}, & |\nu|/t \leq 1, \\\\
///         3 + 2\frac{|\nu|}{t} - \frac{1}{4}\left(3 - \frac{\nu^2}{t^2}\right)^2, & |\nu|/t > 1.
///     \end{cases}
/// $$
pub fn honeycomb(eps: f64, t: f64) -> SpectralFunction {
    SpectralFunction::from_continuous(HoneycombDOS::new(eps, t))
}

//
// Kagome lattice DOS
//

/// Returns the normalized density of states of the Kagome lattice with the hopping constant
/// `t` and the local energy level `eps` (derived from the two-band dispersion law
/// $\varepsilon_\pm(k) = \epsilon + t \left(1 \pm
/// \sqrt{3 + 2[\cos(k_x) + 2\cos(k_x/2)\cos(\sqrt{3}k_y/2)]}\right)$,
/// and a third flat band $\varepsilon_F(k) = \epsilon - 2t$),
/// $$
///     A(\omega) = \frac{1}{3}\delta(\omega-\epsilon+2t) +
///         \frac{2}{3}\frac{|\omega-\epsilon-t|}{\pi^2 t^2\sqrt{z_0(\omega-\epsilon-t)}}
///         K\left(\frac{z_1(\omega-\epsilon-t)}{z_0(\omega-\epsilon-t)}\right)
///         \theta\left((3t)^2 - (\omega-\epsilon-t)^2\right),
/// $$
/// where $K(m)$ is the complete elliptic integral of the first kind and
/// $$
///     z_0(\nu) = \begin{cases}
///         3 + 2\frac{|\nu|}{|t|} - \frac{1}{4}\left(3 - \frac{\nu^2}{t^2}\right)^2, & |\nu|/|t| \leq 1, \\\\
///         4\frac{|\nu|}{|t|}, & |\nu|/|t| > 1,
///     \end{cases} \\\\
///     z_1(\nu) = \begin{cases}
///         4\frac{|\nu|}{|t|}, & |\nu|/|t| \leq 1, \\\\
///         3 + 2\frac{|\nu|}{|t|} - \frac{1}{4}\left(3 - \frac{\nu^2}{t^2}\right)^2, & |\nu|/|t| > 1.
///     \end{cases}
/// $$
/// Positions of the band edges are $\omega_\mathrm{min} = \epsilon+\min(-2t, 4t)$,
/// $\omega_\mathrm{max} = \epsilon+\max(-2t, 4t)$.
pub fn kagome(eps: f64, t: f64) -> SpectralFunction {
    assert!(t != 0.0, "hopping constant must be non-zero");
    SpectralFunction::from_discrete_continuous(
        DiscreteSF::one_resonance(eps - 2.0 * t, 1.0 / 3.0),
        vec![(Rc::new(HoneycombDOS::new(eps + t, t.abs())), 2.0 / 3.0)],
    )
}

//
// Lieb lattice DOS
//

/// Density of states of the Lieb lattice.
///
/// The two dispersive bands $\pm 2t\sqrt{\cos(k_x)^2 + \cos(k_y)^2}$ share the dispersion
/// of the square lattice via $u^2 - 1 = [\cos(2k_x) + \cos(2k_y)]/2$, $u$ being the energy
/// measured from $\epsilon$ in units of $2t$. The DOS is, therefore, that of the square
/// lattice at $x = u^2 - 1$, reweighted by $|u|$, half the Jacobian of that substitution,
/// the other half splitting the weight between the two bands.
struct LiebDOS {
    eps: f64,
    t: f64,
    /// Band edges.
    edges: [f64; 2],
    /// Positions of the logarithmic van Hove singularities.
    singularities: [Singularity; 2],
    prefactor: f64,
}
impl LiebDOS {
    /// Value of $|((\omega-\epsilon)/(2t))^2 - 1|$ below which the regular part is
    /// evaluated using a series expansion. At the threshold both expressions agree to
    /// within a relative error of $10^{-12}$.
    const SERIES_THRESHOLD: f64 = 1e-3;

    fn new(eps: f64, t: f64) -> LiebDOS {
        assert!(t > 0.0, "hopping constant must be positive");
        let prefactor = 1.0 / (PI.powi(2) * t);
        let half_width = 2.0 * SQRT_2 * t;
        LiebDOS {
            eps,
            t,
            edges: [eps - half_width, eps + half_width],
            singularities: [
                Singularity {
                    position: eps - 2.0 * t,
                    law: SingularLaw::Log {
                        c: prefactor,
                        l: prefactor * (2.0 * t).ln(),
                    },
                },
                Singularity {
                    position: eps + 2.0 * t,
                    law: SingularLaw::Log {
                        c: prefactor,
                        l: prefactor * (2.0 * t).ln(),
                    },
                },
            ],
            prefactor,
        }
    }
}
impl ContinuousSF for LiebDOS {
    fn support(&self) -> (f64, f64) {
        self.edges.into()
    }
    fn regular(&self, omega: f64) -> f64 {
        let ax = ((omega - self.eps) / (2.0 * self.t)).abs();
        // Square lattice energy variable, which vanishes at the singular points. The
        // factored form is free of cancellation for ax ≈ 1.
        let y = (ax - 1.0) * (ax + 1.0);
        if y == 0.0 {
            2.0 * self.prefactor * std::f64::consts::LN_2
        } else if y.abs() < Self::SERIES_THRESHOLD {
            // 1 - y^2 rounds to unity for small y, making elliptic_k() overflow. Expand it
            // instead as K(1-y^2) = Λ + (y^2/4)(Λ-1) + O(y^4 ln y), where Λ = ln(4/|y|),
            // and pick up the leftover logarithm through ax - 1 = y/(1+ax).
            let l = (4.0 / y.abs()).ln();
            self.prefactor
                * ((y / (1.0 + ax)) * l
                    + 2.0 * std::f64::consts::LN_2
                    + 0.25 * ax * y.powi(2) * (l - 1.0))
        } else {
            self.prefactor * (ax * (1.0 - y.powi(2)).elliptic_k() + y.abs().ln())
        }
    }
    fn singularities(&self) -> &[Singularity] {
        &self.singularities
    }
    fn asymptotics(&self, p: usize, omega: f64) -> f64 {
        debug_assert!(p <= 1);
        let x = (omega - self.singularities[p].position) / (2.0 * self.t);
        -self.prefactor * x.abs().ln()
    }
    fn asympt_int(&self, _p: usize) -> f64 {
        4.0 * (SQRT_2 - (1.0 + SQRT_2).ln()) / PI.powi(2)
    }
}

/// Returns the normalized density of states of the Lieb lattice with the hopping constant
/// `t` and the local energy level `eps` (derived from the two-band dispersion law
/// $\varepsilon_\pm(k) = \epsilon \pm 2t\sqrt{\cos(k_x)^2 + \cos(k_y)^2}$
/// and a third flat band $\varepsilon_F(k) = \epsilon$),
/// $$
///     A(\omega) = \frac{1}{3} \delta(\omega-\epsilon) +
///         \frac{|\omega-\epsilon|}{3\pi^2 t^2}
///         K\left(1 - \left(\left(\frac{\omega-\epsilon}{2t}\right)^2 - 1\right)^2\right)
///         \theta\left(8t^2 - (\omega-\epsilon)^2\right),
/// $$
/// where $K(m)$ is the complete elliptic integral of the first kind. Positions of the band
/// edges are $\epsilon \pm 2\sqrt{2}t$.
///
/// The weight of the flat band is that of one band out of three, as carried by the DOS
/// averaged over the three sites of the unit cell. It is not the flat band weight of the
/// local DOS, which vanishes on the corner site and equals $1/2$ on the rim sites.
pub fn lieb(eps: f64, t: f64) -> SpectralFunction {
    SpectralFunction::from_discrete_continuous(
        DiscreteSF::one_resonance(eps, 1.0 / 3.0),
        vec![(Rc::new(LiebDOS::new(eps, t)), 2.0 / 3.0)],
    )
}

#[cfg(test)]
mod tests {
    use crate::{SpectralFunction, models};
    use approx::assert_relative_eq;

    fn compute_moment(dos: &SpectralFunction, order: i32) -> f64 {
        dos.integrate(|omega| omega.powi(order), None).unwrap()
    }

    /// Moments $M_n = \sum_{k=0}^n \binom{n}{k} \epsilon^{n-k} \mu_k$ of orders 0..=6 of a
    /// spectral function whose central moments $\mu_1, \ldots, \mu_6$ about `eps` are `mu`.
    fn central_to_moments(eps: f64, mu: [f64; 6]) -> Vec<f64> {
        let mu: Vec<f64> = std::iter::once(1.0).chain(mu).collect();
        (0..mu.len())
            .map(|n| {
                (0..=n)
                    .map(|k| binomial(n, k) * eps.powi((n - k) as i32) * mu[k])
                    .sum()
            })
            .collect()
    }

    fn binomial(n: usize, k: usize) -> f64 {
        (1..=k).map(|i| (n - k + i) as f64 / i as f64).product()
    }

    /// Central moments of a spectral function symmetric about its center, the odd ones
    /// vanishing.
    fn even_central_moments(mu2: f64, mu4: f64, mu6: f64) -> [f64; 6] {
        [0.0, mu2, 0.0, mu4, 0.0, mu6]
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
        let moments_ref = central_to_moments(eps, even_central_moments(mu2, mu4, mu6));

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
        let moments_ref = central_to_moments(
            eps,
            even_central_moments(d.powi(2) / 3.0, d.powi(4) / 5.0, d.powi(6) / 7.0),
        );

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
            even_central_moments(sigma.powi(2), 3.0 * sigma.powi(4), 15.0 * sigma.powi(6)),
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
            even_central_moments(r.powi(2) / 4.0, r.powi(4) / 8.0, 5.0 * r.powi(6) / 64.0),
        );

        let dos = models::semicircle(eps, r);
        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn powerlaw() {
        let eps = 0.5f64;
        let w = 2.0f64;

        // r < 0 makes A(ω) diverge at the band edge, -1 < r < 1 puts the whole of it into
        // the singular part, and r >= 1 into the regular one
        for r in [-0.5f64, 0.0, 0.5, 1.0, 2.5] {
            let dos = models::powerlaw(eps, r, w);
            assert_eq!(dos.support(), Some((eps, eps + w)));

            // A(ω) = (r+1)(ω-ε)^r / w^{r+1}
            let a = |omega: f64| (r + 1.0) * (omega - eps).powf(r) / w.powf(r + 1.0);
            for omega in [eps + 1e-3, eps + 0.7, eps + 1.3, eps + w] {
                assert_relative_eq!(dos.continuous_at(omega), a(omega), max_relative = 1e-12);
            }
            // The band edge itself diverges for r < 0 only
            assert_eq!(
                dos.continuous_at(eps),
                if r < 0.0 { f64::INFINITY } else { a(eps) }
            );
            assert_eq!(dos.continuous_at(eps - 0.5 * w), 0.0);
            assert_eq!(dos.continuous_at(eps + 1.5 * w), 0.0);

            // μ_n = (r+1)w^n/(n+r+1), including the odd orders as A(ω) is one-sided
            let mu = |n: i32| (r + 1.0) * w.powi(n) / (f64::from(n) + r + 1.0);
            let moments_ref = central_to_moments(eps, [mu(1), mu(2), mu(3), mu(4), mu(5), mu(6)]);
            for (order, moment_ref) in moments_ref.iter().enumerate() {
                let moment = compute_moment(&dos, order as i32);
                assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
            }
        }
    }

    #[test]
    fn pseudogap() {
        let eps = 0.5f64;
        let d = 2.0f64;

        // r < 1 puts the whole of A(ω) into the singular part, r >= 1 into the regular one
        for r in [0.5f64, 1.0, 2.5] {
            let dos = models::pseudogap(eps, r, d);
            assert_eq!(dos.support(), Some((eps - d, eps + d)));

            // A(ω) = (1+r)|ω-ε|^r / (2d^{1+r})
            let a = |omega: f64| (1.0 + r) * (omega - eps).abs().powf(r) / (2.0 * d.powf(1.0 + r));
            for omega in [eps - d, eps - 0.7, eps, eps + 1.3, eps + d] {
                assert_relative_eq!(dos.continuous_at(omega), a(omega), max_relative = 1e-12);
            }
            assert_eq!(dos.continuous_at(eps - 1.5 * d), 0.0);

            // μ_n = (1+r)d^n/(n+1+r) for even n, the odd ones vanishing by symmetry
            let mu = |n: i32| (1.0 + r) * d.powi(n) / (f64::from(n) + 1.0 + r);
            let moments_ref = central_to_moments(eps, even_central_moments(mu(2), mu(4), mu(6)));
            for (order, moment_ref) in moments_ref.iter().enumerate() {
                let moment = compute_moment(&dos, order as i32);
                assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
            }
        }
    }

    #[test]
    fn chain() {
        let eps = 0.5f64;
        let t = 2.0f64;
        let moments_ref = central_to_moments(
            eps,
            even_central_moments(2.0 * t.powi(2), 6.0 * t.powi(4), 20.0 * t.powi(6)),
        );

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
        let moments_ref = central_to_moments(
            eps,
            even_central_moments(4.0 * t.powi(2), 36.0 * t.powi(4), 400.0 * t.powi(6)),
        );

        let dos = models::square(eps, t);
        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn triangular() {
        let eps = 0.5f64;
        // The band edges swap places with the sign of t, and μ_n of odd order do not vanish
        for t in [2.0f64, -0.7] {
            let moments_ref = central_to_moments(
                eps,
                [
                    0.0,
                    6.0 * t.powi(2),
                    -12.0 * t.powi(3),
                    90.0 * t.powi(4),
                    -360.0 * t.powi(5),
                    2040.0 * t.powi(6),
                ],
            );

            let dos = models::triangular(eps, t);
            assert_eq!(
                dos.support(),
                Some((eps + (-6.0 * t).min(3.0 * t), eps + (-6.0 * t).max(3.0 * t)))
            );
            for (order, moment_ref) in moments_ref.iter().enumerate() {
                let moment = compute_moment(&dos, order as i32);
                assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
            }
        }
    }

    #[test]
    fn triangular_regular_near_singularity() {
        use crate::ContinuousSF;
        use crate::models::TriangularDOS;

        let (eps, t) = (0.5f64, 2.0f64);
        let dos = TriangularDOS::new(eps, t);
        let prefactor = 1.0 / (std::f64::consts::PI.powi(2) * t);

        // R(ε+(2+d)t) / prefactor at d = \pm 2^-e, evaluated at 200 bits of precision.
        // The series expansion takes over at |d| < 8e-3.
        let ref_values = [
            (1, 0.168390121477016, -0.0938944124091616),
            (4, 0.0408118425415933, -0.0378744776978792),
            (7, 0.00799423604376689, -0.00791868158506950),
            (10, 0.00137611250776708, -0.00137446721459707),
            (20, 2.58250123354371e-6, -2.58249818693858e-6),
        ];

        assert_eq!(dos.regular(eps + 2.0 * t), 0.0);
        for (e, ref_plus, ref_minus) in ref_values {
            for (sign, ref_value) in [(1.0, ref_plus), (-1.0, ref_minus)] {
                let d = sign * 2f64.powi(-e);
                let omega = eps + (2.0 + d) * t;
                assert_relative_eq!(
                    dos.regular(omega),
                    prefactor * ref_value,
                    max_relative = 1e-6
                );
            }
        }
    }

    #[test]
    fn honeycomb() {
        let eps = 0.5f64;
        let t = 2.0f64;
        let moments_ref = central_to_moments(
            eps,
            even_central_moments(3.0 * t.powi(2), 15.0 * t.powi(4), 93.0 * t.powi(6)),
        );

        let dos = models::honeycomb(eps, t);
        assert_eq!(dos.support(), Some((eps - 3.0 * t, eps + 3.0 * t)));

        // Both van Hove singularities diverge, while the Dirac point between them is a
        // zero of A(ω) that only the two logarithms cancelling leaves behind.
        assert_eq!(dos.continuous_at(eps - t), f64::INFINITY);
        assert_eq!(dos.continuous_at(eps + t), f64::INFINITY);
        assert_relative_eq!(dos.continuous_at(eps), 0.0, epsilon = 1e-15);

        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn honeycomb_regular_near_singularity() {
        use crate::ContinuousSF;
        use crate::models::HoneycombDOS;

        let (eps, t) = (0.5f64, 2.0f64);
        let dos = HoneycombDOS::new(eps, t);
        let prefactor = 1.0 / (std::f64::consts::PI.powi(2) * t);

        // Unlike its triangular lattice counterpart, R(ω) does not vanish at the singular
        // points: the two logarithms subtracted there differ in scale by a factor of 2.
        for omega_p in [eps - t, eps + t] {
            assert_relative_eq!(
                dos.regular(omega_p),
                -prefactor * 0.75 * std::f64::consts::LN_2,
                max_relative = 1e-14
            );
        }

        // R(ε \pm (1+δ)t) / prefactor at δ = \pm 2^-e, evaluated at 200 bits of precision.
        // The series expansion takes over at |δ| < 4e-3.
        let ref_values = [
            (1, 0.0943101136887941, -1.29999066584090),
            (4, -0.389118215772244, -0.654568957739481),
            (7, -0.497232083057124, -0.542574386667219),
            (10, -0.516265866766906, -0.523456614923520),
            (20, -0.519854395427088, -0.519866375415643),
        ];

        for (e, ref_plus, ref_minus) in ref_values {
            for (sign, ref_value) in [(1.0f64, ref_plus), (-1.0, ref_minus)] {
                let delta = sign * 2f64.powi(-e);
                // R(ω) is even about ε, both singularities being equivalent
                for band in [-1.0f64, 1.0] {
                    let omega = eps + band * (1.0 + delta) * t;
                    assert_relative_eq!(
                        dos.regular(omega),
                        prefactor * ref_value,
                        max_relative = 1e-8
                    );
                }
            }
        }
    }

    #[test]
    fn kagome() {
        let eps = 0.5f64;
        // The flat band swaps the band edge it is degenerate with as t changes sign, and
        // μ_n of odd order do not vanish.
        for t in [2.0f64, -0.5] {
            let moments_ref = central_to_moments(
                eps,
                [
                    0.0,
                    4.0 * t.powi(2),
                    4.0 * t.powi(3),
                    28.0 * t.powi(4),
                    60.0 * t.powi(5),
                    264.0 * t.powi(6),
                ],
            );

            let dos = models::kagome(eps, t);
            let omega_min = eps + (-2.0 * t).min(4.0 * t);
            let omega_max = eps + (-2.0 * t).max(4.0 * t);
            assert_eq!(dos.support(), Some((omega_min, omega_max)));
            assert_relative_eq!(dos.total_weight(), 1.0, max_relative = 1e-14);

            // The flat band is a δ-function of weight 1/3 at the edge it touches at k = 0
            assert_eq!(dos.discrete().len(), 1);
            assert_eq!(dos.discrete().resonances()[0].eps, eps - 2.0 * t);
            assert_relative_eq!(
                dos.discrete().resonances()[0].weight,
                1.0 / 3.0,
                max_relative = 1e-14
            );

            // The dispersive bands are the honeycomb ones shifted by t, so the van Hove
            // singularities sit at ε and ε+2t with the Dirac point halfway between them
            assert_eq!(dos.continuous_at(eps), f64::INFINITY);
            assert_eq!(dos.continuous_at(eps + 2.0 * t), f64::INFINITY);
            assert_relative_eq!(dos.continuous_at(eps + t), 0.0, epsilon = 1e-15);

            for (order, moment_ref) in moments_ref.iter().enumerate() {
                let moment = compute_moment(&dos, order as i32);
                assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
            }
        }
    }

    #[test]
    fn lieb() {
        use std::f64::consts::{PI, SQRT_2};

        let eps = 0.5f64;
        let t = 2.0f64;
        let moments_ref = central_to_moments(
            eps,
            even_central_moments(
                8.0 * t.powi(2) / 3.0,
                40.0 * t.powi(4) / 3.0,
                224.0 * t.powi(6) / 3.0,
            ),
        );

        let dos = models::lieb(eps, t);
        let half_width = 2.0 * SQRT_2 * t;
        assert_eq!(dos.support(), Some((eps - half_width, eps + half_width)));
        assert_relative_eq!(dos.total_weight(), 1.0, max_relative = 1e-14);

        // The flat band is a δ-function of weight 1/3 at the band center
        assert_eq!(dos.discrete().len(), 1);
        assert_eq!(dos.discrete().resonances()[0].eps, eps);
        assert_relative_eq!(
            dos.discrete().resonances()[0].weight,
            1.0 / 3.0,
            max_relative = 1e-14
        );

        // Both van Hove singularities diverge, while the point where the dispersive bands
        // touch the flat band is a zero of their DOS.
        assert_eq!(dos.continuous_at(eps - 2.0 * t), f64::INFINITY);
        assert_eq!(dos.continuous_at(eps + 2.0 * t), f64::INFINITY);
        assert_relative_eq!(dos.continuous_at(eps), 0.0, epsilon = 1e-15);

        // Being those of a two-dimensional band minimum, the band edges are steps of
        // height K(0)|ω-ε|/(3π^2t^2) = √2/(3πt) rather than zeros
        for omega_e in [eps - half_width, eps + half_width] {
            assert_relative_eq!(
                dos.continuous_at(omega_e),
                SQRT_2 / (3.0 * PI * t),
                max_relative = 1e-12
            );
        }
        assert_eq!(dos.continuous_at(eps + 1.1 * half_width), 0.0);

        for (order, moment_ref) in moments_ref.iter().enumerate() {
            let moment = compute_moment(&dos, order as i32);
            assert_relative_eq!(moment, moment_ref, max_relative = 1e-10);
        }
    }

    #[test]
    fn lieb_regular_near_singularity() {
        use crate::ContinuousSF;
        use crate::models::LiebDOS;

        let (eps, t) = (0.5f64, 2.0f64);
        let dos = LiebDOS::new(eps, t);
        let prefactor = 1.0 / (std::f64::consts::PI.powi(2) * t);

        // At the singular points the subtracted logarithm cancels the divergence of
        // K(1-y^2) up to ln(4)
        for omega_p in [eps - 2.0 * t, eps + 2.0 * t] {
            assert_relative_eq!(
                dos.regular(omega_p),
                2.0 * prefactor * std::f64::consts::LN_2,
                max_relative = 1e-14
            );
        }

        // R(ε \pm (1+δ)2t) / prefactor at δ = \pm 2^-e, evaluated at 200 bits of precision.
        // The series expansion takes over at |δ| < 5e-4; by e = 40 it is the only option
        // left, 1 - y^2 having rounded to unity.
        let ref_values = [
            (2, 1.98946838014527, 0.881213016813116),
            (4, 1.61182181901004, 1.17635167283123),
            (7, 1.42986716906016, 1.34321545818198),
            (10, 1.39374613100041, 1.37885427306285),
            (15, 1.38663282082393, 1.38595591927928),
            (20, 1.3863082429017, 1.38628047936183),
            (40, 1.38629436114574, 1.38629436109404),
        ];

        for (e, ref_plus, ref_minus) in ref_values {
            for (sign, ref_value) in [(1.0f64, ref_plus), (-1.0, ref_minus)] {
                let delta = sign * 2f64.powi(-e);
                // R(ω) is even about ε, both singularities being equivalent
                for band in [-1.0f64, 1.0] {
                    let omega = eps + band * (1.0 + delta) * 2.0 * t;
                    assert_relative_eq!(
                        dos.regular(omega),
                        prefactor * ref_value,
                        max_relative = 1e-8
                    );
                }
            }
        }
    }
}
