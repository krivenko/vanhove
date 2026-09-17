//! Chebyshev interpolation of the regular part of a spectral function.

use crate::ContinuousSF;
use crate::segment::Segment;
use crate::singularity::Singularity;
use crate::util::{chebyshev_coeffs, clenshaw_chebyshev};

/// Chebyshev expansion of a function over one interval.
#[derive(Debug, Clone)]
struct Panel {
    /// Midpoint of the interval, mapping $\omega$ to $x \in [-1, 1]$.
    mid: f64,
    /// Half-width of the interval.
    half_width: f64,
    /// Coefficients $c_k$ of $\sum_k c_k T_k(x)$, the $k = 0$ one already halved.
    coeffs: Box<[f64]>,
    /// Magnitude of the discarded tail relative to the largest coefficient.
    tail: f64,
}

impl Panel {
    /// Number of coefficients the fit starts with.
    const MIN_ORDER: usize = 8;
    /// Number of coefficients the fit refuses to grow beyond. An expansion still short
    /// of `tol` at this point is not about to reach it: the regular part is not smooth
    /// enough at an end of the panel for the coefficients to decay geometrically.
    const MAX_ORDER: usize = 512;

    /// Fit `f` over `segment`, doubling the order until the tail falls below `tol`
    /// relative to the largest coefficient.
    ///
    /// The nodes lie strictly inside the panel, so `f` is never sampled at an end of
    /// it, where a singular point may sit and the regular part be a difference of two
    /// infinities.
    fn fit<F: Fn(f64) -> f64>(segment: Segment, f: &F, tol: f64) -> Panel {
        let (mid, half_width) = (segment.midpoint(), 0.5 * segment.length());
        let mut n = Self::MIN_ORDER;
        loop {
            let coeffs = chebyshev_coeffs(n, |x| f(mid + half_width * x));
            let scale = coeffs.iter().fold(0.0f64, |m, c| m.max(c.abs()));
            // The last two coefficients stand in for everything left out
            let tail = coeffs[n - 2].abs().max(coeffs[n - 1].abs());
            let tail = if scale > 0.0 { tail / scale } else { 0.0 };

            if tail <= tol || n >= Self::MAX_ORDER {
                // Trailing coefficients below the tolerance contribute nothing but work
                let cutoff = tol * scale;
                let kept = coeffs
                    .iter()
                    .rposition(|c| c.abs() > cutoff)
                    .map_or(1, |k| k + 1);
                return Panel {
                    mid,
                    half_width,
                    coeffs: coeffs[..kept].into(),
                    tail,
                };
            }
            n *= 2;
        }
    }

    /// Value of the expansion at `omega`.
    fn eval(&self, omega: f64) -> f64 {
        clenshaw_chebyshev(&self.coeffs, (omega - self.mid) / self.half_width)
    }
}

/// Continuous spectral function whose regular part is stored as a Chebyshev expansion.
///
/// The support is split into panels at the interior singular points, each panel carrying
/// its own expansion: $R(\omega)$ is smooth within a panel but need not be so across
/// a singular point.
#[derive(Debug, Clone)]
pub struct InterpolatedSF {
    support: Segment,
    /// Panels in ascending order of frequency.
    panels: Box<[Panel]>,
    singularities: Box<[Singularity]>,
}

impl InterpolatedSF {
    /// Default relative tolerance of the fit.
    pub(crate) const DEFAULT_TOL: f64 = 1e-12;

    /// Interpolate the regular part of `csf`.
    ///
    /// `tol` is a tolerance on the Chebyshev coefficients relative to the largest of
    /// them, and defaults to $10^{-12}$. Consult [`InterpolatedSF::fit_error()`] for what
    /// the fit actually achieved: a regular part that is not smooth at an end of a
    /// panel converges too slowly to reach any tolerance worth asking for.
    pub fn new(csf: &dyn ContinuousSF, tol: Option<f64>) -> InterpolatedSF {
        InterpolatedSF::from_parts(
            csf.support(),
            csf.singularities().into(),
            |omega| csf.regular(omega),
            tol,
        )
    }

    /// Interpolate a regular part handed over on its own, along with the support and
    /// the singular structure it belongs to.
    ///
    /// `regular` is sampled strictly between consecutive singular points, never at one,
    /// and carries whatever `singularities` describes already subtracted. That matters
    /// where both diverge: their difference has a limit but no value at the point.
    pub fn from_parts<F: Fn(f64) -> f64>(
        support: Segment,
        singularities: Vec<Singularity>,
        regular: F,
        tol: Option<f64>,
    ) -> InterpolatedSF {
        assert!(
            support.is_bounded(),
            "an interpolated spectral function must have a bounded support"
        );
        assert!(
            !support.is_degenerate(),
            "support must be a non-empty segment"
        );
        let tol = tol.unwrap_or(Self::DEFAULT_TOL);
        assert!(tol > 0.0, "fit tolerance must be positive");

        // Panel boundaries: the ends of the support, and the singular points between
        // them, which is every frequency where R(ω) stops being smooth
        let mut breaks: Vec<f64> = singularities
            .iter()
            .map(|s| s.position())
            .filter(|&p| support.strictly_contains(p))
            .collect();
        breaks.push(support.min());
        breaks.push(support.max());
        breaks.sort_unstable_by(f64::total_cmp);
        breaks.dedup();

        InterpolatedSF {
            support,
            panels: breaks
                .windows(2)
                .map(|lr| Panel::fit(Segment::new(lr[0], lr[1]), &regular, tol))
                .collect(),
            singularities: singularities.into_boxed_slice(),
        }
    }

    /// Largest tail left over by the fit, relative to the largest coefficient of its
    /// panel.
    ///
    /// A value above the requested tolerance means the expansion was cut off before
    /// converging, and is an estimate of the relative error of $R(\omega)$.
    pub fn fit_error(&self) -> f64 {
        self.panels.iter().fold(0.0f64, |m, p| m.max(p.tail))
    }
}

impl ContinuousSF for InterpolatedSF {
    fn support(&self) -> Segment {
        self.support
    }
    fn regular(&self, omega: f64) -> f64 {
        // The panels tile the support, so the first one reaching omega owns it; a
        // frequency past the last boundary belongs to the last panel by rounding.
        let panel = self
            .panels
            .iter()
            .find(|p| omega <= p.mid + p.half_width)
            .unwrap_or(&self.panels[self.panels.len() - 1]);
        panel.eval(omega)
    }
    fn singularities(&self) -> &[Singularity] {
        &self.singularities
    }
    fn shifted(&self, by: f64) -> Box<dyn ContinuousSF> {
        Box::new(InterpolatedSF {
            support: self.support.shifted(by),
            panels: self
                .panels
                .iter()
                .map(|p| Panel {
                    mid: p.mid + by,
                    ..p.clone()
                })
                .collect(),
            singularities: self.singularities.iter().map(|s| s.shifted(by)).collect(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::{InterpolatedSF, Panel};
    use crate::ContinuousSF;
    use crate::segment::Segment;
    use crate::singularity::{AsymptTerm, Singularity};
    use approx::assert_abs_diff_eq;

    /// Stand-in spectral function built straight out of its parts.
    struct Model<F> {
        support: Segment,
        singularities: Vec<Singularity>,
        regular: F,
    }

    impl<F: Fn(f64) -> f64 + Send + Sync> ContinuousSF for Model<F> {
        fn support(&self) -> Segment {
            self.support
        }
        fn regular(&self, omega: f64) -> f64 {
            (self.regular)(omega)
        }
        fn singularities(&self) -> &[Singularity] {
            &self.singularities
        }
        fn shifted(&self, _by: f64) -> Box<dyn ContinuousSF> {
            unimplemented!("the test double is never displaced")
        }
    }

    fn model<F>(min: f64, max: f64, regular: F) -> Model<F> {
        Model {
            support: Segment::new(min, max),
            singularities: vec![],
            regular,
        }
    }

    fn with_sing<F>(min: f64, max: f64, singularities: Vec<Singularity>, regular: F) -> Model<F> {
        Model {
            support: Segment::new(min, max),
            singularities,
            regular,
        }
    }

    /// Largest departure of the interpolation from `f` over the support.
    fn worst_error<F: Fn(f64) -> f64>(interp: &InterpolatedSF, f: F) -> f64 {
        let support = interp.support();
        (1..500).fold(0.0f64, |w, i| {
            let omega = support.min() + support.length() * (i as f64) / 500.0;
            w.max((interp.regular(omega) - f(omega)).abs())
        })
    }

    #[test]
    fn smooth_fit() {
        // A function analytic on the panel is caught to machine precision
        let f = |omega: f64| (-(omega - 0.5).powi(2)).exp();
        let interp = InterpolatedSF::new(&model(-2.0, 3.0, f), None);
        assert_eq!(interp.support(), Segment::new(-2.0, 3.0));
        assert!(interp.singularities().is_empty());
        assert!(interp.fit_error() <= 1e-12);
        assert!(worst_error(&interp, f) < 1e-13);

        // The end points of the support are reproduced along with the interior. The
        // tolerance is on the coefficients, hence on the absolute error: f is three
        // orders of magnitude below its peak out here.
        for omega in [-2.0, 3.0] {
            assert_abs_diff_eq!(interp.regular(omega), f(omega), epsilon = 1e-13);
        }
    }

    #[test]
    fn panels_split_at_singular_points() {
        let log = |position| Singularity::new(position, 1.0, vec![AsymptTerm::log(1.0)]);
        // An interior singular point splits the support in two
        let interp = InterpolatedSF::new(&with_sing(-1.0, 2.0, vec![log(0.0)], f64::abs), None);
        assert_eq!(interp.panels.len(), 2);

        // |ω| is a kink each panel resolves exactly, being linear on either side
        assert!(worst_error(&interp, f64::abs) < 1e-14);

        // A singular point sitting at an end of the support adds no panel
        let edge = Singularity::new(-1.0, 1.0, vec![AsymptTerm::power(0.5, 1.0)]);
        let interp = InterpolatedSF::new(&with_sing(-1.0, 2.0, vec![edge], |o: f64| o), None);
        assert_eq!(interp.panels.len(), 1);

        // Two singularities at one point are one break, not two
        let sings = vec![log(0.5), log(0.5)];
        let interp = InterpolatedSF::new(&with_sing(-1.0, 2.0, sings, f64::abs), None);
        assert_eq!(interp.panels.len(), 2);
    }

    #[test]
    fn singularities_carried_over() {
        // The singular part passes through untouched, only R(ω) being approximated
        let sing = Singularity::new(0.5, 2.0, vec![AsymptTerm::power(-0.5, 3.0)]);
        let m = with_sing(-1.0, 2.0, vec![sing.clone()], |omega: f64| omega);
        let interp = InterpolatedSF::new(&m, None);
        assert_eq!(interp.singularities().len(), 1);
        for omega in [-1.0, 0.0, 1.0, 2.0] {
            assert_eq!(interp.singularities()[0].value(omega), sing.value(omega));
        }
    }

    #[test]
    fn singular_point_rescues_a_kink() {
        // Two entire functions meeting at a kink. Left inside one panel the kink caps
        // the fit; reported as a singular point, each side is analytic and converges.
        let f = |omega: f64| {
            if omega < 0.0 {
                omega.exp()
            } else {
                omega.cos()
            }
        };

        let whole = InterpolatedSF::new(&model(-2.0, 2.0, f), None);
        assert!(whole.fit_error() > 1e-13);
        assert!(whole.panels[0].coeffs.len() >= Panel::MAX_ORDER - 1);

        // A singularity with no terms subtracts nothing and only marks the point
        let kink = Singularity::new(0.0, 1.0, vec![]);
        let split = InterpolatedSF::new(&with_sing(-2.0, 2.0, vec![kink], f), None);
        assert_eq!(split.panels.len(), 2);
        assert!(split.fit_error() <= 1e-12);
        assert!(worst_error(&split, f) < 1e-11);

        // Six orders of magnitude of accuracy, off a fortieth of the coefficients
        assert!(worst_error(&whole, f) > 1e6 * worst_error(&split, f));
        assert!(split.panels.iter().all(|p| p.coeffs.len() < 32));
    }

    #[test]
    #[should_panic(expected = "must have a bounded support")]
    fn unbounded_support() {
        let _ = InterpolatedSF::new(&model(f64::NEG_INFINITY, 1.0, |omega: f64| omega), None);
    }
}
