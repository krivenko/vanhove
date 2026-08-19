//! Discrete spectral function

use std::ops::{Add, Mul};

use crate::util;

/// Single discrete resonance $w \delta(\omega - \varepsilon)$.
#[derive(Debug, Clone, Copy)]
pub struct Resonance {
    /// Position of the resonance, $\varepsilon$
    pub eps: f64,
    /// Weight of the resonance, $w$
    pub weight: f64,
}

/// Discrete spectral function with a finite number of resonances,
/// $$
///     A(\omega) = \sum_p w_p \delta(\omega - \varepsilon_p).
/// $$
///
/// The list of resonances is sealed upon construction and is always kept in its
/// canonical form: sorted w.r.t. $\varepsilon_p$, free of duplicate positions and
/// free of negligible weights.
#[derive(Debug, Clone, Default)]
pub struct DiscreteSF {
    /// Resonances $(\varepsilon_p, w_p)$ in the canonical form
    resonances: Vec<Resonance>,
}

/// Build a discrete spectral function out of resonances given in an arbitrary
/// order. Resonances sharing a position are merged, and groups whose total weight
/// cancels out are dropped.
impl FromIterator<Resonance> for DiscreteSF {
    fn from_iter<I: IntoIterator<Item = Resonance>>(resonances: I) -> Self {
        let mut resonances: Vec<Resonance> = resonances
            .into_iter()
            .inspect(|res| debug_assert!(!res.eps.is_nan(), "res.eps must not be NaN"))
            .filter(|res| res.weight != 0.0)
            .map(|res| Resonance {
                // canonicalize -0.0 so that it groups with 0.0 under total_cmp
                eps: if res.eps == 0.0 { 0.0 } else { res.eps },
                weight: res.weight,
            })
            .collect();
        resonances.sort_by(|r1, r2| r1.eps.total_cmp(&r2.eps));
        canonicalize(&mut resonances);
        DiscreteSF { resonances }
    }
}

/// Addition of two discrete spectral functions.
impl Add for DiscreteSF {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let (lhs, rhs) = (self.resonances, rhs.resonances);
        let mut resonances = Vec::with_capacity(lhs.len() + rhs.len());

        let (mut l, mut r) = (lhs.iter().peekable(), rhs.iter().peekable());
        loop {
            let ordering = match (l.peek(), r.peek()) {
                (Some(rl), Some(rr)) => rl.eps.total_cmp(&rr.eps),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => break,
            };
            match ordering {
                std::cmp::Ordering::Less => resonances.push(*l.next().unwrap()),
                std::cmp::Ordering::Greater => resonances.push(*r.next().unwrap()),
                std::cmp::Ordering::Equal => {
                    // Both operands are duplicate-free, so exactly two weights meet here
                    // and their sum needs no compensation to be correctly rounded.
                    let (rl, rr) = (l.next().unwrap(), r.next().unwrap());
                    let weight = rl.weight + rr.weight;
                    // Drop the resonance if the contributions cancel each other out
                    let tol = DiscreteSF::WEIGHT_TOL;
                    if weight.abs() > tol * (rl.weight.abs() + rr.weight.abs()) {
                        resonances.push(Resonance {
                            eps: rl.eps,
                            weight,
                        });
                    }
                }
            }
        }
        DiscreteSF { resonances }
    }
}

/// Multiply all spectral weights by a real number from the right.
impl Mul<f64> for DiscreteSF {
    type Output = Self;

    fn mul(mut self, a: f64) -> Self {
        if a == 0.0 {
            self.resonances.clear();
        } else {
            // Scaling preserves both the ordering w.r.t. ε and the absence of duplicates.
            // A weight can, however, underflow to zero and has to be dropped to keep the
            // list free of negligible weights.
            self.resonances.retain_mut(|res| {
                res.weight *= a;
                res.weight != 0.0
            });
        }
        self
    }
}

/// Multiply all spectral weights by a real number from the left.
impl Mul<DiscreteSF> for f64 {
    type Output = DiscreteSF;
    fn mul(self, sf: DiscreteSF) -> DiscreteSF {
        sf * self
    }
}

/// Iterate over the resonances in the order of increasing $\varepsilon_p$.
impl<'a> IntoIterator for &'a DiscreteSF {
    type Item = &'a Resonance;
    type IntoIter = std::slice::Iter<'a, Resonance>;
    fn into_iter(self) -> Self::IntoIter {
        self.resonances.iter()
    }
}

/// Merge resonances sharing the same position within an ε-sorted list, and remove
/// the groups whose total weight is negligible. Called on a freshly sorted list, this
/// establishes the canonical form.
fn canonicalize(resonances: &mut Vec<Resonance>) {
    let (mut w, mut i) = (0, 0);
    while i < resonances.len() {
        // Find a slice of resonances sharing the same ε
        let start = i;
        let eps = resonances[i].eps;
        while i < resonances.len() && resonances[i].eps == eps {
            i += 1;
        }

        // A lone resonance needs no summation and can only be negligible
        // w.r.t. itself, which `total.abs() > tol * mag` never is
        if i == start + 1 {
            resonances[w] = resonances[start];
            w += 1;
            continue;
        }

        // Sum all weights within the slice, with Neumaier compensation
        let group = &mut resonances[start..i];
        group.sort_unstable_by(|r1, r2| r1.weight.abs().total_cmp(&r2.weight.abs()));
        let total = util::kahan_babushka_neumaier_sum(group.iter().map(|r| r.weight));

        // Keep the resonance if its total weight is not too small compared
        // to the total magnitude of the weights
        let mag: f64 = group.iter().map(|r| r.weight.abs()).sum();
        if total.abs() > DiscreteSF::WEIGHT_TOL * mag {
            resonances[w] = Resonance { eps, weight: total };
            w += 1;
        }
    }
    // Remove the unused elements
    resonances.truncate(w);
}

impl DiscreteSF {
    /// Relative tolerance below which the total weight of a group of resonances sharing
    /// the same position is considered a cancellation artefact and is discarded.
    pub const WEIGHT_TOL: f64 = f64::EPSILON;

    /// Make an empty discrete spectral function (zero resonances)
    pub fn new() -> Self {
        Self::default()
    }

    /// Make a discrete spectral function with one resonance
    pub fn one_resonance(eps: f64, weight: f64) -> DiscreteSF {
        DiscreteSF::from_iter([Resonance { eps, weight }])
    }

    /// Access the resonance list in its canonical form (sorted w.r.t. $\varepsilon_p$,
    /// deduplicated, free of negligible weights).
    pub fn resonances(&self) -> &[Resonance] {
        &self.resonances
    }

    /// Number of resonances in the spectral function
    pub fn len(&self) -> usize {
        self.resonances.len()
    }

    /// Does the spectral function carry no resonances?
    pub fn is_empty(&self) -> bool {
        self.resonances.is_empty()
    }

    /// Iterate over the resonances in the order of increasing $\varepsilon_p$.
    pub fn iter(&self) -> std::slice::Iter<'_, Resonance> {
        self.resonances.iter()
    }

    /// Consume the spectral function and return its resonance list in the canonical form.
    pub fn into_resonances(self) -> Vec<Resonance> {
        self.resonances
    }

    /// Support of the spectral function, i.e. the positions of its lowest and
    /// highest resonances. Returns [`None`] for an empty spectral function.
    pub fn support(&self) -> Option<(f64, f64)> {
        Some((self.resonances.first()?.eps, self.resonances.last()?.eps))
    }

    /// Total spectral weight.
    pub fn total_weight(&self) -> f64 {
        util::kahan_babushka_neumaier_sum(self.iter().map(|r| r.weight))
    }

    /// Find the resonance located at a given position $\varepsilon$, if any.
    pub fn find(&self, eps: f64) -> Option<&Resonance> {
        // -0.0 is stored as 0.0, so the same canonicalization applies to the needle
        let eps = if eps == 0.0 { 0.0 } else { eps };
        let index = self
            .resonances
            .binary_search_by(|res| res.eps.total_cmp(&eps))
            .ok()?;
        Some(&self.resonances[index])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn discrete_dos() {
        let sf = DiscreteSF::new();
        assert!(sf.is_empty());
        assert_eq!(sf.resonances().len(), 0);
        assert_eq!(sf.total_weight(), 0.0);
        assert_eq!(sf.support(), None);

        let sf = DiscreteSF::one_resonance(2.0, 1.0);
        assert_eq!(sf.len(), 1);
        assert_eq!(sf.resonances()[0].eps, 2.0);
        assert_eq!(sf.resonances()[0].weight, 1.0);
        assert_eq!(sf.total_weight(), 1.0);
        assert_eq!(sf.support(), Some((2.0, 2.0)));

        // A resonance of zero weight is not a resonance
        assert!(DiscreteSF::one_resonance(2.0, 0.0).is_empty());

        let sf = sf + DiscreteSF::one_resonance(-1.5, 0.25);
        assert_eq!(sf.len(), 2);
        assert_eq!(sf.total_weight(), 1.25);
        assert_eq!(sf.support(), Some((-1.5, 2.0)));

        // Exactly cancelling contributions annihilate the resonance
        let sf = sf + DiscreteSF::one_resonance(2.0, -1.0);
        assert_eq!(sf.len(), 1);
        assert_eq!(sf.total_weight(), 0.25);
        assert_eq!(sf.support(), Some((-1.5, -1.5)));
    }

    #[test]
    fn discrete_dos_from_iter() {
        let batch: Vec<Resonance> = (0..100)
            .map(|n| Resonance {
                eps: ((n * 37) % 20) as f64 - 10.0,
                weight: 0.01,
            })
            .collect();

        // Building in one go and accumulating pairwise give the same canonical form
        let bulk = DiscreteSF::from_iter(batch.iter().copied());
        let one_by_one = batch.iter().fold(DiscreteSF::new(), |sf, res| {
            sf + DiscreteSF::one_resonance(res.eps, res.weight)
        });

        // 20 distinct levels, each hit 5 times
        assert_eq!(bulk.len(), 20);
        assert_relative_eq!(bulk.total_weight(), 1.0, epsilon = 1e-15);
        assert_eq!(bulk.support(), Some((-10.0, 9.0)));
        for (r1, r2) in bulk.iter().zip(&one_by_one) {
            assert_eq!(r1.eps, r2.eps);
            assert_relative_eq!(r1.weight, r2.weight, epsilon = 1e-15);
        }
        // Sorted w.r.t. ε and free of duplicates
        assert!(bulk.resonances().windows(2).all(|w| w[0].eps < w[1].eps));

        // Exactly cancelling contributions are dropped, keeping the survivors
        let sf = bulk
            + DiscreteSF::from_iter(batch.iter().filter(|r| r.eps < 0.0).map(|r| Resonance {
                eps: r.eps,
                weight: -r.weight,
            }));
        assert_eq!(sf.len(), 10);
        assert!(sf.iter().all(|r| r.eps >= 0.0));
        assert_relative_eq!(sf.total_weight(), 0.5, epsilon = 1e-15);
    }

    #[test]
    fn discrete_dos_find() {
        let sf = DiscreteSF::from_iter((0..64).map(|n| Resonance {
            eps: (n - 32) as f64 / 2.0,
            weight: 1.0 / 64.0,
        }));

        for n in 0..64 {
            let eps = (n - 32) as f64 / 2.0;
            assert_eq!(sf.find(eps).unwrap().eps, eps);
            assert_eq!(sf.find(eps).unwrap().weight, 1.0 / 64.0);
        }
        assert!(sf.find(0.25).is_none());
        assert!(sf.find(100.0).is_none());
        // -0.0 and 0.0 denote the same position
        assert_eq!(sf.find(-0.0).unwrap().eps, 0.0);
        assert_eq!(
            DiscreteSF::one_resonance(-0.0, 1.0).find(0.0).unwrap().eps,
            0.0
        );
    }

    #[test]
    fn discrete_dos_mul() {
        let sf = DiscreteSF::one_resonance(2.0, 1.0) + DiscreteSF::one_resonance(-1.5, 0.25);

        let sf = sf * 4.0;
        assert_eq!(sf.len(), 2);
        assert_eq!(sf.resonances()[0].weight, 1.0);
        assert_eq!(sf.resonances()[1].weight, 4.0);
        assert_eq!(sf.total_weight(), 5.0);

        // Weights underflowing to zero are dropped
        let sf = f64::MIN_POSITIVE * (f64::MIN_POSITIVE * sf);
        assert!(sf.is_empty());
        assert_eq!(sf.total_weight(), 0.0);

        // Scaling by zero empties the spectral function
        let sf = DiscreteSF::one_resonance(2.0, 1.0) * 0.0;
        assert!(sf.is_empty());
        assert_eq!(sf.total_weight(), 0.0);
    }
}
