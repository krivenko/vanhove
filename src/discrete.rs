//! Discrete density of states

use crate::util;

/// Single discrete resonance $w \delta(\omega - \varepsilon)$.
#[derive(Debug, Clone, Copy)]
pub struct Resonance {
    /// Position of the resonance, $\varepsilon$
    pub eps: f64,
    /// Weight of the resonance, $w$
    pub weight: f64,
}

/// Discrete density of states with a finite number of resonances,
/// $$
///     A(\omega) = \sum_p w_p \delta(\omega - \varepsilon_p).
/// $$
#[derive(Debug, Clone)]
pub struct DiscreteDOS {
    /// List of resonances, (\varepsilon_p, w_p)
    resonances: Vec<Resonance>,
    /// resonances[..clean] is sorted w.r.t. \varepsilon_p, deduplicated, zero-weight-free
    clean: usize,
    /// Relative tolerance level for negligible weights
    tol: f64,
}

impl Default for DiscreteDOS {
    fn default() -> Self {
        DiscreteDOS {
            resonances: vec![],
            clean: 0,
            tol: f64::EPSILON,
        }
    }
}

impl DiscreteDOS {
    /// Make an empty discrete DOS (zero resonances)
    pub fn new() -> Self {
        Self::default()
    }

    /// Make a discrete DOS with one resonance
    pub fn one_resonance(eps: f64, weight: f64) -> DiscreteDOS {
        let mut dos = DiscreteDOS::new();
        dos.add(&Resonance { eps, weight });
        dos
    }

    /// Rebuild and access the resonance list in its canonical form (sorted w.r.t.
    /// $\varepsilon_p$, deduplicated, free of negligible weights).
    pub fn resonances(&mut self) -> &[Resonance] {
        if self.clean != self.resonances.len() {
            self.prune();
        }
        &self.resonances
    }

    /// Iterate over the resonances in their current, possibly non-canonical order.
    pub fn iter(&self) -> impl Iterator<Item = &Resonance> {
        self.resonances.iter()
    }

    /// Consume the DOS and return its resonance list in its current,
    /// possibly non-canonical order.
    pub fn into_resonances(self) -> Vec<Resonance> {
        self.resonances
    }

    /// Total spectral weight of the DOS
    pub fn norm(&self) -> f64 {
        util::kahan_babushka_neumaier_sum(self.iter().map(|r| r.weight))
    }

    /// Add a resonance to the DOS
    pub fn add(&mut self, res: &Resonance) {
        debug_assert!(!res.eps.is_nan(), "res.eps must not be NaN");
        self.push(res);
        self.prune_if_needed();
    }

    /// Add a batch of resonances to the DOS, canonicalizing the list at most once
    pub fn extend<I: IntoIterator<Item = Resonance>>(&mut self, resonances: I) {
        let iter = resonances.into_iter();
        self.resonances.reserve(iter.size_hint().0);
        for res in iter {
            debug_assert!(!res.eps.is_nan(), "res.eps must not be NaN");
            self.push(&res);
        }
        self.prune_if_needed();
    }

    /// Multiply all spectral weights by a real number
    pub fn scale(&mut self, a: f64) {
        if a == 0.0 {
            self.resonances.clear();
            self.clean = 0;
            return;
        }

        // Scaling preserves both the ordering w.r.t. ε and the absence of duplicates,
        // so the canonical prefix survives. A weight can, however, underflow to zero
        // and has to be dropped to keep the prefix zero-weight-free.
        let (mut w, mut clean) = (0, self.clean);
        for i in 0..self.resonances.len() {
            let weight = self.resonances[i].weight * a;
            if weight == 0.0 {
                clean -= usize::from(i < self.clean);
            } else {
                self.resonances[w] = Resonance {
                    eps: self.resonances[i].eps,
                    weight,
                };
                w += 1;
            }
        }
        self.resonances.truncate(w);
        self.clean = clean;
    }

    /// Append a resonance to the (possibly unsorted) resonance list
    fn push(&mut self, res: &Resonance) {
        if res.weight != 0.0 {
            // canonicalize -0.0 so it groups with 0.0 under total_cmp
            self.resonances.push(Resonance {
                eps: if res.eps == 0.0 { 0.0 } else { res.eps },
                weight: res.weight,
            });
        }
    }

    /// Canonicalize the resonance list once its unsorted tail has grown big enough.
    ///
    /// The threshold makes `self.clean` grow geometrically between calls, so a
    /// sequence of N additions triggers O(log N) prunes and costs O(N log N) in
    /// total, i.e. O(log N) amortized per addition.
    fn prune_if_needed(&mut self) {
        if self.resonances.len() > 2 * self.clean + 16 {
            self.prune();
        }
    }

    /// Sort, deduplicate and remove small-weight resonances
    fn prune(&mut self) {
        // The singleton fast path below relies on the tolerance being relative
        debug_assert!(self.tol < 1.0, "tol must be a relative tolerance");

        // Sort the resonances w.r.t. ε. `resonances[..clean]` is already sorted,
        // which the stable sort detects as a run and merges rather than re-sorts.
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

            // A lone resonance needs no summation and can only be negligible
            // w.r.t. itself, which `total.abs() > self.tol * mag` never is
            if i == start + 1 {
                res[w] = res[start];
                w += 1;
                continue;
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

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

    #[test]
    fn discrete_dos_extend() {
        // A batch big enough to make `add()` defer the pruning work: the
        // canonical form must be the same either way
        let batch: Vec<Resonance> = (0..100)
            .map(|n| Resonance {
                eps: ((n * 37) % 20) as f64 - 10.0,
                weight: 0.01,
            })
            .collect();

        let mut one_by_one = DiscreteDOS::new();
        for res in &batch {
            one_by_one.add(res);
        }
        let mut bulk = DiscreteDOS::new();
        bulk.extend(batch.iter().copied());

        // 20 distinct levels, each hit 5 times
        assert_eq!(bulk.resonances().len(), 20);
        assert_relative_eq!(bulk.norm(), 1.0, epsilon = 1e-15);
        for (r1, r2) in bulk.resonances().iter().zip(one_by_one.resonances()) {
            assert_eq!(r1.eps, r2.eps);
            assert_eq!(r1.weight, r2.weight);
        }
        // Sorted w.r.t. ε and free of duplicates
        assert!(bulk.resonances().windows(2).all(|w| w[0].eps < w[1].eps));

        // Exactly cancelling contributions are dropped, keeping the survivors
        bulk.extend(batch.iter().map(|r| Resonance {
            eps: r.eps,
            weight: if r.eps < 0.0 { -r.weight } else { 0.0 },
        }));
        assert_eq!(bulk.resonances().len(), 10);
        assert!(bulk.resonances().iter().all(|r| r.eps >= 0.0));
        assert_relative_eq!(bulk.norm(), 0.5, epsilon = 1e-15);
    }

    #[test]
    fn discrete_dos_scale() {
        let mut dos = DiscreteDOS::one_resonance(2.0, 1.0);
        dos.add(&Resonance {
            eps: -1.5,
            weight: 0.25,
        });

        dos.scale(4.0);
        assert_eq!(dos.resonances().len(), 2);
        assert_eq!(dos.resonances()[0].weight, 1.0);
        assert_eq!(dos.resonances()[1].weight, 4.0);
        assert_eq!(dos.norm(), 5.0);

        // Weights underflowing to zero are dropped
        dos.scale(f64::MIN_POSITIVE);
        dos.scale(f64::MIN_POSITIVE);
        assert_eq!(dos.resonances().len(), 0);
        assert_eq!(dos.norm(), 0.0);

        // Scaling by zero empties the DOS
        let mut dos = DiscreteDOS::one_resonance(2.0, 1.0);
        dos.scale(0.0);
        assert_eq!(dos.resonances().len(), 0);
        assert_eq!(dos.norm(), 0.0);
    }
}
