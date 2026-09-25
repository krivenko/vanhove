//! Truncated Laurent series in a small parameter.
//!
//! Arithmetic for series around a pole of finite order, used where two divergences are
//! known to cancel and the finite part is what is wanted. Everything is truncated to a
//! fixed window of powers, so a product loses the powers past it.

use std::ops::{Add, AddAssign, Mul, Neg, Sub};

/// Truncated Laurent series in $\epsilon$, holding the coefficients of $\epsilon^{-w}$
/// through $\epsilon^{w}$.
#[derive(Clone)]
pub struct Laurent {
    /// The powers run from $-w$ to $w$, so the window holds $2w+1$ of them.
    w: usize,
    /// Coefficients from $\epsilon^{-w}$ up, so $\epsilon^k$ sits at `c[k + w]`.
    c: Vec<f64>,
}

impl Laurent {
    pub fn zero(w: usize) -> Laurent {
        Laurent {
            w,
            c: vec![0.0; 2 * w + 1],
        }
    }

    /// A series with no pole, from the coefficients of $\epsilon^0, \epsilon^1, \ldots$.
    ///
    /// Takes the $w+1$ of them from `taylor`, which must contain at least that many.
    pub fn from_taylor(w: usize, taylor: &[f64]) -> Laurent {
        let mut s = Laurent::zero(w);
        s.c[w..].copy_from_slice(&taylor[..w + 1]);
        s
    }

    /// Series for $1/(n + \epsilon)$.
    ///
    /// At $n = 0$ it is the bare pole $1/\epsilon$, and multiplying by that shifts the
    /// window down: the top coefficient of the product would have to come from
    /// $\epsilon^{w+1}$, never in range, so it stays zero. Apply a divisor last, once
    /// the analytic factors have been taken at the full width.
    pub fn reciprocal_linear(w: usize, n: isize) -> Laurent {
        if n == 0 {
            let mut out = Laurent::zero(w);
            out.c[w - 1] = 1.0;
            return out;
        }
        // 1/(n+ε) = (1/n) \sum (-ε/n)^i
        let s = n as f64;
        let mut taylor = Vec::with_capacity(w + 1);
        let mut term = 1.0 / s;
        for _ in 0..=w {
            taylor.push(term);
            term /= -s;
        }
        Laurent::from_taylor(w, &taylor)
    }

    /// Where $\epsilon^k$ sits in `c`, absent if the window does not reach that far.
    fn slot(&self, k: isize) -> Option<usize> {
        let i = k + self.w as isize;
        (i >= 0 && (i as usize) < self.c.len()).then_some(i as usize)
    }

    /// Coefficient of $\epsilon^k$, zero outside the window.
    pub fn at(&self, k: isize) -> f64 {
        self.slot(k).map_or(0.0, |i| self.c[i])
    }

    /// Largest coefficient in absolute value.
    pub fn peak(&self) -> f64 {
        self.c.iter().fold(0.0f64, |x, y| x.max(y.abs()))
    }

    /// The same series with $\epsilon$ replaced by $-\epsilon$.
    pub fn reflected_in_epsilon(&self) -> Laurent {
        let mut out = self.clone();
        for (i, value) in out.c.iter_mut().enumerate() {
            let power = i as isize - self.w as isize;
            if power.rem_euclid(2) == 1 {
                *value = -*value;
            }
        }
        out
    }

    /// $d/d\epsilon$, term by term.
    pub fn diff(&self) -> Laurent {
        let mut r = Laurent::zero(self.w);
        for (i, &v) in self.c.iter().enumerate() {
            let k = i as isize - self.w as isize;
            if k == 0 || v == 0.0 {
                continue;
            }
            if let Some(slot) = r.slot(k - 1) {
                r.c[slot] += v * k as f64;
            }
        }
        r
    }
}

/// Arithmetic on two series assumes they share a window.
macro_rules! same_window {
    ($a:expr, $b:expr) => {
        debug_assert_eq!($a.w, $b.w, "Laurent series of unequal windows");
    };
}

impl Neg for Laurent {
    type Output = Laurent;
    fn neg(mut self) -> Laurent {
        self.c.iter_mut().for_each(|v| *v = -*v);
        self
    }
}

impl Neg for &Laurent {
    type Output = Laurent;
    fn neg(self) -> Laurent {
        self.clone().neg()
    }
}

impl AddAssign<&Laurent> for Laurent {
    fn add_assign(&mut self, other: &Laurent) {
        same_window!(self, other);
        for (a, b) in self.c.iter_mut().zip(&other.c) {
            *a += b;
        }
    }
}

impl AddAssign<Laurent> for Laurent {
    fn add_assign(&mut self, other: Laurent) {
        *self += &other;
    }
}

impl Add<&Laurent> for &Laurent {
    type Output = Laurent;
    fn add(self, other: &Laurent) -> Laurent {
        let mut out = self.clone();
        out += other;
        out
    }
}

impl Sub<&Laurent> for &Laurent {
    type Output = Laurent;
    fn sub(self, other: &Laurent) -> Laurent {
        same_window!(self, other);
        let mut out = self.clone();
        for (a, b) in out.c.iter_mut().zip(&other.c) {
            *a -= b;
        }
        out
    }
}

/// The product, truncated back to the power window the two operands share.
impl Mul<&Laurent> for &Laurent {
    type Output = Laurent;
    fn mul(self, other: &Laurent) -> Laurent {
        same_window!(self, other);
        let mut r = Laurent::zero(self.w);
        for (i, &a) in self.c.iter().enumerate() {
            if a == 0.0 {
                continue;
            }
            for (j, &b) in other.c.iter().enumerate() {
                let k = (i + j) as isize - self.w as isize - other.w as isize;
                if let Some(slot) = r.slot(k) {
                    r.c[slot] += a * b;
                }
            }
        }
        r
    }
}

impl Mul<f64> for Laurent {
    type Output = Laurent;
    fn mul(mut self, f: f64) -> Laurent {
        self.c.iter_mut().for_each(|v| *v *= f);
        self
    }
}

impl Mul<f64> for &Laurent {
    type Output = Laurent;
    fn mul(self, f: f64) -> Laurent {
        self.clone() * f
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic() {
        let w = 2;
        // ε^-1 + 2 + 3ε
        let mut a = Laurent::zero(w);
        (a.c[w - 1], a.c[w], a.c[w + 1]) = (1.0, 2.0, 3.0);
        // 1 + ε
        let b = Laurent::from_taylor(w, &[1.0, 1.0, 0.0]);

        // Negation turns every power
        assert_eq!((-&a).at(-1), -1.0);
        assert_eq!((-&a).at(1), -3.0);
        // Subtraction is a composition of negation and addition
        for k in -2..=2 {
            assert_eq!((&a - &b).at(k), (&a + &(-&b)).at(k));
        }

        // (ε^-1 + 2 + 3ε)(1 + ε) = ε^-1 + 3 + 5ε + 3ε^2
        let p = &a * &b;
        assert_eq!([p.at(-1), p.at(0), p.at(1), p.at(2)], [1.0, 3.0, 5.0, 3.0]);

        // A product reaching past the window is dropped
        let (c, d) = (
            Laurent::from_taylor(w, &[0.0, 0.0, 1.0]),
            Laurent::from_taylor(w, &[0.0, 1.0, 0.0]),
        );
        assert!((-2..=2).all(|k| (&c * &d).at(k) == 0.0));

        // Accumulating a series and its negative leaves zero
        let mut acc = Laurent::zero(w);
        acc += &a;
        acc += &a * -1.0;
        assert!((-2..=2).all(|k| acc.at(k) == 0.0));
    }

    #[test]
    fn peak_is_the_largest_coefficient() {
        let w = 2;
        let mut a = Laurent::zero(w);
        (a.c[w - 1], a.c[w], a.c[w + 1]) = (1.0, -7.5, 3.0);
        assert_eq!(a.peak(), 7.5);
        assert_eq!(Laurent::zero(w).peak(), 0.0);
    }

    /// A pole applied last leaves the top coefficient at zero rather than filling it
    /// from data the window never held.
    #[test]
    fn a_pole_cannot_fill_the_top_slot() {
        let w = 2;
        let analytic = Laurent::from_taylor(w, &[1.0, 1.0, 1.0]);
        let shifted = &analytic * &Laurent::reciprocal_linear(w, 0);
        assert_eq!(shifted.at(-1), 1.0);
        assert_eq!(shifted.at(0), 1.0);
        assert_eq!(shifted.at(1), 1.0);
        assert_eq!(shifted.at(w as isize), 0.0);
    }
}
