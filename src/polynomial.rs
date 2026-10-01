//! Polynomial in one variable, held by its coefficients.

use std::ops::{AddAssign, Index, IndexMut, MulAssign};

use crate::util::binomials;

/// Polynomial of degree $n$, $\sum_{i=0}^n c_i x^i$.
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    // Coefficients $c_i$ from $c_0$ up. Never empty.
    coeffs: Vec<f64>,
}

impl Polynomial {
    /// A polynomial of degree `n` with all coefficients set to zero.
    pub fn zeros(n: usize) -> Polynomial {
        Polynomial {
            coeffs: vec![0.0; n + 1],
        }
    }

    /// Degree $n$, which is one less than the number of coefficients held.
    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    /// Coefficient of $x^i$, zero for $i > n$.
    pub fn coeff(&self, i: usize) -> f64 {
        self.coeffs.get(i).copied().unwrap_or(0.0)
    }

    /// Become a polynomial of degree `n`, dropping any coefficient past it and filling
    /// with zeros to reach it.
    pub fn resize(&mut self, n: usize) {
        self.coeffs.resize(n + 1, 0.0);
    }

    /// Add `a` times `other`, reaching this one's degree up to the greater of the two
    /// where it has to.
    pub fn add_scaled(&mut self, other: &Polynomial, a: f64) {
        if other.degree() > self.degree() {
            self.resize(other.degree());
        }
        for (entry, add) in self.coeffs.iter_mut().zip(&other.coeffs) {
            *entry += a * add;
        }
    }

    /// Whether every coefficient is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|c| *c == 0.0)
    }

    /// Drop trailing zero coefficients, so that the degree is the highest power that
    /// weighs anything. A polynomial of nothing but zeros keeps its constant term.
    pub fn trim(&mut self) {
        while self.coeffs.len() > 1 && self.coeffs.last() == Some(&0.0) {
            self.coeffs.pop();
        }
    }

    /// The same polynomial read about a displaced origin, $Q(x) = P(x - b)$.
    pub fn shifted(&self, b: f64) -> Polynomial {
        // $P(x-b) = \sum_k x^k \sum_{i \ge k} c_i \binom{i}{k} (-b)^{i-k}$, a triangular
        // map onto a polynomial of the same degree.
        let mut out = Polynomial::zeros(self.degree());
        for (i, &c) in self.coeffs.iter().enumerate() {
            let c_row = binomials(i);
            let mut power = 1.0;
            for j in 0..=i {
                out[i - j] += c * c_row[j] * power;
                power *= -b;
            }
        }
        out
    }

    /// Value at `x`.
    pub fn eval(&self, x: f64) -> f64 {
        // Horner algorithm from the highest power with a non-zero coefficient.
        let Some(top) = self.coeffs.iter().rposition(|c| *c != 0.0) else {
            return 0.0;
        };
        let mut acc = self.coeffs[top];
        for i in (0..top).rev() {
            acc = acc * x + self.coeffs[i];
        }
        acc
    }
}

/// Add a polynomial of any degree.
impl AddAssign<&Polynomial> for Polynomial {
    fn add_assign(&mut self, other: &Polynomial) {
        self.add_scaled(other, 1.0);
    }
}

impl MulAssign<f64> for Polynomial {
    fn mul_assign(&mut self, a: f64) {
        for c in &mut self.coeffs {
            *c *= a;
        }
    }
}

/// Panics past the last coefficient, as for a slice.
impl Index<usize> for Polynomial {
    type Output = f64;
    fn index(&self, i: usize) -> &f64 {
        &self.coeffs[i]
    }
}

/// Panics past the last coefficient, as for a slice.
impl IndexMut<usize> for Polynomial {
    fn index_mut(&mut self, i: usize) -> &mut f64 {
        &mut self.coeffs[i]
    }
}

#[cfg(test)]
mod tests {
    use super::Polynomial;
    use approx::assert_relative_eq;

    fn poly(coeffs: &[f64]) -> Polynomial {
        let mut p = Polynomial::zeros(coeffs.len() - 1);
        for (i, &c) in coeffs.iter().enumerate() {
            p[i] = c;
        }
        p
    }

    #[test]
    fn eval() {
        // 2 - 3x + x^3, read at a few points against the written-out form
        let p = poly(&[2.0, -3.0, 0.0, 1.0]);
        for x in [-2.5f64, -1.0, 0.0, 0.5, 4.0] {
            assert_relative_eq!(p.eval(x), 2.0 - 3.0 * x + x.powi(3), max_relative = 1e-14);
        }
        assert_eq!(Polynomial::zeros(3).eval(3.0), 0.0);
        assert_eq!(Polynomial::zeros(0).eval(3.0), 0.0);
    }

    #[test]
    fn eval_at_an_infinite_argument() {
        // Horner started from a trailing zero coefficient would turn an infinite
        // argument into a NaN, so `eval` starts below one.
        let constant = poly(&[2.0, 0.0, 0.0]);
        assert_eq!(constant.eval(f64::NEG_INFINITY), 2.0);
        let linear = poly(&[7.0, -1.0]);
        assert_eq!(linear.eval(f64::NEG_INFINITY), f64::INFINITY);
        assert_eq!(Polynomial::zeros(2).eval(f64::NEG_INFINITY), 0.0);
    }

    #[test]
    fn shifted() {
        // $P(x-b)$ against the written-out form, read at a few points
        let p = poly(&[2.0, -3.0, 0.0, 1.0]);
        for b in [-1.5f64, 0.0, 0.25, 4.0] {
            let q = p.shifted(b);
            assert_eq!(q.degree(), p.degree());
            for x in [-2.0f64, 0.0, 1.3] {
                assert_relative_eq!(q.eval(x), p.eval(x - b), max_relative = 1e-13);
            }
        }
        // Shifting back returns the polynomial it came from
        let there_and_back = p.shifted(3.0).shifted(-3.0);
        for i in 0..=p.degree() {
            assert_relative_eq!(there_and_back[i], p[i], max_relative = 1e-13);
        }
    }

    #[test]
    fn add_and_scale() {
        let mut p = poly(&[1.0, 2.0]);
        p += &poly(&[10.0, 20.0, 30.0]);
        assert_eq!((p.degree(), p[0], p[1], p[2]), (2, 11.0, 22.0, 30.0));
        p *= 0.5;
        assert_eq!((p[0], p[1], p[2]), (5.5, 11.0, 15.0));
        // Adding a shorter one leaves the degree alone
        p += &poly(&[1.0]);
        assert_eq!((p.degree(), p[0]), (2, 6.5));
    }

    #[test]
    fn add_scaled() {
        let mut p = poly(&[1.0, 2.0]);
        p.add_scaled(&poly(&[10.0, 20.0, 30.0]), -0.5);
        assert_eq!((p.degree(), p[0], p[1], p[2]), (2, -4.0, -8.0, -15.0));
        // A weight of zero leaves the coefficients alone but still reaches the degree
        let mut q = poly(&[1.0]);
        q.add_scaled(&poly(&[0.0, 0.0, 7.0]), 0.0);
        assert_eq!((q.degree(), q[0], q[2]), (2, 1.0, 0.0));
    }

    #[test]
    fn trim() {
        let mut p = poly(&[1.0, 2.0, 0.0, 0.0]);
        p.trim();
        assert_eq!(p.degree(), 1);
        // Cancellation can empty the top, and nothing but the constant term survives
        let mut zero = poly(&[0.0, 0.0, 0.0]);
        zero.trim();
        assert_eq!((zero.degree(), zero[0]), (0, 0.0));
    }

    #[test]
    fn coefficients_past_the_end_are_zero() {
        let p = poly(&[1.0, 2.0]);
        assert_eq!(
            (p.coeff(0), p.coeff(1), p.coeff(2), p.coeff(99)),
            (1.0, 2.0, 0.0, 0.0)
        );
        assert_eq!(p.degree(), 1);
    }
}
