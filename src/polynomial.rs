//! Polynomial in one variable, held by its coefficients.

use std::ops::{Index, IndexMut};

/// Polynomial of degree $n$, $\sum_{i=0}^n c_i x^i$.
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    // Coefficients $c_i$ from $c_0$ up.
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
    fn coefficients_past_the_end_are_zero() {
        let p = poly(&[1.0, 2.0]);
        assert_eq!(
            (p.coeff(0), p.coeff(1), p.coeff(2), p.coeff(99)),
            (1.0, 2.0, 0.0, 0.0)
        );
        assert_eq!(p.degree(), 1);
    }
}
