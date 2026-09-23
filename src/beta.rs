//! The beta function and its derivatives with respect to its parameters.
//!
//! A convolution of two singular parts is an incomplete beta function: a power against
//! a power over a finite range is $B_z(a,b)$, and a logarithm on either factor is
//! a derivative with respect to one of the two exponents. The derivatives here are
//! therefore with respect to the *parameters*, not the argument, which is what no library
//! offers.

use std::ops::{AddAssign, Mul};

use special::Gamma;

use crate::laurent::Laurent;
use crate::util::{Table, alternating_sign, binomials, polygamma};

//
// General functions
//

/// $n$-th derivative of a product from those of its factors,
/// $(fg)^{(n)} = \sum_j \binom{n}{j} f^{(n-j)} g^{(j)}$.
fn leibniz(n: usize, f: &[f64], g: &[f64]) -> f64 {
    let c_row = binomials(n);
    (0..=n).map(|j| c_row[j] * f[n - j] * g[j]).sum()
}

//
// Gamma function
//

/// Derivatives of $\Gamma$ at `z`, up to order `n`.
///
/// $\Gamma' = \Gamma\psi$ differentiated by Leibniz, which needs `z` to be no pole.
fn gamma_derivatives(z: f64, n: usize) -> Vec<f64> {
    let psi: Vec<f64> = (0..=n).map(|j| polygamma(j as u32, z)).collect();
    let mut d = vec![Gamma::gamma(z)];
    for k in 0..n {
        let next = leibniz(k, &d, &psi);
        d.push(next);
    }
    d
}

/// Derivatives of $1/\Gamma$ at `z`, up to order `n`.
///
/// $1/\Gamma$ is entire, so these stay finite at the poles of $\Gamma$, which is the
/// whole reason the Beta derivatives are taken through it: $B(a, b)$ is a
/// perfectly good number where $a+b$ is a non-positive integer, and only the
/// way of writing it down falls over there.
/// The recurrence $z\\,\Gamma(z) = \Gamma(z+1)$ walks the argument up to where
/// $1/\Gamma = \exp(-\ln\Gamma)$ may be differentiated directly, leaving a polynomial
/// factor behind.
fn recip_gamma_derivatives(z: f64, n: usize) -> Vec<f64> {
    // 1/Γ(z) = z(z+1)...(z+N-1) / Γ(z+N), with N chosen to put the argument past 1
    let shift = if z >= 1.0 {
        0
    } else {
        (2.0 - z).ceil() as usize
    };
    let raised = z + shift as f64;

    // Q = 1/Γ(raised), Q' = -Qψ(raised)
    let psi: Vec<f64> = (0..=n).map(|j| polygamma(j as u32, raised)).collect();
    let mut q = vec![Gamma::gamma(raised).recip()];
    for k in 0..n {
        let next = -leibniz(k, &q, &psi);
        q.push(next);
    }
    if shift == 0 {
        return q;
    }

    // Coefficients of the polynomial z(z+1)...(z+N-1), lowest power first
    let mut poly = vec![1.0];
    for m in 0..shift {
        let mut next = vec![0.0; poly.len() + 1];
        for (i, &p) in poly.iter().enumerate() {
            next[i] += p * m as f64;
            next[i + 1] += p;
        }
        poly = next;
    }
    // Its derivatives at z, the j-th being Σ_i i!/(i-j)! poly_i z^{i-j}
    let p: Vec<f64> = (0..=n)
        .map(|j| {
            poly.iter()
                .enumerate()
                .skip(j)
                .map(|(i, &c)| {
                    c * ((i - j + 1)..=i).map(|f| f as f64).product::<f64>()
                        * z.powi((i - j) as i32)
                })
                .sum()
        })
        .collect();

    (0..=n).map(|k| leibniz(k, &p, &q)).collect()
}

/// $\Gamma(-n-\epsilon)$ as a Laurent series.
///
/// $z\\,\Gamma(z) = \Gamma(z+1)$ walked $n+1$ times isolates the pole as a bare
/// $1/\epsilon$ against a factor that is analytic there:
/// $$
///     \Gamma(-n-\epsilon) = \frac{(-1)^{n+1}}{\epsilon}
///         \frac{\Gamma(1-\epsilon)}{\prod_{k=1}^{n}(k+\epsilon)}.
/// $$
fn gamma_laurent(n: usize, w: usize) -> Laurent {
    let order = 2 * w;
    let gd = gamma_derivatives(1.0, order);
    let numerator: Vec<f64> = (0..=order)
        .map(|j| gd[j] * alternating_sign(j) / Gamma::gamma(j as f64 + 1.0))
        .collect();

    // The denominator is a product of linear factors, so it divides out one exact
    // reciprocal at a time. The pole goes on last, once every analytic factor has been
    // taken at the full width: it shifts the window down, and the top coefficient it
    // would need was never in range to begin with.
    let mut analytic = Laurent::from_taylor(w, &numerator);
    for k in 1..=n {
        analytic = &analytic * &Laurent::reciprocal_linear(w, k as isize);
    }
    &(&analytic * &Laurent::reciprocal_linear(w, 0)) * -alternating_sign(n)
}

//
// Beta function
//

/// Table of $\partial_a^j \partial_b^k B(a, b)$.
///
/// Taken through $B = \Gamma(a)\Gamma(b)\cdot\frac{1}{\Gamma(a+b)}$ rather than through
/// $\ln B$, so that a pole of $\Gamma(a+b)$ is the zero of an entire function rather than
/// an infinity to cancel against another.
pub fn beta_derivatives(a: f64, b: f64, j_max: usize, k_max: usize) -> Table {
    let g_a = gamma_derivatives(a, j_max);
    let g_b = gamma_derivatives(b, k_max);
    let r = recip_gamma_derivatives(a + b, j_max + k_max);
    leibniz_table(&g_a, &g_b, |t| r[t], j_max, k_max, &0.0)
}

/// The Leibniz sum for $\Gamma(a)\,\Gamma(b)\cdot\frac{1}{\Gamma(a+b)}$, from the
/// derivatives of the three factors, at every $\partial_a^j \partial_b^k$ up to
/// `j_max` and `k_max`.
///
/// `g_a` are plain numbers, $a$ sitting nowhere near a pole. The other two factors carry
/// whatever $\Gamma(b)$ does - numbers where it is finite, series in $\epsilon$ where it
/// is not - and `zero` says what an empty sum of those is. `r` is reached at $p+q$, that
/// factor being the only one to see both parameters.
fn leibniz_table<T>(
    g_a: &[f64],
    g_b: &[T],
    r: impl Fn(usize) -> T,
    j_max: usize,
    k_max: usize,
    zero: &T,
) -> Vec<Vec<T>>
where
    T: Clone + AddAssign<T>,
    for<'a> &'a T: Mul<&'a T, Output = T> + Mul<f64, Output = T>,
{
    // Pascal's triangle down to the deeper of the two orders, so that a row is built
    // once rather than once per row of the other index
    let c_rows: Vec<Vec<f64>> = (0..=j_max.max(k_max)).map(binomials).collect();
    let mut table = vec![vec![zero.clone(); k_max + 1]; j_max + 1];
    for j in 0..=j_max {
        let c_row_j = &c_rows[j];
        for k in 0..=k_max {
            let c_row_k = &c_rows[k];
            for p in 0..=j {
                for q in 0..=k {
                    let scale = c_row_j[p] * c_row_k[q] * g_a[j - p];
                    let product = &g_b[k - q] * &r(p + q);
                    table[j][k] += &product * scale;
                }
            }
        }
    }
    table
}

/// Table of $\partial_a^j \partial_b^k B(a,b)$ at $b = -m + \epsilon$, as Laurent series.
///
/// $\Gamma(b)$ is the only factor with a pole there, and $B$ is taken as
/// $\Gamma(a)\Gamma(b)\cdot 1/\Gamma(a+b)$ so that $\Gamma(a+b)$ enters as the zero of
/// an entire function rather than as an infinity of its own.
pub fn beta_derivatives_laurent(
    a: f64,
    m: usize,
    j_max: usize,
    k_max: usize,
    w: usize,
) -> Vec<Vec<Laurent>> {
    let g_a = gamma_derivatives(a, j_max);

    // Γ^(k)(-m+ε): the pole series read with ε for -ε, then differentiated in ε
    let mut g_b = vec![gamma_laurent(m, w).reflected_in_epsilon()];
    for k in 0..k_max {
        g_b.push(g_b[k].diff());
    }

    // R = 1/Γ at a+b = a-m+ε, entire and so an ordinary Taylor series
    let rg = recip_gamma_derivatives(a - m as f64, j_max + k_max + 2 * w + 2);
    let r_series = |t: usize| {
        let taylor: Vec<f64> = (0..=2 * w)
            .map(|i| rg[t + i] / Gamma::gamma(i as f64 + 1.0))
            .collect();
        Laurent::from_taylor(w, &taylor)
    };

    leibniz_table(&g_a, &g_b, r_series, j_max, k_max, &Laurent::zero(w))
}

//
// Incomplete beta function
//

/// Advance $P_n = (1-b)_n/n!$ and its derivatives in $b$ from $n$ to $n+1$, in place.
///
/// $P^{(k)}_{n+1} = [(1-b+n)P^{(k)}_n - k\\,P^{(k-1)}_n]/(n+1)$, walked from the top so
/// that each entry is read before it is overwritten.
fn advance_pochhammer(poch_der: &mut [f64], b: f64, n: usize) {
    let (shift, step) = (1.0 - b + n as f64, (n + 1) as f64);
    for k in (0..poch_der.len()).rev() {
        let lower = if k > 0 {
            k as f64 * poch_der[k - 1]
        } else {
            0.0
        };
        poch_der[k] = (shift * poch_der[k] - lower) / step;
    }
}

/// Largest number of terms the series is allowed before it is declared not to converge.
const MAX_SERIES_TERMS: usize = 4096;

/// $\partial_a^j \partial_b^k \sum_{n\ge0} \frac{(1-b)_n}{n!}\frac{z^{a+n}}{a+n}$.
///
/// The series for $B_z(a,b)$, differentiated term by term. Its two factors depend on
/// one parameter each, so there is no Leibniz rule to apply between them: $\partial_a$
/// reaches only $z^{a+n}/(a+n)$ and $\partial_b$ only the Pochhammer.
///
/// The Pochhammer is carried from one term to the next by [`advance_pochhammer()`]
/// rather than recomputed, which is what keeps the summation linear in the number of
/// terms.
fn inc_beta_derivatives_series(
    a: f64,
    b: f64,
    z: f64,
    j_max: usize,
    k_max: usize,
    tol: f64,
) -> Table {
    let mut out = vec![vec![0.0; k_max + 1]; j_max + 1];
    if z <= 0.0 {
        return out;
    }
    let ln_z = z.ln();
    let c_rows: Vec<Vec<f64>> = (0..=j_max).map(binomials).collect();

    // P_n = (1-b)_n/n! and its derivatives in b, starting from P_0 = 1
    let mut poch_der = vec![0.0; k_max + 1];
    poch_der[0] = 1.0;
    let mut power = z.powf(a);
    let mut scale = 0.0f64;

    for n in 0..MAX_SERIES_TERMS {
        let an = a + n as f64;
        let mut largest = 0.0f64;
        for j in 0..=j_max {
            // ∂_a^j [z^{a+n}/(a+n)] = z^{a+n} Σ_i C(j,i) ln^{j-i}(z) (-1)^i i!/(a+n)^{i+1}
            let mut d = 0.0;
            let mut falling = 1.0;
            for (i, &c) in c_rows[j].iter().enumerate() {
                d += c * ln_z.powi((j - i) as i32) * falling / an.powi(i as i32 + 1);
                falling *= -((i + 1) as f64);
            }
            d *= power;
            for k in 0..=k_max {
                let term = poch_der[k] * d;
                out[j][k] += term;
                largest = largest.max(term.abs());
                scale = scale.max(out[j][k].abs());
            }
        }
        if n > 0 && largest <= tol * scale {
            break;
        }

        advance_pochhammer(&mut poch_der, b, n);
        power *= z;
    }
    out
}

/// The series of [`inc_beta_derivatives_series()`] with its first parameter at $-m + \epsilon$.
///
/// Every term is finite but the one at $n = m$, whose $1/(a+n)$ is a bare $1/\epsilon$.
/// Carrying the whole sum as a Laurent series keeps that pole in hand until it meets
/// the one in the complete Beta and cancels.
fn inc_beta_derivatives_laurent(
    m: usize,
    b: f64,
    z: f64,
    j_max: usize,
    k_max: usize,
    w: usize,
    tol: f64,
) -> Vec<Vec<Laurent>> {
    let mut out = vec![vec![Laurent::zero(w); k_max + 1]; j_max + 1];
    if z <= 0.0 {
        return out;
    }
    let ln_z = z.ln();

    // P_n = (1-b)_n/n! and its derivatives in b, starting from P_0 = 1
    let mut poch_der = vec![0.0; k_max + 1];
    poch_der[0] = 1.0;
    let mut scale = 0.0f64;

    for n in 0..MAX_SERIES_TERMS {
        let shift = n as isize - m as isize;

        // z^{a+n} = z^shift e^{ε ln z}, as a Taylor series in ε
        let base = if shift >= 0 {
            z.powi(shift as i32)
        } else {
            1.0 / z.powi((-shift) as i32)
        };
        let mut taylor = Vec::with_capacity(2 * w + 1);
        let mut term = base;
        for i in 0..=2 * w {
            taylor.push(term);
            term *= ln_z / (i + 1) as f64;
        }
        let d = &Laurent::from_taylor(w, &taylor) * &Laurent::reciprocal_linear(w, shift);

        let mut largest = 0.0f64;
        let mut dk = d.clone();
        for row in out.iter_mut() {
            for (value, &weight) in row.iter_mut().zip(&poch_der) {
                let piece = &dk * weight;
                largest = largest.max(piece.peak());
                *value += &piece;
                scale = scale.max(value.peak());
            }
            dk = dk.diff();
        }
        if n > m && largest <= tol * scale {
            break;
        }

        advance_pochhammer(&mut poch_der, b, n);
    }
    out
}

/// Tolerance the series are summed to.
const TOL: f64 = 1e-17;

/// How near a non-positive integer `b` has to be for the expansion about it to be taken
/// instead of the plain reflection.
///
/// The reflection subtracts two halves that each grow as $1/\epsilon$ there, and loses
/// digits to it long before $\epsilon$ reaches zero - five of them at
/// $\epsilon = 10^{-12}$. The other branch subtracts the two as series, where the poles
/// cancel term by term, and is limited only by where the expansion is cut: its error
/// goes as $\epsilon^{w+1}$. The two are worth about the same here, and the worst either
/// does is a part in $10^{12}$.
const POLE_WIDTH: f64 = 1e-4;

/// Whether `b` is near enough to a pole of $\Gamma$ that only the reflection as a whole
/// is a number.
fn near_a_pole(b: f64) -> bool {
    let nearest = (-b).round();
    nearest >= 0.0 && (-b - nearest).abs() < POLE_WIDTH
}

/// The complete beta's derivative table and the deficit's, as
/// [`inc_beta_derivatives_split()`] hands them back.
pub type SplitTables = (Table, Table);

/// $\partial_a^j \partial_b^k$ of $B(a,b)$, and of the deficit $B(a,b) - B_z(a,b)$ that
/// the argument falling short of one leaves behind.
///
/// Their difference is [`inc_beta_derivatives()`], but a caller separating a singular
/// part from a regular one wants them apart: against a $\Delta^{-b}$ out front the first
/// is what survives as $\Delta \to 0$ and the second is what stays analytic there, the
/// two prefactors cancelling in the latter.
///
/// [`None`] where `b` sits near a non-positive integer. $\Gamma(b)$ has a pole there, so
/// both halves are infinite and only their difference is a number.
pub fn inc_beta_derivatives_split(
    a: f64,
    b: f64,
    z: f64,
    j_max: usize,
    k_max: usize,
) -> Option<SplitTables> {
    assert!(
        a > 0.0,
        "the first parameter of an incomplete beta must be positive"
    );
    assert!((0.0..=1.0).contains(&z), "the argument must lie in [0, 1]");
    if near_a_pole(b) {
        return None;
    }
    let complete = beta_derivatives(a, b, j_max, k_max);
    // Past z = 1/2 the reflection already is the deficit, read transposed because
    // exchanging the parameters exchanges which derivative is which. Below it the series
    // gives $B_z$ whole and the deficit is taken from the complete beta, far enough from
    // the singular point that the subtraction costs nothing.
    let deficit = if z > 0.5 {
        let reflected = inc_beta_derivatives_series(b, a, 1.0 - z, k_max, j_max, TOL);
        (0..=j_max)
            .map(|j| (0..=k_max).map(|k| reflected[k][j]).collect())
            .collect()
    } else {
        let whole = inc_beta_derivatives_series(a, b, z, j_max, k_max, TOL);
        (0..=j_max)
            .map(|j| (0..=k_max).map(|k| complete[j][k] - whole[j][k]).collect())
            .collect()
    };
    Some((complete, deficit))
}

/// Table of $\partial_a^j \partial_b^k B_z(a,b)$ for $j \le$ `j_max` and $k \le$ `k_max`.
///
/// $B_z(a,b) = \int_0^z t^{a-1}(1-t)^{b-1} dt$, with `a` positive and `b` free to be
/// negative or a non-integer. The parameters are what is differentiated, not the
/// argument: a logarithm under the integral is what one of these derivatives is.
///
/// The series of the section above converges for $z < 1$ and any real $b$, since $a > 0$
/// leaves no denominator to vanish. Past $z = 1/2$ the reflection
/// $B_z(a,b) = B(a,b) - B_{1-z}(b,a)$ turns it into the same series in $1-z$ with the
/// two parameters exchanged, and converges the faster the nearer $z$ is to one.
pub fn inc_beta_derivatives(a: f64, b: f64, z: f64, j_max: usize, k_max: usize) -> Table {
    assert!(
        a > 0.0,
        "the first parameter of an incomplete beta must be positive"
    );
    assert!((0.0..=1.0).contains(&z), "the argument must lie in [0, 1]");
    // Below the crossover the series gives $B_z$ whole, and no pole of $\Gamma$ enters
    // to be cancelled: $a > 0$ leaves no denominator to vanish
    if z <= 0.5 {
        return inc_beta_derivatives_series(a, b, z, j_max, k_max, TOL);
    }
    if let Some((complete, deficit)) = inc_beta_derivatives_split(a, b, z, j_max, k_max) {
        return (0..=j_max)
            .map(|j| {
                (0..=k_max)
                    .map(|k| complete[j][k] - deficit[j][k])
                    .collect()
            })
            .collect();
    }
    // Only the pole is left: Γ(b) has one at every non-positive integer, and so does the
    // one term of the series whose denominator b+n vanishes. Carrying each as a Laurent
    // series in ε about b = -m lets the two cancel before anything is evaluated, which
    // is why this branch cannot hand the halves back one at a time.
    let m = (-b).round() as usize;
    {
        let width = j_max + k_max + 3;
        let complete = beta_derivatives_laurent(a, m, j_max, k_max, width);
        let reflected = inc_beta_derivatives_laurent(m, a, 1.0 - z, k_max, j_max, width, TOL);
        // The poles cancel, so the difference is analytic and its non-negative part is
        // the Taylor expansion of B_z about b = -m. Summing it at the ε we were given
        // answers for that b rather than for the integer beside it.
        let eps = b + m as f64;
        (0..=j_max)
            .map(|j| {
                (0..=k_max)
                    .map(|k| {
                        let mut total = 0.0;
                        let mut power = 1.0;
                        for p in 0..=(width as isize) {
                            total += (complete[j][k].at(p) - reflected[k][j].at(p)) * power;
                            power *= eps;
                        }
                        total
                    })
                    .collect()
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::segment::Segment;
    use crate::util::bilby_integrate;
    use approx::assert_relative_eq;

    /// $\int_0^z t^{a-1}(1-t)^{b-1}\ln^j t \ln^k(1-t)\\,dt$ by quadrature, which is what
    /// the derivative table claims to be.
    fn by_quadrature(a: f64, b: f64, z: f64, j: usize, k: usize) -> f64 {
        let f = |t: f64| {
            if t <= 0.0 || t >= 1.0 {
                return 0.0;
            }
            t.powf(a - 1.0)
                * (1.0 - t).powf(b - 1.0)
                * t.ln().powi(j as i32)
                * (1.0 - t).ln().powi(k as i32)
        };
        bilby_integrate(f, Segment::new(0.0, z), 1e-13)
            .unwrap()
            .value
    }

    #[test]
    fn against_quadrature() {
        // b negative, b a non-integer, b a positive integer where the series terminates
        // but its derivatives do not, and z on either side of the branch crossover
        for &(a, b) in &[
            (0.5f64, -0.2f64),
            (0.5, -1.7),
            (1.5, 0.3),
            (0.25, 2.0),
            (1.0, 1.0),
            (2.5, -0.5),
            // b a non-positive integer, where the two halves of the reflection are
            // each infinite: r = 0 is what two inverse square roots come to
            (0.5, 0.0),
            (1.5, 0.0),
            (0.5, -1.0),
            (0.25, -2.0),
        ] {
            for &z in &[0.1f64, 0.4, 0.49, 0.51, 0.75, 0.95] {
                let table = inc_beta_derivatives(a, b, z, 2, 2);
                for (j, row) in table.iter().enumerate() {
                    for (k, &got) in row.iter().enumerate() {
                        let want = by_quadrature(a, b, z, j, k);
                        assert_relative_eq!(got, want, max_relative = 1e-11, epsilon = 1e-13);
                    }
                }
            }
        }
    }

    #[test]
    fn near_a_pole() {
        // Neither branch is comfortable just off a non-positive integer: the generic
        // one cancels two halves that are each growing without bound, and the pole one
        // answers a question about the integer rather than the argument given. Between
        // them they have to stay honest anyway.
        let (a, z) = (0.5f64, 0.9f64);
        for eps in [1e-12f64, 1e-10, 1e-9, 5e-9, 1e-8, 1e-7, 1e-6, 1e-4] {
            for side in [1.0f64, -1.0] {
                for m in [0.0f64, 1.0, 2.0] {
                    let b = -m + side * eps;
                    let got = inc_beta_derivatives(a, b, z, 0, 0)[0][0];
                    let want = by_quadrature(a, b, z, 0, 0);
                    assert_relative_eq!(got, want, max_relative = 1e-7);
                }
            }
        }
    }

    #[test]
    fn the_two_branches_agree() {
        // Either side of the crossover the answer must not jump
        for &(a, b) in &[(0.5f64, -0.2f64), (1.5, 0.3), (0.25, 2.0)] {
            let below = inc_beta_derivatives(a, b, 0.5, 2, 2);
            let above = inc_beta_derivatives(a, b, 0.5 + 1e-11, 2, 2);
            for (lower, upper) in below.iter().flatten().zip(above.iter().flatten()) {
                assert_relative_eq!(lower, upper, max_relative = 1e-9);
            }
        }
    }

    #[test]
    fn the_poles_cancel() {
        // Both halves of the reflection are infinite at a non-positive integer b, and
        // the answer is finite only because their pole parts agree exactly
        for &(a, m) in &[(0.5f64, 0usize), (1.5, 0), (0.5, 1), (0.25, 2)] {
            let (j_max, k_max) = (2usize, 2usize);
            let width = j_max + k_max + 3;
            let complete = beta_derivatives_laurent(a, m, j_max, k_max, width);
            let reflected = inc_beta_derivatives_laurent(m, a, 0.25, k_max, j_max, width, 1e-17);
            let mut seen = 0;
            for (j, row) in complete.iter().enumerate() {
                for (k, entry) in row.iter().enumerate() {
                    for p in 1..=(k + 1) {
                        let (x, y) = (entry.at(-(p as isize)), reflected[k][j].at(-(p as isize)));
                        if x.abs() > 1e-12 {
                            seen += 1;
                            assert_relative_eq!(x, y, max_relative = 1e-10);
                        }
                    }
                }
            }
            assert!(seen > 0, "a = {a}, m = {m} has no pole to cancel");
        }
    }

    #[test]
    fn the_ends_are_the_complete_beta() {
        // z = 1 is the complete beta, derivatives and all; z = 0 is nothing
        for &(a, b) in &[(0.5f64, 0.75f64), (1.5, 2.5), (2.0, 0.25)] {
            let whole = inc_beta_derivatives(a, b, 1.0, 2, 2);
            let complete = beta_derivatives(a, b, 2, 2);
            for (got, want) in whole.iter().flatten().zip(complete.iter().flatten()) {
                assert_relative_eq!(got, want, max_relative = 1e-12);
            }
            for row in inc_beta_derivatives(a, b, 0.0, 2, 2) {
                assert!(row.iter().all(|v| *v == 0.0));
            }
        }
    }
}
