//! The beta function and its derivatives with respect to its parameters.
//!
//! A convolution of two singular parts is a beta function: a power against a power over
//! a finite range is (a,b)$, and a logarithm on either factor is a derivative with
//! respect to one of the two exponents. The derivatives here are therefore with respect
//! to the *parameters*, not the argument, which is what no library offers.

use special::Gamma;

use crate::util::{binomials, is_natural, polygamma};

/// Derivatives of $\Gamma$ at `x`, up to order `n`.
///
/// $\Gamma' = \Gamma\psi$ differentiated by Leibniz, which needs `x` to be no pole.
fn gamma_derivatives(x: f64, n: usize) -> Vec<f64> {
    let psi: Vec<f64> = (0..=n).map(|j| polygamma(j as u32, x)).collect();
    let mut d = vec![Gamma::gamma(x)];
    for k in 0..n {
        let c = binomials(k);
        d.push((0..=k).map(|j| c[j] * d[k - j] * psi[j]).sum());
    }
    d
}

/// Derivatives of $1/\Gamma$ at `z`, up to order `n`.
///
/// $1/\Gamma$ is entire, so these stay finite at the poles of $\Gamma$, which is the
/// whole reason the Beta derivatives are taken through it: an outer stretch lands on
/// $\alpha+\gamma = -r_2$, a pole whenever the other exponent is a non-negative integer.
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
    let w = z + shift as f64;

    // Q = 1/Γ(w), Q' = -Qψ(w)
    let psi: Vec<f64> = (0..=n).map(|j| polygamma(j as u32, w)).collect();
    let mut q = vec![Gamma::gamma(w).recip()];
    for k in 0..n {
        let c = binomials(k);
        q.push(-(0..=k).map(|j| c[j] * q[k - j] * psi[j]).sum::<f64>());
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

    (0..=n)
        .map(|k| {
            let c = binomials(k);
            (0..=k).map(|i| c[i] * p[k - i] * q[i]).sum()
        })
        .collect()
}

/// Table of $\partial_\alpha^a \partial_\gamma^b B(\alpha, \gamma)$.
///
/// Taken through $B = \Gamma(\alpha)\Gamma(\gamma)\cdot\frac{1}{\Gamma(\alpha+\gamma)}$
/// rather than through $\ln B$, so that a pole of $\Gamma(\alpha+\gamma)$ is the zero of
/// an entire function rather than an infinity to cancel against another.
fn beta_derivatives(alpha: f64, gamma: f64, a_max: usize, b_max: usize) -> Vec<Vec<f64>> {
    let ga = gamma_derivatives(alpha, a_max);
    let gc = gamma_derivatives(gamma, b_max);
    let r = recip_gamma_derivatives(alpha + gamma, a_max + b_max);
    (0..=a_max)
        .map(|a| {
            let ca = binomials(a);
            (0..=b_max)
                .map(|b| {
                    let cb = binomials(b);
                    (0..=a)
                        .flat_map(|i| (0..=b).map(move |j| (i, j)))
                        .map(|(i, j)| ca[i] * cb[j] * ga[a - i] * gc[b - j] * r[i + j])
                        .sum()
                })
                .collect()
        })
        .collect()
}

/// Truncated Laurent series in $\epsilon$, holding the coefficients of $\epsilon^{-w}$
/// through $\epsilon^{w}$.
#[derive(Clone)]
struct Laurent {
    w: usize,
    c: Vec<f64>,
}

impl Laurent {
    fn zero(w: usize) -> Laurent {
        Laurent {
            w,
            c: vec![0.0; 2 * w + 1],
        }
    }

    /// A series with no pole, from the coefficients of $\epsilon^0, \epsilon^1, \ldots$
    fn from_taylor(w: usize, taylor: &[f64]) -> Laurent {
        let mut s = Laurent::zero(w);
        for (j, &v) in taylor.iter().enumerate() {
            if w + j < s.c.len() {
                s.c[w + j] = v;
            }
        }
        s
    }

    /// Coefficient of $\epsilon^k$, zero outside the window.
    fn at(&self, k: isize) -> f64 {
        let i = k + self.w as isize;
        if i < 0 || i as usize >= self.c.len() {
            0.0
        } else {
            self.c[i as usize]
        }
    }

    fn mul(&self, other: &Laurent) -> Laurent {
        let mut r = Laurent::zero(self.w);
        for (i, &a) in self.c.iter().enumerate() {
            if a == 0.0 {
                continue;
            }
            for (j, &b) in other.c.iter().enumerate() {
                let k = (i + j) as isize - other.w as isize;
                if k >= 0 && (k as usize) < r.c.len() {
                    r.c[k as usize] += a * b;
                }
            }
        }
        r
    }

    /// $d/d\epsilon$, term by term.
    fn diff(&self) -> Laurent {
        let mut r = Laurent::zero(self.w);
        for (i, &v) in self.c.iter().enumerate() {
            let k = i as isize - self.w as isize;
            if k != 0 && v != 0.0 && i > 0 {
                r.c[i - 1] += v * k as f64;
            }
        }
        r
    }

    fn scaled(&self, f: f64) -> Laurent {
        Laurent {
            w: self.w,
            c: self.c.iter().map(|v| v * f).collect(),
        }
    }

    fn add(&mut self, other: &Laurent) {
        for (a, b) in self.c.iter_mut().zip(&other.c) {
            *a += b;
        }
    }
}

/// $\Gamma(-n-\epsilon)$ as a Laurent series.
///
/// $z\\,\Gamma(z) = \Gamma(z+1)$ walked $n+1$ times isolates the pole as a bare
/// $1/\epsilon$ against a factor that is analytic there:
/// $$
///     \Gamma(-n-\epsilon) = \frac{(-1)^{n+1}}{\epsilon}
///         \frac{\Gamma(1-\epsilon)}{\prod_{k=1}^{n}(k+\epsilon)}.
/// $$
fn gamma_pole_series(n: usize, w: usize) -> Laurent {
    let order = 2 * w;
    let gd = gamma_derivatives(1.0, order);
    let factorial = |k: usize| (1..=k).map(|i| i as f64).product::<f64>();
    let numerator: Vec<f64> = (0..=order)
        .map(|j| {
            let sign = if j.is_multiple_of(2) { 1.0 } else { -1.0 };
            gd[j] * sign / factorial(j)
        })
        .collect();

    // Coefficients of Π_{k=1}^{n}(k+ε), lowest power first, and its reciprocal series
    let mut den = vec![1.0];
    for k in 1..=n {
        let mut next = vec![0.0; den.len() + 1];
        for (i, &c) in den.iter().enumerate() {
            next[i] += c * k as f64;
            next[i + 1] += c;
        }
        den = next;
    }
    let mut inv = vec![1.0 / den[0]];
    for j in 1..=order {
        let acc: f64 = (1..=j.min(den.len() - 1))
            .map(|i| den[i] * inv[j - i])
            .sum();
        inv.push(-acc / den[0]);
    }

    let analytic = Laurent::from_taylor(w, &numerator).mul(&Laurent::from_taylor(w, &inv));
    let mut over_eps = Laurent::zero(w);
    over_eps.c[w - 1] = 1.0;
    let sign = if n.is_multiple_of(2) { -1.0 } else { 1.0 };
    analytic.mul(&over_eps).scaled(sign)
}

/// Largest number of terms the series is allowed before it is declared not to converge.
const MAX_TERMS: usize = 4096;

/// $\partial_p^j \partial_q^k \sum_{n\ge0} \frac{(1-q)_n}{n!}\frac{w^{p+n}}{p+n}$.
///
/// The series for $B_w(p,q)$, differentiated term by term. Its two factors depend on
/// one parameter each, so there is no Leibniz rule to apply between them: $\partial_p$
/// reaches only $w^{p+n}/(p+n)$ and $\partial_q$ only the Pochhammer.
///
/// The Pochhammer is carried by the recurrence
/// $C^{(k)}_{n+1} = [(1-q+n)C^{(k)}_n - k\\,C^{(k-1)}_n]/(n+1)$ on $C_n = (1-q)_n/n!$
/// rather than through $\ln(1-q)_n$. That divides out the factorial before it can
/// overflow, and it stays right where a factor vanishes — $q$ a positive integer
/// terminates the series without terminating its derivatives.
fn series(p: f64, q: f64, w: f64, j_max: usize, k_max: usize, tol: f64) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0; k_max + 1]; j_max + 1];
    if w <= 0.0 {
        return out;
    }
    let ln_w = w.ln();
    let rows: Vec<Vec<f64>> = (0..=j_max).map(binomials).collect();

    // C_n = (1-q)_n/n! and its derivatives in q, starting from C_0 = 1
    let mut c = vec![0.0; k_max + 1];
    c[0] = 1.0;
    let mut power = w.powf(p);
    let mut scale = 0.0f64;

    for n in 0..MAX_TERMS {
        let pn = p + n as f64;
        let mut largest = 0.0f64;
        for j in 0..=j_max {
            // ∂_p^j [w^{p+n}/(p+n)] = w^{p+n} Σ_i C(j,i) ln^{j-i}(w) (-1)^i i!/(p+n)^{i+1}
            let mut d = 0.0;
            let mut falling = 1.0;
            for (i, &coeff) in rows[j].iter().enumerate() {
                d += coeff * ln_w.powi((j - i) as i32) * falling / pn.powi(i as i32 + 1);
                falling *= -((i + 1) as f64);
            }
            d *= power;
            for k in 0..=k_max {
                let term = c[k] * d;
                out[j][k] += term;
                largest = largest.max(term.abs());
                scale = scale.max(out[j][k].abs());
            }
        }
        if n > 0 && largest <= tol * scale {
            break;
        }

        let shift = 1.0 - q + n as f64;
        let step = (n + 1) as f64;
        for k in (0..=k_max).rev() {
            let lower = if k > 0 { k as f64 * c[k - 1] } else { 0.0 };
            c[k] = (shift * c[k] - lower) / step;
        }
        power *= w;
    }
    out
}

/// $1/(\text{shift} + \epsilon)$ as a Laurent series, a pole where `shift` vanishes.
fn reciprocal_linear(shift: isize, w: usize) -> Laurent {
    if shift == 0 {
        let mut out = Laurent::zero(w);
        out.c[w - 1] = 1.0;
        return out;
    }
    // 1/(s+ε) = (1/s) Σ (-ε/s)^i
    let s = shift as f64;
    let mut taylor = Vec::with_capacity(w + 1);
    let mut term = 1.0 / s;
    for _ in 0..=w {
        taylor.push(term);
        term /= -s;
    }
    Laurent::from_taylor(w, &taylor)
}

/// The same series with `ε` replaced by `-ε`.
fn reflected_in_epsilon(series: &Laurent) -> Laurent {
    let mut out = series.clone();
    for (i, value) in out.c.iter_mut().enumerate() {
        let power = i as isize - series.w as isize;
        if power.rem_euclid(2) == 1 {
            *value = -*value;
        }
    }
    out
}

/// $\partial_b^c B(a,b)$ at $b = -m + \epsilon$, as Laurent series.
///
/// $\Gamma(b)$ is the only factor with a pole there, and $B$ is taken as
/// $\Gamma(a)\Gamma(b)\cdot 1/\Gamma(a+b)$ so that $\Gamma(a+b)$ enters as the zero of
/// an entire function rather than as an infinity of its own.
fn complete_at_pole(a: f64, m: usize, j_max: usize, k_max: usize, w: usize) -> Vec<Vec<Laurent>> {
    let factorial = |k: usize| (1..=k).map(|i| i as f64).product::<f64>();
    let g_a = gamma_derivatives(a, j_max);

    // Γ^(c)(-m+ε): the pole series read with ε for -ε, then differentiated in ε
    let mut g_b = vec![reflected_in_epsilon(&gamma_pole_series(m, w))];
    for c in 0..k_max {
        g_b.push(g_b[c].diff());
    }

    // R = 1/Γ at a+b = a-m+ε, entire and so an ordinary Taylor series
    let rg = recip_gamma_derivatives(a - m as f64, j_max + k_max + 2 * w + 2);
    let r_series = |t: usize| {
        let taylor: Vec<f64> = (0..=2 * w).map(|i| rg[t + i] / factorial(i)).collect();
        Laurent::from_taylor(w, &taylor)
    };

    let mut table = vec![vec![Laurent::zero(w); k_max + 1]; j_max + 1];
    for j in 0..=j_max {
        let cj = binomials(j);
        for k in 0..=k_max {
            let ck = binomials(k);
            for i in 0..=j {
                for l in 0..=k {
                    let scale = cj[i] * ck[l] * g_a[j - i];
                    let term = g_b[k - l].mul(&r_series(i + l)).scaled(scale);
                    table[j][k].add(&term);
                }
            }
        }
    }
    table
}

/// The series of [`series()`] with its first parameter at $-m + \epsilon$.
///
/// Every term is finite but the one at $n = m$, whose $1/(p+n)$ is a bare $1/\epsilon$.
/// Carrying the whole sum as a Laurent series keeps that pole in hand until it meets
/// the one in the complete Beta and cancels.
fn series_at_pole(
    m: usize,
    q: f64,
    v: f64,
    p_max: usize,
    q_max: usize,
    w: usize,
    tol: f64,
) -> Vec<Vec<Laurent>> {
    let mut out = vec![vec![Laurent::zero(w); q_max + 1]; p_max + 1];
    if v <= 0.0 {
        return out;
    }
    let ln_v = v.ln();

    let mut c = vec![0.0; q_max + 1];
    c[0] = 1.0;
    let mut scale = 0.0f64;

    for n in 0..MAX_TERMS {
        let shift = n as isize - m as isize;

        // v^{p+n} = v^shift e^{ε ln v}, as a Taylor series in ε
        let base = if shift >= 0 {
            v.powi(shift as i32)
        } else {
            1.0 / v.powi((-shift) as i32)
        };
        let mut taylor = Vec::with_capacity(2 * w + 1);
        let mut term = base;
        for i in 0..=2 * w {
            taylor.push(term);
            term *= ln_v / (i + 1) as f64;
        }
        let d = Laurent::from_taylor(w, &taylor).mul(&reciprocal_linear(shift, w));

        let mut largest = 0.0f64;
        let mut dk = d.clone();
        let peak = |l: &Laurent| l.c.iter().fold(0.0f64, |x, y| x.max(y.abs()));
        for row in out.iter_mut() {
            for (value, &weight) in row.iter_mut().zip(&c) {
                let piece = dk.scaled(weight);
                largest = largest.max(peak(&piece));
                value.add(&piece);
                scale = scale.max(peak(value));
            }
            dk = dk.diff();
        }
        if n > m && largest <= tol * scale {
            break;
        }

        let shift_q = 1.0 - q + n as f64;
        let step = (n + 1) as f64;
        for k in (0..=q_max).rev() {
            let lower = if k > 0 { k as f64 * c[k - 1] } else { 0.0 };
            c[k] = (shift_q * c[k] - lower) / step;
        }
    }
    out
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
/// two parameters exchanged, which is the branch a convolution spends its time on: the
/// argument approaches one exactly as the frequency approaches a singular point.
// Taken up by the convolution of two singular parts, which is the next step
#[allow(dead_code)]
pub fn incomplete_beta_derivatives(
    a: f64,
    b: f64,
    z: f64,
    j_max: usize,
    k_max: usize,
) -> Vec<Vec<f64>> {
    assert!(
        a > 0.0,
        "the first parameter of an incomplete beta must be positive"
    );
    assert!((0.0..=1.0).contains(&z), "the argument must lie in [0, 1]");
    const TOL: f64 = 1e-17;

    if z <= 0.5 {
        return series(a, b, z, j_max, k_max, TOL);
    }
    // Exchanging the parameters exchanges which derivative is which, so the reflected
    // table is read transposed throughout.
    //
    // Where `b` is a non-positive integer both halves are infinite and their difference
    // is not: Γ(b) has a pole, and so does the one term of the series whose denominator
    // b+n vanishes. Carrying each as a Laurent series in ε about b = -m keeps the two
    // poles until they meet.
    if is_natural(-b) {
        let m = (-b) as usize;
        let width = j_max + k_max + 3;
        let complete = complete_at_pole(a, m, j_max, k_max, width);
        let reflected = series_at_pole(m, a, 1.0 - z, k_max, j_max, width, TOL);
        return (0..=j_max)
            .map(|j| {
                (0..=k_max)
                    .map(|k| complete[j][k].at(0) - reflected[k][j].at(0))
                    .collect()
            })
            .collect();
    }
    let complete = beta_derivatives(a, b, j_max, k_max);
    let reflected = series(b, a, 1.0 - z, k_max, j_max, TOL);
    (0..=j_max)
        .map(|j| {
            (0..=k_max)
                .map(|k| complete[j][k] - reflected[k][j])
                .collect()
        })
        .collect()
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
                let table = incomplete_beta_derivatives(a, b, z, 2, 2);
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
    fn the_two_branches_agree() {
        // Either side of the crossover the answer must not jump
        for &(a, b) in &[(0.5f64, -0.2f64), (1.5, 0.3), (0.25, 2.0)] {
            let below = incomplete_beta_derivatives(a, b, 0.5, 2, 2);
            let above = incomplete_beta_derivatives(a, b, 0.5 + 1e-11, 2, 2);
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
            let complete = complete_at_pole(a, m, j_max, k_max, width);
            let reflected = series_at_pole(m, a, 0.25, k_max, j_max, width, 1e-17);
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
            let whole = incomplete_beta_derivatives(a, b, 1.0, 2, 2);
            let complete = beta_derivatives(a, b, 2, 2);
            for (got, want) in whole.iter().flatten().zip(complete.iter().flatten()) {
                assert_relative_eq!(got, want, max_relative = 1e-12);
            }
            for row in incomplete_beta_derivatives(a, b, 0.0, 2, 2) {
                assert!(row.iter().all(|v| *v == 0.0));
            }
        }
    }
}
