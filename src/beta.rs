//! The beta function and its derivatives with respect to its parameters.
//!
//! A convolution of two singular parts is an incomplete beta function: a power against
//! a power over a finite range is $B_z(a,b)$, and a logarithm on either factor is
//! a derivative with respect to one of the two exponents. The derivatives here are
//! therefore with respect to the *parameters*, not the argument, which is what no library
//! provides.

use crate::singularity::power_log_integral;
use crate::util::{Table, alternating_sign, binomials, polygamma, subtract_tables};
use special::Gamma;

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

/// Derivatives of $\Gamma(z)$, up to order `n`, where `z` is no pole.
fn gamma_derivatives(z: f64, n: usize) -> Vec<f64> {
    // $\Gamma' = \Gamma\psi$ differentiated by the general Leibniz rule
    let psi: Vec<f64> = (0..=n).map(|j| polygamma(j as u32, z)).collect();
    let mut d = vec![Gamma::gamma(z)];
    for k in 0..n {
        let next = leibniz(k, &d, &psi);
        d.push(next);
    }
    d
}

/// Derivatives of $1/\Gamma(z)$ up to order `n`.
///
/// $1/\Gamma$ is entire, so these stay finite at the poles of $\Gamma$, which is the
/// whole reason the beta-function derivatives are taken through it: $B(a, b)$ is a
/// perfectly good number where $a+b$ is a non-positive integer, and only the
/// way of writing it down falls over there.
fn recip_gamma_derivatives(z: f64, n: usize) -> Vec<f64> {
    // The recurrence $z\,\Gamma(z) = \Gamma(z+1)$ walks the argument up to where
    // $1/\Gamma = \exp(-\ln\Gamma)$ may be differentiated directly, leaving a
    // polynomial factor behind:
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

//
// Beta function
//

/// Table of $\partial_a^j \partial_b^k B(a, b)$, finite where $a+b$ is a non-positive
/// integer.
pub fn beta_derivatives(a: f64, b: f64, j_max: usize, k_max: usize) -> Table {
    // Taken through $B = \Gamma(a)\Gamma(b)\cdot\frac{1}{\Gamma(a+b)}$ rather than
    // through $\ln B$, so that a pole of $\Gamma(a+b)$ is the zero of an entire function
    // rather than an infinity to cancel against another
    let g_a = gamma_derivatives(a, j_max);
    let g_b = gamma_derivatives(b, k_max);
    let inv_g_ab = recip_gamma_derivatives(a + b, j_max + k_max);
    leibniz_table(&g_a, &g_b, &inv_g_ab, j_max, k_max)
}

/// The Leibniz sum for $\Gamma(a)\,\Gamma(b)\cdot\frac{1}{\Gamma(a+b)}$, from the
/// derivatives of the three factors, at every $\partial_a^j \partial_b^k$ up to
/// `j_max` and `k_max`.
///
/// `inv_g_ab` is indexed at $p+q$, that factor being the only one to see both parameters,
/// so it runs to `j_max + k_max`.
fn leibniz_table(g_a: &[f64], g_b: &[f64], inv_g_ab: &[f64], j_max: usize, k_max: usize) -> Table {
    // Pascal's triangle down to the deeper of the two orders, so that a row is built once
    let c_rows: Vec<Vec<f64>> = (0..=j_max.max(k_max)).map(binomials).collect();
    let mut table = vec![vec![0.0; k_max + 1]; j_max + 1];
    for j in 0..=j_max {
        let c_row_j = &c_rows[j];
        for k in 0..=k_max {
            let c_row_k = &c_rows[k];
            for p in 0..=j {
                for q in 0..=k {
                    let scale = c_row_j[p] * c_row_k[q] * g_a[j - p];
                    table[j][k] += g_b[k - q] * inv_g_ab[p + q] * scale;
                }
            }
        }
    }
    table
}

//
// Incomplete beta function
//

/// Advance $P_n = (1-b)_n/n!$ and its derivatives in $b$ from $n$ to $n+1$, in place.
fn advance_pochhammer(poch_der: &mut [f64], b: f64, n: usize) {
    // $P^{(k)}_{n+1} = [(1-b+n)P^{(k)}_n - k\,P^{(k-1)}_n]/(n+1)$, walked from the top
    // so that each entry is read before it is overwritten
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

/// $\partial_a^j \partial_b^k B_z(a,b)$ from the ascending series
/// $\sum_{n\ge0} \frac{(1-b)_n}{n!}\frac{z^{a+n}}{a+n}$, summed until the largest term
/// falls below `tol` times the largest entry so far.
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
    // Differentiated term by term. The two factors take one parameter each, so no
    // Leibniz rule comes between them: $\partial_a$ reaches only $z^{a+n}/(a+n)$ and
    // $\partial_b$ only the Pochhammer symbol, which `advance_pochhammer()` propagates
    // from one term to the next and so keeps the summation linear in the term count.
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

/// $\partial_\delta^k \int_{l_1}^{l_2} e^{\delta u} du = \int_{l_1}^{l_2} u^k e^{\delta u} du$,
/// for $k$ up to `k_max`.
///
/// An `l1` of $-\infty$ is [`exp_moments_to()`], which unlike a finite range has a pole
/// at $\delta = 0$.
fn exp_moments(delta: f64, k_max: usize, l1: f64, l2: f64) -> Vec<f64> {
    if l1 == f64::NEG_INFINITY {
        return exp_moments_to(delta, k_max, l2);
    }
    if (delta * l1).abs() <= 2.0 {
        // ∫ u^k Σ_i (δu)^i/i! du = Σ_i δ^i/i! (l2^{k+i+1} - l1^{k+i+1})/(k+i+1)
        return (0..=k_max)
            .map(|k| {
                let mut total = 0.0;
                let mut weight = 1.0;
                for i in 0..MAX_EXP_MOMENT_TERMS {
                    let p = (k + i + 1) as i32;
                    let term = weight * (l2.powi(p) - l1.powi(p)) / f64::from(p);
                    total += term;
                    if i > 2 && term.abs() <= f64::EPSILON * total.abs() {
                        break;
                    }
                    weight *= delta / (i + 1) as f64;
                }
                total
            })
            .collect();
    }
    let (upper, lower) = (
        exp_moments_to(delta, k_max, l2),
        exp_moments_to(delta, k_max, l1),
    );
    let mut out = Vec::with_capacity(k_max + 1);
    for k in 0..=k_max {
        out.push(upper[k] - lower[k]);
    }
    out
}

/// Terms the Taylor branch of [`exp_moments()`] is allowed, which at $|\delta l_1| \le 2$
/// it never approaches.
const MAX_EXP_MOMENT_TERMS: usize = 64;

/// $\int_{-\infty}^{l} u^k e^{\delta u} du$ for $k$ up to `k_max`.
///
/// The integral converges only for $\delta > 0$ and this is its continuation elsewhere.
/// It has a pole at $\delta = 0$, of order $k+1$, which a finite range does not, so
/// [`exp_moments()`] is what to ask for near zero.
fn exp_moments_to(delta: f64, k_max: usize, l: f64) -> Vec<f64> {
    // $e^{\delta l}\sum_{i\le k}\frac{(-1)^{k-i}k!}{i!\,\delta^{k-i+1}}l^i$
    let e = (delta * l).exp();
    (0..=k_max)
        .map(|k| {
            let c_row = binomials(k);
            let mut total = 0.0;
            // C(k,i)(k-i)! = k!/i!, walked down from i = k where (-1)^{k-i}(k-i)! is one
            let mut falling = 1.0;
            for i in (0..=k).rev() {
                let weight = c_row[i] * falling / delta.powi((k - i + 1) as i32);
                total += weight * l.powi(i as i32) * e;
                falling *= -((k - i + 1) as f64);
            }
            total
        })
        .collect()
}

/// $\partial_a^j \partial_b^k B_z(a,b)$ for $z$ past $1/2$, with no pole anywhere.
///
/// Nothing here forms the complete beta, which at a non-positive integer $b$ is no
/// number at all, and the term count does not grow as $z$ approaches one.
fn inc_beta_derivatives_near_one(
    a: f64,
    b: f64,
    z: f64,
    j_max: usize,
    k_max: usize,
    tol: f64,
) -> Table {
    // Splitting $\int_y^1 s^{b-1}(1-s)^{a-1} ds$ at $s = 1/2$, where $y = 1 - z$, and
    // expanding $(1-s)^{a-1}$ over the lower piece,
    // $$
    //     B_z(a,b) = B_{1/2}(a,b)
    //         + \sum_{n\ge0}\frac{(1-a)_n}{n!}\frac{(1/2)^{b+n} - y^{b+n}}{b+n}.
    // $$
    // Each summand is $\int_{\ln y}^{-\ln 2} e^{(b+n)u} du$ and so is entire in $b$,
    // worth $\ln\frac{1}{2y}$ where the denominator vanishes. Both pieces fall off at a
    // ratio of $1/2$ or better, whatever $y$. As in `inc_beta_derivatives_series()`, the
    // two factors take one parameter each, $\partial_a$ reaching only the Pochhammer
    // symbol and $\partial_b$ only the integral, so no Leibniz rule comes between them.
    let y = 1.0 - z;
    let (l1, l2) = (y.ln(), -std::f64::consts::LN_2);
    let mut out = inc_beta_derivatives_series(a, b, 0.5, j_max, k_max, tol);

    // P_n = (1-a)_n/n! and its derivatives in a, starting from P_0 = 1
    let mut poch_der = vec![0.0; j_max + 1];
    poch_der[0] = 1.0;

    // Each entry is weighed against its own running size, not against the largest of
    // them. At a negative `b` the table spans orders of magnitude - the $y^{b+n}$ of the
    // first few terms is enormous, while $\ln^j t$ suppresses the same end of the range
    // and leaves $\partial_a^j B_z$ small - and one shared scale would cut the small
    // entries off while the large ones are still converging.
    let mut scale: Table = out
        .iter()
        .map(|row| row.iter().map(|v| v.abs()).collect())
        .collect();

    for n in 0..MAX_SERIES_TERMS {
        let moments = exp_moments(b + n as f64, k_max, l1, l2);
        let mut converged = n > 0;
        for (j, row) in out.iter_mut().enumerate() {
            for (k, value) in row.iter_mut().enumerate() {
                let term = poch_der[j] * moments[k];
                *value += term;
                scale[j][k] = scale[j][k].max(value.abs());
                converged &= term.abs() <= tol * scale[j][k];
            }
        }
        if converged {
            break;
        }

        advance_pochhammer(&mut poch_der, a, n);
    }
    out
}

/// Tolerance the series are summed to.
const TOL: f64 = 1e-17;

/// How near a whole number $\rho$ may come before the pole of $\Gamma(-\rho)$ swamps
/// the finite part beside it.
///
/// Splitting a beta function about such a $\rho$ puts that pole into both parts, the
/// complete value and the remainder alike, and only their difference is finite.
/// $\Gamma(-\rho)$ grows as $1/\epsilon$, so a caller that forms both and subtracts
/// them back loses roughly $\log_{10}(1/\epsilon)$ digits to the cancellation: four at
/// this width, ten by $\epsilon = 10^{-12}$. [`inc_beta_derivatives()`] is under no such
/// restriction, reaching the difference without forming either part.
pub const WHOLE_EXPONENT_WIDTH: f64 = 1e-4;

/// Whether $\rho$ lies within [`WHOLE_EXPONENT_WIDTH`] of a non-negative integer.
pub fn near_a_whole_exponent(rho: f64) -> bool {
    let nearest = rho.round();
    nearest >= 0.0 && (rho - nearest).abs() < WHOLE_EXPONENT_WIDTH
}

/// The complete beta's derivative table and the deficit's, as
/// [`inc_beta_derivatives_split()`] returns them.
pub type SplitTables = (Table, Table);

/// $\partial_a^j \partial_b^k$ of $B(a,b)$, and of the deficit $B(a,b) - B_z(a,b)$ that
/// the argument falling short of one leaves behind.
///
/// Their difference is [`inc_beta_derivatives()`], but a caller separating a singular
/// part from a regular one wants them apart.
///
/// [`None`] where `b` sits near a non-positive integer. $\Gamma(b)$ has a pole there, so
/// both pieces are infinite and only their difference is a number.
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
    if near_a_whole_exponent(-b) {
        return None;
    }
    let complete = beta_derivatives(a, b, j_max, k_max);
    // Past z = 1/2 the series taken at $1-z$ with the parameters exchanged already is
    // the deficit, read transposed because exchanging them exchanges which derivative is
    // which. Below it the series gives $B_z$ whole and the deficit is taken from the
    // complete beta.
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
/// negative or a non-integer. Each derivative introduces a power of a logarithm under
/// the integral.
///
/// No pole of $\Gamma$ is met at any $z < 1$, whatever real $b$ comes to. Only $z = 1$
/// itself is the complete beta, which at a non-positive integer `b` is no number at all.
pub fn inc_beta_derivatives(a: f64, b: f64, z: f64, j_max: usize, k_max: usize) -> Table {
    assert!(
        a > 0.0,
        "the first parameter of an incomplete beta must be positive"
    );
    assert!((0.0..=1.0).contains(&z), "the argument must lie in [0, 1]");
    // Below the crossover the expansion is in $z$ about the lower end of the range,
    // above it in $1-z$ about the upper end, and each falls off at a ratio of $1/2$ or
    // better
    if z <= 0.5 {
        return inc_beta_derivatives_series(a, b, z, j_max, k_max, TOL);
    }
    if z == 1.0 {
        return beta_derivatives(a, b, j_max, k_max);
    }
    inc_beta_derivatives_near_one(a, b, z, j_max, k_max, TOL)
}

//
// Outer beta function
//

/// Taylor coefficients of $(1-x)^{\alpha}$, up to $x^{n_{max}}$.
fn binomial_series(alpha: f64, n_max: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(n_max + 1);
    let mut c = 1.0;
    for n in 0..=n_max {
        out.push(c);
        c *= (n as f64 - alpha) / (n + 1) as f64;
    }
    out
}

/// Taylor coefficients of $\ln^m(1-x)$ for every $m$ up to `m_max`, each one taken to
/// $x^{n_{max}}$.
fn log_powers(m_max: usize, n_max: usize) -> Vec<Vec<f64>> {
    let mut base = vec![0.0; n_max + 1];
    for (i, c) in base.iter_mut().enumerate().skip(1) {
        *c = -1.0 / i as f64;
    }
    let mut out = vec![vec![0.0; n_max + 1]];
    out[0][0] = 1.0;
    for m in 1..=m_max {
        let mut next = vec![0.0; n_max + 1];
        for i in 0..=n_max {
            let p = out[m - 1][i];
            if p == 0.0 {
                continue;
            }
            for l in 0..=(n_max - i) {
                next[i + l] += p * base[l];
            }
        }
        out.push(next);
    }
    out
}

/// Terms a series of the given `ratio` needs to reach [`TOL`], with room to spare for
/// the polynomial in $i$ the series here multiply $q^i$ by.
fn series_length(ratio: f64) -> usize {
    // A bare $q^i$ needs $\ln(\text{TOL})/\ln q$ terms, 57 at $q = \tfrac12$. The
    // polynomial's degree grows with the derivative order asked for, and $i^4 q^i$ needs
    // 82, which the 32 covers.
    let n = (TOL.ln() / ratio.ln()).ceil() as usize;
    (n + 32).min(MAX_SERIES_TERMS)
}

/// Table of $\partial_a^j \partial_b^k B^{\mathrm{out}}_U(a,b)$, where
/// $$
///     B^{\mathrm{out}}_U(a,b) = \int_0^U (1+u)^{a-1} u^{b-1} du.
/// $$
///
/// This is the integral an outer stretch of a convolution comes to, at $a = r_2+1$ and
/// $b = r_1+1$ once the distance to the singular point is scaled out. Both parameters
/// are then positive, and the exponent the pair produces enters only as
/// $\rho = a+b-1$.
///
/// $\partial_a$ puts one $\ln(1+u)$ under the integral and $\partial_b$ one $\ln u$, so
/// index $j$ counts the first and index $k$ the second.
pub fn outer_beta_derivatives(a: f64, b: f64, u_max: f64, j_max: usize, k_max: usize) -> Table {
    assert!(b > 0.0, "the outer beta needs b > 0 to converge at u = 0");
    assert!(
        u_max >= 0.0 && !u_max.is_nan(),
        "the outer beta runs over a non-negative range"
    );
    // The range is cut at $u = 1$ so that neither expansion is asked to run past a
    // ratio of $\tfrac12$, however large $U$ is
    let mut out = outer_beta_below_one(a, b, u_max.min(1.0), j_max, k_max);
    if u_max > 1.0 {
        let far = outer_beta_above_one(a, b, u_max, j_max, k_max);
        for j in 0..=j_max {
            for k in 0..=k_max {
                out[j][k] += far[j][k];
            }
        }
    }
    out
}

/// Table of $\partial_a^j \partial_b^k B(b, -\rho)$ at $\rho = a+b-1$: the complete beta,
/// written in the outer beta's own parameters.
///
/// This is $B^{\mathrm{out}}_U$ at an infinite $U$ only for $\rho < 0$. The integral
/// diverges for every $\rho \ge 0$, logarithmically already at $\rho = 0$, and what the
/// table holds past that is its continuation in $\rho$, which need not be positive.
///
/// At a whole $\rho \ge 0$ it is no number at all as $B(b, -\rho)$ has a pole there.
pub fn outer_beta_complete(a: f64, b: f64, j_max: usize, k_max: usize) -> Table {
    assert!(b > 0.0, "the outer beta needs b > 0 to converge at u = 0");
    // At a whole rho the i = n term divides by zero and the table contains inf and NaN.
    debug_assert!(
        !near_a_whole_exponent(a + b - 1.0),
        "the complete outer beta has a pole at a whole rho"
    );
    let mut out = outer_beta_below_one(a, b, 1.0, j_max, k_max);
    let far = outer_beta_far_half(
        a,
        b,
        j_max,
        k_max,
        f64::NEG_INFINITY,
        -std::f64::consts::LN_2,
    );
    for (j, row) in out.iter_mut().enumerate() {
        for (k, entry) in row.iter_mut().enumerate() {
            *entry += far[j][k];
        }
    }
    out
}

/// [`outer_beta_complete()`] and [`outer_beta_tail()`], the two pieces
/// $B^{\mathrm{out}}_U$ is the difference of.
///
/// [`None`] where $\rho$ sits within [`WHOLE_EXPONENT_WIDTH`] of a whole number, neither piece
/// being a number there and only their difference, and [`None`] for a `u_max` that is
/// not finite.
pub fn outer_beta_derivatives_split(
    a: f64,
    b: f64,
    u_max: f64,
    j_max: usize,
    k_max: usize,
) -> Option<SplitTables> {
    assert!(b > 0.0, "the outer beta needs b > 0 to converge at u = 0");
    assert!(
        u_max >= 0.0,
        "the outer beta runs over a non-negative range"
    );
    let rho = a + b - 1.0;
    if !u_max.is_finite() || near_a_whole_exponent(rho) {
        return None;
    }
    let complete = outer_beta_complete(a, b, j_max, k_max);
    let tail = outer_beta_tail(a, b, u_max, &complete, j_max, k_max);
    Some((complete, tail))
}

/// Table of $\partial_a^j \partial_b^k [B(b, -\rho) - B^{\mathrm{out}}_U(a,b)]$ at
/// $U =$ `u_min`.
///
/// `complete` is that $B(b, -\rho)$, [`outer_beta_complete()`] at the same parameters.
/// It does not depend on the range, so a caller taking several tails of one outer beta
/// computes it once.
pub fn outer_beta_tail(
    a: f64,
    b: f64,
    u_min: f64,
    complete: &Table,
    j_max: usize,
    k_max: usize,
) -> Table {
    // Past the cut the tail is the expansion in $y$ over $(0, 1/(1+U)]$. Short of the
    // cut that expansion does not reach, and the tail is the complete value less the
    // integral over the range instead.
    if u_min >= 1.0 {
        let l1 = (1.0 / (1.0 + u_min)).ln();
        return outer_beta_far_half(a, b, j_max, k_max, f64::NEG_INFINITY, l1);
    }
    let whole = outer_beta_derivatives(a, b, u_min, j_max, k_max);
    subtract_tables(complete, &whole)
}

/// The part of [`outer_beta_derivatives()`] over $u \le 1$, where `u_max` must lie.
fn outer_beta_below_one(a: f64, b: f64, u_max: f64, j_max: usize, k_max: usize) -> Table {
    // $t = u/(1+u)$ turns the integrand into $t^{b-1}(1-t)^{-a-b}$ over $[0, w]$ for
    // $w = U/(1+U) \le \tfrac12$, with $\ln(1+u) = -\ln(1-t)$ and
    // $\ln u = \ln t - \ln(1-t)$. Expanding $(1-t)^{-a-b}$ against its logarithms
    // leaves $\int_0^w t^{b+n-1}\ln^M t\,dt$, which `power_log_integral()` gives in
    // closed form, and $b + n > 0$ keeps every one of them finite.
    let mut out = vec![vec![0.0; k_max + 1]; j_max + 1];
    if u_max <= 0.0 {
        return out;
    }
    let w = u_max / (1.0 + u_max);
    let n_max = series_length(w);
    let degree = j_max + k_max;

    // $\lambda^{(m)} = (1-t)^{-a-b}(-\ln(1-t))^m$, the sign of the logarithm folded in
    let base = binomial_series(-(a + b), n_max);
    let logs = log_powers(degree, n_max);
    let mut lambda = vec![vec![0.0; n_max + 1]; degree + 1];
    for m in 0..=degree {
        let sign = alternating_sign(m);
        for i in 0..=n_max {
            let p = base[i];
            if p == 0.0 {
                continue;
            }
            for l in 0..=(n_max - i) {
                lambda[m][i + l] += sign * p * logs[m][l];
            }
        }
    }

    // $\int_0^w t^{b+n-1}\ln^M t\\,dt$, wanted at every $n$ and every $M \le k_{max}$
    let mut moment = vec![vec![0.0; k_max + 1]; n_max + 1];
    for (n, row) in moment.iter_mut().enumerate() {
        for (m, entry) in row.iter_mut().enumerate() {
            *entry = power_log_integral(b + n as f64 - 1.0, m as u8, w);
        }
    }

    let c_rows: Vec<Vec<f64>> = (0..=k_max).map(binomials).collect();
    for j in 0..=j_max {
        for k in 0..=k_max {
            let mut total = 0.0;
            for p in 0..=k {
                let weight = c_rows[k][p];
                for n in 0..=n_max {
                    total += weight * lambda[j + p][n] * moment[n][k - p];
                }
            }
            out[j][k] = total;
        }
    }
    out
}

/// The part of [`outer_beta_derivatives()`] over $1 < u \le U$.
fn outer_beta_above_one(a: f64, b: f64, u_max: f64, j_max: usize, k_max: usize) -> Table {
    let l1 = (1.0 / (1.0 + u_max)).ln();
    outer_beta_far_half(a, b, j_max, k_max, l1, -std::f64::consts::LN_2)
}

/// The part of the outer beta beyond $u = 1$, over the range in $y$ between $e^{l_1}$
/// and $e^{l_2}$.
///
/// An `l1` of $-\infty$ is the range running all the way out to $u = \infty$.
fn outer_beta_far_half(a: f64, b: f64, j_max: usize, k_max: usize, l1: f64, l2: f64) -> Table {
    // $y = 1/(1+u)$ turns the integrand into $(1-y)^{b-1}y^{-\rho-1}$, with
    // $\ln(1+u) = -\ln y$ and $\ln u = \ln(1-y) - \ln y$. Expanding $(1-y)^{b-1}$
    // against its logarithms leaves $\int y^{i-\rho-1}\ln^M y\,dy$, which is
    // `exp_moments()` at $d = i - \rho$. The expansion runs at a ratio of $\tfrac12$
    // whatever the range, $u = 1$ being $y = \tfrac12$.
    let rho = a + b - 1.0;
    let n_max = series_length(0.5);
    let degree = j_max + k_max;

    // $\xi^{(q)} = (1-y)^{b-1}\ln^q(1-y)$
    let base = binomial_series(b - 1.0, n_max);
    let logs = log_powers(k_max, n_max);
    let mut xi = vec![vec![0.0; k_max + 1]; n_max + 1];
    for q in 0..=k_max {
        for (i, &p) in base.iter().enumerate() {
            if p == 0.0 {
                continue;
            }
            for l in 0..=(n_max - i) {
                xi[i + l][q] += p * logs[q][l];
            }
        }
    }

    let c_rows: Vec<Vec<f64>> = (0..=k_max).map(binomials).collect();
    let mut out = vec![vec![0.0; k_max + 1]; j_max + 1];
    for (i, xi_row) in xi.iter().enumerate() {
        let e = exp_moments(i as f64 - rho, degree, l1, l2);
        for (j, row) in out.iter_mut().enumerate() {
            for (k, entry) in row.iter_mut().enumerate() {
                let mut total = 0.0;
                for q in 0..=k {
                    let m = j + k - q;
                    total += c_rows[k][q] * alternating_sign(m) * xi_row[q] * e[m];
                }
                *entry += total;
            }
        }
    }
    out
}

/// The whole of [`outer_beta_complete()`] about a whole $\rho$ but for its pole:
/// `out[j][l]` is $\partial_b^j f_l$, where $B^{\mathrm{out}}_\infty = F(\delta)/\delta$
/// at $\rho = n - \delta$ and $F(\delta) = \sum_l f_l\delta^l$.
///
/// The pole there is of order one and no worse, so $F(\delta)$ is analytic.
///
/// A caller reads the pole of $\partial_b^j\partial_\delta^k$ off $f_0$ alone and its
/// $\delta^0$ part off $f_{k+1}$, and reaches the outer beta's own two directions
/// through $\partial_a = -\partial_\delta$ and
/// $\partial_b\vert_a = \partial_b - \partial_\delta$. The parameters are held apart
/// as $b$ and $\delta$ because $\rho$ and not $a$ is what the pole turns on.
pub fn outer_beta_entire_factor(b: f64, n: usize, j_max: usize, l_max: usize) -> Table {
    // Cut the range at $u = 1$. Below it is `outer_beta_below_one()`, finite in
    // $\delta$, and above it $y = 1/(1+u)$ gives
    // $\int_0^{1/2}(1-y)^{b-1}y^{-\rho-1}dy$, whose expansion
    // $\sum_{i\ge0}\frac{(1-b)_i}{i!}\frac{2^{\rho-i}}{i-\rho}$ runs at a ratio of
    // $\tfrac12$. Exactly one of its terms is singular at $\rho = n$, so
    // $$
    //     F(\delta) = 2^{-\delta}\left[\frac{(1-b)_n}{n!}
    //         + \delta\sum_{i\ne n}\frac{(1-b)_i}{i!}\frac{2^{n-i}}{i-n+\delta}\right]
    //         + \delta B^{\mathrm{out}}_1,
    // $$
    // and the residue $f_0 = (1-b)_n/n!$ is a term of that series rather than a $\Gamma$
    // product. Where $b$ is a positive integer no greater than $n$ the Pochhammer symbol
    // is exactly zero, which is the pair of locally constant factors.
    let factorial = |l: usize| Gamma::gamma(l as f64 + 1.0);
    let n_max = series_length(0.5).max(n);

    // g[l][j] = ∂_b^j of the δ^l coefficient of the bracket, whose i-th summand takes
    // (1-b)_i/i! and δ/(d+δ) = Σ_{l≥1} (-1)^{l-1} δ^l / d^l
    let mut g = vec![vec![0.0; j_max + 1]; l_max + 1];
    let mut poch = vec![0.0; j_max + 1];
    poch[0] = 1.0;
    for i in 0..=n_max {
        if i == n {
            for (j, &p) in poch.iter().enumerate() {
                g[0][j] += p;
            }
        } else {
            let d = i as f64 - n as f64;
            let mut c = 2.0f64.powi(n as i32 - i as i32) / d;
            for g_row in g.iter_mut().skip(1) {
                for (j, &p) in poch.iter().enumerate() {
                    g_row[j] += p * c;
                }
                c /= -d;
            }
        }
        advance_pochhammer(&mut poch, b, i);
    }

    let mut out = vec![vec![0.0; l_max + 1]; j_max + 1];
    for (j, row) in out.iter_mut().enumerate() {
        for (l, entry) in row.iter_mut().enumerate() {
            for (v, g_row) in g.iter().enumerate().take(l + 1) {
                let ln2_power = (-std::f64::consts::LN_2).powi((l - v) as i32);
                *entry += g_row[j] * ln2_power / factorial(l - v);
            }
        }
    }
    if l_max == 0 {
        return out;
    }

    // $B^{out}_1$ is taken in its own two parameters, so the table it returns is
    // rotated by $\partial_\delta = -\partial_a$ and $\partial_b|_\delta = \partial_b -
    // \partial_a$ before it is read. $\delta B^{out}_1$ reaches $f_l$ through
    // $\partial_\delta^{l-1} B^{out}_1 / (l-1)!$
    let near = outer_beta_below_one(n as f64 + 1.0 - b, b, 1.0, j_max + l_max - 1, j_max);
    let c_rows: Vec<Vec<f64>> = (0..=j_max).map(binomials).collect();
    for (j, row) in out.iter_mut().enumerate() {
        for (l, entry) in row.iter_mut().enumerate().skip(1) {
            let m = l - 1;
            let mut total = 0.0;
            for (i, &c) in c_rows[j].iter().enumerate() {
                total += c * alternating_sign(m + i) * near[i + m][j - i];
            }
            *entry += total / factorial(m);
        }
    }
    out
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

    /// $\int_0^U (1+u)^{a-1}u^{b-1}\ln^j(1+u)\ln^k u\,du$ by quadrature.
    fn outer_by_quadrature(a: f64, b: f64, u_max: f64, j: usize, k: usize) -> f64 {
        let f = |u: f64| {
            if u <= 0.0 {
                return 0.0;
            }
            (1.0 + u).powf(a - 1.0)
                * u.powf(b - 1.0)
                * (1.0 + u).ln().powi(j as i32)
                * u.ln().powi(k as i32)
        };
        bilby_integrate(f, Segment::new(0.0, u_max), 1e-13)
            .unwrap()
            .value
    }

    #[test]
    fn the_outer_beta_against_quadrature() {
        for &(a, b) in &[(0.5f64, 0.5f64), (1.5, 0.5), (0.25, 1.25), (2.0, 1.0)] {
            for &u_max in &[0.01f64, 0.5, 1.0, 3.0, 40.0] {
                let table = outer_beta_derivatives(a, b, u_max, 2, 2);
                for (j, row) in table.iter().enumerate() {
                    for (k, &got) in row.iter().enumerate() {
                        let want = outer_by_quadrature(a, b, u_max, j, k);
                        assert_relative_eq!(got, want, max_relative = 1e-9, epsilon = 1e-13);
                    }
                }
            }
        }
    }

    #[test]
    fn the_complete_outer_beta_is_the_beta_reparameterised() {
        // $B^{out}_\infty(a,b) = B(b, 1-a-b)$, under the same rotation of the
        // derivatives as the incomplete pair. Taken above and below $\rho = 0$, where
        // the integral stops converging and the continuation takes over, and clear of
        // every whole $\rho$, which is the one place neither side is a number.
        let (j_max, k_max) = (2usize, 2usize);
        for &(a, b) in &[
            (0.5f64, 0.25f64),
            (0.1, 0.3),
            (1.5, 0.7),
            (1.25, 0.9),
            (2.4, 1.3),
        ] {
            let want = beta_derivatives(b, 1.0 - a - b, k_max, j_max + k_max);
            let mine = outer_beta_complete(a, b, j_max, k_max);
            for j in 0..=j_max {
                for k in 0..=k_max {
                    let c_row = binomials(k);
                    let mut rotated = 0.0;
                    for i in 0..=k {
                        rotated += c_row[i] * alternating_sign(j + i) * want[k - i][j + i];
                    }
                    assert_relative_eq!(mine[j][k], rotated, max_relative = 1e-9, epsilon = 1e-12);
                }
            }
        }
    }

    #[test]
    fn the_outer_split_adds_back_up() {
        // The two endpoints of the same exp_moments(): the complete value less the tail
        // beyond $U$ is the integral over $[0, U]$, however far out $U$ runs, and on
        // either side of the cut at $u = 1$ that the tail is read about.
        let (j_max, k_max) = (2usize, 2usize);
        for &(a, b) in &[(0.5f64, 0.25f64), (0.1, 0.3), (1.5, 0.7), (2.4, 1.3)] {
            for &u_max in &[0.02f64, 0.6, 1.0, 7.0, 1e3, 1e8] {
                let (complete, tail) =
                    outer_beta_derivatives_split(a, b, u_max, j_max, k_max).unwrap();
                let whole = outer_beta_derivatives(a, b, u_max, j_max, k_max);
                for j in 0..=j_max {
                    for k in 0..=k_max {
                        let got = complete[j][k] - tail[j][k];
                        let scale = complete[j][k].abs().max(whole[j][k].abs());
                        assert_relative_eq!(
                            got,
                            whole[j][k],
                            max_relative = 1e-9,
                            epsilon = 1e-12 * scale
                        );
                    }
                }
            }
        }

        // Nothing to divide where both pieces are infinite, whatever the range
        assert!(outer_beta_derivatives_split(0.5, 0.5, 4.0, 1, 1).is_none());
        assert!(outer_beta_derivatives_split(0.5, 0.5, 0.2, 1, 1).is_none());
        assert!(outer_beta_derivatives_split(1.5, 0.5, f64::INFINITY, 1, 1).is_none());
    }

    #[test]
    fn the_outer_beta_is_the_incomplete_one_reparameterised() {
        // $t = u/(1+u)$ gives $B^{out}_U(a,b) = B_w(b, 1-a-b)$ at $w = U/(1+U)$, so
        // $\partial_a = -\partial_B$ and $\partial_b = \partial_A - \partial_B$. The
        // two routes share no code, and this one reaches $U$ where quadrature cannot.
        let (j_max, k_max) = (2usize, 2usize);
        for &(a, b) in &[(0.5f64, 0.5f64), (1.5, 0.5), (0.25, 1.25), (1.1, 0.9)] {
            for &u_max in &[0.3f64, 1.0, 7.0, 1e3, 1e8] {
                let w = u_max / (1.0 + u_max);
                let inc = inc_beta_derivatives(b, 1.0 - a - b, w, k_max, j_max + k_max);
                let mine = outer_beta_derivatives(a, b, u_max, j_max, k_max);
                for j in 0..=j_max {
                    for k in 0..=k_max {
                        let c_row = binomials(k);
                        let mut want = 0.0;
                        for i in 0..=k {
                            want += c_row[i] * alternating_sign(j + i) * inc[k - i][j + i];
                        }
                        assert_relative_eq!(mine[j][k], want, max_relative = 1e-8, epsilon = 1e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn a_whole_rho_is_no_special_case_for_the_outer_beta() {
        // $\rho = a+b-1$ whole is where `outer_beta_complete()` has a pole of
        // $\Gamma(-\rho)$. Over a finite range it is one parameter value among others,
        // and the table walks through it: the pole belongs to the range running all the
        // way out, not to $\rho$ being whole.
        let b = 0.5f64;
        for n in [0.0f64, 1.0, 2.0] {
            for eps in [0.0f64, 1e-12, -1e-12, 1e-6, -1e-6] {
                let a = n + eps + 1.0 - b;
                let table = outer_beta_derivatives(a, b, 50.0, 1, 1);
                for (j, row) in table.iter().enumerate() {
                    for (k, &got) in row.iter().enumerate() {
                        let want = outer_by_quadrature(a, b, 50.0, j, k);
                        assert_relative_eq!(got, want, max_relative = 1e-9, epsilon = 1e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn the_entire_factor_rebuilds_the_beta() {
        // F(δ)/δ is B(a, -n+δ) wherever B is a number at all, and the same row by row:
        // $\partial_a$ is taken at a fixed δ, so it commutes with the expansion
        let j_max = 2;
        for &(a, n) in &[(0.5f64, 0usize), (0.5, 1), (1.5, 2), (0.25, 3), (2.5, 1)] {
            // Truncated at $\delta^8$, so the comparison stays near the expansion point
            let f = outer_beta_entire_factor(a, n, j_max, 8);
            for &delta in &[1e-2f64, -1e-2, 0.05, -0.05] {
                let want = beta_derivatives(a, -(n as f64) + delta, j_max, 0);
                for (j, row) in f.iter().enumerate() {
                    let mut series = 0.0;
                    for (l, &c) in row.iter().enumerate() {
                        series += c * delta.powi(l as i32);
                    }
                    assert_relative_eq!(series / delta, want[j][0], max_relative = 1e-9);
                }
            }
        }
    }

    #[test]
    fn the_residue_is_a_polynomial() {
        // f_0 is $(-1)^n/n! \prod_{i=1}^{n}(a-i)$, so the whole pole of B about a
        // non-positive integer is elementary - no Γ survives into it
        for &(a, n) in &[(0.5f64, 0usize), (0.5, 1), (1.5, 2), (0.25, 3), (3.0, 1)] {
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            let want = (1..=n).fold(sign / Gamma::gamma(n as f64 + 1.0), |p, i| {
                p * (a - i as f64)
            });
            assert_relative_eq!(
                outer_beta_entire_factor(a, n, 0, 0)[0][0],
                want,
                epsilon = 1e-14
            );
        }

        // At a = n the product has a zero factor and the pole goes.
        assert_eq!(outer_beta_entire_factor(1.0, 1, 0, 0)[0][0], 0.0);
        assert_eq!(outer_beta_entire_factor(2.0, 2, 0, 0)[0][0], 0.0);
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
            // b a non-positive integer, where the complete beta is no number at all:
            // rho = 0 is what two inverse square roots come to
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
    fn a_pole_is_no_special_case() {
        // A non-positive integer `b` is where the complete beta has a pole, so any
        // route through it returns an infinity. Nothing singles it out here: the answer
        // walks through, from either side and at the integer itself.
        let (a, z) = (0.5f64, 0.9f64);
        for eps in [0.0f64, 1e-12, 1e-10, 1e-9, 5e-9, 1e-8, 1e-7, 1e-6, 1e-4] {
            for side in [1.0f64, -1.0] {
                for m in [0.0f64, 1.0, 2.0] {
                    let b = -m + side * eps;
                    let got = inc_beta_derivatives(a, b, z, 0, 0)[0][0];
                    let want = by_quadrature(a, b, z, 0, 0);
                    assert_relative_eq!(got, want, max_relative = 1e-11);
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
    fn a_whole_exponent_is_no_special_case() {
        // At a non-positive integer b the complete beta is infinite, while the
        // expansion about z = 1 is an ordinary sum of finite terms. Every derivative
        // comes out, at every z the branch covers.
        for &(a, b) in &[
            (0.5f64, 0.0f64),
            (0.5, -1.0),
            (0.5, -2.0),
            (0.5, -3.0),
            (1.5, -1.0),
            (0.25, -2.0),
        ] {
            // Nearer to one than this the quadrature is no longer the more accurate of
            // the two; `a_whole_exponent_against_a_closed_form` takes over from here.
            for &z in &[0.51f64, 0.75, 0.9, 1.0 - 1e-6] {
                let table = inc_beta_derivatives(a, b, z, 2, 2);
                for (j, row) in table.iter().enumerate() {
                    for (k, &got) in row.iter().enumerate() {
                        let want = by_quadrature(a, b, z, j, k);
                        assert_relative_eq!(got, want, max_relative = 1e-9, epsilon = 1e-13);
                    }
                }
            }
        }
    }

    #[test]
    fn a_whole_exponent_against_a_closed_form() {
        // b = 0 and a = 1/2 integrate in elementary terms, which pins the branch where
        // the quadrature cannot follow it:
        // $\int_0^z t^{-1/2}(1-t)^{-1} dt = \ln\frac{1+\sqrt z}{1-\sqrt z}$, read as
        // $\ln\frac{(1+\sqrt z)^2}{1-z}$ so that no subtraction cancels near $z = 1$.
        //
        // The comparison is at $1-z$ rather than at a nominal small number, because
        // `1.0 - (1.0 - y) != y` below about 1e-3 and $B_z$ magnifies the difference by
        // $1/(1-z)$: the two would part company at 8e-7 by $1-z = 10^{-12}$, saying
        // nothing about the branch.
        for &z in &[
            0.51f64,
            0.9,
            1.0 - 1e-3,
            1.0 - 1e-6,
            1.0 - 1e-9,
            1.0 - 1e-12,
        ] {
            let y = 1.0 - z;
            let want = ((1.0 + z.sqrt()).powi(2) / y).ln();
            let got = inc_beta_derivatives(0.5, 0.0, z, 0, 0)[0][0];
            assert_relative_eq!(got, want, max_relative = 1e-14);
        }
    }

    #[test]
    fn the_ends_are_the_complete_beta() {
        // z = 1 is the complete beta, derivatives and all; z = 0 gives zero
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
