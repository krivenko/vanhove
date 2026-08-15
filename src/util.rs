/// Utility functions
use bilby::{
    QuadratureError, QuadratureResult, adaptive_integrate, integrate_infinite,
    integrate_semi_infinite_lower, integrate_semi_infinite_upper,
};

/// Kahan-Babuška-Neumaier summation algorithm
pub fn kahan_babushka_neumaier_sum<I: Iterator<Item = f64>>(input: I) -> f64 {
    let (mut sum, mut c) = (0.0f64, 0.0f64);
    for x in input {
        let t = sum + x;
        c += if sum.abs() >= x.abs() {
            (sum - t) + x
        } else {
            (x - t) + sum
        };
        sum = t;
    }
    sum + c
}

/// Call a bilby adaptive integration function depending on the integration limits.
pub fn bilby_integrate<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    tol: f64,
) -> Result<QuadratureResult<f64>, QuadratureError> {
    let a_inf = a.is_infinite();
    let b_inf = b.is_infinite();
    // Choose a bilby call depending on the integration limits
    match (a_inf, b_inf) {
        (false, false) => adaptive_integrate(f, a, b, tol),
        (false, true) => {
            assert!(b > 0.0);
            integrate_semi_infinite_upper(f, a, tol)
        }
        (true, false) => {
            assert!(a < 0.0);
            integrate_semi_infinite_lower(f, b, tol)
        }
        (true, true) => {
            assert!(a < 0.0 && b > 0.0);
            integrate_infinite(f, tol)
        }
    }
}

/// Fermi step function
pub fn fermi(x: f64) -> f64 {
    if x <= 0.0 {
        1.0 / (1.0 + x.exp())
    } else {
        1.0 - fermi(-x)
    }
}

#[cfg(test)]
mod tests {
    use crate::util;
    use approx::assert_abs_diff_eq;

    #[test]
    fn kahan_babushka_neumaier_sum() {
        let v = vec![10000.0f64, 3.14159f64, 2.71828f64];
        assert_abs_diff_eq!(
            util::kahan_babushka_neumaier_sum(v.into_iter()),
            10005.85987,
            epsilon = 1e-12
        );
    }

    #[test]
    fn bilby_integrate() {
        let pi = std::f64::consts::PI;
        let inf = std::f64::INFINITY;
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| x.cos() * x.cos(), -pi, pi, 1e-12)
                .unwrap()
                .value,
            pi,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| (-x * x / 2.0).exp(), -inf, inf, 1e-12)
                .unwrap()
                .value,
            (2.0 * pi).sqrt(),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| (2.0 * x).exp(), -inf, 0.0, 1e-12)
                .unwrap()
                .value,
            0.5,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            util::bilby_integrate(move |x| (-2.0 * x).exp(), 0.0, inf, 1e-12)
                .unwrap()
                .value,
            0.5,
            epsilon = 1e-12
        );
    }
}
