/// Utility functions

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

#[cfg(test)]
mod tests {
    use crate::util;

    #[test]
    fn kahan_babushka_neumaier_sum() {
        let v = vec![10000.0f64, 3.14159f64, 2.71828f64];
        assert_eq!(
            util::kahan_babushka_neumaier_sum(v.into_iter()),
            10005.85987
        );
    }
}
