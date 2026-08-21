# vanhove

[![crates.io](https://img.shields.io/crates/v/vanhove.svg)](https://crates.io/crates/vanhove)
[![docs](https://img.shields.io/badge/docs-vanhove-blue)](https://krivenko.github.io/vanhove/)
[![license](https://img.shields.io/badge/license-MIT%20or%20Apache--2.0-blue.svg)](#license)

Model densities of states with integrable van Hove singularities, and accurate
integration of spectral functions.

Many quantities in condensed matter physics are spectral integrals
$\int A(\omega) f(\omega) d\omega$ of a spectral function $A(\omega)$ against
some function $f(\omega)$. Doing this numerically is awkward whenever
$A(\omega)$ has van Hove singularities: they defeat generic adaptive quadrature
even though the integrals themselves converge.

This crate handles them by splitting each singular model into a smooth part
$R(\omega)$ and analytically known asymptotics $S_p(\omega)$ around every
singular point $\Omega_p$,

```math
\int A(\omega) f(\omega) d\omega = \int R(\omega) f(\omega) d\omega
    + \sum_p \int S_p(\omega) [f(\omega) - f(\Omega_p)] d\omega
    + \sum_p f(\Omega_p) \int S_p(\omega) d\omega.
```

Only the first two terms are integrated numerically, and both integrands are
regular; the last term uses a closed-form value of $\int S_p(\omega) d\omega$.
The numerical work is delegated to the
[bilby](https://crates.io/crates/bilby) quadrature.

## Models

| Function                    | Density of states                                            |
| --------------------------- | ------------------------------------------------------------ |
| `discrete(levels, weights)` | discrete levels, $\sum_p w_p \delta(\omega - \varepsilon_p)$ |
| `flat(eps, d, delta)`       | flat band with Fermi-like edges, sharp for `delta = 0`        |
| `gaussian(eps, sigma)`      | Gaussian                                                     |
| `semicircle(eps, r)`        | semicircle (Wigner), square-root band edges                  |
| `chain(eps, t)`             | linear chain, square-root divergent edges                    |
| `square(eps, t)`            | square lattice, logarithmic van Hove peak                    |

Every model is normalized to unit spectral weight. They can be scaled by real
numbers and added together, so mixed discrete/continuous spectra are built by
simple arithmetic.

## Usage

```toml
[dependencies]
vanhove = "0.1"
```

```rust
use num_complex::Complex64;
use vanhove::models::{discrete, semicircle, square};

// 80% square-lattice band plus 20% of a single discrete level at ω = 3
let dos = 0.8 * square(0.0, 1.0) + 0.2 * discrete(&[3.0], &[1.0]);
assert_eq!(dos.total_weight(), 1.0);

// First spectral moment ∫A(ω) ω dω
let m1 = dos.integrate(|omega| omega, None).unwrap();

// Retarded Green's function of a semicircular band at ω + i0⁺
let z = Complex64::new(0.5, 1e-3);
let g = semicircle(0.0, 2.0)
    .integrate_complex(|omega| 1.0 / (z - omega))
    .unwrap();
```

`integrate()` takes an optional absolute tolerance (`1e-10` by default) and
returns a `Result`, since the underlying adaptive quadrature can fail to
converge for a badly behaved `f`.

## Documentation

API documentation is available at
<https://krivenko.github.io/vanhove/>.

Formulas in the docs are rendered with KaTeX, so building them locally requires
the extra header:

```sh
RUSTDOCFLAGS="--html-in-header katex-header.html" cargo doc --no-deps --open
```

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or
[MIT license](LICENSE-MIT) at your option.
