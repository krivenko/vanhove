# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - Unreleased

### Added

- `SpectralFunction::conv()` accepts two continuous spectral functions. The
  singular structure of the result is derived rather than supplied: its van Hove
  points are the pairwise sums of the frequencies where either factor stops being
  smooth, and the asymptotics at each follows from the asymptotic forms of the
  two.
- `models::simple_cubic()`, density of states of a simple cubic lattice, computed
  as the convolution of a linear chain with a square lattice.
- The continuous spectral function interface is public: the `ContinuousSF` trait,
  `singularity::Singularity` and `singularity::AsymptTerm` describing a singular
  part in closed form, `interp::InterpolatedSF`, and
  `SpectralFunction::from_continuous()`. A model defined outside the crate can
  now be turned into a spectral function.
- `segment::Segment`, a segment of the frequency axis.
- The `theory` module, a prose account of the splitting of a spectral function
  into discrete, regular and singular parts, what `integrate()` does with it, and
  the mathematics behind `conv()`.
- `SpectralFunction` is `Send` and `Sync`, so that frequency scans can be spread
  over threads.

### Changed

- `SpectralFunction::conv()` takes the tolerance the result is interpolated to.
- `SpectralFunction::support()` and `DiscreteSF::support()` return a `Segment` in
  place of a pair of frequencies.
- Adding a spectral function to one it shares a continuous contribution with sums
  the weights of that contribution rather than listing it twice, and drops it where
  the weights cancel.

## [0.1.1] - 2026-09-10

### Added

- Unary minus and subtraction for `SpectralFunction` and `DiscreteSF`.
- `models::bethe()`, density of states of a Bethe lattice with a finite coordination
  number.

## [0.1.0] - 2026-08-28

Initial public release.

[0.2.0]: https://github.com/krivenko/vanhove/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/krivenko/vanhove/releases/tag/v0.1.1
[0.1.0]: https://github.com/krivenko/vanhove/releases/tag/v0.1.0
