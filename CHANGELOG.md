# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - Unreleased

### Added

- `segment::Segment`, a segment of the frequency axis.

### Changed

- `SpectralFunction::support()` and `DiscreteSF::support()` return a `Segment` in
  place of a pair of frequencies.

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
