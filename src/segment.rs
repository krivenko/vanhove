//! Segments of the frequency axis.

use std::ops::Add;

/// Segment of the frequency axis, $[\omega_{min}, \omega_{max}]$.
///
/// A segment of zero length is valid, since it is the support of a discrete spectral
/// function of a single resonance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    min: f64,
    max: f64,
}

impl Segment {
    /// Make a segment spanning `min` to `max`.
    pub fn new(min: f64, max: f64) -> Segment {
        assert!(min <= max, "segment bounds must satisfy min <= max");
        assert!(
            min < max || min.is_finite(),
            "a segment of zero length must sit at a finite frequency"
        );
        Segment { min, max }
    }

    /// Lower end of the segment, $\omega_{min}$.
    pub fn min(&self) -> f64 {
        self.min
    }

    /// Upper end of the segment, $\omega_{max}$.
    pub fn max(&self) -> f64 {
        self.max
    }

    /// Whether `omega` lies within the segment, its ends included.
    pub fn contains(&self, omega: f64) -> bool {
        self.min <= omega && omega <= self.max
    }

    /// Whether `omega` lies within the segment, its ends excluded.
    pub fn strictly_contains(&self, omega: f64) -> bool {
        self.min < omega && omega < self.max
    }

    /// Whether the segment is of a finite length.
    pub fn is_bounded(&self) -> bool {
        self.min.is_finite() && self.max.is_finite()
    }

    /// Whether the segment is of zero length.
    pub fn is_degenerate(&self) -> bool {
        self.min == self.max
    }

    /// Length of the segment, infinite where it is unbounded.
    pub fn length(&self) -> f64 {
        self.max - self.min
    }

    /// Frequency halfway between the ends.
    pub fn midpoint(&self) -> f64 {
        // Halving each end in turn keeps the sum from overflowing
        0.5 * self.min + 0.5 * self.max
    }

    /// Image of the segment under the reflection $\nu \mapsto \omega - \nu$.
    pub fn mirrored(&self, omega: f64) -> Segment {
        Segment {
            min: omega - self.max,
            max: omega - self.min,
        }
    }

    /// The same segment displaced in frequency by `by`.
    pub fn shifted(&self, by: f64) -> Segment {
        Segment {
            min: self.min + by,
            max: self.max + by,
        }
    }

    /// The parts of the segment below and above `omega`.
    ///
    /// `omega` must lie within the segment. The two parts share the frequency they are
    /// split at, so a split at an end of the segment leaves a degenerate part.
    pub fn split_at(&self, omega: f64) -> (Segment, Segment) {
        assert!(
            self.contains(omega),
            "a segment can only be split at a frequency within it"
        );
        (
            Segment {
                min: self.min,
                max: omega,
            },
            Segment {
                min: omega,
                max: self.max,
            },
        )
    }

    /// Smallest segment containing all of `segments`, or [`None`] where there are none.
    pub fn hull<I: IntoIterator<Item = Segment>>(segments: I) -> Option<Segment> {
        segments.into_iter().reduce(|hull, segment| Segment {
            min: f64::min(hull.min, segment.min),
            max: f64::max(hull.max, segment.max),
        })
    }

    /// Overlap of two segments, or [`None`] where they are disjoint.
    ///
    /// Segments meeting at a single frequency overlap in the degenerate segment there.
    pub fn intersection(&self, other: Segment) -> Option<Segment> {
        let min = f64::max(self.min, other.min);
        let max = f64::min(self.max, other.max);
        (min <= max).then_some(Segment { min, max })
    }
}

/// Sum of two segments, $[\omega^1_{min} + \omega^2_{min}, \omega^1_{max} + \omega^2_{max}]$.
impl Add for Segment {
    type Output = Segment;

    fn add(self, other: Segment) -> Segment {
        Segment {
            min: self.min + other.min,
            max: self.max + other.max,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Segment;

    #[test]
    fn new() {
        let segment = Segment::new(-1.5, 2.0);
        assert_eq!(segment.min(), -1.5);
        assert_eq!(segment.max(), 2.0);

        // A single resonance is supported on a segment of zero length
        let point = Segment::new(2.0, 2.0);
        assert_eq!(point.min(), point.max());
    }

    #[test]
    #[should_panic(expected = "min <= max")]
    fn new_backwards() {
        let _ = Segment::new(2.0, -1.5);
    }

    #[test]
    #[should_panic(expected = "min <= max")]
    fn new_nan() {
        let _ = Segment::new(f64::NAN, 1.0);
    }

    #[test]
    #[should_panic(expected = "must sit at a finite frequency")]
    fn new_degenerate_at_infinity() {
        // Its length would be the difference of two infinities, which is no number
        let _ = Segment::new(f64::INFINITY, f64::INFINITY);
    }

    #[test]
    #[should_panic(expected = "must sit at a finite frequency")]
    fn new_degenerate_at_negative_infinity() {
        let _ = Segment::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    }

    #[test]
    fn contains() {
        let segment = Segment::new(-1.5, 2.0);
        for omega in [-1.5, -0.4, 0.0, 1.9, 2.0] {
            assert!(segment.contains(omega), "{omega} is within the segment");
        }
        for omega in [-1.6, 2.1, f64::NEG_INFINITY, f64::INFINITY, f64::NAN] {
            assert!(!segment.contains(omega), "{omega} is not");
        }

        // A segment of zero length holds the one frequency it sits at
        let point = Segment::new(2.0, 2.0);
        assert!(point.contains(2.0));
        assert!(!point.contains(2.0f64.next_up()));

        // An unbounded segment holds every finite frequency and the end it reaches
        let unbounded = Segment::new(f64::NEG_INFINITY, f64::INFINITY);
        assert!(unbounded.contains(1e300));
        assert!(unbounded.contains(f64::INFINITY));
        assert!(!unbounded.contains(f64::NAN));
    }

    #[test]
    fn strictly_contains() {
        let segment = Segment::new(-1.5, 2.0);
        for omega in [-1.4, 0.0, 1.9] {
            assert!(segment.strictly_contains(omega), "{omega} is inside");
        }
        for omega in [-1.5, 2.0, -1.6, 2.1, f64::NAN] {
            assert!(!segment.strictly_contains(omega), "{omega} is not");
        }

        // One ulp in from either end is already inside
        assert!(segment.strictly_contains((-1.5f64).next_up()));
        assert!(segment.strictly_contains(2.0f64.next_down()));

        // A segment of zero length has no inside at all
        assert!(!Segment::new(2.0, 2.0).strictly_contains(2.0));
    }

    #[test]
    fn is_bounded() {
        assert!(Segment::new(-1.5, 2.0).is_bounded());
        assert!(Segment::new(2.0, 2.0).is_bounded());
        assert!(!Segment::new(f64::NEG_INFINITY, 2.0).is_bounded());
        assert!(!Segment::new(-1.5, f64::INFINITY).is_bounded());
        assert!(!Segment::new(f64::NEG_INFINITY, f64::INFINITY).is_bounded());
    }

    #[test]
    fn is_degenerate() {
        assert!(Segment::new(2.0, 2.0).is_degenerate());
        assert!(!Segment::new(-1.5, 2.0).is_degenerate());

        // The shortest segment that is not degenerate spans one ulp
        assert!(!Segment::new(2.0, 2.0f64.next_up()).is_degenerate());

        // Zero of either sign is one frequency
        assert!(Segment::new(-0.0, 0.0).is_degenerate());
    }

    #[test]
    fn length() {
        assert_eq!(Segment::new(-1.5, 2.0).length(), 3.5);
        assert_eq!(Segment::new(2.0, 2.0).length(), 0.0);
        assert_eq!(Segment::new(-1.5, f64::INFINITY).length(), f64::INFINITY);
    }

    #[test]
    fn midpoint() {
        assert_eq!(Segment::new(-1.5, 2.0).midpoint(), 0.25);
        assert_eq!(Segment::new(2.0, 2.0).midpoint(), 2.0);

        // The ends are halved before they are added, leaving room for the widest
        // segment of finite frequencies
        let huge = Segment::new(-f64::MAX, f64::MAX);
        assert_eq!(huge.midpoint(), 0.0);
    }

    #[test]
    fn mirrored() {
        let segment = Segment::new(-1.5, 2.0);
        assert_eq!(segment.mirrored(0.0), Segment::new(-2.0, 1.5));
        assert_eq!(segment.mirrored(1.0), Segment::new(-1.0, 2.5));

        // The reflection is its own inverse, and turns the segment around without
        // stretching it
        assert_eq!(segment.mirrored(0.5).mirrored(0.5), segment);
        let mirrored = segment.mirrored(3.0);
        assert_eq!(mirrored.length(), segment.length());
        assert_eq!(mirrored.midpoint(), 3.0 - segment.midpoint());
    }

    #[test]
    fn split_at() {
        let segment = Segment::new(-1.5, 2.0);
        let (lower, upper) = segment.split_at(0.5);
        assert_eq!(lower, Segment::new(-1.5, 0.5));
        assert_eq!(upper, Segment::new(0.5, 2.0));
        assert_eq!(lower.length() + upper.length(), segment.length());

        // A split at an end leaves a degenerate part beside the whole segment
        let (lower, upper) = segment.split_at(-1.5);
        assert!(lower.is_degenerate());
        assert_eq!(upper, segment);
    }

    #[test]
    #[should_panic(expected = "within it")]
    fn split_at_outside() {
        let _ = Segment::new(-1.5, 2.0).split_at(2.5);
    }

    #[test]
    fn add() {
        assert_eq!(
            Segment::new(-1.5, 2.0) + Segment::new(0.5, 3.0),
            Segment::new(-1.0, 5.0)
        );

        // A segment of zero length displaces the other without widening it
        let point = Segment::new(2.0, 2.0);
        assert_eq!(point + point, Segment::new(4.0, 4.0));
        let segment = Segment::new(-1.5, 2.0);
        assert_eq!(segment + point, segment.shifted(2.0));

        // Lengths add, so the sum is as wide as the two together
        let (x, y) = (Segment::new(-1.5, 2.0), Segment::new(0.5, 3.0));
        assert_eq!((x + y).length(), x.length() + y.length());

        // One unbounded operand reaches as far in the sum
        let upper = Segment::new(0.0, f64::INFINITY);
        assert_eq!(segment + upper, Segment::new(-1.5, f64::INFINITY));
        let whole = Segment::new(f64::NEG_INFINITY, f64::INFINITY);
        assert_eq!(segment + whole, whole);
    }

    #[test]
    fn hull() {
        // Nothing to take the hull of
        assert_eq!(Segment::hull([]), None);

        let segment = Segment::new(-1.5, 2.0);
        assert_eq!(Segment::hull([segment]), Some(segment));

        // The hull reaches around segments that are disjoint, and is not stretched by
        // one that another already contains
        let spread = [
            Segment::new(0.0, 1.0),
            Segment::new(-1.5, -0.5),
            Segment::new(0.25, 0.75),
            Segment::new(3.0, 4.0),
        ];
        assert_eq!(Segment::hull(spread), Some(Segment::new(-1.5, 4.0)));

        // One unbounded segment makes the hull unbounded
        let mixed = [Segment::new(0.0, 1.0), Segment::new(-1.5, f64::INFINITY)];
        assert_eq!(
            Segment::hull(mixed),
            Some(Segment::new(-1.5, f64::INFINITY))
        );
    }

    #[test]
    fn intersection() {
        let segment = Segment::new(-1.5, 2.0);

        // Overlapping in part, and one lying inside the other
        assert_eq!(
            segment.intersection(Segment::new(0.0, 4.0)),
            Some(Segment::new(0.0, 2.0))
        );
        assert_eq!(
            segment.intersection(Segment::new(-0.5, 0.5)),
            Some(Segment::new(-0.5, 0.5))
        );
        assert_eq!(segment.intersection(segment), Some(segment));

        // Segments meeting at one frequency overlap in the degenerate segment there
        let touching = segment.intersection(Segment::new(2.0, 5.0)).unwrap();
        assert_eq!(touching, Segment::new(2.0, 2.0));
        assert!(touching.is_degenerate());

        // Disjoint segments do not overlap at all
        assert_eq!(segment.intersection(Segment::new(2.5, 5.0)), None);
        assert_eq!(Segment::new(2.5, 5.0).intersection(segment), None);

        // An unbounded segment cuts down to whatever it is met with
        let unbounded = Segment::new(f64::NEG_INFINITY, f64::INFINITY);
        assert_eq!(unbounded.intersection(segment), Some(segment));
        assert_eq!(segment.intersection(unbounded), Some(segment));
    }
}
