//! Segments of the frequency axis.

/// Segment of the frequency axis, $[\omega_{min}, \omega_{max}]$.
///
/// A segment of zero length is a valid one, a discrete spectral function of a single
/// resonance being supported on exactly that.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    min: f64,
    max: f64,
}

impl Segment {
    /// Make a segment spanning `min` to `max`.
    pub fn new(min: f64, max: f64) -> Segment {
        assert!(min <= max, "segment bounds must satisfy min <= max");
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

    /// The same segment displaced in frequency by `by`.
    pub fn shifted(&self, by: f64) -> Segment {
        Segment {
            min: self.min + by,
            max: self.max + by,
        }
    }

    /// The parts of the segment below and above `omega`, which must lie within it.
    ///
    /// The two share the frequency they are split at, so a split at an end of the
    /// segment leaves a degenerate part rather than an empty one.
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
