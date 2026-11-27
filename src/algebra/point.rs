//! Point jets: jets of one sample whose components are stored inline.
//!
//! The jet structs are generic over their payload. The array payload of the
//! public jets keeps every component on the heap, even for one sample, so
//! each jet operation allocates. [`Point`] holds one complex sample by value:
//! `Jet0`, `Jet1` and `Jet2` over it are `Copy`-sized values whose operations
//! never allocate.
//!
//! The payload implements the same capability traits as the arrays, so the
//! product, quotient and chain rules are the generic ones of the jet structs,
//! and every operation performs the floating-point operations of the array
//! jet at one sample in the same order: results agree bit for bit.

use nalgebra::ComplexField;

use super::{
    JetAdditive, JetBilinear, JetConjugate, JetConstant, JetOneLike, JetReciprocal, JetScaleBy,
    JetZeroLike, jet_one::Jet1, jet_two::Jet2, jet_zero::Jet0,
};

/// One complex sample, stored inline: the payload of point jets.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point<C>(pub(crate) C);

pub(crate) type PointJet0<C, P> = Jet0<Point<C>, P>;
pub(crate) type PointJet1<C, P> = Jet1<Point<C>, P>;
pub(crate) type PointJet2<C, P> = Jet2<Point<C>, P>;

impl<C: ComplexField + Copy> JetZeroLike for Point<C> {
    fn jet_zeros_like(_: &Self) -> Self {
        Self(C::zero())
    }
}

impl<C: ComplexField + Copy> JetOneLike for Point<C> {
    fn jet_ones_like(_: &Self) -> Self {
        Self(C::one())
    }
}

impl<C: ComplexField + Copy> JetAdditive for Point<C> {
    fn jet_add(&self, rhs: &Self) -> Self {
        Self(self.0 + rhs.0)
    }

    fn jet_subtract(&self, rhs: &Self) -> Self {
        Self(self.0 - rhs.0)
    }

    fn jet_negate(&self) -> Self {
        Self(-self.0)
    }
}

impl<C: ComplexField + Copy> JetBilinear for Point<C> {
    fn jet_multiply(&self, rhs: &Self) -> Self {
        Self(self.0 * rhs.0)
    }
}

impl<C: ComplexField + Copy> JetScaleBy for Point<C> {
    type Scalar = C;

    fn jet_scale_by(&self, value: C) -> Self {
        Self(self.0 * value)
    }
}

impl<C: ComplexField + Copy> JetConstant for Point<C> {
    type Scalar = C;

    fn jet_constant_like(&self, value: C) -> Self {
        Self(value)
    }
}

impl<C: ComplexField + Copy> JetConjugate for Point<C> {
    fn jet_conjugate(&self) -> Self {
        Self(self.0.conjugate())
    }
}

impl<C: ComplexField + Copy> JetReciprocal for Point<C> {
    fn jet_elementwise_reciprocal(&self) -> Self {
        Self(C::one() / self.0)
    }
}

fn finite<C: ComplexField + Copy>(point: &Point<C>) -> bool {
    point.0.real().is_finite() && point.0.imaginary().is_finite()
}

impl<C: ComplexField + Copy, P> PointJet0<C, P> {
    /// The principal square root, as `ArrayJet0::sqrt`.
    pub(crate) fn sqrt(&self) -> Self {
        Self::new(Point(self.value().0.sqrt()))
    }

    pub(crate) fn all_finite(&self) -> bool {
        finite(self.value())
    }
}

impl<C: ComplexField + Copy, P> PointJet1<C, P> {
    /// The principal square root and its derivative, as `ArrayJet1::sqrt`.
    pub(crate) fn sqrt(&self) -> Self {
        let result = self.value().0.sqrt();
        let two = C::one() + C::one();
        let derivative = self.first().0 / (two * result);

        Self::from_parts(Point(result), Point(derivative))
    }

    pub(crate) fn all_finite(&self) -> bool {
        finite(self.value()) && finite(self.first())
    }
}

impl<C: ComplexField + Copy, P> PointJet2<C, P> {
    /// The principal square root and its first two derivatives, as
    /// `ArrayJet2::sqrt`.
    pub(crate) fn sqrt(&self) -> Self {
        let y = self.value().0.sqrt();
        let two = C::one() + C::one();
        let four = two + two;
        let (x1, x2) = (self.first().0, self.second().0);

        let first = x1 / (two * y);
        let second = x2 / (two * y) - (x1 * x1) / (four * y * y * y);

        Self::from_parts(Point(y), Point(first), Point(second))
    }

    pub(crate) fn all_finite(&self) -> bool {
        finite(self.value()) && finite(self.first()) && finite(self.second())
    }
}
