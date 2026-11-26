//! The algebra of one entry of a dense jet matrix, and the inline form of a
//! jet.
//!
//! The 4×4 backends work on small dense matrices of jets and need only
//! arithmetic and read access to the sampled values of each entry, not the
//! arrays that hold them. [`EntryAlgebra`] is that narrower interface. Every
//! [`ScalarAlgebra`] jet provides it, and so do the point jets of
//! [`super::point`], whose operations do not allocate.
//!
//! [`InlineJet`] names the form a backend computes in: a point jet for a
//! one-sample (`Ix0`) `Jet0`, `Jet1` or `Jet2`, the jet itself otherwise.

use std::fmt::Debug;

use ndarray::{Array, Ix0, Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn, arr0};

use crate::ComplexScalar;

use super::{
    ArrayJet0, ArrayJet1, ArrayJet2, ArrayJetBivariate1, ArrayJetBivariate2, ScalarAlgebra,
    point::{Point, PointJet0, PointJet1, PointJet2},
};

mod sealed {
    use std::fmt::Debug;

    use crate::ComplexScalar;

    /// Arithmetic on one entry of a dense jet matrix.
    ///
    /// The operations mean what the [`ScalarAlgebra`](crate::ScalarAlgebra)
    /// operations of the same names mean, and for a `ScalarAlgebra` jet they
    /// are those operations.
    pub trait EntryAlgebra: Clone + Debug {
        type Scalar: ComplexScalar;

        /// Zero with the sampled shape of this entry.
        fn zero_like(&self) -> Self;

        /// Unity with the sampled shape of this entry.
        fn one_like(&self) -> Self;

        /// A constant with the sampled shape of this entry.
        fn constant(&self, value: Self::Scalar) -> Self;

        fn add(&self, rhs: &Self) -> Self;

        fn subtract(&self, rhs: &Self) -> Self;

        fn negate(&self) -> Self;

        fn multiply(&self, rhs: &Self) -> Self;

        fn reciprocal(&self) -> Self;

        fn scale(&self, coefficient: Self::Scalar) -> Self;

        /// The principal square root.
        fn sqrt(&self) -> Self;

        /// Whether every value and derivative component is finite.
        fn all_finite(&self) -> bool;

        fn square(&self) -> Self {
            self.multiply(self)
        }

        fn divide(&self, rhs: &Self) -> Self {
            self.multiply(&rhs.reciprocal())
        }

        /// The number of samples.
        fn sample_count(&self) -> usize;

        /// The value of every sample, in the logical order of the samples.
        fn sample_values(&self) -> impl Iterator<Item = Self::Scalar> + '_;

        /// The constant whose samples are `values`, shaped like `like`;
        /// `None` when the count does not match.
        fn from_sample_values(like: &Self, values: Vec<Self::Scalar>) -> Option<Self>;
    }
}

pub(crate) use sealed::EntryAlgebra;

impl<A> EntryAlgebra for A
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
{
    type Scalar = A::Scalar;

    fn zero_like(&self) -> Self {
        ScalarAlgebra::zero_like(self)
    }

    fn one_like(&self) -> Self {
        ScalarAlgebra::one_like(self)
    }

    fn constant(&self, value: A::Scalar) -> Self {
        ScalarAlgebra::constant(self, value)
    }

    fn add(&self, rhs: &Self) -> Self {
        ScalarAlgebra::add(self, rhs)
    }

    fn subtract(&self, rhs: &Self) -> Self {
        ScalarAlgebra::subtract(self, rhs)
    }

    fn negate(&self) -> Self {
        ScalarAlgebra::negate(self)
    }

    fn multiply(&self, rhs: &Self) -> Self {
        ScalarAlgebra::multiply(self, rhs)
    }

    fn reciprocal(&self) -> Self {
        ScalarAlgebra::reciprocal(self)
    }

    fn scale(&self, coefficient: A::Scalar) -> Self {
        ScalarAlgebra::scale(self, coefficient)
    }

    fn sqrt(&self) -> Self {
        ScalarAlgebra::sqrt(self)
    }

    fn all_finite(&self) -> bool {
        ScalarAlgebra::all_finite(self)
    }

    fn square(&self) -> Self {
        ScalarAlgebra::square(self)
    }

    fn divide(&self, rhs: &Self) -> Self {
        ScalarAlgebra::divide(self, rhs)
    }

    fn sample_count(&self) -> usize {
        self.value().len()
    }

    fn sample_values(&self) -> impl Iterator<Item = A::Scalar> + '_ {
        self.value().iter().copied()
    }

    fn from_sample_values(like: &Self, values: Vec<A::Scalar>) -> Option<Self> {
        Array::from_shape_vec(like.value().raw_dim(), values)
            .ok()
            .map(A::lift_constant)
    }
}

/// [`EntryAlgebra`] for a point jet, whose payload operations are the
/// generic jet rules over [`Point`].
macro_rules! point_entry {
    ($jet:ident) => {
        impl<C, P> EntryAlgebra for $jet<C, P>
        where
            C: ComplexScalar,
            P: Clone + Debug,
        {
            type Scalar = C;

            fn zero_like(&self) -> Self {
                Self::constant_like(self.value(), C::zero())
            }

            fn one_like(&self) -> Self {
                Self::constant_like(self.value(), C::one())
            }

            fn constant(&self, value: C) -> Self {
                Self::constant_like(self.value(), value)
            }

            fn add(&self, rhs: &Self) -> Self {
                $jet::add(self, rhs)
            }

            fn subtract(&self, rhs: &Self) -> Self {
                $jet::subtract(self, rhs)
            }

            fn negate(&self) -> Self {
                $jet::negate(self)
            }

            fn multiply(&self, rhs: &Self) -> Self {
                $jet::multiply(self, rhs)
            }

            fn reciprocal(&self) -> Self {
                $jet::reciprocal(self)
            }

            fn scale(&self, coefficient: C) -> Self {
                $jet::scale_by(self, coefficient)
            }

            fn sqrt(&self) -> Self {
                $jet::sqrt(self)
            }

            fn all_finite(&self) -> bool {
                $jet::all_finite(self)
            }

            fn sample_count(&self) -> usize {
                1
            }

            fn sample_values(&self) -> impl Iterator<Item = C> + '_ {
                std::iter::once(self.value().0)
            }

            fn from_sample_values(_: &Self, values: Vec<C>) -> Option<Self> {
                match values[..] {
                    [value] => Some(Self::constant(Point(value))),
                    _ => None,
                }
            }
        }
    };
}

point_entry!(PointJet0);
point_entry!(PointJet1);
point_entry!(PointJet2);

/// The form of a jet that the 4×4 backends compute in.
///
/// One-sample (`Ix0`) `Jet0`, `Jet1` and `Jet2` jets map to point jets,
/// whose operations do not allocate; every other jet maps to itself. The
/// conversions copy values and derivatives exactly.
///
/// This trait is implemented for every jet family and sample dimension of
/// lamina-core and is not intended to be implemented downstream.
#[doc(hidden)]
pub trait InlineJet: ScalarAlgebra {
    #[doc(hidden)]
    type Inline: EntryAlgebra<Scalar = Self::Scalar>;

    #[doc(hidden)]
    fn to_inline(&self) -> Self::Inline;

    #[doc(hidden)]
    fn from_inline(inline: Self::Inline) -> Self;
}

impl<C: ComplexScalar, P: Clone + Debug> InlineJet for ArrayJet0<C, Ix0, P> {
    type Inline = PointJet0<C, P>;

    fn to_inline(&self) -> Self::Inline {
        PointJet0::new(Point(self.value()[()]))
    }

    fn from_inline(inline: Self::Inline) -> Self {
        Self::new(arr0(inline.value().0))
    }
}

impl<C: ComplexScalar, P: Clone + Debug> InlineJet for ArrayJet1<C, Ix0, P> {
    type Inline = PointJet1<C, P>;

    fn to_inline(&self) -> Self::Inline {
        PointJet1::from_parts(Point(self.value()[()]), Point(self.first()[()]))
    }

    fn from_inline(inline: Self::Inline) -> Self {
        Self::from_parts(arr0(inline.value().0), arr0(inline.first().0))
    }
}

impl<C: ComplexScalar, P: Clone + Debug> InlineJet for ArrayJet2<C, Ix0, P> {
    type Inline = PointJet2<C, P>;

    fn to_inline(&self) -> Self::Inline {
        PointJet2::from_parts(
            Point(self.value()[()]),
            Point(self.first()[()]),
            Point(self.second()[()]),
        )
    }

    fn from_inline(inline: Self::Inline) -> Self {
        Self::from_parts(
            arr0(inline.value().0),
            arr0(inline.first().0),
            arr0(inline.second().0),
        )
    }
}

/// Jets that compute in their own form.
macro_rules! inline_as_itself {
    ($jet:ident: $($dimension:ty),+) => {$(
        impl<C: ComplexScalar, P: Clone + Debug> InlineJet for $jet<C, $dimension, P> {
            type Inline = Self;

            fn to_inline(&self) -> Self {
                self.clone()
            }

            fn from_inline(inline: Self) -> Self {
                inline
            }
        }
    )+};
}

inline_as_itself!(ArrayJet0: Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);
inline_as_itself!(ArrayJet1: Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);
inline_as_itself!(ArrayJet2: Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);
inline_as_itself!(ArrayJetBivariate1: Ix0, Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);
inline_as_itself!(ArrayJetBivariate2: Ix0, Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);
