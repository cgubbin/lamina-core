//! Lifting an oriented tensor material's stack-frame entries.
//!
//! The stack-frame tensor `R T(k₀) Rᵀ` of an oriented material is evaluated
//! and rotated once per Taylor coefficient and `k₀` sample; each needed entry
//! is read from the rotated tensor and composed with `k₀` by the chain rule.
//! `R` is real and constant, so it commutes with the derivatives.
//!
//! Each entry is the same element of the same rotated tensor that lifting the
//! entry on its own would compute, and the composition performs the same
//! floating-point operations as `compose_sampled_function`, so the results
//! agree bit for bit with the per-entry lift.

use std::fmt::Debug;

use nalgebra::ComplexField;
use ndarray::{Array, Dimension, Ix0, Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn};
use num_traits::Zero;

use crate::{
    ComplexScalar,
    algebra::{
        ArrayJet0, ArrayJet1, ArrayJet2, ArrayJetBivariate1, ArrayJetBivariate2,
        FirstOrderExpansion, InlineJet, Point, PointJet0, PointJet1, PointJet2,
        SecondOrderExpansion,
    },
    domain::{ComplexPlane, RealAxis},
    material::{
        DerivativeOrder, Rotation, Tensor3,
        tensor::{TensorHandle, TensorMaterialPoint, TensorMaterialValuePoint},
    },
};

/// The real scalar type of a tensor handle.
type Real<T> = <<T as TensorHandle>::Scalar as ComplexField>::RealField;

/// Which constitutive tensor to evaluate.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Quantity {
    Permittivity,
    Permeability,
}

const QUANTITIES: [Quantity; 2] = [Quantity::Permittivity, Quantity::Permeability];

/// Evaluates a tensor material's (crystal-frame) tensor at one `k₀` sample in
/// a spectral domain.
#[doc(hidden)]
pub trait StackFrameEvaluator<T: TensorHandle> {
    fn tensor(material: &T, quantity: Quantity, k0: T::Scalar) -> Tensor3<T::Scalar>;
}

/// [`StackFrameEvaluator`] with derivatives with respect to `k₀`.
#[doc(hidden)]
pub trait StackFrameDerivativeEvaluator<T: TensorHandle>: StackFrameEvaluator<T> {
    fn tensor_derivative(
        material: &T,
        quantity: Quantity,
        k0: T::Scalar,
        order: DerivativeOrder,
    ) -> Tensor3<T::Scalar>;
}

impl<T> StackFrameEvaluator<T> for RealAxis
where
    T: TensorHandle,
    Real<T>: Copy,
{
    fn tensor(material: &T, quantity: Quantity, k0: T::Scalar) -> Tensor3<T::Scalar> {
        let point = material.point();
        match quantity {
            Quantity::Permittivity => point.permittivity_at(k0.real()),
            Quantity::Permeability => point.permeability_at(k0.real()),
        }
    }
}

impl<T> StackFrameDerivativeEvaluator<T> for RealAxis
where
    T: TensorHandle,
    T::Point: TensorMaterialPoint<T::Scalar>,
    Real<T>: Copy,
{
    fn tensor_derivative(
        material: &T,
        quantity: Quantity,
        k0: T::Scalar,
        order: DerivativeOrder,
    ) -> Tensor3<T::Scalar> {
        debug_assert!(k0.imaginary() == T::Scalar::zero().real());
        let point = material.point();
        match quantity {
            Quantity::Permittivity => point.permittivity_derivative_at(k0.real(), order),
            Quantity::Permeability => point.permeability_derivative_at(k0.real(), order),
        }
    }
}

impl<T> StackFrameEvaluator<T> for ComplexPlane
where
    T: TensorHandle,
    T::Point: TensorMaterialPoint<T::Scalar>,
    Real<T>: Copy,
{
    fn tensor(material: &T, quantity: Quantity, k0: T::Scalar) -> Tensor3<T::Scalar> {
        let point = material.point();
        match quantity {
            Quantity::Permittivity => point.permittivity_complex_at(k0),
            Quantity::Permeability => point.permeability_complex_at(k0),
        }
    }
}

impl<T> StackFrameDerivativeEvaluator<T> for ComplexPlane
where
    T: TensorHandle,
    T::Point: TensorMaterialPoint<T::Scalar>,
    Real<T>: Copy,
{
    fn tensor_derivative(
        material: &T,
        quantity: Quantity,
        k0: T::Scalar,
        order: DerivativeOrder,
    ) -> Tensor3<T::Scalar> {
        let point = material.point();
        match quantity {
            Quantity::Permittivity => point.permittivity_complex_derivative_at(k0, order),
            Quantity::Permeability => point.permeability_complex_derivative_at(k0, order),
        }
    }
}

/// An oriented material, its rotation into the stack frame, and the
/// stack-frame entries a lift needs.
#[doc(hidden)]
pub struct StackFrame<'a, T: TensorHandle, const N: usize> {
    pub(crate) material: &'a T,
    pub(crate) rotation: &'a Rotation<Real<T>>,
    pub(crate) positions: [(usize, usize); N],
}

impl<T, const N: usize> StackFrame<'_, T, N>
where
    T: TensorHandle,
    Real<T>: Copy,
{
    /// The needed entries of one rotated coefficient at one sample.
    fn entries(&self, tensor: Tensor3<T::Scalar>) -> [T::Scalar; N] {
        let rotated = self.rotation.rotate(&tensor);
        self.positions.map(|(row, column)| rotated[row][column])
    }

    /// The needed entries of one rotated coefficient at every sample.
    fn sampled<D: Dimension>(
        &self,
        k0: &Array<T::Scalar, D>,
        coefficient: impl Fn(T::Scalar) -> Tensor3<T::Scalar>,
    ) -> [Array<T::Scalar, D>; N] {
        let rotated = k0.map(|&k0| self.rotation.rotate(&coefficient(k0)));
        self.positions
            .map(|(row, column)| rotated.map(|tensor| tensor[row][column]))
    }
}

/// The stack-frame ε and μ entries of an oriented tensor material, lifted
/// through a jet family in domain `E`.
///
/// Implemented for every jet family and sample dimension of lamina-core, under
/// the same domain capabilities as [`ConstitutiveLift`](super::ConstitutiveLift):
/// value-only jets need tensor values, every other jet their derivatives.
#[doc(hidden)]
pub trait StackFrameLift<E, T: TensorHandle>: Sized {
    /// The `[ε, μ]` entries at `frame.positions`.
    fn stack_frame_entries<const N: usize>(
        frame: &StackFrame<'_, T, N>,
        vacuum_angular_wavenumber: &Self,
    ) -> [[Self; N]; 2];

    /// [`Self::stack_frame_entries`] in the inline form of the jet.
    fn stack_frame_entries_inline<const N: usize>(
        frame: &StackFrame<'_, T, N>,
        vacuum_angular_wavenumber: &Self,
    ) -> [[Self::Inline; N]; 2]
    where
        Self: InlineJet,
    {
        Self::stack_frame_entries(frame, vacuum_angular_wavenumber)
            .map(|entries| entries.map(|entry| entry.to_inline()))
    }
}

fn array_jet0<C, D, P, E, T, const N: usize>(
    frame: &StackFrame<'_, T, N>,
    k0: &ArrayJet0<C, D, P>,
) -> [[ArrayJet0<C, D, P>; N]; 2]
where
    C: ComplexScalar + Copy,
    D: Dimension,
    P: Clone + Debug,
    T: TensorHandle<Scalar = C>,
    Real<T>: Copy,
    E: StackFrameEvaluator<T>,
{
    QUANTITIES.map(|quantity| {
        frame
            .sampled(k0.value(), |k| E::tensor(frame.material, quantity, k))
            .map(ArrayJet0::new)
    })
}

/// First-order array jets: value and first coefficient, composed with `k₀`.
macro_rules! array_first {
    ($name:ident, $jet:ident) => {
        fn $name<C, D, P, E, T, const N: usize>(
            frame: &StackFrame<'_, T, N>,
            k0: &$jet<C, D, P>,
        ) -> [[$jet<C, D, P>; N]; 2]
        where
            C: ComplexScalar + Copy,
            D: Dimension,
            P: Clone + Debug,
            T: TensorHandle<Scalar = C>,
            Real<T>: Copy,
            E: StackFrameDerivativeEvaluator<T>,
        {
            QUANTITIES.map(|quantity| {
                let values = frame.sampled(k0.value(), |k| E::tensor(frame.material, quantity, k));
                let mut firsts = frame
                    .sampled(k0.value(), |k| {
                        E::tensor_derivative(frame.material, quantity, k, DerivativeOrder::First)
                    })
                    .into_iter();
                values.map(|value| {
                    let first = firsts.next().expect("one coefficient per entry");
                    $jet::compose_sampled_function(k0, FirstOrderExpansion::new(value, first))
                })
            })
        }
    };
}

/// Second-order array jets: value, first and second coefficients, composed
/// with `k₀`.
macro_rules! array_second {
    ($name:ident, $jet:ident) => {
        fn $name<C, D, P, E, T, const N: usize>(
            frame: &StackFrame<'_, T, N>,
            k0: &$jet<C, D, P>,
        ) -> [[$jet<C, D, P>; N]; 2]
        where
            C: ComplexScalar + Copy,
            D: Dimension,
            P: Clone + Debug,
            T: TensorHandle<Scalar = C>,
            Real<T>: Copy,
            E: StackFrameDerivativeEvaluator<T>,
        {
            QUANTITIES.map(|quantity| {
                let values = frame.sampled(k0.value(), |k| E::tensor(frame.material, quantity, k));
                let derivative = |order| {
                    frame
                        .sampled(k0.value(), |k| {
                            E::tensor_derivative(frame.material, quantity, k, order)
                        })
                        .into_iter()
                };
                let mut firsts = derivative(DerivativeOrder::First);
                let mut seconds = derivative(DerivativeOrder::Second);
                values.map(|value| {
                    let first = firsts.next().expect("one coefficient per entry");
                    let second = seconds.next().expect("one coefficient per entry");
                    $jet::compose_sampled_function(
                        k0,
                        SecondOrderExpansion::new(value, first, second),
                    )
                })
            })
        }
    };
}

array_first!(array_jet1, ArrayJet1);
array_first!(array_bivariate1, ArrayJetBivariate1);
array_second!(array_jet2, ArrayJet2);
array_second!(array_bivariate2, ArrayJetBivariate2);

/// [`StackFrameLift`] for array jets of the given dimensions, lifting on the
/// arrays (the inline form, where it differs, is converted).
macro_rules! array_lift {
    ($jet:ident, $lift:ident, $evaluator:ident: $($dimension:ty),+) => {$(
        impl<C, P, E, T> StackFrameLift<E, T> for $jet<C, $dimension, P>
        where
            C: ComplexScalar + Copy,
            P: Clone + Debug,
            T: TensorHandle<Scalar = C>,
            Real<T>: Copy,
            E: $evaluator<T>,
        {
            fn stack_frame_entries<const N: usize>(
                frame: &StackFrame<'_, T, N>,
                k0: &Self,
            ) -> [[Self; N]; 2] {
                $lift::<C, $dimension, P, E, T, N>(frame, k0)
            }
        }
    )+};
}

array_lift!(ArrayJet0, array_jet0, StackFrameEvaluator: Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);
array_lift!(ArrayJet1, array_jet1, StackFrameDerivativeEvaluator: Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);
array_lift!(ArrayJet2, array_jet2, StackFrameDerivativeEvaluator: Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);
array_lift!(ArrayJetBivariate1, array_bivariate1, StackFrameDerivativeEvaluator: Ix0, Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);
array_lift!(ArrayJetBivariate2, array_bivariate2, StackFrameDerivativeEvaluator: Ix0, Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, IxDyn);

impl<C, P, E, T> StackFrameLift<E, T> for ArrayJet0<C, Ix0, P>
where
    C: ComplexScalar + Copy,
    P: Clone + Debug,
    T: TensorHandle<Scalar = C>,
    Real<T>: Copy,
    E: StackFrameEvaluator<T>,
{
    fn stack_frame_entries<const N: usize>(
        frame: &StackFrame<'_, T, N>,
        k0: &Self,
    ) -> [[Self; N]; 2] {
        array_jet0::<C, Ix0, P, E, T, N>(frame, k0)
    }

    fn stack_frame_entries_inline<const N: usize>(
        frame: &StackFrame<'_, T, N>,
        k0: &Self,
    ) -> [[<Self as InlineJet>::Inline; N]; 2] {
        let k = k0.value()[()];
        QUANTITIES.map(|quantity| {
            frame
                .entries(E::tensor(frame.material, quantity, k))
                .map(|value| PointJet0::new(Point(value)))
        })
    }
}

impl<C, P, E, T> StackFrameLift<E, T> for ArrayJet1<C, Ix0, P>
where
    C: ComplexScalar + Copy,
    P: Clone + Debug,
    T: TensorHandle<Scalar = C>,
    Real<T>: Copy,
    E: StackFrameDerivativeEvaluator<T>,
{
    fn stack_frame_entries<const N: usize>(
        frame: &StackFrame<'_, T, N>,
        k0: &Self,
    ) -> [[Self; N]; 2] {
        array_jet1::<C, Ix0, P, E, T, N>(frame, k0)
    }

    /// `f′ k₀′`, as `ArrayJet1::compose_sampled_function` at one sample.
    fn stack_frame_entries_inline<const N: usize>(
        frame: &StackFrame<'_, T, N>,
        k0: &Self,
    ) -> [[<Self as InlineJet>::Inline; N]; 2] {
        let (k, k_first) = (k0.value()[()], k0.first()[()]);
        QUANTITIES.map(|quantity| {
            let values = frame.entries(E::tensor(frame.material, quantity, k));
            let firsts = frame.entries(E::tensor_derivative(
                frame.material,
                quantity,
                k,
                DerivativeOrder::First,
            ));
            std::array::from_fn(|index| {
                PointJet1::from_parts(Point(values[index]), Point(firsts[index] * k_first))
            })
        })
    }
}

impl<C, P, E, T> StackFrameLift<E, T> for ArrayJet2<C, Ix0, P>
where
    C: ComplexScalar + Copy,
    P: Clone + Debug,
    T: TensorHandle<Scalar = C>,
    Real<T>: Copy,
    E: StackFrameDerivativeEvaluator<T>,
{
    fn stack_frame_entries<const N: usize>(
        frame: &StackFrame<'_, T, N>,
        k0: &Self,
    ) -> [[Self; N]; 2] {
        array_jet2::<C, Ix0, P, E, T, N>(frame, k0)
    }

    /// `f′ k₀′` and `f″ k₀′² + f′ k₀″`, as
    /// `ArrayJet2::compose_sampled_function` at one sample.
    fn stack_frame_entries_inline<const N: usize>(
        frame: &StackFrame<'_, T, N>,
        k0: &Self,
    ) -> [[<Self as InlineJet>::Inline; N]; 2] {
        let (k, k_first, k_second) = (k0.value()[()], k0.first()[()], k0.second()[()]);
        let k_first_squared = k_first * k_first;
        QUANTITIES.map(|quantity| {
            let values = frame.entries(E::tensor(frame.material, quantity, k));
            let derivative =
                |order| frame.entries(E::tensor_derivative(frame.material, quantity, k, order));
            let firsts = derivative(DerivativeOrder::First);
            let seconds = derivative(DerivativeOrder::Second);
            std::array::from_fn(|index| {
                PointJet2::from_parts(
                    Point(values[index]),
                    Point(firsts[index] * k_first),
                    Point(seconds[index] * k_first_squared + firsts[index] * k_second),
                )
            })
        })
    }
}
