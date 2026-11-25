//! Constitutive regions: isotropic, diagonal or full ε and μ in the stack
//! frame.

use nalgebra::ComplexField;
use ndarray::Dimension;
use num_traits::Float;
use thiserror::Error;

use super::{
    AnalyticalMaterialHandle, AnalyticalTensorMaterialHandle, ConstitutiveEvaluator,
    ConstitutiveLift, DerivativeOrder, EvaluateDifferentiableMaterial,
    EvaluateDifferentiableMeromorphicMaterial, EvaluateMaterial, EvaluateMeromorphicMaterial,
    Rotation, Sampled, Tensor3,
};
use crate::{ComplexScalar, algebra::ScalarAlgebra, stack::StackRegion};

/// Structural class of a region's stack-frame ε and μ.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConstitutiveClass {
    /// Scalar ε and μ.
    Isotropic,
    /// ε and μ diagonal in the stack frame.
    Diagonal,
    /// Any 3×3 ε and μ.
    Full,
}

/// The constitutive description of one stack region: a scalar material, or
/// a tensor material with its crystal-to-stack [`Rotation`].
///
/// The class is decided when the region is built, from the material's
/// [`PRINCIPAL_AXES`](super::TensorMaterial::PRINCIPAL_AXES) declaration and
/// the rotation's exact entries, never from tensor values:
///
/// | material | rotation | class | separable |
/// | --- | --- | --- | --- |
/// | scalar | none | isotropic | yes |
/// | principal axes | signed axis permutation | diagonal | yes |
/// | principal axes | a crystal axis exactly on `y` | full | yes |
/// | principal axes | other | full | no |
/// | general | any | full | no |
///
/// A region is *separable* when its stack-frame ε and μ have no `xy`, `yx`,
/// `yz` or `zy` component, so that TE and TM decouple in it. The angle
/// constructors of [`Rotation`] give exact entries at multiples of `π/2`, so
/// an optic axis along a stack axis is diagonal and one in the `xz` plane is
/// separable.
#[derive(Clone, Debug, PartialEq)]
pub struct Constitutive<S, T, R> {
    kind: Kind<S, T, R>,
}

#[derive(Clone, Debug, PartialEq)]
enum Kind<S, T, R> {
    Isotropic(S),
    Oriented {
        material: T,
        rotation: Rotation<R>,
        diagonal: bool,
        separable: bool,
    },
}

/// A region of materials with a complex continuation and derivatives.
pub type AnalyticalMaterialRegion<C> = Constitutive<
    AnalyticalMaterialHandle<C>,
    AnalyticalTensorMaterialHandle<C>,
    <C as ComplexField>::RealField,
>;

impl<S, T, R> From<S> for Constitutive<S, T, R> {
    /// An isotropic region of a scalar material.
    fn from(material: S) -> Self {
        Self {
            kind: Kind::Isotropic(material),
        }
    }
}

impl<S, C> Constitutive<S, AnalyticalTensorMaterialHandle<C>, C::RealField>
where
    C: ComplexScalar,
    C::RealField: Float,
{
    /// A tensor material oriented in the stack frame by `rotation`.
    pub fn oriented(
        material: AnalyticalTensorMaterialHandle<C>,
        rotation: Rotation<C::RealField>,
    ) -> Self {
        let principal = material.principal_axes();
        Self {
            kind: Kind::Oriented {
                diagonal: principal && rotation.is_signed_permutation(),
                separable: principal && rotation.has_crystal_axis_on_y(),
                material,
                rotation,
            },
        }
    }
}

impl<S, T, R> Constitutive<S, T, R> {
    /// The structural class.
    pub fn class(&self) -> ConstitutiveClass {
        match self.kind {
            Kind::Isotropic(_) => ConstitutiveClass::Isotropic,
            Kind::Oriented { diagonal: true, .. } => ConstitutiveClass::Diagonal,
            Kind::Oriented { .. } => ConstitutiveClass::Full,
        }
    }

    /// `true` when the class guarantees ε and μ have no `xy`, `yx`, `yz` or
    /// `zy` component.
    pub fn is_separable(&self) -> bool {
        match self.kind {
            Kind::Isotropic(_) => true,
            Kind::Oriented { separable, .. } => separable,
        }
    }

    /// The scalar material of an isotropic region.
    pub fn scalar(&self) -> Option<&S> {
        match &self.kind {
            Kind::Isotropic(material) => Some(material),
            Kind::Oriented { .. } => None,
        }
    }

    /// The tensor material and its recorded rotation, for a diagonal or full
    /// region.
    pub fn tensor(&self) -> Option<(&T, &Rotation<R>)> {
        match &self.kind {
            Kind::Isotropic(_) => None,
            Kind::Oriented {
                material, rotation, ..
            } => Some((material, rotation)),
        }
    }
}

/// Stack-frame ε or μ of one region, in the region's class.
///
/// Isotropic regions carry the scalar, diagonal regions the stack-frame
/// diagonal `(xx, yy, zz)`, full regions all nine entries, row-major.
#[derive(Clone, Debug, PartialEq)]
pub enum ConstitutiveTensor<J> {
    Isotropic(J),
    Diagonal([J; 3]),
    Full([[J; 3]; 3]),
}

/// A region's stack-frame tensors are unusable at an evaluated point.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum ConstitutiveError {
    #[error("{region:?}: a stack-frame permittivity or permeability entry is not finite")]
    NonFinite { region: StackRegion },

    #[error(
        "{region:?}: the normal permittivity ε_zz or permeability μ_zz vanishes; \
         the normal fields are eliminated through it"
    )]
    NormalComponentVanishes { region: StackRegion },
}

/// A region whose class the formulation does not support.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("{region:?}: a {class:?} constitutive class is not supported by this formulation")]
pub struct UnsupportedConstitutive {
    pub region: StackRegion,
    pub class: ConstitutiveClass,
}

/// One stack-frame entry `(R T Rᵀ)_{row, column}` of an oriented tensor
/// material, seen as a scalar material whose permittivity is `ε_{row,column}`
/// and whose permeability is `μ_{row,column}`.
///
/// Lifting an entry through the scalar jet algebra gives its values and
/// spectral derivatives: `R` is real and constant, so it commutes with both.
#[doc(hidden)]
pub struct StackFrameEntry<C: ComplexScalar> {
    material: AnalyticalTensorMaterialHandle<C>,
    rotation: Rotation<C::RealField>,
    row: usize,
    column: usize,
}

impl<C> StackFrameEntry<C>
where
    C: ComplexScalar,
    C::RealField: Copy,
{
    fn rotate(&self, tensor: Tensor3<C>) -> C {
        let r = self.rotation.matrix();
        let mut entry = C::zero();
        for (k, tensor_row) in tensor.iter().enumerate() {
            for (l, &value) in tensor_row.iter().enumerate() {
                entry += value * C::from_real(r[self.row][k]) * C::from_real(r[self.column][l]);
            }
        }
        entry
    }
}

/// The scalar evaluation traits for [`StackFrameEntry`], each method rotating
/// the matching tensor evaluation.
macro_rules! rotated_entry {
    ($trait:ident, $k0:ty { $($item:item)* } $($method:ident => $point:ident $(($order:ident))?;)+) => {
        impl<C> $trait<C> for StackFrameEntry<C>
        where
            C: ComplexScalar,
            C::RealField: Copy,
        {
            $($item)*
            $(
                fn $method<I>(&self, k0: I $(, $order: DerivativeOrder)?) -> I::Mapped<C>
                where
                    I: Sampled<Elem = $k0>,
                {
                    k0.map(|k0| self.rotate(self.material.point().$point(k0 $(, $order)?)))
                }
            )+
        }
    };
}

rotated_entry!(EvaluateMaterial, C::RealField { type Real = C::RealField; }
    evaluate_relative_permittivity => permittivity_at;
    evaluate_relative_permeability => permeability_at;
);
rotated_entry!(EvaluateDifferentiableMaterial, C::RealField {}
    evaluate_relative_permittivity_derivative => permittivity_derivative_at(order);
    evaluate_relative_permeability_derivative => permeability_derivative_at(order);
);
rotated_entry!(EvaluateMeromorphicMaterial, C {}
    evaluate_relative_permittivity_complex => permittivity_complex_at;
    evaluate_relative_permeability_complex => permeability_complex_at;
);
rotated_entry!(EvaluateDifferentiableMeromorphicMaterial, C {}
    evaluate_relative_permittivity_complex_derivative => permittivity_complex_derivative_at(order);
    evaluate_relative_permeability_complex_derivative => permeability_complex_derivative_at(order);
);

/// Class-preserving evaluation of a region's stack-frame ε and μ through the
/// jet algebra.
///
/// Isotropic regions evaluate exactly as through [`ConstitutiveLift`].
/// Diagonal and full regions are validated at every sampled point: every
/// entry finite, and `ε_zz ≠ 0`, `μ_zz ≠ 0`.
#[doc(hidden)]
pub trait TensorConstitutiveLift<E, M>: ScalarAlgebra {
    /// Return `[ε, μ]` in the stack frame; `region` names the region in an
    /// error.
    fn constitutive_tensors(
        region: StackRegion,
        material: &M,
        vacuum_angular_wavenumber: &Self,
    ) -> Result<[ConstitutiveTensor<Self>; 2], ConstitutiveError>;
}

impl<J, E, S, C>
    TensorConstitutiveLift<E, Constitutive<S, AnalyticalTensorMaterialHandle<C>, C::RealField>>
    for J
where
    J: ConstitutiveLift<E, S> + ConstitutiveLift<E, StackFrameEntry<C>>,
    J: ScalarAlgebra<Scalar = C>,
    J::Dimension: Dimension,
    C: ComplexScalar,
    C::RealField: Copy,
    E: ConstitutiveEvaluator<C, J::Dimension, S>
        + ConstitutiveEvaluator<C, J::Dimension, StackFrameEntry<C>>,
{
    fn constitutive_tensors(
        region: StackRegion,
        material: &Constitutive<S, AnalyticalTensorMaterialHandle<C>, C::RealField>,
        k0: &Self,
    ) -> Result<[ConstitutiveTensor<Self>; 2], ConstitutiveError> {
        let (tensor, rotation, diagonal) = match &material.kind {
            Kind::Isotropic(scalar) => {
                return Ok([
                    ConstitutiveTensor::Isotropic(
                        <J as ConstitutiveLift<E, S>>::relative_permittivity(scalar, k0),
                    ),
                    ConstitutiveTensor::Isotropic(
                        <J as ConstitutiveLift<E, S>>::relative_permeability(scalar, k0),
                    ),
                ]);
            }
            Kind::Oriented {
                material,
                rotation,
                diagonal,
                ..
            } => (material, rotation, *diagonal),
        };

        let entry = |row: usize, column: usize| StackFrameEntry {
            material: tensor.clone(),
            rotation: *rotation,
            row,
            column,
        };
        let epsilon = |entry: &StackFrameEntry<C>| {
            <J as ConstitutiveLift<E, _>>::relative_permittivity(entry, k0)
        };
        let mu = |entry: &StackFrameEntry<C>| {
            <J as ConstitutiveLift<E, _>>::relative_permeability(entry, k0)
        };

        // Every entry finite; the last, row-major, is the normal component.
        let validate = |entries: &[&J]| {
            if !entries
                .iter()
                .all(|entry| entry.value().iter().all(|value| value.is_finite()))
            {
                return Err(ConstitutiveError::NonFinite { region });
            }
            match entries.last() {
                Some(zz) if zz.value().iter().any(|value| *value == C::zero()) => {
                    Err(ConstitutiveError::NormalComponentVanishes { region })
                }
                _ => Ok(()),
            }
        };

        if diagonal {
            let cells = [0, 1, 2].map(|axis| entry(axis, axis));
            let [epsilon, mu] = [cells.each_ref().map(epsilon), cells.each_ref().map(mu)];
            validate(&epsilon.each_ref())?;
            validate(&mu.each_ref())?;
            Ok([
                ConstitutiveTensor::Diagonal(epsilon),
                ConstitutiveTensor::Diagonal(mu),
            ])
        } else {
            let cells = [0, 1, 2].map(|row| [0, 1, 2].map(|column| entry(row, column)));
            let epsilon = cells.each_ref().map(|row| row.each_ref().map(epsilon));
            let mu = cells.each_ref().map(|row| row.each_ref().map(mu));
            validate(&epsilon.iter().flatten().collect::<Vec<_>>())?;
            validate(&mu.iter().flatten().collect::<Vec<_>>())?;
            Ok([
                ConstitutiveTensor::Full(epsilon),
                ConstitutiveTensor::Full(mu),
            ])
        }
    }
}

#[cfg(test)]
mod tests {
    use ndarray::arr0;
    use num_complex::Complex64;

    use num_traits::Zero;

    use super::*;
    use crate::{
        ComplexPlane, DifferentiableMaterial, DifferentiableMeromorphicMaterial,
        DifferentiableMeromorphicTensorMaterial, DifferentiableTensorMaterial, FiniteLayerIndex,
        Material, MeromorphicMaterial, MeromorphicTensorMaterial, RealAxis, Scalar, SignedAxis,
        TensorMaterial,
        algebra::{ArrayJet0, ArrayJet1, ArrayJet2},
        test_support::{
            jet::{H, HoloJ0, HoloJ1, HoloJ2, P, RealJ0, RealJ1, RealJ2},
            materials::{QuadraticDispersion, quadratic},
        },
    };

    type C = Complex64;
    type Region = AnalyticalMaterialRegion<C>;

    /// diag(o, o, e) for ε and μ, from two scalar models; `PRINCIPAL`
    /// selects the structural declaration.
    struct Uniaxial<const PRINCIPAL: bool> {
        ordinary: QuadraticDispersion,
        extraordinary: QuadraticDispersion,
    }

    fn diagonal<X: Copy + Zero>(ordinary: X, extraordinary: X) -> Tensor3<X> {
        let zero = X::zero();
        [
            [ordinary, zero, zero],
            [zero, ordinary, zero],
            [zero, zero, extraordinary],
        ]
    }

    macro_rules! axes {
        ($self:ident, $method:ident, $k:ident $(, $order:ident)?) => {
            $k.map(|k| diagonal(
                $self.ordinary.$method(Scalar(k) $(, $order)?),
                $self.extraordinary.$method(Scalar(k) $(, $order)?),
            ))
        };
    }

    impl<const PRINCIPAL: bool> TensorMaterial for Uniaxial<PRINCIPAL> {
        type Real = f64;
        const PRINCIPAL_AXES: bool = PRINCIPAL;

        fn relative_permittivity_tensor<I, X>(&self, k: I) -> I::Mapped<Tensor3<X>>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = f64>,
        {
            axes!(self, relative_permittivity, k)
        }

        fn relative_permeability_tensor<I, X>(&self, k: I) -> I::Mapped<Tensor3<X>>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = f64>,
        {
            axes!(self, relative_permeability, k)
        }
    }

    impl<const PRINCIPAL: bool> DifferentiableTensorMaterial for Uniaxial<PRINCIPAL> {
        fn relative_permittivity_tensor_derivative<I, X>(
            &self,
            k: I,
            order: DerivativeOrder,
        ) -> I::Mapped<Tensor3<X>>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = f64>,
        {
            axes!(self, relative_permittivity_derivative, k, order)
        }

        fn relative_permeability_tensor_derivative<I, X>(
            &self,
            k: I,
            order: DerivativeOrder,
        ) -> I::Mapped<Tensor3<X>>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = f64>,
        {
            axes!(self, relative_permeability_derivative, k, order)
        }
    }

    impl<const PRINCIPAL: bool> MeromorphicTensorMaterial for Uniaxial<PRINCIPAL> {
        fn relative_permittivity_tensor_complex<I, X>(&self, k: I) -> I::Mapped<Tensor3<X>>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = X>,
        {
            axes!(self, relative_permittivity_complex, k)
        }

        fn relative_permeability_tensor_complex<I, X>(&self, k: I) -> I::Mapped<Tensor3<X>>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = X>,
        {
            axes!(self, relative_permeability_complex, k)
        }
    }

    impl<const PRINCIPAL: bool> DifferentiableMeromorphicTensorMaterial for Uniaxial<PRINCIPAL> {
        fn relative_permittivity_tensor_complex_derivative<I, X>(
            &self,
            k: I,
            order: DerivativeOrder,
        ) -> I::Mapped<Tensor3<X>>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = X>,
        {
            axes!(self, relative_permittivity_complex_derivative, k, order)
        }

        fn relative_permeability_tensor_complex_derivative<I, X>(
            &self,
            k: I,
            order: DerivativeOrder,
        ) -> I::Mapped<Tensor3<X>>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = X>,
        {
            axes!(self, relative_permeability_complex_derivative, k, order)
        }
    }

    /// Lossy, dispersive and magnetic on both axes.
    fn uniaxial<const PRINCIPAL: bool>() -> Uniaxial<PRINCIPAL> {
        Uniaxial {
            ordinary: quadratic(2.2, 0.3, -0.05, 1.1, 0.02, 0.01),
            extraordinary: quadratic(4.6, -0.4, 0.07, 0.9, -0.03, 0.02),
        }
    }

    fn oriented<const PRINCIPAL: bool>(rotation: Rotation<f64>) -> Region {
        Region::oriented(
            AnalyticalTensorMaterialHandle::new(uniaxial::<PRINCIPAL>()),
            rotation,
        )
    }

    fn tensors<J, E>(region: &Region, k0: &J) -> [ConstitutiveTensor<J>; 2]
    where
        J: TensorConstitutiveLift<E, Region>,
    {
        J::constitutive_tensors(StackRegion::LeftExterior, region, k0).unwrap()
    }

    fn full<J: Clone>(tensor: &ConstitutiveTensor<J>) -> [[J; 3]; 3] {
        match tensor {
            ConstitutiveTensor::Full(full) => full.clone(),
            other => panic!(
                "expected a full tensor, got {:?}",
                std::mem::discriminant(other)
            ),
        }
    }

    #[test]
    fn class_and_separability_follow_the_structural_table() {
        let isotropic = Region::from(AnalyticalMaterialHandle::new(quadratic(
            2.0, 0.1, 0.0, 1.0, 0.0, 0.0,
        )));
        assert_eq!(isotropic.class(), ConstitutiveClass::Isotropic);
        assert!(isotropic.is_separable() && isotropic.scalar().is_some());
        assert!(isotropic.tensor().is_none());

        let permutation =
            Rotation::signed_permutation(SignedAxis::PlusZ, SignedAxis::MinusY, SignedAxis::PlusX)
                .unwrap();
        let rows = [
            (true, permutation, ConstitutiveClass::Diagonal, true),
            (
                true,
                Rotation::from_euler(0.0, 0.8, 0.0),
                ConstitutiveClass::Full,
                true,
            ),
            (
                true,
                Rotation::from_optic_axis(0.8, 0.0),
                ConstitutiveClass::Full,
                true,
            ),
            (
                true,
                Rotation::from_euler(0.3, 0.8, 0.0),
                ConstitutiveClass::Full,
                false,
            ),
            (false, Rotation::identity(), ConstitutiveClass::Full, false),
        ];

        for (index, (principal, rotation, class, separable)) in rows.into_iter().enumerate() {
            let region = if principal {
                oriented::<true>(rotation)
            } else {
                oriented::<false>(rotation)
            };
            assert_eq!(region.class(), class, "row {index}");
            assert_eq!(region.is_separable(), separable, "row {index}");
            assert!(region.scalar().is_none());
            assert_eq!(region.tensor().unwrap().1, &rotation, "row {index}");
        }
    }

    #[test]
    fn optic_axis_along_a_stack_axis_by_angles_is_diagonal() {
        use std::f64::consts::FRAC_PI_2;

        let k0 = RealJ0::new(arr0(C::from(1.3)));
        let model = uniaxial::<true>();
        let o: C = model.ordinary.relative_permittivity(Scalar(1.3));
        let e: C = model.extraordinary.relative_permittivity(Scalar(1.3));

        for (theta, phi, expected) in [
            (FRAC_PI_2, 0.0, [e, o, o]),
            (FRAC_PI_2, FRAC_PI_2, [o, e, o]),
            (0.0, 0.0, [o, o, e]),
        ] {
            let region = oriented::<true>(Rotation::from_optic_axis(theta, phi));
            assert_eq!(region.class(), ConstitutiveClass::Diagonal);
            let [epsilon, _] = tensors::<_, RealAxis>(&region, &k0);
            let ConstitutiveTensor::Diagonal(diagonal) = epsilon else {
                panic!("expected a diagonal tensor")
            };
            assert_eq!(diagonal.map(|entry| entry.value()[()]), expected);
        }
    }

    #[test]
    fn stack_frame_tensor_is_r_t_r_transpose() {
        // For diag(o, o, e) with c mapped onto â: T = o I + (e − o) â âᵀ.
        let k = 0.9;
        let k0 = RealJ0::new(arr0(C::from(k)));
        let model = uniaxial::<true>();

        for rotation in [
            Rotation::from_optic_axis(0.7, 0.0),
            Rotation::from_optic_axis(1.1, 2.3),
            Rotation::from_euler(0.4, -1.3, 2.0),
        ] {
            let m = rotation.matrix();
            let axis = [m[0][2], m[1][2], m[2][2]];
            let [epsilon, mu] = tensors::<_, RealAxis>(&oriented::<true>(rotation), &k0);

            for (tensor, o, e) in [
                (
                    full(&epsilon),
                    model.ordinary.relative_permittivity::<_, C>(Scalar(k)),
                    model.extraordinary.relative_permittivity::<_, C>(Scalar(k)),
                ),
                (
                    full(&mu),
                    model.ordinary.relative_permeability::<_, C>(Scalar(k)),
                    model.extraordinary.relative_permeability::<_, C>(Scalar(k)),
                ),
            ] {
                for i in 0..3 {
                    for j in 0..3 {
                        let delta = if i == j { o } else { C::zero() };
                        let expected = delta + (e - o) * axis[i] * axis[j];
                        let actual = tensor[i][j].value()[()];
                        assert!((actual - expected).norm() < 1e-14, "({i}, {j})");
                    }
                }
            }
        }
    }

    #[test]
    fn optic_axis_in_the_xz_plane_couples_only_x_and_z() {
        let (theta, k) = (0.7_f64, 0.9);
        let model = uniaxial::<true>();
        let o: C = model.ordinary.relative_permittivity(Scalar(k));
        let e: C = model.extraordinary.relative_permittivity(Scalar(k));

        let region = oriented::<true>(Rotation::from_optic_axis(theta, 0.0));
        assert_eq!(region.class(), ConstitutiveClass::Full);
        assert!(region.is_separable());
        let [epsilon, _] = tensors::<_, RealAxis>(&region, &RealJ0::new(arr0(C::from(k))));
        let epsilon = full(&epsilon).map(|row| row.map(|entry| entry.value()[()]));

        let coupling = (e - o) * theta.sin() * theta.cos();
        assert!((epsilon[0][2] - coupling).norm() < 1e-14);
        assert!((epsilon[2][0] - coupling).norm() < 1e-14);
        for (i, j) in [(0, 1), (1, 0), (1, 2), (2, 1)] {
            assert_eq!(epsilon[i][j], C::zero(), "({i}, {j})");
        }
    }

    #[test]
    fn isotropic_region_lifts_bit_identically_to_its_scalar_material() {
        let handle = AnalyticalMaterialHandle::new(quadratic(2.2, 0.3, -0.05, 1.1, 0.02, 0.01));
        let region = Region::from(handle.clone());

        macro_rules! check {
            ($jet:ty, $domain:ty, $k0:expr) => {{
                let k0 = $k0;
                let [epsilon, mu] = tensors::<$jet, $domain>(&region, &k0);
                let scalar_epsilon =
                    <$jet as ConstitutiveLift<$domain, _>>::relative_permittivity(&handle, &k0);
                let scalar_mu =
                    <$jet as ConstitutiveLift<$domain, _>>::relative_permeability(&handle, &k0);
                assert_eq!(epsilon, ConstitutiveTensor::Isotropic(scalar_epsilon));
                assert_eq!(mu, ConstitutiveTensor::Isotropic(scalar_mu));
            }};
        }

        let real = C::from(1.3);
        let complex = C::new(1.3, -0.2);
        let first = |k: C| ArrayJet1::<C, _, P>::from_parts(arr0(k), arr0(C::from(1.0)));
        let second =
            |k: C| ArrayJet2::<C, _, P>::from_parts(arr0(k), arr0(C::from(1.0)), arr0(C::zero()));
        check!(RealJ1, RealAxis, first(real));
        check!(RealJ2, RealAxis, second(real));

        let first = |k: C| ArrayJet1::<C, _, H>::from_parts(arr0(k), arr0(C::from(1.0)));
        let second =
            |k: C| ArrayJet2::<C, _, H>::from_parts(arr0(k), arr0(C::from(1.0)), arr0(C::zero()));
        check!(HoloJ1, ComplexPlane, first(complex));
        check!(HoloJ2, ComplexPlane, second(complex));
    }

    #[test]
    fn rotated_dispersive_tensor_derivatives_match_central_differences() {
        let region = oriented::<true>(Rotation::from_euler(0.4, -1.3, 2.0));
        let h = 1e-4;

        for k in [C::from(1.3), C::new(1.3, -0.2)] {
            let at = |k: C| -> [[[C; 3]; 3]; 2] {
                let k0 = ArrayJet0::<C, _, H>::new(arr0(k));
                tensors::<HoloJ0, ComplexPlane>(&region, &k0)
                    .map(|tensor| full(&tensor).map(|row| row.map(|entry| entry.value()[()])))
            };
            let (minus, centre, plus) = (at(k - h), at(k), at(k + h));

            let k0 = ArrayJet2::<C, _, H>::from_parts(arr0(k), arr0(C::from(1.0)), arr0(C::zero()));
            let lifted = tensors::<HoloJ2, ComplexPlane>(&region, &k0).map(|tensor| full(&tensor));

            for quantity in 0..2 {
                for i in 0..3 {
                    for j in 0..3 {
                        let jet = &lifted[quantity][i][j];
                        let first = (plus[quantity][i][j] - minus[quantity][i][j]) / (2.0 * h);
                        let second = (plus[quantity][i][j] - 2.0 * centre[quantity][i][j]
                            + minus[quantity][i][j])
                            / (h * h);
                        assert_eq!(jet.value()[()], centre[quantity][i][j]);
                        assert!((jet.first()[()] - first).norm() < 1e-9, "first ({i}, {j})");
                        assert!(
                            (jet.second()[()] - second).norm() < 1e-6,
                            "second ({i}, {j})"
                        );
                    }
                }
            }
        }

        // The real axis agrees with the complex plane there.
        let k0 = ArrayJet2::<C, _, P>::from_parts(
            arr0(C::from(1.3)),
            arr0(C::from(1.0)),
            arr0(C::zero()),
        );
        let real = tensors::<RealJ2, RealAxis>(&region, &k0).map(|tensor| full(&tensor));
        let k0 = ArrayJet2::<C, _, H>::from_parts(
            arr0(C::from(1.3)),
            arr0(C::from(1.0)),
            arr0(C::zero()),
        );
        let complex = tensors::<HoloJ2, ComplexPlane>(&region, &k0).map(|tensor| full(&tensor));
        for quantity in 0..2 {
            for i in 0..3 {
                for j in 0..3 {
                    let (r, c) = (&real[quantity][i][j], &complex[quantity][i][j]);
                    assert_eq!(r.value()[()], c.value()[()]);
                    assert_eq!(r.first()[()], c.first()[()]);
                    assert_eq!(r.second()[()], c.second()[()]);
                }
            }
        }
    }

    #[test]
    fn vanishing_normal_component_and_non_finite_entries_are_refused() {
        // ε_e(k) = −1 + k² vanishes at k = 1; with c along z it is ε_zz.
        let enz = Uniaxial::<true> {
            ordinary: quadratic(2.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            extraordinary: quadratic(-1.0, 0.0, 1.0, 1.0, 0.0, 0.0),
        };
        let region = Region::oriented(
            AnalyticalTensorMaterialHandle::new(enz),
            Rotation::identity(),
        );
        let layer = StackRegion::Layer(FiniteLayerIndex::new(2));

        let at = |k: f64| {
            <RealJ0 as TensorConstitutiveLift<RealAxis, _>>::constitutive_tensors(
                layer,
                &region,
                &RealJ0::new(arr0(C::from(k))),
            )
        };
        assert_eq!(
            at(1.0).unwrap_err(),
            ConstitutiveError::NormalComponentVanishes { region: layer }
        );
        assert!(at(1.5).is_ok());
        assert_eq!(
            at(f64::INFINITY).unwrap_err(),
            ConstitutiveError::NonFinite { region: layer }
        );
    }
}
