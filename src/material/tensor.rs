//! Tensor-valued material capabilities and their erased handles.
//!
//! The four traits mirror the scalar capability levels. Each returns the
//! relative permittivity and permeability as 3×3 tensors, row-major
//! `[[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]]`, in the material's own
//! (crystal) frame. Every scalar material is also a tensor material, with
//! `ε I` and `μ I`.

use std::{any::type_name, fmt, sync::Arc};

use num_traits::{One, Zero};

use super::{
    DerivativeOrder, DifferentiableMaterial, DifferentiableMeromorphicMaterial, Material,
    MeromorphicMaterial, Sampled, Scalar,
};
use crate::ComplexScalar;

/// A 3×3 tensor, row-major.
pub type Tensor3<C> = [[C; 3]; 3];

fn isotropic<C: Copy + Zero>(value: C) -> Tensor3<C> {
    let zero = C::zero();
    [
        [value, zero, zero],
        [zero, value, zero],
        [zero, zero, value],
    ]
}

/// Tensor material model evaluated on the real spectral axis.
///
/// The tensor counterpart of [`Material`]; see there for the spectral
/// variable and its units.
pub trait TensorMaterial {
    /// Real scalar type used by the material.
    type Real: One + Zero;

    /// `true` when ε and μ are diagonal in the material's own frame at every
    /// wavenumber, by construction.
    ///
    /// The structural class of an oriented region is decided from this
    /// declaration and its rotation, never from tensor values.
    const PRINCIPAL_AXES: bool = false;

    /// Return the relative permittivity tensor in the material frame.
    fn relative_permittivity_tensor<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = Self::Real>;

    /// Return the relative permeability tensor in the material frame.
    ///
    /// The default implementation represents a non-magnetic material.
    fn relative_permeability_tensor<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = Self::Real>,
    {
        vacuum_angular_wavenumber.map(|_| isotropic(C::one()))
    }
}

/// Tensor material model with derivatives on the real spectral axis.
///
/// The tensor counterpart of [`DifferentiableMaterial`].
pub trait DifferentiableTensorMaterial: TensorMaterial {
    /// Return a derivative of the relative permittivity tensor with respect
    /// to `k₀`.
    fn relative_permittivity_tensor_derivative<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = Self::Real>;

    /// Return a derivative of the relative permeability tensor with respect
    /// to `k₀`.
    ///
    /// The default implementation represents a non-magnetic material.
    fn relative_permeability_tensor_derivative<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
        _order: DerivativeOrder,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = Self::Real>,
    {
        vacuum_angular_wavenumber.map(|_| [[C::zero(); 3]; 3])
    }
}

/// Tensor material model with a meromorphic continuation into the complex
/// spectral plane.
///
/// The tensor counterpart of [`MeromorphicMaterial`], whose conditions on the
/// continuation apply entry by entry.
pub trait MeromorphicTensorMaterial: TensorMaterial {
    /// Return the complex relative permittivity tensor.
    fn relative_permittivity_tensor_complex<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = C>;

    /// Return the complex relative permeability tensor.
    fn relative_permeability_tensor_complex<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = C>;
}

/// Meromorphic tensor material model with complex derivatives.
///
/// The tensor counterpart of [`DifferentiableMeromorphicMaterial`].
pub trait DifferentiableMeromorphicTensorMaterial:
    MeromorphicTensorMaterial + DifferentiableTensorMaterial
{
    /// Return a complex derivative of the relative permittivity tensor.
    fn relative_permittivity_tensor_complex_derivative<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = C>;

    /// Return a complex derivative of the relative permeability tensor.
    ///
    /// The default implementation represents a non-magnetic material.
    fn relative_permeability_tensor_complex_derivative<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
        _order: DerivativeOrder,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = C>,
    {
        vacuum_angular_wavenumber.map(|_| [[C::zero(); 3]; 3])
    }
}

/// Every scalar material is a tensor material with `ε I` and `μ I`, diagonal
/// in its own frame.
macro_rules! scalar_as_isotropic {
    ($tensor:ident for $scalar:ident, $k0:ty { $($item:item)* } $($method:ident => $scalar_method:ident $(($order:ident))?;)+) => {
        impl<M> $tensor for M
        where
            M: $scalar,
            M::Real: Copy,
        {
            $($item)*
            $(
                fn $method<I, C>(&self, k0: I $(, $order: DerivativeOrder)?) -> I::Mapped<Tensor3<C>>
                where
                    C: ComplexScalar<RealField = Self::Real>,
                    I: Sampled<Elem = $k0>,
                {
                    k0.map(|k0| isotropic(self.$scalar_method(Scalar(k0) $(, $order)?)))
                }
            )+
        }
    };
}

scalar_as_isotropic!(TensorMaterial for Material, Self::Real {
    type Real = M::Real;
    const PRINCIPAL_AXES: bool = true;
}
    relative_permittivity_tensor => relative_permittivity;
    relative_permeability_tensor => relative_permeability;
);
scalar_as_isotropic!(DifferentiableTensorMaterial for DifferentiableMaterial, Self::Real {}
    relative_permittivity_tensor_derivative => relative_permittivity_derivative(order);
    relative_permeability_tensor_derivative => relative_permeability_derivative(order);
);
scalar_as_isotropic!(MeromorphicTensorMaterial for MeromorphicMaterial, C {}
    relative_permittivity_tensor_complex => relative_permittivity_complex;
    relative_permeability_tensor_complex => relative_permeability_complex;
);
scalar_as_isotropic!(DifferentiableMeromorphicTensorMaterial for DifferentiableMeromorphicMaterial, C {}
    relative_permittivity_tensor_complex_derivative => relative_permittivity_complex_derivative(order);
    relative_permeability_tensor_complex_derivative => relative_permeability_complex_derivative(order);
);

/// Object-safe pointwise form of [`DifferentiableMeromorphicTensorMaterial`],
/// implemented for every such material and erased by
/// [`AnalyticalTensorMaterialHandle`].
pub(crate) trait TensorMaterialPoint<C: ComplexScalar>: Send + Sync {
    fn permittivity_at(&self, k0: C::RealField) -> Tensor3<C>;
    fn permeability_at(&self, k0: C::RealField) -> Tensor3<C>;
    fn permittivity_derivative_at(&self, k0: C::RealField, order: DerivativeOrder) -> Tensor3<C>;
    fn permeability_derivative_at(&self, k0: C::RealField, order: DerivativeOrder) -> Tensor3<C>;
    fn permittivity_complex_at(&self, k0: C) -> Tensor3<C>;
    fn permeability_complex_at(&self, k0: C) -> Tensor3<C>;
    fn permittivity_complex_derivative_at(&self, k0: C, order: DerivativeOrder) -> Tensor3<C>;
    fn permeability_complex_derivative_at(&self, k0: C, order: DerivativeOrder) -> Tensor3<C>;
}

impl<M, C> TensorMaterialPoint<C> for M
where
    M: DifferentiableMeromorphicTensorMaterial<Real = C::RealField> + Send + Sync,
    C: ComplexScalar,
    C::RealField: Copy,
{
    fn permittivity_at(&self, k0: C::RealField) -> Tensor3<C> {
        self.relative_permittivity_tensor(Scalar(k0))
    }

    fn permeability_at(&self, k0: C::RealField) -> Tensor3<C> {
        self.relative_permeability_tensor(Scalar(k0))
    }

    fn permittivity_derivative_at(&self, k0: C::RealField, order: DerivativeOrder) -> Tensor3<C> {
        self.relative_permittivity_tensor_derivative(Scalar(k0), order)
    }

    fn permeability_derivative_at(&self, k0: C::RealField, order: DerivativeOrder) -> Tensor3<C> {
        self.relative_permeability_tensor_derivative(Scalar(k0), order)
    }

    fn permittivity_complex_at(&self, k0: C) -> Tensor3<C> {
        self.relative_permittivity_tensor_complex(Scalar(k0))
    }

    fn permeability_complex_at(&self, k0: C) -> Tensor3<C> {
        self.relative_permeability_tensor_complex(Scalar(k0))
    }

    fn permittivity_complex_derivative_at(&self, k0: C, order: DerivativeOrder) -> Tensor3<C> {
        self.relative_permittivity_tensor_complex_derivative(Scalar(k0), order)
    }

    fn permeability_complex_derivative_at(&self, k0: C, order: DerivativeOrder) -> Tensor3<C> {
        self.relative_permeability_tensor_complex_derivative(Scalar(k0), order)
    }
}

/// Type-erased tensor material with a complex continuation and derivatives.
///
/// The handle records whether the material declared
/// [`TensorMaterial::PRINCIPAL_AXES`]. Equality compares handle identity, as
/// for the scalar handles.
pub struct AnalyticalTensorMaterialHandle<C: ComplexScalar> {
    type_name: &'static str,
    principal_axes: bool,
    inner: Arc<dyn TensorMaterialPoint<C>>,
}

impl<C> AnalyticalTensorMaterialHandle<C>
where
    C: ComplexScalar + 'static,
    C::RealField: Copy + 'static,
{
    /// Erase a concrete tensor material.
    pub fn new<M>(material: M) -> Self
    where
        M: DifferentiableMeromorphicTensorMaterial<Real = C::RealField> + Send + Sync + 'static,
    {
        Self {
            type_name: type_name::<M>(),
            principal_axes: M::PRINCIPAL_AXES,
            inner: Arc::new(material),
        }
    }
}

impl<C: ComplexScalar> AnalyticalTensorMaterialHandle<C> {
    pub(crate) fn principal_axes(&self) -> bool {
        self.principal_axes
    }

    pub(crate) fn point(&self) -> &dyn TensorMaterialPoint<C> {
        &*self.inner
    }
}

impl<C: ComplexScalar> Clone for AnalyticalTensorMaterialHandle<C> {
    fn clone(&self) -> Self {
        Self {
            type_name: self.type_name,
            principal_axes: self.principal_axes,
            inner: Arc::clone(&self.inner),
        }
    }
}

impl<C: ComplexScalar> PartialEq for AnalyticalTensorMaterialHandle<C> {
    fn eq(&self, other: &Self) -> bool {
        self.type_name == other.type_name && Arc::ptr_eq(&self.inner, &other.inner)
    }
}

impl<C: ComplexScalar> fmt::Debug for AnalyticalTensorMaterialHandle<C> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("AnalyticalTensorMaterialHandle")
            .field("type", &self.type_name)
            .field("principal_axes", &self.principal_axes)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::*;
    use crate::Constant;

    type C = Complex64;

    /// ε(k) = 2 + k², μ(k) = 1 + k: a dispersive scalar model.
    struct Quadratic;

    impl Material for Quadratic {
        type Real = f64;

        fn relative_permittivity<I, X>(&self, k: I) -> I::Mapped<X>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = f64>,
        {
            k.map(|k| X::from_real(2.0 + k * k))
        }

        fn relative_permeability<I, X>(&self, k: I) -> I::Mapped<X>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = f64>,
        {
            k.map(|k| X::from_real(1.0 + k))
        }
    }

    impl DifferentiableMaterial for Quadratic {
        fn relative_permittivity_derivative<I, X>(
            &self,
            k: I,
            order: DerivativeOrder,
        ) -> I::Mapped<X>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = f64>,
        {
            k.map(|k| {
                X::from_real(match order {
                    DerivativeOrder::First => 2.0 * k,
                    DerivativeOrder::Second => 2.0,
                    DerivativeOrder::Third => 0.0,
                })
            })
        }

        fn relative_permeability_derivative<I, X>(
            &self,
            k: I,
            order: DerivativeOrder,
        ) -> I::Mapped<X>
        where
            X: ComplexScalar<RealField = f64>,
            I: Sampled<Elem = f64>,
        {
            k.map(|_| {
                X::from_real(match order {
                    DerivativeOrder::First => 1.0,
                    _ => 0.0,
                })
            })
        }
    }

    #[test]
    fn scalar_material_is_epsilon_and_mu_times_identity() {
        let k = 1.5;
        let epsilon: Tensor3<C> = Quadratic.relative_permittivity_tensor(Scalar(k));
        let mu: Tensor3<C> = Quadratic.relative_permeability_tensor(Scalar(k));
        assert_eq!(epsilon, isotropic(C::from(2.0 + k * k)));
        assert_eq!(mu, isotropic(C::from(1.0 + k)));

        for (order, d_epsilon, d_mu) in [
            (DerivativeOrder::First, 2.0 * k, 1.0),
            (DerivativeOrder::Second, 2.0, 0.0),
        ] {
            let epsilon: Tensor3<C> =
                Quadratic.relative_permittivity_tensor_derivative(Scalar(k), order);
            let mu: Tensor3<C> =
                Quadratic.relative_permeability_tensor_derivative(Scalar(k), order);
            assert_eq!(epsilon, isotropic(C::from(d_epsilon)));
            assert_eq!(mu, isotropic(C::from(d_mu)));
        }
        const { assert!(<Quadratic as TensorMaterial>::PRINCIPAL_AXES) };
    }

    #[test]
    fn complex_continuation_of_a_scalar_material_is_isotropic() {
        let constant = Constant::magnetodielectric(4.0, 2.0);
        let k = C::new(1.0, 0.5);
        let epsilon: Tensor3<C> = constant.relative_permittivity_tensor_complex(Scalar(k));
        let mu: Tensor3<C> = constant.relative_permeability_tensor_complex(Scalar(k));
        let d_epsilon: Tensor3<C> = constant
            .relative_permittivity_tensor_complex_derivative(Scalar(k), DerivativeOrder::First);
        assert_eq!(epsilon, isotropic(C::from(4.0)));
        assert_eq!(mu, isotropic(C::from(2.0)));
        assert_eq!(d_epsilon, [[C::from(0.0); 3]; 3]);
    }

    #[test]
    fn handles_record_the_declaration_and_compare_by_identity() {
        let handle = AnalyticalTensorMaterialHandle::<C>::new(Constant::dielectric(2.0));
        assert!(handle.principal_axes());
        assert_eq!(handle.clone(), handle);
        assert_ne!(
            handle,
            AnalyticalTensorMaterialHandle::<C>::new(Constant::dielectric(2.0))
        );
        assert_eq!(handle.point().permittivity_at(1.0), isotropic(C::from(2.0)));
    }
}
