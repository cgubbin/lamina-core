mod builder;
mod layer;
mod termination;

pub use builder::StackBuilder;
pub use layer::Layer;
pub use termination::{Pml, StackSide, StretchProfile, Termination, Wall, WallKind};

use nalgebra::ComplexField;
use num_traits::Float;

use crate::{
    ComplexScalar, DifferentiableMaterial, DifferentiableMeromorphicMaterial, FiniteLayerIndex,
    IncidentSide, MeromorphicMaterial,
    material::{
        AnalyticalMaterialHandle, AnalyticalMaterialRegion, Constitutive,
        DifferentiableMaterialHandle, Material, MaterialHandle, MeromorphicMaterialHandle,
    },
};

/// Heterogeneous stack supporting real-axis constitutive evaluation.
pub type MaterialStack<C> = Stack<MaterialHandle<C>, <C as ComplexField>::RealField>;

/// Heterogeneous stack supporting real-axis derivatives.
pub type DifferentiableMaterialStack<C> =
    Stack<DifferentiableMaterialHandle<C>, <C as ComplexField>::RealField>;

/// Heterogeneous stack supporting complex continuation.
pub type MeromorphicMaterialStack<C> =
    Stack<MeromorphicMaterialHandle<C>, <C as ComplexField>::RealField>;

/// Heterogeneous stack supporting complex continuation and derivatives.
pub type AnalyticalMaterialStack<C> =
    Stack<AnalyticalMaterialHandle<C>, <C as ComplexField>::RealField>;

/// Stack of regions, each isotropic or tensor-valued, supporting complex
/// continuation and derivatives.
pub type AnalyticalTensorStack<C> =
    Stack<AnalyticalMaterialRegion<C>, <C as ComplexField>::RealField>;

/// One region of a stack: an exterior medium or a finite layer.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StackRegion {
    LeftExterior,
    Layer(FiniteLayerIndex),
    RightExterior,
}

/// A planar stack of finite layers between two exterior media.
///
/// Finite layers are stored in geometric left-to-right order. The exterior
/// media are not included in `len()` or layer iteration.
///
/// Each side carries a [`Termination`]: by default the exterior medium is
/// semi-infinite ([`Termination::Open`]); a side may instead be closed by a
/// [`Wall`] or a wall-backed [`Pml`] at an absolute position in the stack
/// frame, with the exterior medium filling the gap.
#[derive(Clone, Debug, PartialEq)]
pub struct Stack<M, F> {
    left_exterior: M,
    right_exterior: M,
    layers_left_to_right: Vec<Layer<M, F>>,
    left_termination: Termination<F>,
    right_termination: Termination<F>,
}

impl<M, F> Stack<M, F> {
    #[cfg(test)]
    pub(crate) fn new(
        left_exterior: M,
        layers_left_to_right: Vec<Layer<M, F>>,
        right_exterior: M,
    ) -> Self {
        Self {
            left_exterior,
            right_exterior,
            layers_left_to_right,
            left_termination: Termination::Open,
            right_termination: Termination::Open,
        }
    }

    pub fn builder(left_exterior: M, right_exterior: M) -> StackBuilder<M, F> {
        StackBuilder::new(left_exterior, right_exterior)
    }
}

impl<C> Stack<MaterialHandle<C>, C::RealField>
where
    C: ComplexScalar + Copy + 'static,
    C::RealField: Copy + 'static,
{
    pub fn from_materials<L, U>(
        left_exterior: L,
        right_exterior: U,
    ) -> StackBuilder<MaterialHandle<C>, C::RealField>
    where
        L: Material<Real = C::RealField> + Send + Sync + 'static,
        U: Material<Real = C::RealField> + Send + Sync + 'static,
        C::RealField: Float,
    {
        StackBuilder::from_materials(left_exterior, right_exterior)
    }
}

impl<C> Stack<AnalyticalMaterialHandle<C>, C::RealField>
where
    C: ComplexScalar + Copy + 'static,
    C::RealField: Copy + 'static,
{
    pub fn from_analytical_materials<L, U>(
        left_exterior: L,
        right_exterior: U,
    ) -> StackBuilder<AnalyticalMaterialHandle<C>, C::RealField>
    where
        L: DifferentiableMeromorphicMaterial<Real = C::RealField> + Send + Sync + 'static,
        U: DifferentiableMeromorphicMaterial<Real = C::RealField> + Send + Sync + 'static,
        C::RealField: Float,
    {
        StackBuilder::from_analytical_materials(left_exterior, right_exterior)
    }
}

impl<C> Stack<DifferentiableMaterialHandle<C>, C::RealField>
where
    C: ComplexScalar + Copy + 'static,
    C::RealField: Copy + 'static,
{
    pub fn from_differentiable_materials<L, U>(
        left_exterior: L,
        right_exterior: U,
    ) -> StackBuilder<DifferentiableMaterialHandle<C>, C::RealField>
    where
        L: DifferentiableMaterial<Real = C::RealField> + Send + Sync + 'static,
        U: DifferentiableMaterial<Real = C::RealField> + Send + Sync + 'static,
        C::RealField: Float,
    {
        StackBuilder::from_differentiable_materials(left_exterior, right_exterior)
    }
}

impl<C> Stack<MeromorphicMaterialHandle<C>, C::RealField>
where
    C: ComplexScalar + Copy + 'static,
    C::RealField: Copy + 'static,
{
    pub fn from_meromorphic_materials<L, U>(
        left_exterior: L,
        right_exterior: U,
    ) -> StackBuilder<MeromorphicMaterialHandle<C>, C::RealField>
    where
        L: MeromorphicMaterial<Real = C::RealField> + Send + Sync + 'static,
        U: MeromorphicMaterial<Real = C::RealField> + Send + Sync + 'static,
        C::RealField: Float,
    {
        StackBuilder::from_meromorphic_materials(left_exterior, right_exterior)
    }
}

impl<M, F> Stack<M, F> {
    /// Return the left semi-infinite exterior material.
    pub fn left_exterior(&self) -> &M {
        &self.left_exterior
    }

    /// Return the right semi-infinite exterior material.
    pub fn right_exterior(&self) -> &M {
        &self.right_exterior
    }

    /// Return the finite layers in geometric left-to-right order.
    pub fn layers_left_to_right(&self) -> &[Layer<M, F>] {
        &self.layers_left_to_right
    }

    /// Return the termination on the left side (open by default).
    pub fn left_termination(&self) -> &Termination<F> {
        &self.left_termination
    }

    /// Return the termination on the right side (open by default).
    pub fn right_termination(&self) -> &Termination<F> {
        &self.right_termination
    }

    /// Return `true` when both sides are open.
    pub fn is_open(&self) -> bool {
        self.left_termination.is_open() && self.right_termination.is_open()
    }

    /// Replace the left termination.
    pub fn with_left_termination(mut self, termination: impl Into<Termination<F>>) -> Self {
        self.left_termination = termination.into();
        self
    }

    /// Replace the right termination.
    pub fn with_right_termination(mut self, termination: impl Into<Termination<F>>) -> Self {
        self.right_termination = termination.into();
        self
    }

    /// Return the number of finite layers.
    ///
    /// The two exterior media are not included.
    pub fn len(&self) -> usize {
        self.layers_left_to_right.len()
    }

    /// Return `true` when the stack contains no finite layers.
    ///
    /// The exterior media do not affect this result.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate over finite layers in geometric left-to-right order.
    pub fn iter(&self) -> impl DoubleEndedIterator<Item = &Layer<M, F>> {
        self.layers_left_to_right.iter()
    }

    pub(crate) fn incident_exterior(&self, incident_side: IncidentSide) -> &M {
        match incident_side {
            IncidentSide::Left => self.left_exterior(),

            IncidentSide::Right => self.right_exterior(),
        }
    }

    #[cfg(test)]
    pub(crate) fn into_parts(self) -> (M, Vec<Layer<M, F>>, M) {
        (
            self.left_exterior,
            self.layers_left_to_right,
            self.right_exterior,
        )
    }
}

impl<S, T, R, F> From<Stack<S, F>> for Stack<Constitutive<S, T, R>, F> {
    /// Lift a stack of scalar materials to isotropic regions, unchanged.
    fn from(stack: Stack<S, F>) -> Self {
        Self {
            left_exterior: stack.left_exterior.into(),
            right_exterior: stack.right_exterior.into(),
            layers_left_to_right: stack
                .layers_left_to_right
                .into_iter()
                .map(|layer| {
                    let (material, thickness) = layer.into_parts();
                    Layer::new(material.into(), thickness)
                })
                .collect(),
            left_termination: stack.left_termination,
            right_termination: stack.right_termination,
        }
    }
}

impl<S, T, R, F> Stack<Constitutive<S, T, R>, F> {
    /// Return `true` when every region, exteriors included, is separable:
    /// TE and TM then decouple throughout the stack.
    pub fn is_separable(&self) -> bool {
        self.left_exterior.is_separable()
            && self.right_exterior.is_separable()
            && self.iter().all(|layer| layer.material().is_separable())
    }
}

#[cfg(test)]
mod tests {
    use lamina_units::Length;
    use num_complex::Complex64;

    use super::*;
    use crate::{ComplexScalar, Constant, Material, Sampled, material::MaterialHandle};

    #[derive(Clone, Debug, PartialEq)]
    struct TestMaterial {
        id: usize,
        relative_permittivity: f64,
    }

    impl Material for TestMaterial {
        type Real = f64;

        fn relative_permeability<I, C>(&self, vacuum_angular_wavenumber: I) -> I::Mapped<C>
        where
            I: Sampled<Elem = Self::Real>,
            C: ComplexScalar<RealField = Self::Real>,
        {
            vacuum_angular_wavenumber.map(|_| C::one())
        }

        fn relative_permittivity<I, C>(&self, vacuum_angular_wavenumber: I) -> I::Mapped<C>
        where
            I: Sampled<Elem = Self::Real>,
            C: ComplexScalar<RealField = Self::Real>,
        {
            let epsilon = C::from_real(self.relative_permittivity);

            vacuum_angular_wavenumber.map(|_| epsilon)
        }
    }

    fn test_stack() -> Stack<TestMaterial, f64> {
        Stack::builder(
            TestMaterial {
                id: 0,
                relative_permittivity: 1.0,
            },
            TestMaterial {
                id: 3,
                relative_permittivity: 4.0,
            },
        )
        .layer(
            TestMaterial {
                id: 1,
                relative_permittivity: 2.0,
            },
            Length::nanometres(100.0),
        )
        .layer(
            TestMaterial {
                id: 2,
                relative_permittivity: 3.0,
            },
            Length::micrometres(2.0),
        )
        .finalise()
    }

    #[test]
    fn builder_preserves_exterior_media() {
        let stack = test_stack();

        assert_eq!(stack.left_exterior().id, 0);
        assert_eq!(stack.right_exterior().id, 3);
    }

    #[test]
    fn builder_preserves_layer_order() {
        let stack = test_stack();

        let ids = stack
            .layers_left_to_right()
            .iter()
            .map(|layer| layer.material().id)
            .collect::<Vec<_>>();

        assert_eq!(ids, vec![1, 2]);
    }

    #[test]
    fn len_and_is_empty_describe_only_finite_layers() {
        let stack = test_stack();

        assert_eq!(stack.len(), 2);
        assert!(!stack.is_empty());

        let empty = Stack::<TestMaterial, f64>::builder(
            TestMaterial {
                id: 0,
                relative_permittivity: 1.0,
            },
            TestMaterial {
                id: 1,
                relative_permittivity: 1.0,
            },
        )
        .finalise();

        assert_eq!(empty.len(), 0);
        assert!(empty.is_empty());
    }

    #[test]
    fn iter_follows_left_to_right_order() {
        let stack = test_stack();

        let ids = stack
            .iter()
            .map(|layer| layer.material().id)
            .collect::<Vec<_>>();

        assert_eq!(ids, vec![1, 2]);
    }

    #[test]
    fn incident_exterior_selects_requested_side() {
        let stack = test_stack();

        assert_eq!(stack.incident_exterior(IncidentSide::Left).id, 0,);

        assert_eq!(stack.incident_exterior(IncidentSide::Right).id, 3,);
    }

    #[test]
    fn layer_preserves_material_and_thickness() {
        let layer = Layer::new(
            TestMaterial {
                id: 7,
                relative_permittivity: 2.5,
            },
            Length::nanometres(350.0),
        );

        assert_eq!(layer.material().id, 7);
        assert_eq!(layer.thickness(), Length::nanometres(350.0),);

        let (material, thickness) = layer.into_parts();

        assert_eq!(material.id, 7);
        assert_eq!(thickness, Length::nanometres(350.0),);
    }

    #[test]
    fn material_handle_supports_heterogeneous_core_materials() {
        type C = Complex64;
        type Handle = MaterialHandle<C>;

        let stack = StackBuilder::<Handle, f64>::from_materials(
            Constant::dielectric(1.0),
            TestMaterial {
                id: 10,
                relative_permittivity: 2.25,
            },
        )
        .material_layer(Constant::dielectric(4.0), Length::centimetres(100.0))
        .material_layer(
            TestMaterial {
                id: 11,
                relative_permittivity: 6.25,
            },
            Length::centimetres(200.0),
        )
        .finalise();

        assert_eq!(stack.len(), 2);
    }

    #[test]
    fn into_parts_preserves_stack_order() {
        let stack = test_stack();

        let (left, layers, right) = stack.into_parts();

        assert_eq!(left.id, 0);
        assert_eq!(right.id, 3);

        assert_eq!(
            layers
                .iter()
                .map(|layer| layer.material().id)
                .collect::<Vec<_>>(),
            vec![1, 2],
        );
    }

    mod regions {
        use ndarray::Ix0;

        use super::*;
        use crate::{
            AnalyticalTensorMaterialHandle, ConstitutiveClass, Rotation, SignedAxis,
            ValidationConfig, compile_canonical_constant_stack, test_support::jet::HoloJ2,
        };

        type Region = AnalyticalMaterialRegion<Complex64>;

        /// A principal-axis model (a scalar seen as a tensor) under `rotation`.
        fn oriented(rotation: Rotation<f64>) -> Region {
            Region::oriented(
                AnalyticalTensorMaterialHandle::new(Constant::dielectric(2.5)),
                rotation,
            )
        }

        fn diagonal() -> Region {
            oriented(
                Rotation::signed_permutation(
                    SignedAxis::PlusY,
                    SignedAxis::MinusX,
                    SignedAxis::PlusZ,
                )
                .unwrap(),
            )
        }

        fn isotropic(epsilon: f64) -> Region {
            AnalyticalMaterialHandle::new(Constant::dielectric(epsilon)).into()
        }

        fn stack_with(left: Region, layer: Region) -> AnalyticalTensorStack<Complex64> {
            Stack::builder(left, isotropic(2.0))
                .layer(isotropic(3.0), Length::centimetres(0.1))
                .layer(layer, Length::centimetres(0.2))
                .finalise()
        }

        #[test]
        fn scalar_stack_lifts_and_returns_unchanged() {
            let scalar: AnalyticalMaterialStack<Complex64> =
                Stack::from_analytical_materials(Constant::vacuum(), Constant::new(2.1, 1.0))
                    .analytical_layer(Constant::new(3.0, 1.2), Length::centimetres(0.12))
                    .finalise()
                    .with_right_termination(Wall::new(
                        WallKind::PerfectMagnetic,
                        Length::centimetres(0.2),
                    ));
            let lifted = AnalyticalTensorStack::<Complex64>::from(scalar.clone());

            assert!(lifted.is_separable());
            assert_eq!(
                lifted.left_exterior().scalar(),
                Some(scalar.left_exterior())
            );
            assert_eq!(
                lifted.right_exterior().scalar(),
                Some(scalar.right_exterior())
            );
            assert_eq!(
                lifted.layers_left_to_right()[0].material().scalar(),
                Some(scalar.layers_left_to_right()[0].material())
            );
            assert_eq!(lifted.right_termination(), scalar.right_termination());
        }

        #[test]
        fn classes_and_separability_of_a_mixed_stack() {
            let about_y = oriented(Rotation::from_euler(0.0, 0.6, 0.0));
            let separable = stack_with(diagonal(), about_y);
            assert!(separable.is_separable());
            assert_eq!(
                separable.left_exterior().class(),
                ConstitutiveClass::Diagonal
            );
            assert_eq!(
                separable.layers_left_to_right()[1].material().class(),
                ConstitutiveClass::Full
            );

            let general = oriented(Rotation::from_euler(0.3, 0.6, 0.0));
            assert!(!stack_with(isotropic(1.0), general).is_separable());
        }

        #[test]
        fn wall_in_a_uniaxial_exterior_records_the_gap_class() {
            let stack =
                stack_with(oriented(Rotation::identity()), isotropic(2.0)).with_left_termination(
                    Wall::new(WallKind::PerfectElectric, Length::centimetres(-0.05)),
                );

            let compiled = compile_canonical_constant_stack::<_, HoloJ2>(
                &stack,
                Ix0(),
                &ValidationConfig::permissive(),
            )
            .unwrap();
            assert_eq!(
                compiled.left_exterior().class(),
                ConstitutiveClass::Diagonal
            );
            assert_eq!(
                compiled.left_exterior().tensor().unwrap().1,
                &Rotation::identity()
            );
            assert!(!compiled.left_termination().is_open());
        }
    }
}
