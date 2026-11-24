mod builder;
mod layer;
mod termination;

pub use builder::StackBuilder;
pub use layer::Layer;
pub use termination::{Pml, StackSide, StretchProfile, Termination, Wall, WallKind};

use nalgebra::ComplexField;
use num_traits::Float;

use crate::{
    ComplexScalar, DifferentiableMaterial, DifferentiableMeromorphicMaterial, IncidentSide,
    MeromorphicMaterial,
    material::{
        AnalyticalMaterialHandle, DifferentiableMaterialHandle, Material, MaterialHandle,
        MeromorphicMaterialHandle,
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
}
