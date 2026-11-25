use ndarray::Dimension;
use num_traits::{FromPrimitive, One};
use thiserror::Error;

use crate::{
    ComplexScalar, FiniteLayerIndex, StackRegion,
    algebra::{Jet, RealScalarAlgebra, ScalarAlgebra},
    backend::{ChannelEntries, LayerQuantities, PrincipalConstitutive},
    material::{ConstitutiveError, TensorConstitutiveSpectralFirstLift},
    observable::{
        LayerAggregateError,
        layer::{
            IntegratedHermitianCrossStateProducts,
            integration::{
                IntegratedBilinearCrossStateProducts, project_integrated_bilinear_field_overlap,
            },
            overlap::BilinearLayerNormalization,
        },
    },
};

use super::{
    LayerProjectionError, Layers, integration::project_integrated_field_norms,
    project::IntegratedLayerData,
};

#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum LayerEnergyError {
    #[error(transparent)]
    Aggregation(#[from] LayerAggregateError),

    #[error(transparent)]
    Projection(#[from] LayerProjectionError),

    #[error(transparent)]
    Constitutive(#[from] ConstitutiveError),

    #[error(
        "integrated layer count {layer_count} does not match differentiable \
     material count {material_count}"
    )]
    MaterialCountMismatch {
        layer_count: usize,
        material_count: usize,
    },
}

/// Integrated electromagnetic energy associated with one finite layer.
///
/// The electric and magnetic components are calculated using either:
///
/// - nondispersive constitutive weights, `Re(epsilon)` and `Re(mu)`; or
/// - Brillouin weights,
///   `Re[d(k0 epsilon)/d k0]` and `Re[d(k0 mu)/d k0]`.
///
/// Evaluator methods determine which definition is used. The default
/// differentiable-material path uses Brillouin energy, while the explicit
/// fallback uses nondispersive energy.
///
/// Results are normalized per unit incident-wave power flux.
///
/// ```text
/// total = electric + magnetic
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct LayerEnergy<R> {
    electric: R,
    magnetic: R,
    total: R,
}

impl<R> LayerEnergy<R> {
    pub(crate) const fn new(electric: R, magnetic: R, total: R) -> Self {
        Self {
            electric,
            magnetic,
            total,
        }
    }

    /// Return the integrated electric-energy contribution.
    pub fn electric(&self) -> &R {
        &self.electric
    }

    /// Return the integrated magnetic-energy contribution.
    pub fn magnetic(&self) -> &R {
        &self.magnetic
    }

    /// Return the total integrated layer energy.
    pub fn total(&self) -> &R {
        &self.total
    }

    /// Consume the value and return `(electric, magnetic, total)`.
    pub fn into_parts(self) -> (R, R, R) {
        (self.electric, self.magnetic, self.total)
    }

    /// Transform all energy components.
    pub fn map<U>(self, mut map: impl FnMut(R) -> U) -> LayerEnergy<U> {
        LayerEnergy {
            electric: map(self.electric),
            magnetic: map(self.magnetic),
            total: map(self.total),
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct BrillouinLayerInput<A> {
    integrated: IntegratedLayerData<IntegratedHermitianCrossStateProducts<A>, A>,
    derivative: PrincipalConstitutive<A>,
}

impl<A> BrillouinLayerInput<A> {
    pub(crate) const fn new(
        integrated: IntegratedLayerData<IntegratedHermitianCrossStateProducts<A>, A>,
        derivative: PrincipalConstitutive<A>,
    ) -> Self {
        Self {
            integrated,
            derivative,
        }
    }

    #[cfg(test)]
    pub(crate) fn derivative(&self) -> &PrincipalConstitutive<A> {
        &self.derivative
    }

    pub(crate) fn into_parts(
        self,
    ) -> (
        IntegratedLayerData<IntegratedHermitianCrossStateProducts<A>, A>,
        PrincipalConstitutive<A>,
    ) {
        (self.integrated, self.derivative)
    }
}

impl<A> Layers<IntegratedLayerData<IntegratedHermitianCrossStateProducts<A>, A>> {
    pub(crate) fn into_brillouin_layers<'a, E, M>(
        self,
        materials: impl ExactSizeIterator<Item = &'a M>,
        vacuum_angular_wavenumber: &A,
    ) -> Result<Layers<BrillouinLayerInput<A>>, LayerEnergyError>
    where
        M: 'a,
        A: ScalarAlgebra + TensorConstitutiveSpectralFirstLift<E, M> + Clone,
        A::Scalar: ComplexScalar,
        A::Dimension: Dimension,
    {
        if self.len() != materials.len() {
            return Err(LayerEnergyError::MaterialCountMismatch {
                layer_count: self.len(),
                material_count: materials.len(),
            });
        }

        Ok(Layers::new(
            self.into_inner()
                .into_iter()
                .zip(materials)
                .enumerate()
                .map(|(index, (layer, material))| {
                    let derivative = PrincipalConstitutive::new(
                        StackRegion::Layer(FiniteLayerIndex::new(index)),
                        A::constitutive_tensors_spectral_first(material, vacuum_angular_wavenumber),
                    )?;
                    Ok(BrillouinLayerInput::new(layer, derivative))
                })
                .collect::<Result<_, ConstitutiveError>>()?,
        ))
    }
}

#[derive(Clone, Debug)]
pub(crate) struct BrillouinBilinearLayerInput<A> {
    integrated: IntegratedLayerData<IntegratedBilinearCrossStateProducts<A>, A>,
    derivative: PrincipalConstitutive<A>,
}

impl<A> BrillouinBilinearLayerInput<A> {
    pub(crate) const fn new(
        integrated: IntegratedLayerData<IntegratedBilinearCrossStateProducts<A>, A>,
        derivative: PrincipalConstitutive<A>,
    ) -> Self {
        Self {
            integrated,
            derivative,
        }
    }

    pub(crate) fn into_parts(
        self,
    ) -> (
        IntegratedLayerData<IntegratedBilinearCrossStateProducts<A>, A>,
        PrincipalConstitutive<A>,
    ) {
        (self.integrated, self.derivative)
    }
}

impl<A> Layers<IntegratedLayerData<IntegratedBilinearCrossStateProducts<A>, A>> {
    pub(crate) fn into_brillouin_layers<'a, E, M>(
        self,
        materials: impl ExactSizeIterator<Item = &'a M>,
        vacuum_angular_wavenumber: &A,
    ) -> Result<Layers<BrillouinBilinearLayerInput<A>>, LayerEnergyError>
    where
        M: 'a,
        A: ScalarAlgebra + TensorConstitutiveSpectralFirstLift<E, M> + Clone,
        A::Scalar: ComplexScalar,
        A::Dimension: Dimension,
    {
        if self.len() != materials.len() {
            return Err(LayerEnergyError::MaterialCountMismatch {
                layer_count: self.len(),
                material_count: materials.len(),
            });
        }

        Ok(Layers::new(
            self.into_inner()
                .into_iter()
                .zip(materials)
                .enumerate()
                .map(|(index, (layer, material))| {
                    let derivative = PrincipalConstitutive::new(
                        StackRegion::Layer(FiniteLayerIndex::new(index)),
                        A::constitutive_tensors_spectral_first(material, vacuum_angular_wavenumber),
                    )?;
                    Ok(BrillouinBilinearLayerInput::new(layer, derivative))
                })
                .collect::<Result<_, ConstitutiveError>>()?,
        ))
    }
}

pub(crate) fn energy_density_prefactor<A>(vacuum_angular_wavenumber: &A) -> A::RealJet
where
    A: RealScalarAlgebra,
    A::RealJet: ScalarAlgebra,
    <A::RealJet as Jet>::Scalar: FromPrimitive,
{
    let vacuum = vacuum_angular_wavenumber.real();

    let quarter =
        <A::RealJet as Jet>::Scalar::from_f64(0.25).expect("one quarter must be representable");

    A::RealJet::filled_constant_like(vacuum.value(), quarter)
}

impl<A> IntegratedLayerData<IntegratedHermitianCrossStateProducts<A>, A> {
    fn into_nondispersive_energy(
        self,
        vacuum_angular_wavenumber: &A,
        parallel_angular_wavenumber: &A,
        normalization: &A::RealJet,
    ) -> LayerEnergy<A::RealJet>
    where
        A: RealScalarAlgebra,
        A::RealJet: ScalarAlgebra,
        <A::RealJet as Jet>::Scalar: One,
    {
        let (state_products, quantities) = self.into_parts();

        let (electric, magnetic) = project_integrated_field_norms(
            &state_products,
            &quantities,
            vacuum_angular_wavenumber,
            parallel_angular_wavenumber,
        )
        .weighted(
            quantities
                .channel()
                .map(|entry| entry.real().multiply(normalization)),
        );

        project_energy(electric, magnetic)
    }
}

impl<A> Layers<IntegratedLayerData<IntegratedHermitianCrossStateProducts<A>, A>> {
    pub(crate) fn into_nondispersive_energy(
        self,
        vacuum_angular_wavenumber: &A,
        parallel_angular_wavenumber: &A,
    ) -> Layers<LayerEnergy<A::RealJet>>
    where
        A: RealScalarAlgebra,
        A::RealJet: ScalarAlgebra,
        <A::RealJet as Jet>::Scalar: FromPrimitive + One,
    {
        let normalization = energy_density_prefactor(vacuum_angular_wavenumber);

        self.map(|layer| {
            layer.into_nondispersive_energy(
                vacuum_angular_wavenumber,
                parallel_angular_wavenumber,
                &normalization,
            )
        })
    }
}

/// Brillouin weights `∂(k0 x)/∂k0 = x + k0 ∂x/∂k0` of the entries
/// `quantities` sees, from their intrinsic spectral derivatives
/// `derivative`, which are distinct from the outer derivative coordinates of
/// `A`.
fn brillouin_weights<A>(
    vacuum_angular_wavenumber: &A,
    quantities: &LayerQuantities<A>,
    derivative: &PrincipalConstitutive<A>,
) -> ChannelEntries<A>
where
    A: ScalarAlgebra,
{
    quantities.channel().zip(
        derivative.channel(quantities.polarisation()),
        |value, first| value.add(&vacuum_angular_wavenumber.multiply(first)),
    )
}

impl<A> BrillouinLayerInput<A> {
    fn into_brillouin_energy(
        self,
        vacuum_angular_wavenumber: &A,
        parallel_angular_wavenumber: &A,
        normalization: &A::RealJet,
    ) -> LayerEnergy<A::RealJet>
    where
        A: RealScalarAlgebra,
        A::RealJet: ScalarAlgebra,
        <A::RealJet as Jet>::Scalar: One,
    {
        let (layer, derivative) = self.into_parts();
        let (state_products, quantities) = layer.into_parts();

        let weights = brillouin_weights(vacuum_angular_wavenumber, &quantities, &derivative)
            .map(|weight| weight.real().multiply(normalization));

        let (electric, magnetic) = project_integrated_field_norms(
            &state_products,
            &quantities,
            vacuum_angular_wavenumber,
            parallel_angular_wavenumber,
        )
        .weighted(weights);

        project_energy(electric, magnetic)
    }
}

fn project_energy<R>(electric: R, magnetic: R) -> LayerEnergy<R>
where
    R: ScalarAlgebra,
{
    let total = electric.add(&magnetic);

    LayerEnergy::new(electric, magnetic, total)
}

impl<A> BrillouinBilinearLayerInput<A> {
    fn into_qnm_normalisation(
        self,
        vacuum_angular_wavenumber: &A,
        parallel_angular_wavenumber: &A,
    ) -> BilinearLayerNormalization<A>
    where
        A: ScalarAlgebra,
    {
        let (layer, derivative) = self.into_parts();

        let (state_products, quantities) = layer.into_parts();

        /*
         * These are the unweighted, unconjugated field products:
         *
         * electric = ∫ E · E dz
         * magnetic = ∫ H · H dz
         */
        let field_overlap = project_integrated_bilinear_field_overlap(
            &state_products,
            &quantities,
            &quantities,
            vacuum_angular_wavenumber,
            vacuum_angular_wavenumber,
            parallel_angular_wavenumber,
            parallel_angular_wavenumber,
        );

        let (electric, magnetic) = field_overlap.weighted(brillouin_weights(
            vacuum_angular_wavenumber,
            &quantities,
            &derivative,
        ));

        BilinearLayerNormalization::new(electric, magnetic)
    }
}

impl<A> Layers<BrillouinLayerInput<A>> {
    pub(crate) fn into_brillouin_energy(
        self,
        vacuum_angular_wavenumber: &A,
        parallel_angular_wavenumber: &A,
    ) -> Layers<LayerEnergy<A::RealJet>>
    where
        A: RealScalarAlgebra,
        A::RealJet: ScalarAlgebra,
        <A::RealJet as Jet>::Scalar: FromPrimitive + One,
    {
        let normalization = energy_density_prefactor(vacuum_angular_wavenumber);

        self.map(|each| {
            each.into_brillouin_energy(
                vacuum_angular_wavenumber,
                parallel_angular_wavenumber,
                &normalization,
            )
        })
    }
}

impl<A> Layers<BrillouinBilinearLayerInput<A>> {
    pub(crate) fn into_qnm_normalisation(
        self,
        vacuum_angular_wavenumber: &A,
        parallel_angular_wavenumber: &A,
    ) -> Layers<BilinearLayerNormalization<A>>
    where
        A: ScalarAlgebra,
    {
        self.map(|each| {
            each.into_qnm_normalisation(vacuum_angular_wavenumber, parallel_angular_wavenumber)
        })
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use ndarray::{Ix0, arr0};
    use num_complex::Complex64;

    use super::*;

    use crate::{
        Polarisation,
        algebra::{ArrayJet0, ArrayJet1, ComplexJet, Jet0, RealParameter},
        backend::LayerQuantities,
        material::DifferentiableMaterialHandle,
        observable::layer::{
            IntegratedHermitianCrossStateProducts, Layers, project::IntegratedLayerData,
        },
    };

    type C = Complex64;

    type A0 = ArrayJet0<C, Ix0, RealParameter>;
    type R0 = <A0 as ComplexJet>::RealJet;

    type A1 = ArrayJet1<C, Ix0, RealParameter>;
    type R1 = <A1 as ComplexJet>::RealJet;

    const TOLERANCE: f64 = 1.0e-12;

    fn c(real: f64, imaginary: f64) -> C {
        C::new(real, imaginary)
    }

    fn jet(value: C) -> A0 {
        Jet0::new(arr0(value))
    }

    fn real_jet(value: f64) -> R0 {
        Jet0::new(arr0(value))
    }

    fn scalar(value: &R0) -> f64 {
        value.value()[()]
    }

    fn jet1(value: C, first: C) -> A1 {
        A1::from_parts(arr0(value), arr0(first))
    }

    fn scalar1_value(value: &R1) -> f64 {
        value.value()[()]
    }

    fn scalar1_first(value: &R1) -> f64 {
        value.first()[()]
    }

    fn state_products(
        field_field: f64,
        secondary_secondary: f64,
    ) -> IntegratedHermitianCrossStateProducts<A0> {
        IntegratedHermitianCrossStateProducts::new(
            jet(c(field_field, 0.0)),
            jet(c(secondary_secondary, 0.0)),
            jet(c(0.0, 0.0)),
            jet(c(0.0, 0.0)),
        )
    }

    fn quantities(polarisation: Polarisation, epsilon: C, mu: C) -> LayerQuantities<A0> {
        LayerQuantities::test_fixture(jet(c(3.0, 0.0)), jet(epsilon), jet(mu), polarisation)
    }

    fn integrated_layer(
        polarisation: Polarisation,
        epsilon: C,
        mu: C,
        field_field: f64,
        secondary_secondary: f64,
    ) -> IntegratedLayerData<IntegratedHermitianCrossStateProducts<A0>, A0> {
        IntegratedLayerData::new(
            state_products(field_field, secondary_secondary),
            quantities(polarisation, epsilon, mu),
        )
    }

    fn derivatives(epsilon_first: C, mu_first: C) -> PrincipalConstitutive<A0> {
        PrincipalConstitutive::Isotropic {
            epsilon: jet(epsilon_first),
            mu: jet(mu_first),
        }
    }

    fn brillouin_layer(
        polarisation: Polarisation,
        epsilon: C,
        mu: C,
        epsilon_first: C,
        mu_first: C,
        field_field: f64,
        secondary_secondary: f64,
    ) -> BrillouinLayerInput<A0> {
        BrillouinLayerInput::new(
            integrated_layer(polarisation, epsilon, mu, field_field, secondary_secondary),
            derivatives(epsilon_first, mu_first),
        )
    }

    #[test]
    fn layer_energy_preserves_component_order() {
        let energy = LayerEnergy::new(1, 2, 3);

        assert_eq!(energy.electric(), &1);
        assert_eq!(energy.magnetic(), &2);
        assert_eq!(energy.total(), &3);

        assert_eq!(energy.into_parts(), (1, 2, 3),);
    }

    #[test]
    fn layer_energy_map_transforms_all_components() {
        let energy = LayerEnergy::new(1, 2, 3);

        let mapped = energy.map(|value| value * 10);

        assert_eq!(mapped.electric(), &10);
        assert_eq!(mapped.magnetic(), &20);
        assert_eq!(mapped.total(), &30);
    }

    #[test]
    fn layer_energy_map_supports_non_clone_storage() {
        #[derive(Debug, PartialEq)]
        struct NonClone(i32);

        let energy = LayerEnergy::new(NonClone(1), NonClone(2), NonClone(3));

        let mapped = energy.map(|value| value.0 * 10);

        assert_eq!(mapped.electric(), &10);
        assert_eq!(mapped.magnetic(), &20);
        assert_eq!(mapped.total(), &30);
    }

    #[test]
    fn energy_density_prefactor_is_one_quarter() {
        let prefactor = energy_density_prefactor(&jet(c(4.0, 0.0)));

        assert_relative_eq!(
            scalar(&prefactor),
            0.25,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn energy_density_prefactor_has_zero_outer_derivative() {
        let vacuum = jet1(c(4.0, 0.0), c(3.0, 0.0));

        let prefactor = energy_density_prefactor(&vacuum);

        assert_relative_eq!(
            scalar1_value(&prefactor),
            0.25,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );

        assert_relative_eq!(
            scalar1_first(&prefactor),
            0.0,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn project_energy_total_is_exact_component_sum() {
        let energy = project_energy(real_jet(2.5), real_jet(3.5));

        assert_relative_eq!(
            scalar(energy.total()),
            scalar(energy.electric()) + scalar(energy.magnetic()),
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn nondispersive_energy_uses_real_constitutive_weights() {
        let layer = integrated_layer(
            Polarisation::TransverseElectric,
            c(2.0, 100.0),
            c(3.0, -200.0),
            5.0,
            7.0,
        );

        let prefactor = real_jet(0.25);

        let energy =
            layer.into_nondispersive_energy(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)), &prefactor);

        let electric_norm = 5.0;

        let magnetic_norm = 7.0 / 4.0 + 5.0 * (0.6 * 0.6 / (4.0 * (3.0 * 3.0 + 200.0 * 200.0)));

        let expected_electric = 0.25 * electric_norm * 2.0;

        let expected_magnetic = 0.25 * magnetic_norm * 3.0;

        assert_relative_eq!(
            scalar(energy.electric()),
            expected_electric,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );

        assert_relative_eq!(
            scalar(energy.magnetic()),
            expected_magnetic,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );

        assert_relative_eq!(
            scalar(energy.total()),
            expected_electric + expected_magnetic,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn nondispersive_energy_coefficients_use_only_real_constitutive_parts() {
        let prefactor = real_jet(0.25);

        let first = integrated_layer(
            Polarisation::TransverseElectric,
            c(2.0, 100.0),
            c(3.0, -200.0),
            5.0,
            7.0,
        )
        .into_nondispersive_energy(&jet(c(2.0, 0.0)), &jet(c(0.0, 0.0)), &prefactor);

        let second = integrated_layer(
            Polarisation::TransverseElectric,
            c(2.0, -17.0),
            c(3.0, 29.0),
            5.0,
            7.0,
        )
        .into_nondispersive_energy(&jet(c(2.0, 0.0)), &jet(c(0.0, 0.0)), &prefactor);

        assert_eq!(first, second);
    }

    #[test]
    fn nondispersive_sequence_preserves_count_and_order() {
        let layers = Layers::new(vec![
            integrated_layer(
                Polarisation::TransverseElectric,
                c(2.0, 0.0),
                c(3.0, 0.0),
                5.0,
                7.0,
            ),
            integrated_layer(
                Polarisation::TransverseElectric,
                c(4.0, 0.0),
                c(3.0, 0.0),
                5.0,
                7.0,
            ),
        ]);

        let energy = layers.into_nondispersive_energy(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)));

        assert_eq!(energy.len(), 2);

        assert!(
            scalar(energy.first().unwrap().electric(),)
                < scalar(energy.last().unwrap().electric(),),
            "larger epsilon must remain associated with the second layer",
        );
    }

    #[test]
    fn empty_nondispersive_sequence_remains_empty() {
        let layers: Layers<IntegratedLayerData<IntegratedHermitianCrossStateProducts<A0>, A0>> =
            Layers::new(Vec::new());

        let energy = layers.into_nondispersive_energy(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)));

        assert!(energy.is_empty());
    }

    #[test]
    fn brillouin_input_preserves_layer_and_derivative_data() {
        let integrated = integrated_layer(
            Polarisation::TransverseElectric,
            c(2.0, 0.0),
            c(3.0, 0.0),
            5.0,
            7.0,
        );

        let derivatives = derivatives(c(11.0, 0.0), c(13.0, 0.0));

        let input = BrillouinLayerInput::new(integrated, derivatives);

        assert_eq!(
            input.derivative().diagonal()[0][0].value()[()],
            c(11.0, 0.0),
        );

        assert_eq!(
            input.derivative().diagonal()[1][0].value()[()],
            c(13.0, 0.0),
        );
    }

    #[test]
    fn brillouin_layer_projects_total_as_component_sum() {
        let input = brillouin_layer(
            Polarisation::TransverseElectric,
            c(2.0, 0.0),
            c(3.0, 0.0),
            c(5.0, 0.0),
            c(7.0, 0.0),
            5.0,
            7.0,
        );

        let prefactor = real_jet(0.25);

        let energy = input.into_brillouin_energy(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)), &prefactor);

        assert_relative_eq!(
            scalar(energy.total()),
            scalar(energy.electric()) + scalar(energy.magnetic()),
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn zero_intrinsic_derivatives_match_nondispersive_energy() {
        let prefactor = real_jet(0.25);

        let nondispersive = integrated_layer(
            Polarisation::TransverseMagnetic,
            c(2.0, 0.0),
            c(3.0, 0.0),
            5.0,
            7.0,
        )
        .into_nondispersive_energy(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)), &prefactor);

        let brillouin = brillouin_layer(
            Polarisation::TransverseMagnetic,
            c(2.0, 0.0),
            c(3.0, 0.0),
            c(0.0, 0.0),
            c(0.0, 0.0),
            5.0,
            7.0,
        )
        .into_brillouin_energy(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)), &prefactor);

        assert_eq!(brillouin, nondispersive,);
    }

    #[test]
    fn brillouin_sequence_preserves_count_and_order() {
        let layers = Layers::new(vec![
            brillouin_layer(
                Polarisation::TransverseElectric,
                c(2.0, 0.0),
                c(3.0, 0.0),
                c(1.0, 0.0),
                c(0.0, 0.0),
                5.0,
                7.0,
            ),
            brillouin_layer(
                Polarisation::TransverseElectric,
                c(2.0, 0.0),
                c(3.0, 0.0),
                c(4.0, 0.0),
                c(0.0, 0.0),
                5.0,
                7.0,
            ),
        ]);

        let energy = layers.into_brillouin_energy(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)));

        assert_eq!(energy.len(), 2);

        assert!(
            scalar(energy.first().unwrap().electric(),)
                < scalar(energy.last().unwrap().electric(),),
            "larger epsilon derivative must remain in the second layer",
        );
    }

    #[test]
    fn empty_brillouin_sequence_remains_empty() {
        let layers: Layers<BrillouinLayerInput<A0>> = Layers::new(Vec::new());

        let energy = layers.into_brillouin_energy(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)));

        assert!(energy.is_empty());
    }

    use crate::material::{DerivativeOrder, DifferentiableMaterial, Material, Sampled};

    #[derive(Clone, Copy, Debug)]
    struct LinearDifferentiableMaterial {
        epsilon_slope: f64,
        mu_slope: f64,
    }

    impl Material for LinearDifferentiableMaterial {
        type Real = f64;

        fn relative_permittivity<I, X>(&self, vacuum_wavenumber: I) -> I::Mapped<X>
        where
            I: Sampled<Elem = Self::Real>,
            X: ComplexScalar<RealField = f64>,
        {
            vacuum_wavenumber.map(|k| X::from_real(1.0 + self.epsilon_slope * k))
        }

        fn relative_permeability<I, X>(&self, vacuum_wavenumber: I) -> I::Mapped<X>
        where
            I: Sampled<Elem = Self::Real>,
            X: ComplexScalar<RealField = f64>,
        {
            vacuum_wavenumber.map(|k| X::from_real(1.0 + self.mu_slope * k))
        }
    }

    impl DifferentiableMaterial for LinearDifferentiableMaterial {
        fn relative_permittivity_derivative<I, X>(
            &self,
            vacuum_wavenumber: I,
            order: DerivativeOrder,
        ) -> I::Mapped<X>
        where
            I: Sampled<Elem = Self::Real>,
            X: ComplexScalar<RealField = f64>,
        {
            match order {
                DerivativeOrder::First => {
                    vacuum_wavenumber.map(|_| X::from_real(self.epsilon_slope))
                }

                DerivativeOrder::Second | DerivativeOrder::Third => {
                    vacuum_wavenumber.map(|_| X::from_real(0.0))
                }
            }
        }

        fn relative_permeability_derivative<I, X>(
            &self,
            vacuum_wavenumber: I,
            order: DerivativeOrder,
        ) -> I::Mapped<X>
        where
            I: Sampled<Elem = Self::Real>,
            X: ComplexScalar<RealField = f64>,
        {
            match order {
                DerivativeOrder::First => vacuum_wavenumber.map(|_| X::from_real(self.mu_slope)),

                DerivativeOrder::Second | DerivativeOrder::Third => {
                    vacuum_wavenumber.map(|_| X::from_real(0.0))
                }
            }
        }
    }

    fn material_handle(epsilon_slope: f64, mu_slope: f64) -> DifferentiableMaterialHandle<C> {
        DifferentiableMaterialHandle::new(LinearDifferentiableMaterial {
            epsilon_slope,
            mu_slope,
        })
    }

    #[test]
    fn brillouin_pairing_rejects_too_few_materials() {
        let layers = Layers::new(vec![
            integrated_layer(
                Polarisation::TransverseElectric,
                c(2.0, 0.0),
                c(3.0, 0.0),
                5.0,
                7.0,
            ),
            integrated_layer(
                Polarisation::TransverseElectric,
                c(2.0, 0.0),
                c(3.0, 0.0),
                5.0,
                7.0,
            ),
        ]);

        let materials = [material_handle(2.0, 3.0)];

        let error = layers
            .into_brillouin_layers::<crate::domain::RealAxis, _>(
                materials.iter(),
                &jet(c(2.0, 0.0)),
            )
            .expect_err("each integrated layer requires one material");

        assert_eq!(
            error,
            LayerEnergyError::MaterialCountMismatch {
                layer_count: 2,
                material_count: 1,
            },
        );
    }

    #[test]
    fn brillouin_pairing_rejects_too_many_materials() {
        let layers = Layers::new(vec![integrated_layer(
            Polarisation::TransverseElectric,
            c(2.0, 0.0),
            c(3.0, 0.0),
            5.0,
            7.0,
        )]);

        let materials = [material_handle(2.0, 3.0), material_handle(5.0, 7.0)];

        let error = layers
            .into_brillouin_layers::<crate::domain::RealAxis, _>(
                materials.iter(),
                &jet(c(2.0, 0.0)),
            )
            .expect_err("each integrated layer requires one material");

        assert_eq!(
            error,
            LayerEnergyError::MaterialCountMismatch {
                layer_count: 1,
                material_count: 2,
            },
        );
    }

    #[test]
    fn brillouin_pairing_preserves_material_order() {
        let layers = Layers::new(vec![
            integrated_layer(
                Polarisation::TransverseElectric,
                c(2.0, 0.0),
                c(3.0, 0.0),
                5.0,
                7.0,
            ),
            integrated_layer(
                Polarisation::TransverseElectric,
                c(2.0, 0.0),
                c(3.0, 0.0),
                5.0,
                7.0,
            ),
        ]);

        let materials = [material_handle(2.0, 3.0), material_handle(5.0, 7.0)];

        let paired = layers
            .into_brillouin_layers::<crate::domain::RealAxis, _>(
                materials.iter(),
                &jet(c(11.0, 0.0)),
            )
            .unwrap();

        assert_eq!(paired.len(), 2);

        assert_eq!(
            paired.first().unwrap().derivative().diagonal()[0][0].value()[()],
            c(2.0, 0.0),
        );

        assert_eq!(
            paired.first().unwrap().derivative().diagonal()[1][0].value()[()],
            c(3.0, 0.0),
        );

        assert_eq!(
            paired.last().unwrap().derivative().diagonal()[0][0].value()[()],
            c(5.0, 0.0),
        );

        assert_eq!(
            paired.last().unwrap().derivative().diagonal()[1][0].value()[()],
            c(7.0, 0.0),
        );
    }

    #[test]
    fn empty_layers_pair_with_empty_materials() {
        let layers: Layers<IntegratedLayerData<IntegratedHermitianCrossStateProducts<A0>, A0>> =
            Layers::new(Vec::new());

        let materials: [DifferentiableMaterialHandle<C>; 0] = [];

        let paired = layers
            .into_brillouin_layers::<crate::domain::RealAxis, _>(
                materials.iter(),
                &jet(c(2.0, 0.0)),
            )
            .unwrap();

        assert!(paired.is_empty());
    }
}
