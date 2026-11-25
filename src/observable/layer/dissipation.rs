//! Analytically integrated dissipation in homogeneous finite layers.

use num_traits::{FromPrimitive, One};

use crate::{
    algebra::{Jet, RealScalarAlgebra, ScalarAlgebra},
    observable::layer::IntegratedHermitianCrossStateProducts,
};

use super::{Layers, integration::project_integrated_field_norms, project::IntegratedLayerData};

/// Normalized power dissipated within one finite layer.
///
/// The components are spatial integrals over the complete layer and are
/// normalized by the magnitude of the incident-wave power flux.
///
/// For isotropic media:
///
/// ```text
/// total = electric + magnetic
/// ```
///
/// Passive media should produce a non-negative `total`, subject to numerical
/// error. Active media may produce a negative value.
#[derive(Clone, Debug, PartialEq)]
pub struct LayerDissipation<R> {
    electric: R,
    magnetic: R,
    total: R,
}

impl<R> LayerDissipation<R> {
    pub(crate) const fn new(electric: R, magnetic: R, total: R) -> Self {
        Self {
            electric,
            magnetic,
            total,
        }
    }

    /// Return the integrated electric dissipation.
    pub fn electric(&self) -> &R {
        &self.electric
    }

    /// Return the integrated magnetic dissipation.
    pub fn magnetic(&self) -> &R {
        &self.magnetic
    }

    /// Return the total integrated dissipation.
    pub fn total(&self) -> &R {
        &self.total
    }

    /// Consume the result and return `(electric, magnetic, total)`.
    pub fn into_parts(self) -> (R, R, R) {
        (self.electric, self.magnetic, self.total)
    }

    /// Transform every component.
    pub fn map<U>(self, mut map: impl FnMut(R) -> U) -> LayerDissipation<U> {
        LayerDissipation {
            electric: map(self.electric),
            magnetic: map(self.magnetic),
            total: map(self.total),
        }
    }
}

impl<A> IntegratedLayerData<IntegratedHermitianCrossStateProducts<A>, A> {
    /// Project this integrated layer into normalized electric, magnetic, and
    /// total dissipation.
    fn into_dissipation(
        self,
        vacuum_angular_wavenumber: &A,
        parallel_angular_wavenumber: &A,
    ) -> LayerDissipation<A::RealJet>
    where
        A: RealScalarAlgebra,
        A::RealJet: ScalarAlgebra,
        <A::RealJet as Jet>::Scalar: FromPrimitive + One,
    {
        let (state, quantities) = self.into_parts();

        let half =
            <A::RealJet as Jet>::Scalar::from_f64(0.5).expect("one half must be representable");

        let common = vacuum_angular_wavenumber.real().scale(half);

        let (electric, magnetic) = project_integrated_field_norms(
            &state,
            &quantities,
            vacuum_angular_wavenumber,
            parallel_angular_wavenumber,
        )
        .weighted(
            quantities
                .channel()
                .map(|entry| common.multiply(&entry.imaginary())),
        );

        let total = electric.add(&magnetic);

        LayerDissipation::new(electric, magnetic, total)
    }
}

impl<A> Layers<IntegratedLayerData<IntegratedHermitianCrossStateProducts<A>, A>> {
    /// Project every integrated finite layer into normalized dissipation.
    ///
    /// Results preserve physical left-to-right finite-layer order.
    pub(crate) fn into_dissipation(
        self,
        vacuum_angular_wavenumber: &A,
        parallel_angular_wavenumber: &A,
    ) -> Layers<LayerDissipation<A::RealJet>>
    where
        A: RealScalarAlgebra,
        A::RealJet: ScalarAlgebra,
        <A::RealJet as Jet>::Scalar: FromPrimitive + One,
    {
        self.map(|each| {
            each.into_dissipation(vacuum_angular_wavenumber, parallel_angular_wavenumber)
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
        algebra::{ArrayJet0, ComplexJet, Jet0, RealParameter},
        backend::LayerQuantities,
        observable::layer::{
            IntegratedHermitianCrossStateProducts, Layers, project::IntegratedLayerData,
        },
    };

    type C = Complex64;

    type A0 = ArrayJet0<C, Ix0, RealParameter>;
    type R0 = <A0 as ComplexJet>::RealJet;

    const TOLERANCE: f64 = 1.0e-12;

    fn c(real: f64, imaginary: f64) -> C {
        C::new(real, imaginary)
    }

    fn jet(value: C) -> A0 {
        Jet0::new(arr0(value))
    }

    fn scalar(value: &R0) -> f64 {
        value.value()[()]
    }

    fn quantities(polarisation: Polarisation, epsilon: C, mu: C) -> LayerQuantities<A0> {
        LayerQuantities::test_fixture(jet(c(3.0, 0.0)), jet(epsilon), jet(mu), polarisation)
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

    #[test]
    fn layer_dissipation_preserves_component_order() {
        let dissipation = LayerDissipation::new(1, 2, 3);

        assert_eq!(dissipation.electric(), &1);
        assert_eq!(dissipation.magnetic(), &2);
        assert_eq!(dissipation.total(), &3);

        assert_eq!(dissipation.into_parts(), (1, 2, 3),);
    }

    #[test]
    fn layer_dissipation_map_transforms_every_component() {
        let dissipation = LayerDissipation::new(1, 2, 3);

        let mapped = dissipation.map(|value| value * 10);

        assert_eq!(mapped.electric(), &10);
        assert_eq!(mapped.magnetic(), &20);
        assert_eq!(mapped.total(), &30);
    }

    #[test]
    fn integrated_layer_projects_total_as_component_sum() {
        let layer = integrated_layer(
            Polarisation::TransverseElectric,
            c(2.0, 0.4),
            c(3.0, 0.2),
            5.0,
            7.0,
        );

        let dissipation = layer.into_dissipation(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)));

        assert_relative_eq!(
            scalar(dissipation.total()),
            scalar(dissipation.electric()) + scalar(dissipation.magnetic()),
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn zero_material_loss_produces_zero_dissipation() {
        let layer = integrated_layer(
            Polarisation::TransverseMagnetic,
            c(2.0, 0.0),
            c(3.0, 0.0),
            5.0,
            7.0,
        );

        let dissipation = layer.into_dissipation(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)));

        assert_relative_eq!(
            scalar(dissipation.electric()),
            0.0,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );

        assert_relative_eq!(
            scalar(dissipation.magnetic()),
            0.0,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );

        assert_relative_eq!(
            scalar(dissipation.total()),
            0.0,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn layer_sequence_preserves_count_and_order() {
        let layers = Layers::new(vec![
            integrated_layer(
                Polarisation::TransverseElectric,
                c(2.0, 0.1),
                c(3.0, 0.0),
                5.0,
                7.0,
            ),
            integrated_layer(
                Polarisation::TransverseElectric,
                c(2.0, 0.3),
                c(3.0, 0.0),
                5.0,
                7.0,
            ),
        ]);

        let projected = layers.into_dissipation(&jet(c(2.0, 0.0)), &jet(c(0.6, 0.0)));

        assert_eq!(projected.len(), 2);

        let first = projected.first().unwrap();

        let second = projected.last().unwrap();

        assert!(
            scalar(first.total()) < scalar(second.total()),
            "larger electric loss should remain in the second layer",
        );
    }
}
