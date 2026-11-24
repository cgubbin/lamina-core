//! Narrow supported-backend characteristic evaluation for local derivatives.
use super::{
    ComplexPlaneEvaluator, ComplexPlaneState, OwnedModeFields,
    mode::{ExteriorWaves, reconstruct_mode_fields},
};
use crate::{
    AnalyticalMaterialHandle, AnalyticalMaterialStack, CanonicalCoordinates,
    EvaluateMeromorphicMaterial, ExteriorWavevectors, ModeJet1, PlaneWaveDeterminant, Polarisation,
    Scalar, ScalarAlgebra,
    backend::{Backend, Scatter2, Transfer2},
    input::{CanonicalStack, canonical::CanonicalLayer},
    material::{ConstitutiveEvaluator, ConstitutiveLift},
};
use crate::{ElectromagneticFields, FieldSampling, SpatialResponse, backend::ModalSolutionSource};
use ndarray::{Array0, Ix0, arr0};
use num_complex::Complex64 as C;
use std::{error::Error, sync::Arc};

type J = ModeJet1<C, Ix0>;

/// Independently seeded physical inputs for a first-order modal direction.
#[derive(Clone, Copy)]
pub struct ModalCharacteristicDirection<'a> {
    pub stack: &'a AnalyticalMaterialStack<C>,
    pub coordinates: &'a CanonicalCoordinates<J>,
    pub exterior: &'a ExteriorWavevectors<J>,
    pub polarisation: Polarisation,
    /// Layer thickness rates, in centimetres per parameter unit.
    pub thickness_rates: &'a [f64],
    /// Total relative permittivity and permeability rates in left-to-right order.
    pub constitutive_rates: &'a [(C, C)],
}

#[derive(Clone)]
struct SeededMaterial {
    epsilon: J,
    mu: J,
}

struct SeededConstitutive;

impl ConstitutiveEvaluator<C, Ix0, SeededMaterial> for SeededConstitutive {
    fn relative_permittivity(material: &SeededMaterial, _: &Array0<C>) -> Array0<C> {
        material.epsilon.value().clone()
    }

    fn relative_permeability(material: &SeededMaterial, _: &Array0<C>) -> Array0<C> {
        material.mu.value().clone()
    }
}

impl ConstitutiveLift<SeededConstitutive, SeededMaterial> for J {
    fn relative_permittivity(material: &SeededMaterial, _: &Self) -> Self {
        material.epsilon.clone()
    }

    fn relative_permeability(material: &SeededMaterial, _: &Self) -> Self {
        material.mu.clone()
    }
}

impl ModalCharacteristicDirection<'_> {
    fn canonical_stack(&self) -> Result<CanonicalStack<SeededMaterial, J>, std::io::Error> {
        let invalid = |message| std::io::Error::new(std::io::ErrorKind::InvalidInput, message);
        if !self.coordinates.vacuum_angular_wavenumber().all_finite()
            || !self.coordinates.parallel_angular_wavenumber().all_finite()
            || !self.exterior.left().all_finite()
            || !self.exterior.right().all_finite()
        {
            return Err(invalid("nonfinite coordinate jet"));
        }
        let layers = self.stack.layers_left_to_right();
        if self.thickness_rates.len() != layers.len()
            || self.constitutive_rates.len() != layers.len() + 2
        {
            return Err(invalid("direction does not match stack topology"));
        }
        let k0 = Scalar(self.coordinates.vacuum_angular_wavenumber().value()[()]);
        let seed = |material: &AnalyticalMaterialHandle<C>, (epsilon_rate, mu_rate)| {
            let epsilon = J::from_parts(
                arr0(material.evaluate_relative_permittivity_complex(k0)),
                arr0(epsilon_rate),
            );
            let mu = J::from_parts(
                arr0(material.evaluate_relative_permeability_complex(k0)),
                arr0(mu_rate),
            );
            if !epsilon.all_finite() || !mu.all_finite() {
                return Err(invalid("nonfinite constitutive jet"));
            }
            Ok(SeededMaterial { epsilon, mu })
        };
        let left = seed(self.stack.left_exterior(), self.constitutive_rates[0])?;
        let right = seed(
            self.stack.right_exterior(),
            self.constitutive_rates[layers.len() + 1],
        )?;
        let layers = layers
            .iter()
            .zip(self.thickness_rates)
            .zip(&self.constitutive_rates[1..])
            .map(|((layer, &rate), &constitutive_rate)| {
                let thickness = layer.thickness().as_centimetres();
                if !thickness.is_finite() || thickness <= 0.0 || !rate.is_finite() {
                    return Err(invalid("invalid layer thickness jet"));
                }
                let material = seed(layer.material(), constitutive_rate)?;
                Ok(CanonicalLayer::new(
                    material,
                    J::from_parts(arr0(C::from(thickness)), arr0(C::from(rate))),
                ))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(CanonicalStack::new(left, right, layers))
    }
}

/// Unscaled homogeneous fields and their first derivatives at the same samples.
#[derive(Clone, Debug)]
pub struct DirectionalModeFields {
    pub fields: OwnedModeFields,
    pub tangent: OwnedModeFields,
}

/// Checked characteristic evaluation on a caller-supplied canonical stack.
///
/// Implemented only by the built-in isotropic modal backends. This keeps
/// derivative postprocessing explicit about its numerical backend while
/// allowing a selected finite-layer thickness jet to reuse core compilation.
pub trait ModalDerivativeBackend: Clone {
    /// Evaluate the characteristic with independently seeded material and coordinate jets.
    fn directional_characteristic(
        &self,
        direction: ModalCharacteristicDirection<'_>,
    ) -> Result<J, Arc<dyn Error + Send + Sync>>;

    /// Recover field jets at material-attached coordinates. The backend's
    /// homogeneous seed is used as the gauge, multiplied by a constant amplitude.
    fn directional_fields(
        &self,
        direction: ModalCharacteristicDirection<'_>,
        sampling: &FieldSampling<f64>,
        amplitude: C,
    ) -> Result<DirectionalModeFields, Arc<dyn Error + Send + Sync>>;

    fn checked_characteristic(
        &self,
        stack: crate::CanonicalStack<AnalyticalMaterialHandle<C>, J>,
        coordinates: &CanonicalCoordinates<J>,
        exterior: &ExteriorWavevectors<J>,
        polarisation: Polarisation,
    ) -> Result<J, Arc<dyn Error + Send + Sync>>;
}

macro_rules! implementation {
    ($backend:ty) => {
        impl ModalDerivativeBackend for $backend {
            fn directional_characteristic(
                &self,
                direction: ModalCharacteristicDirection<'_>,
            ) -> Result<J, Arc<dyn Error + Send + Sync>> {
                let stack = direction.canonical_stack().map_err(super::shared_error)?;
                let determinant = <Self as Backend<J, SeededConstitutive>>::solve(
                    self,
                    direction.coordinates,
                    &stack,
                    direction.exterior,
                    direction.polarisation,
                )
                .map_err(super::shared_error)?
                .determinant()
                .into_inner();
                if !determinant.all_finite() {
                    return Err(Arc::new(std::io::Error::other(
                        "nonfinite directional characteristic",
                    )));
                }
                Ok(determinant)
            }

            fn directional_fields(
                &self,
                direction: ModalCharacteristicDirection<'_>,
                sampling: &FieldSampling<f64>,
                amplitude: C,
            ) -> Result<DirectionalModeFields, Arc<dyn Error + Send + Sync>> {
                if !amplitude.re.is_finite() || !amplitude.im.is_finite() || amplitude == C::ZERO {
                    return Err(Arc::new(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "invalid field amplitude",
                    )));
                }
                let stack = direction.canonical_stack().map_err(super::shared_error)?;
                let workspace = <Self as Backend<J, SeededConstitutive>>::retain(
                    self,
                    direction.coordinates,
                    &stack,
                    direction.exterior,
                    direction.polarisation,
                )
                .map_err(super::shared_error)?;
                let state = ComplexPlaneState::new(
                    direction.coordinates.clone(),
                    direction.exterior.clone(),
                    &stack,
                    workspace,
                    direction.polarisation,
                );
                let seed = state
                    .workspace()
                    .modal_boundary_solution()
                    .map_err(super::shared_error)?
                    .scaled(&J::constant(arr0(amplitude)));
                let resolved = sampling
                    .resolve(direction.stack)
                    .map_err(super::shared_error)?;
                let (electric, magnetic) =
                    reconstruct_mode_fields(&state, &seed, &resolved, ExteriorWaves::Outgoing)
                        .map_err(super::shared_error)?
                        .into_parts();
                let primal =
                    ElectromagneticFields::new(electric.value().clone(), magnetic.value().clone());
                let tangent =
                    ElectromagneticFields::new(electric.first().clone(), magnetic.first().clone());
                for fields in [&primal, &tangent] {
                    let electric = fields.electric();
                    let magnetic = fields.magnetic();
                    if [
                        electric.x(),
                        electric.y(),
                        electric.z(),
                        magnetic.x(),
                        magnetic.y(),
                        magnetic.z(),
                    ]
                    .into_iter()
                    .flatten()
                    .any(|v| !v.re.is_finite() || !v.im.is_finite())
                    {
                        return Err(Arc::new(std::io::Error::other(
                            "nonfinite directional fields",
                        )));
                    }
                }
                Ok(DirectionalModeFields {
                    fields: SpatialResponse::new(primal, resolved.clone()),
                    tangent: SpatialResponse::new(tangent, resolved),
                })
            }

            fn checked_characteristic(
                &self,
                stack: CanonicalStack<AnalyticalMaterialHandle<C>, J>,
                coordinates: &CanonicalCoordinates<J>,
                exterior: &ExteriorWavevectors<J>,
                polarisation: Polarisation,
            ) -> Result<J, Arc<dyn Error + Send + Sync>> {
                let evaluator = ComplexPlaneEvaluator::from_canonical_stack(stack, self.clone());
                evaluator
                    .checked_determinant(coordinates, exterior, polarisation)
                    .map(PlaneWaveDeterminant::into_inner)
                    .map_err(|error| Arc::new(error) as Arc<dyn Error + Send + Sync>)
            }
        }
    };
}

implementation!(Scatter2);
implementation!(Transfer2);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Constant, Stack};
    use lamina_units::Length;

    #[test]
    fn independent_permittivity_and_permeability_jets_match_perturbed_characteristics() {
        let constant = |value| J::from_parts(arr0(C::from(value)), arr0(C::ZERO));
        let coordinates = CanonicalCoordinates::new(constant(2.0), constant(1.0));
        let exterior = ExteriorWavevectors::new(constant(3.0f64.sqrt()), constant(3.0f64.sqrt()));
        let scatter = Scatter2::new();
        for polarisation in [
            Polarisation::TransverseElectric,
            Polarisation::TransverseMagnetic,
        ] {
            let evaluate = |step: f64, rates: (C, C)| {
                let stack =
                    Stack::from_analytical_materials(Constant::vacuum(), Constant::vacuum())
                        .analytical_layer(
                            Constant::new(2.0 + step * rates.0.re, 1.2 + step * rates.1.re),
                            Length::centimetres(0.3 + step * 0.07),
                        )
                        .finalise();
                let input = ModalCharacteristicDirection {
                    stack: &stack,
                    coordinates: &coordinates,
                    exterior: &exterior,
                    polarisation,
                    thickness_rates: &[0.07],
                    constitutive_rates: &[(C::ZERO, C::ZERO), rates, (C::ZERO, C::ZERO)],
                };
                let sampling = FieldSampling::new().layer(
                    0,
                    crate::LayerSampling::uniform_with_boundaries(3, true, true),
                );
                (
                    scatter.directional_characteristic(input).unwrap(),
                    scatter
                        .directional_fields(input, &sampling, C::ONE)
                        .unwrap(),
                )
            };
            for rates in [(C::ONE, C::ZERO), (C::ZERO, C::ONE)] {
                let h = 1e-5;
                let (plus, plus_fields) = evaluate(h, rates);
                let (minus, minus_fields) = evaluate(-h, rates);
                let (base, base_fields) = evaluate(0.0, rates);
                let finite_difference = (plus.value()[()] - minus.value()[()]) / (2.0 * h);
                let derivative = base.first()[()];
                assert!(derivative.norm() > 1e-4);
                assert!((derivative - finite_difference).norm() < 1e-9);
                fn components(fields: &OwnedModeFields) -> [&ndarray::Array1<C>; 6] {
                    let e = fields.quantity().electric();
                    let h = fields.quantity().magnetic();
                    [e.x(), e.y(), e.z(), h.x(), h.y(), h.z()]
                }
                for ((plus, minus), tangent) in components(&plus_fields.fields)
                    .into_iter()
                    .zip(components(&minus_fields.fields))
                    .zip(components(&base_fields.tangent))
                {
                    for ((plus, minus), tangent) in plus.iter().zip(minus).zip(tangent) {
                        let fd = (*plus - *minus) / (2.0 * h);
                        assert!(
                            (fd - *tangent).norm() < 1e-8,
                            "{polarisation:?} {rates:?} fd={fd}, jet={tangent}"
                        );
                    }
                }
            }
        }
    }
}
