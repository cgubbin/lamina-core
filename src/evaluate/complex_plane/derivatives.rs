//! Narrow supported-backend characteristic evaluation for local derivatives.
use super::{
    ComplexPlaneEvaluator, ComplexPlaneState, OwnedModeFields,
    mode::{ExteriorWaves, reconstruct_mode_fields},
};
use crate::{
    AnalyticalMaterialHandle, CanonicalCoordinates, ComplexPlane, ConstitutiveClass,
    ConstitutiveError, ConstitutiveTensor, ExteriorWavevectors, FiniteLayerIndex, ModeJet1,
    PlaneWaveDeterminant, Polarisation, ScalarAlgebra, Stack, StackRegion, TensorConstitutiveLift,
    algebra::{ArrayJet0, HolomorphicParameter},
    backend::{Backend, Scatter2, Transfer2},
    input::{CanonicalStack, canonical::CanonicalLayer, compile_terminations},
};
use crate::{ElectromagneticFields, FieldSampling, SpatialResponse, backend::ModalSolutionSource};
use ndarray::{Ix0, arr0};
use num_complex::Complex64 as C;
use std::{error::Error, sync::Arc};

type J = ModeJet1<C, Ix0>;
type ValueJet = ArrayJet0<C, Ix0, HolomorphicParameter>;

/// Independently seeded physical inputs for a first-order modal direction.
pub struct ModalCharacteristicDirection<'a, M = AnalyticalMaterialHandle<C>, R = (C, C)> {
    pub stack: &'a Stack<M, f64>,
    pub coordinates: &'a CanonicalCoordinates<J>,
    pub exterior: &'a ExteriorWavevectors<J>,
    pub polarisation: Polarisation,
    /// Layer thickness rates, in centimetres per parameter unit.
    pub thickness_rates: &'a [f64],
    /// Rates of each region's relative permittivity and permeability, left
    /// exterior, layers left to right, then right exterior.
    pub constitutive_rates: &'a [R],
}

impl<M, R> Clone for ModalCharacteristicDirection<'_, M, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<M, R> Copy for ModalCharacteristicDirection<'_, M, R> {}

/// The rates of one region's relative permittivity and permeability along a
/// modal direction, in the stack frame: total `(ε', μ')` of an isotropic
/// region, or `[ε', μ']` entry by entry.
///
/// A material-frame rate `T'` enters as `R T' Rᵀ`
/// ([`Rotation::rotate`](crate::Rotation::rotate)).
pub trait ConstitutiveRates {
    fn stack_frame(&self) -> [ConstitutiveTensor<C>; 2];
}

impl ConstitutiveRates for (C, C) {
    fn stack_frame(&self) -> [ConstitutiveTensor<C>; 2] {
        [
            ConstitutiveTensor::Isotropic(self.0),
            ConstitutiveTensor::Isotropic(self.1),
        ]
    }
}

impl ConstitutiveRates for [ConstitutiveTensor<C>; 2] {
    fn stack_frame(&self) -> [ConstitutiveTensor<C>; 2] {
        self.clone()
    }
}

/// A region's stack-frame `[ε, μ]` seeded with their rates.
#[derive(Clone)]
struct SeededMaterial([ConstitutiveTensor<J>; 2]);

struct SeededConstitutive;

impl TensorConstitutiveLift<SeededConstitutive, SeededMaterial> for J {
    fn constitutive_tensors(
        _: StackRegion,
        material: &SeededMaterial,
        _: &Self,
    ) -> Result<[ConstitutiveTensor<Self>; 2], ConstitutiveError> {
        Ok(material.0.clone())
    }

    fn is_separable(_: &SeededMaterial) -> bool {
        true
    }
}

/// Seed one stack-frame tensor `value` with its stack-frame `rate`.
fn seed(
    region: StackRegion,
    value: ConstitutiveTensor<ValueJet>,
    rate: ConstitutiveTensor<C>,
) -> Result<ConstitutiveTensor<J>, ConstitutiveError> {
    let jet = |value: ValueJet, rate| J::from_parts(value.into_inner(), arr0(rate));
    if let (ConstitutiveTensor::Isotropic(value), ConstitutiveTensor::Isotropic(rate)) =
        (&value, &rate)
    {
        return Ok(ConstitutiveTensor::Isotropic(jet(value.clone(), *rate)));
    }
    let unsupported = ConstitutiveError::Unsupported {
        region,
        class: ConstitutiveClass::Full,
    };
    let [x, y, z] = value.into_diagonal().ok_or(unsupported)?;
    let [rate_x, rate_y, rate_z] = rate.into_diagonal().ok_or(unsupported)?;
    Ok(ConstitutiveTensor::Diagonal([
        jet(x, rate_x),
        jet(y, rate_y),
        jet(z, rate_z),
    ]))
}

impl<M, R> ModalCharacteristicDirection<'_, M, R>
where
    ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
    R: ConstitutiveRates,
{
    /// The seeded canonical stack, with walls and PMLs fixed in the stack
    /// frame: the right gap changes by `−Σ thickness_rates`. A direction that
    /// thickens a stack closed on its last interface fails with
    /// [`crate::StackCompileError::OutwardThicknessDerivative`].
    fn canonical_stack(
        &self,
    ) -> Result<CanonicalStack<SeededMaterial, J>, Arc<dyn Error + Send + Sync>> {
        let stack = self.open_canonical_stack()?;
        let thickening = self.thickness_rates.iter().sum::<f64>() > 0.0;
        let outward_layer = self
            .thickness_rates
            .iter()
            .position(|rate| *rate > 0.0)
            .filter(|_| thickening);
        let (left, right) = compile_terminations(self.stack, stack.layers(), Ix0(), outward_layer)
            .map_err(super::shared_error)?;
        Ok(stack.with_terminations(left, right))
    }

    fn open_canonical_stack(
        &self,
    ) -> Result<CanonicalStack<SeededMaterial, J>, Arc<dyn Error + Send + Sync>> {
        let invalid = |message| {
            super::shared_error(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                message,
            ))
        };
        if !self.coordinates.vacuum_angular_wavenumber().all_finite()
            || !self.coordinates.parallel_angular_wavenumber().all_finite()
            || !self.exterior.all(|kappa| kappa.all_finite())
        {
            return Err(invalid("nonfinite coordinate jet"));
        }
        let layers = self.stack.layers_left_to_right();
        if self.thickness_rates.len() != layers.len()
            || self.constitutive_rates.len() != layers.len() + 2
        {
            return Err(invalid("direction does not match stack topology"));
        }
        let k0 = ValueJet::new(arr0(
            self.coordinates.vacuum_angular_wavenumber().value()[()],
        ));
        let seeded = |region, material: &M, rates: &R| {
            let [epsilon, mu] =
                <ValueJet as TensorConstitutiveLift<ComplexPlane, M>>::constitutive_tensors(
                    region, material, &k0,
                )
                .map_err(super::shared_error)?;
            let [epsilon_rate, mu_rate] = rates.stack_frame();
            let tensors = [
                seed(region, epsilon, epsilon_rate).map_err(super::shared_error)?,
                seed(region, mu, mu_rate).map_err(super::shared_error)?,
            ];
            if !tensors
                .iter()
                .all(|tensor| tensor.entries().all(J::all_finite))
            {
                return Err(invalid("nonfinite constitutive jet"));
            }
            Ok(SeededMaterial(tensors))
        };
        let left = seeded(
            StackRegion::LeftExterior,
            self.stack.left_exterior(),
            &self.constitutive_rates[0],
        )?;
        let right = seeded(
            StackRegion::RightExterior,
            self.stack.right_exterior(),
            &self.constitutive_rates[layers.len() + 1],
        )?;
        let layers = layers
            .iter()
            .zip(self.thickness_rates)
            .zip(&self.constitutive_rates[1..])
            .enumerate()
            .map(|(index, ((layer, &rate), constitutive_rates))| {
                let thickness = layer.thickness().as_centimetres();
                if !thickness.is_finite() || thickness <= 0.0 || !rate.is_finite() {
                    return Err(invalid("invalid layer thickness jet"));
                }
                let material = seeded(
                    StackRegion::Layer(FiniteLayerIndex::new(index)),
                    layer.material(),
                    constitutive_rates,
                )?;
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
/// Implemented only by the built-in 2×2 modal backends. This keeps
/// derivative postprocessing explicit about its numerical backend while
/// allowing a selected finite-layer thickness jet to reuse core compilation.
pub trait ModalDerivativeBackend: Clone {
    /// Evaluate the characteristic with independently seeded material and coordinate jets.
    fn directional_characteristic<M, R>(
        &self,
        direction: ModalCharacteristicDirection<'_, M, R>,
    ) -> Result<J, Arc<dyn Error + Send + Sync>>
    where
        ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
        R: ConstitutiveRates;

    /// Recover field jets at material-attached coordinates. The backend's
    /// homogeneous seed is used as the gauge, multiplied by a constant amplitude.
    fn directional_fields<M, R>(
        &self,
        direction: ModalCharacteristicDirection<'_, M, R>,
        sampling: &FieldSampling<f64>,
        amplitude: C,
    ) -> Result<DirectionalModeFields, Arc<dyn Error + Send + Sync>>
    where
        ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
        R: ConstitutiveRates;

    fn checked_characteristic<M>(
        &self,
        stack: crate::CanonicalStack<M, J>,
        coordinates: &CanonicalCoordinates<J>,
        exterior: &ExteriorWavevectors<J>,
        polarisation: Polarisation,
    ) -> Result<J, Arc<dyn Error + Send + Sync>>
    where
        M: Clone,
        J: TensorConstitutiveLift<ComplexPlane, M>;
}

macro_rules! implementation {
    ($backend:ty) => {
        impl ModalDerivativeBackend for $backend {
            fn directional_characteristic<M, R>(
                &self,
                direction: ModalCharacteristicDirection<'_, M, R>,
            ) -> Result<J, Arc<dyn Error + Send + Sync>>
            where
                ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
                R: ConstitutiveRates,
            {
                let stack = direction.canonical_stack()?;
                let determinant = <Self as Backend<J, SeededConstitutive>>::solve(
                    self,
                    direction.coordinates,
                    &stack,
                    direction.exterior,
                    direction.polarisation.into(),
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

            fn directional_fields<M, R>(
                &self,
                direction: ModalCharacteristicDirection<'_, M, R>,
                sampling: &FieldSampling<f64>,
                amplitude: C,
            ) -> Result<DirectionalModeFields, Arc<dyn Error + Send + Sync>>
            where
                ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
                R: ConstitutiveRates,
            {
                if !amplitude.re.is_finite() || !amplitude.im.is_finite() || amplitude == C::ZERO {
                    return Err(Arc::new(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "invalid field amplitude",
                    )));
                }
                let stack = direction.canonical_stack()?;
                let workspace = <Self as Backend<J, SeededConstitutive>>::retain(
                    self,
                    direction.coordinates,
                    &stack,
                    direction.exterior,
                    direction.polarisation.into(),
                )
                .map_err(super::shared_error)?;
                let state = ComplexPlaneState::new(
                    direction.coordinates.clone(),
                    direction.exterior.clone(),
                    &stack,
                    workspace,
                    direction.polarisation.into(),
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

            fn checked_characteristic<M>(
                &self,
                stack: CanonicalStack<M, J>,
                coordinates: &CanonicalCoordinates<J>,
                exterior: &ExteriorWavevectors<J>,
                polarisation: Polarisation,
            ) -> Result<J, Arc<dyn Error + Send + Sync>>
            where
                M: Clone,
                J: TensorConstitutiveLift<ComplexPlane, M>,
            {
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

    #[test]
    fn tensor_entry_rates_match_perturbed_characteristics() {
        use crate::{
            AnalyticalTensorStack, Rotation,
            test_support::materials::{ComplexQuadratic, principal_region},
        };
        let constant = |value| J::from_parts(arr0(C::from(value)), arr0(C::ZERO));
        let coordinates = CanonicalCoordinates::new(constant(2.0), constant(1.0));
        let exterior = ExteriorWavevectors::new(constant(3.0f64.sqrt()), constant(3.0f64.sqrt()));
        let material = |epsilon: f64| ComplexQuadratic {
            epsilon: [C::new(epsilon, 0.05), C::ZERO, C::ZERO],
            mu: [C::ONE, C::ZERO, C::ZERO],
        };
        let zero = || {
            [
                ConstitutiveTensor::Isotropic(C::ZERO),
                ConstitutiveTensor::Isotropic(C::ZERO),
            ]
        };
        let sampling = FieldSampling::new().layer(
            0,
            crate::LayerSampling::uniform_with_boundaries(3, true, true),
        );
        // The optic axis along z and along y; the rate of the extraordinary
        // permittivity is a material-frame rate, rotated into the stack frame.
        for rotation in [
            Rotation::identity(),
            Rotation::from_optic_axis(std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2),
        ] {
            let stack = |step: f64| -> AnalyticalTensorStack<C> {
                Stack::builder(
                    AnalyticalMaterialHandle::new(Constant::vacuum()).into(),
                    AnalyticalMaterialHandle::new(Constant::vacuum()).into(),
                )
                .layer(
                    principal_region(
                        [material(2.2), material(2.2), material(3.5 + step)],
                        rotation,
                    ),
                    Length::centimetres(0.3),
                )
                .finalise()
            };
            let mut rate = [[C::ZERO; 3]; 3];
            rate[2][2] = C::ONE;
            let rate = rotation.rotate(&rate);
            let rates = [
                zero(),
                [
                    ConstitutiveTensor::Diagonal([rate[0][0], rate[1][1], rate[2][2]]),
                    ConstitutiveTensor::Isotropic(C::ZERO),
                ],
                zero(),
            ];
            for polarisation in [
                Polarisation::TransverseElectric,
                Polarisation::TransverseMagnetic,
            ] {
                let evaluate = |step: f64| {
                    let stack = stack(step);
                    let direction = ModalCharacteristicDirection {
                        stack: &stack,
                        coordinates: &coordinates,
                        exterior: &exterior,
                        polarisation,
                        thickness_rates: &[0.0],
                        constitutive_rates: &rates,
                    };
                    (
                        Scatter2::new()
                            .directional_characteristic(direction)
                            .unwrap(),
                        Transfer2::new()
                            .directional_fields(direction, &sampling, C::ONE)
                            .unwrap(),
                    )
                };
                let h = 1e-5;
                let (plus, plus_fields) = evaluate(h);
                let (minus, minus_fields) = evaluate(-h);
                let (base, base_fields) = evaluate(0.0);
                let finite_difference = (plus.value()[()] - minus.value()[()]) / (2.0 * h);
                let derivative = base.first()[()];
                // TE sees ε_yy only, TM ε_xx and ε_zz.
                let seen = (rotation == Rotation::identity())
                    == (polarisation == Polarisation::TransverseMagnetic);
                assert_eq!(
                    derivative.norm() > 1e-4,
                    seen,
                    "{polarisation:?} {derivative}"
                );
                assert!(
                    (derivative - finite_difference).norm() < 1e-9,
                    "{derivative} {finite_difference}"
                );
                let tangent = base_fields.tangent.quantity();
                let (plus, minus) = (
                    plus_fields.fields.quantity(),
                    minus_fields.fields.quantity(),
                );
                for (plus, minus, tangent) in [
                    (plus.electric(), minus.electric(), tangent.electric()),
                    (plus.magnetic(), minus.magnetic(), tangent.magnetic()),
                ] {
                    for (plus, minus, tangent) in [
                        (plus.x(), minus.x(), tangent.x()),
                        (plus.y(), minus.y(), tangent.y()),
                        (plus.z(), minus.z(), tangent.z()),
                    ] {
                        for ((plus, minus), tangent) in plus.iter().zip(minus).zip(tangent) {
                            assert!(((*plus - *minus) / (2.0 * h) - *tangent).norm() < 1e-8);
                        }
                    }
                }
            }
        }
    }
}
