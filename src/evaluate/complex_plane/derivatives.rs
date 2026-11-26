//! Narrow supported-backend characteristic evaluation for local derivatives.
use super::{
    ComplexPlaneEvaluator, ComplexPlaneState, OwnedModeFields,
    mode::{ExteriorWaves, reconstruct_mode_fields},
};
use crate::{
    AnalyticalMaterialHandle, CanonicalCoordinates, Channel, ComplexPlane, ConstitutiveError,
    ConstitutiveTensor, ExteriorWavevectors, FiniteLayerIndex, ModeJet1, PlaneWaveDeterminant,
    Polarisation, ScalarAlgebra, Stack, StackRegion, TensorConstitutiveLift, VectorField,
    algebra::{ArrayJet0, HolomorphicParameter},
    backend::{Backend, Scatter2, Scatter4, Transfer2, Transfer4},
    input::{CanonicalStack, canonical::CanonicalLayer, compile_terminations},
};
use crate::{
    ElectromagneticFields, FieldSampling, SpatialResponse,
    backend::{CoupledSolution, ModalSolutionSource},
};
use ndarray::{Ix0, arr0};
use num_complex::Complex64 as C;
use std::{error::Error, sync::Arc};

type J = ModeJet1<C, Ix0>;
type ValueJet = ArrayJet0<C, Ix0, HolomorphicParameter>;

/// Independently seeded physical inputs for a first-order modal direction.
///
/// `polarisation` names the channel: a [`Polarisation`] for TE or TM, or a
/// [`Channel`] (the coupled channel needs a 4×4 backend).
pub struct ModalCharacteristicDirection<
    'a,
    M = AnalyticalMaterialHandle<C>,
    R = (C, C),
    P = Polarisation,
> {
    pub stack: &'a Stack<M, f64>,
    pub coordinates: &'a CanonicalCoordinates<J>,
    pub exterior: &'a ExteriorWavevectors<J>,
    pub polarisation: P,
    /// Layer thickness rates, in centimetres per parameter unit.
    pub thickness_rates: &'a [f64],
    /// Rates of each region's relative permittivity and permeability, left
    /// exterior, layers left to right, then right exterior.
    pub constitutive_rates: &'a [R],
}

impl<M, R, P: Copy> Clone for ModalCharacteristicDirection<'_, M, R, P> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<M, R, P: Copy> Copy for ModalCharacteristicDirection<'_, M, R, P> {}

/// The rates of one region's relative permittivity and permeability along a
/// modal direction, in the stack frame: total `(ε', μ')` of an isotropic
/// region, or `[ε', μ']` entry by entry, diagonal for the 2×2 backends and
/// any for the 4×4 backends.
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

    fn is_separable(material: &SeededMaterial) -> bool {
        // TE and TM decouple when the xy and yz entries vanish.
        material.0.iter().all(|tensor| match tensor {
            ConstitutiveTensor::Full(entries) => {
                [(0, 1), (1, 0), (1, 2), (2, 1)]
                    .into_iter()
                    .all(|(row, column)| {
                        let entry = &entries[row][column];
                        entry.value()[()] == C::ZERO && entry.first()[()] == C::ZERO
                    })
            }
            _ => true,
        })
    }
}

/// Seed one stack-frame tensor `value` with its stack-frame `rate`, in the
/// smaller class of the two.
fn seed(value: ConstitutiveTensor<ValueJet>, rate: ConstitutiveTensor<C>) -> ConstitutiveTensor<J> {
    let jet = |value: &ValueJet, rate: C| J::from_parts(value.value().clone(), arr0(rate));
    match (value, rate) {
        (ConstitutiveTensor::Isotropic(value), ConstitutiveTensor::Isotropic(rate)) => {
            ConstitutiveTensor::Isotropic(jet(&value, rate))
        }
        (value, rate) => match (value.clone().into_diagonal(), rate.clone().into_diagonal()) {
            (Some(value), Some(rate)) => {
                ConstitutiveTensor::Diagonal(std::array::from_fn(|k| jet(&value[k], rate[k])))
            }
            _ => {
                let rate = match rate {
                    ConstitutiveTensor::Full(rate) => rate,
                    rate => {
                        let diagonal = rate.into_diagonal().unwrap_or([C::ZERO; 3]);
                        std::array::from_fn(|row| {
                            std::array::from_fn(|column| {
                                if row == column {
                                    diagonal[row]
                                } else {
                                    C::ZERO
                                }
                            })
                        })
                    }
                };
                let value = value.full();
                ConstitutiveTensor::Full(std::array::from_fn(|row| {
                    std::array::from_fn(|column| jet(&value[row][column], rate[row][column]))
                }))
            }
        },
    }
}

impl<M, R, P> ModalCharacteristicDirection<'_, M, R, P>
where
    ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
    R: ConstitutiveRates,
    P: Into<Channel> + Copy,
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
            let tensors = [seed(epsilon, epsilon_rate), seed(mu, mu_rate)];
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

/// One mode of [`ModalDerivativeBackend::directional_forms`]: the direction
/// it is evaluated along and the constant amplitude multiplying the
/// backend's homogeneous seed, as for
/// [`ModalDerivativeBackend::directional_fields`].
pub struct ModalFormOperand<'a, M = AnalyticalMaterialHandle<C>, R = (C, C), P = Polarisation> {
    pub direction: ModalCharacteristicDirection<'a, M, R, P>,
    pub amplitude: C,
}

impl<M, R, P: Copy> Clone for ModalFormOperand<'_, M, R, P> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<M, R, P: Copy> Copy for ModalFormOperand<'_, M, R, P> {}

/// A quadratic form `∫ f_a^{†,ᵀ} W f_b dz` of the six components
/// `f = (E_x, E_y, E_z, Z₀H_x, Z₀H_y, Z₀H_z)` of two modes.
#[derive(Clone, Copy, Debug)]
pub struct ModalForm<'a> {
    /// Conjugate `f_a` (a Hermitian form) or not (a bilinear form).
    pub conjugate: bool,
    /// The `6 × 6` weight `W`, rows indexing `f_a`: one for every region, or
    /// one per region (left exterior, layers left to right, right exterior).
    pub weights: &'a [[[C; 6]; 6]],
}

/// Integrals of [`ModalForm`]s, form by form, each region by region: left
/// exterior, layers left to right, right exterior, in cm times the
/// fields' units squared, with their first derivatives along the
/// operands' common direction.
pub type ModalFormIntegrals = Vec<Vec<J>>;

/// The regions of a stack with `layers` finite layers, left to right.
fn regions(layers: usize) -> impl Iterator<Item = StackRegion> {
    std::iter::once(StackRegion::LeftExterior)
        .chain((0..layers).map(|index| StackRegion::Layer(FiniteLayerIndex::new(index))))
        .chain(std::iter::once(StackRegion::RightExterior))
}

/// Whether two operands are the same mode by reference, evaluated once.
fn same_operand<M, R, P>(modes: &[ModalFormOperand<'_, M, R, P>; 2]) -> bool
where
    P: Into<Channel> + Copy,
{
    let [a, b] = modes.map(|mode| mode.direction);
    std::ptr::eq(a.stack, b.stack)
        && std::ptr::eq(a.coordinates, b.coordinates)
        && std::ptr::eq(a.exterior, b.exterior)
        && std::ptr::eq(a.thickness_rates, b.thickness_rates)
        && std::ptr::eq(a.constitutive_rates, b.constitutive_rates)
        && a.polarisation.into() == b.polarisation.into()
        && modes[0].amplitude == modes[1].amplitude
}

/// A jet conjugated as a function of a real parameter: value and derivative.
fn conjugated(jet: &J) -> J {
    J::from_parts(
        jet.value().mapv(|v| v.conj()),
        jet.first().mapv(|v| v.conj()),
    )
}

/// The forms of two mode solutions over every region of a stack with
/// `layers` finite layers.
fn integrate_forms(
    solutions: [&CoupledSolution<'_, J>; 2],
    layers: usize,
    forms: &[ModalForm<'_>],
) -> Result<ModalFormIntegrals, Arc<dyn Error + Send + Sync>> {
    let count = layers + 2;
    if forms
        .iter()
        .any(|form| form.weights.len() != 1 && form.weights.len() != count)
    {
        return Err(Arc::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "form weights do not match the stack's regions",
        )));
    }
    let lift = |weight: &[[C; 6]; 6]| -> [[J; 6]; 6] {
        weight.map(|row| row.map(|entry| J::constant(arr0(entry))))
    };
    let mut integrals = vec![Vec::with_capacity(count); forms.len()];
    for (index, region) in regions(layers).enumerate() {
        // The forms of one pairing over one region in one call.
        for conjugate in [true, false] {
            let selected: Vec<usize> = (0..forms.len())
                .filter(|&form| forms[form].conjugate == conjugate)
                .collect();
            if selected.is_empty() {
                continue;
            }
            let weights: Vec<[[J; 6]; 6]> = selected
                .iter()
                .map(|&form| {
                    let weights = forms[form].weights;
                    lift(&weights[if weights.len() == 1 { 0 } else { index }])
                })
                .collect();
            let values = solutions[0]
                .forms(
                    solutions[1],
                    region,
                    conjugate.then_some(conjugated),
                    &weights,
                )
                .ok_or_else(|| {
                    Arc::new(std::io::Error::other(format!(
                        "no closed-form integral over {region:?}: a PML side or incoming waves"
                    ))) as Arc<dyn Error + Send + Sync>
                })?;
            for (form, value) in selected.into_iter().zip(values) {
                if !value.all_finite() {
                    return Err(Arc::new(std::io::Error::other("nonfinite form integral")));
                }
                integrals[form].push(value);
            }
        }
    }
    Ok(integrals)
}

/// Checked characteristic evaluation on a caller-supplied canonical stack.
///
/// Implemented only by the built-in 2×2 modal backends. This keeps
/// derivative postprocessing explicit about its numerical backend while
/// allowing a selected finite-layer thickness jet to reuse core compilation.
pub trait ModalDerivativeBackend: Clone {
    /// Evaluate the characteristic with independently seeded material and coordinate jets.
    fn directional_characteristic<M, R, P>(
        &self,
        direction: ModalCharacteristicDirection<'_, M, R, P>,
    ) -> Result<J, Arc<dyn Error + Send + Sync>>
    where
        ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
        R: ConstitutiveRates,
        P: Into<Channel> + Copy;

    /// Recover field jets at material-attached coordinates. The backend's
    /// homogeneous seed is used as the gauge, multiplied by a constant amplitude.
    fn directional_fields<M, R, P>(
        &self,
        direction: ModalCharacteristicDirection<'_, M, R, P>,
        sampling: &FieldSampling<f64>,
        amplitude: C,
    ) -> Result<DirectionalModeFields, Arc<dyn Error + Send + Sync>>
    where
        ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
        R: ConstitutiveRates,
        P: Into<Channel> + Copy;

    /// Closed-form integrals over every region of quadratic forms of two
    /// modes' fields (see [`ModalForm`]), each mode the backend's
    /// homogeneous seed along its direction times its amplitude, the gauge
    /// of [`Self::directional_fields`]. Exact over finite layers, wall gaps
    /// and open exteriors (whose fields hold only their outgoing waves),
    /// including at degenerate and defective layer eigenvalues; a side
    /// closed by a PML is refused. The two modes may belong to different
    /// stacks with the same layer thicknesses; operands that are the same
    /// by reference are evaluated once.
    fn directional_forms<M, R, P>(
        &self,
        modes: [ModalFormOperand<'_, M, R, P>; 2],
        forms: &[ModalForm<'_>],
    ) -> Result<ModalFormIntegrals, Arc<dyn Error + Send + Sync>>
    where
        ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
        R: ConstitutiveRates,
        P: Into<Channel> + Copy;

    fn checked_characteristic<M>(
        &self,
        stack: crate::CanonicalStack<M, J>,
        coordinates: &CanonicalCoordinates<J>,
        exterior: &ExteriorWavevectors<J>,
        polarisation: impl Into<Channel>,
    ) -> Result<J, Arc<dyn Error + Send + Sync>>
    where
        M: Clone,
        J: TensorConstitutiveLift<ComplexPlane, M>;
}

fn invalid_amplitude(amplitude: C) -> Result<(), Arc<dyn Error + Send + Sync>> {
    if !amplitude.re.is_finite() || !amplitude.im.is_finite() || amplitude == C::ZERO {
        return Err(Arc::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "invalid field amplitude",
        )));
    }
    Ok(())
}

/// Primal and tangent fields, refusing non-finite samples.
fn directional_mode_fields(
    primal: ElectromagneticFields<VectorField<C, ndarray::Ix1>>,
    tangent: ElectromagneticFields<VectorField<C, ndarray::Ix1>>,
    resolved: crate::ResolvedFieldSampling<f64>,
) -> Result<DirectionalModeFields, Arc<dyn Error + Send + Sync>> {
    for fields in [&primal, &tangent] {
        let (electric, magnetic) = (fields.electric(), fields.magnetic());
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

/// Fields of a 2×2 backend's mode along a direction.
macro_rules! scalar_fields {
    ($backend:expr, $direction:expr, $sampling:expr, $amplitude:expr) => {{
        let (direction, sampling, amplitude) = ($direction, $sampling, $amplitude);
        invalid_amplitude(amplitude)?;
        let stack = direction.canonical_stack()?;
        let channel = direction.polarisation.into();
        let workspace = <Self as Backend<J, SeededConstitutive>>::retain(
            $backend,
            direction.coordinates,
            &stack,
            direction.exterior,
            channel,
        )
        .map_err(super::shared_error)?;
        let state = ComplexPlaneState::new(
            direction.coordinates.clone(),
            direction.exterior.clone(),
            &stack,
            workspace,
            channel,
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
        directional_mode_fields(
            ElectromagneticFields::new(electric.value().clone(), magnetic.value().clone()),
            ElectromagneticFields::new(electric.first().clone(), magnetic.first().clone()),
            resolved,
        )
    }};
}

/// Fields of a 4×4 backend's mode along a direction.
fn coupled_fields<B, M, R, P>(
    backend: &B,
    direction: ModalCharacteristicDirection<'_, M, R, P>,
    sampling: &FieldSampling<f64>,
    amplitude: C,
) -> Result<DirectionalModeFields, Arc<dyn Error + Send + Sync>>
where
    B: Backend<J, SeededConstitutive, Workspace = crate::backend::CoupledWorkspace<J>>,
    ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
    R: ConstitutiveRates,
    P: Into<Channel> + Copy,
{
    invalid_amplitude(amplitude)?;
    let stack = direction.canonical_stack()?;
    let workspace = backend
        .retain(
            direction.coordinates,
            &stack,
            direction.exterior,
            direction.polarisation.into(),
        )
        .map_err(super::shared_error)?;
    let mode = workspace.mode().map_err(super::shared_error)?;
    let solution = workspace
        .mode_solution(&mode)
        .scaled(&J::constant(arr0(amplitude)));
    let resolved = sampling
        .resolve(direction.stack)
        .map_err(super::shared_error)?;
    let centimetres =
        |length: &lamina_units::Length<f64>| J::constant(arr0(C::from(length.as_centimetres())));
    let mut samples: [Vec<J>; 6] = Default::default();
    for position in resolved.positions() {
        let (region, depth) = match position {
            crate::FieldPosition::LeftExterior { distance } => {
                (StackRegion::LeftExterior, centimetres(distance))
            }
            crate::FieldPosition::RightExterior { distance } => {
                (StackRegion::RightExterior, centimetres(distance))
            }
            crate::FieldPosition::Layer { index, position } => {
                let thickness = stack.layers()[index.get()].thickness_cm();
                let depth = match position {
                    crate::ResolvedLayerPosition::FromLeft(distance) => centimetres(distance),
                    crate::ResolvedLayerPosition::FromRight(distance) => {
                        thickness.subtract(&centimetres(distance))
                    }
                    crate::ResolvedLayerPosition::Fraction(fraction) => {
                        thickness.scale(C::from(*fraction))
                    }
                };
                (StackRegion::Layer(*index), depth)
            }
        };
        let [electric, magnetic] = solution.fields(region, &depth).ok_or_else(|| {
            Arc::new(std::io::Error::other("fields unavailable")) as Arc<dyn Error + Send + Sync>
        })?;
        for (slot, value) in samples.iter_mut().zip(electric.into_iter().chain(magnetic)) {
            slot.push(value);
        }
    }
    let part = |range: std::ops::Range<usize>, first: bool| {
        let [x, y, z]: [ndarray::Array1<C>; 3] = std::array::from_fn(|k| {
            samples[range.start + k]
                .iter()
                .map(|jet| {
                    if first {
                        jet.first()[()]
                    } else {
                        jet.value()[()]
                    }
                })
                .collect()
        });
        VectorField::new(x, y, z).map_err(super::shared_error)
    };
    directional_mode_fields(
        ElectromagneticFields::new(part(0..3, false)?, part(3..6, false)?),
        ElectromagneticFields::new(part(0..3, true)?, part(3..6, true)?),
        resolved,
    )
}

/// The retained 4×4 workspace of a direction's seeded stack, with its
/// number of finite layers.
fn coupled_workspace<B, M, R, P>(
    backend: &B,
    direction: ModalCharacteristicDirection<'_, M, R, P>,
) -> Result<(crate::backend::CoupledWorkspace<J>, usize), Arc<dyn Error + Send + Sync>>
where
    B: Backend<J, SeededConstitutive, Workspace = crate::backend::CoupledWorkspace<J>>,
    ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
    R: ConstitutiveRates,
    P: Into<Channel> + Copy,
{
    let stack = direction.canonical_stack()?;
    let workspace = backend
        .retain(
            direction.coordinates,
            &stack,
            direction.exterior,
            direction.polarisation.into(),
        )
        .map_err(super::shared_error)?;
    Ok((workspace, stack.layers().len()))
}

/// Forms of two modes of a 4×4 backend: its own modal solutions.
fn coupled_forms<B, M, R, P>(
    backend: &B,
    modes: [ModalFormOperand<'_, M, R, P>; 2],
    forms: &[ModalForm<'_>],
) -> Result<ModalFormIntegrals, Arc<dyn Error + Send + Sync>>
where
    B: Backend<J, SeededConstitutive, Workspace = crate::backend::CoupledWorkspace<J>>,
    ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
    R: ConstitutiveRates,
    P: Into<Channel> + Copy,
{
    for mode in &modes {
        invalid_amplitude(mode.amplitude)?;
    }
    let same = same_operand(&modes);
    let (first, layers) = coupled_workspace(backend, modes[0].direction)?;
    let second = match same {
        true => None,
        false => Some(coupled_workspace(backend, modes[1].direction)?.0),
    };
    fn solution(
        workspace: &crate::backend::CoupledWorkspace<J>,
        amplitude: C,
    ) -> Result<CoupledSolution<'_, J>, Arc<dyn Error + Send + Sync>> {
        let mode = workspace.mode().map_err(super::shared_error)?;
        Ok(workspace
            .mode_solution(&mode)
            .scaled(&J::constant(arr0(amplitude))))
    }
    let a = solution(&first, modes[0].amplitude)?;
    let b = match &second {
        Some(workspace) => Some(solution(workspace, modes[1].amplitude)?),
        None => None,
    };
    integrate_forms([&a, b.as_ref().unwrap_or(&a)], layers, forms)
}

/// The six field components of a 2×2 backend's mode at every interface,
/// left to right, and the 4×4 workspace of its channel that integrates
/// them; the mode keeps the 2×2 backend's gauge.
macro_rules! scalar_interface_fields {
    ($backend:expr, $direction:expr, $amplitude:expr) => {{
        let (direction, amplitude) = ($direction, $amplitude);
        invalid_amplitude(amplitude)?;
        let stack = direction.canonical_stack()?;
        let channel = direction.polarisation.into();
        let workspace = <Self as Backend<J, SeededConstitutive>>::retain(
            $backend,
            direction.coordinates,
            &stack,
            direction.exterior,
            channel,
        )
        .map_err(super::shared_error)?;
        let state = ComplexPlaneState::new(
            direction.coordinates.clone(),
            direction.exterior.clone(),
            &stack,
            workspace,
            channel,
        );
        let seed = state
            .workspace()
            .modal_boundary_solution()
            .map_err(super::shared_error)?
            .scaled(&J::constant(arr0(amplitude)));
        let layers = stack.layers().len();
        // The left interface of every layer, then the last interface.
        let zero = || vec![0.0];
        let sampling = (0..layers)
            .fold(FieldSampling::new(), |sampling, index| {
                sampling.layer(
                    index,
                    crate::LayerSampling::offsets_in(zero(), lamina_units::LengthUnit::Centimetre),
                )
            })
            .right_exterior(crate::ExteriorSampling::point(
                lamina_units::Length::centimetres(0.0),
            ));
        let resolved = sampling
            .resolve(direction.stack)
            .map_err(super::shared_error)?;
        let (electric, magnetic) =
            reconstruct_mode_fields(&state, &seed, &resolved, ExteriorWaves::Outgoing)
                .map_err(super::shared_error)?
                .into_parts();
        let (value, first) = (
            [electric.value(), magnetic.value()],
            [electric.first(), magnetic.first()],
        );
        let component = |field: usize, index: usize, sample: usize| {
            let pick = |part: &VectorField<C, ndarray::Ix1>| match index {
                0 => part.x()[sample],
                1 => part.y()[sample],
                _ => part.z()[sample],
            };
            J::from_parts(arr0(pick(value[field])), arr0(pick(first[field])))
        };
        let fields: Vec<[J; 6]> = (0..=layers)
            .map(|sample| std::array::from_fn(|k| component(k / 3, k % 3, sample)))
            .collect();
        let coupled = <crate::backend::Transfer4 as Backend<J, SeededConstitutive>>::retain(
            &crate::backend::Transfer4::new(),
            direction.coordinates,
            &stack,
            direction.exterior,
            channel,
        )
        .map_err(super::shared_error)?;
        (fields, coupled, layers)
    }};
}

/// Forms of two modes of a 2×2 backend: their fields at the interfaces,
/// carried through each region by the 4×4 generators of their channel.
macro_rules! scalar_forms {
    ($backend:expr, $modes:expr, $forms:expr) => {{
        let (modes, forms) = ($modes, $forms);
        let same = same_operand(&modes);
        let (fields_a, workspace_a, layers) =
            scalar_interface_fields!($backend, modes[0].direction, modes[0].amplitude);
        let second = match same {
            true => None,
            false => {
                let (fields, workspace, _) =
                    scalar_interface_fields!($backend, modes[1].direction, modes[1].amplitude);
                Some((fields, workspace))
            }
        };
        let unavailable = || {
            Arc::new(std::io::Error::other(
                "interface fields do not match the stack",
            )) as Arc<dyn Error + Send + Sync>
        };
        let a = workspace_a
            .interface_solution(&fields_a)
            .ok_or_else(unavailable)?;
        let b = match &second {
            Some((fields, workspace)) => Some(
                workspace
                    .interface_solution(fields)
                    .ok_or_else(unavailable)?,
            ),
            None => None,
        };
        integrate_forms([&a, b.as_ref().unwrap_or(&a)], layers, forms)
    }};
}

macro_rules! implementation {
    ($backend:ty, scalar) => {
        implementation!(
            $backend,
            |backend, direction, sampling, amplitude| scalar_fields!(
                backend, direction, sampling, amplitude
            ),
            |backend, modes, forms| scalar_forms!(backend, modes, forms)
        );
    };
    ($backend:ty, coupled) => {
        implementation!(
            $backend,
            |backend, direction, sampling, amplitude| {
                coupled_fields(backend, direction, sampling, amplitude)
            },
            |backend, modes, forms| coupled_forms(backend, modes, forms)
        );
    };
    (
        $backend:ty,
        |$b:ident, $d:ident, $s:ident, $a:ident| $fields:expr,
        |$fb:ident, $fm:ident, $ff:ident| $forms:expr
    ) => {
        impl ModalDerivativeBackend for $backend {
            fn directional_characteristic<M, R, P>(
                &self,
                direction: ModalCharacteristicDirection<'_, M, R, P>,
            ) -> Result<J, Arc<dyn Error + Send + Sync>>
            where
                ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
                R: ConstitutiveRates,
                P: Into<Channel> + Copy,
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

            fn directional_fields<M, R, P>(
                &self,
                direction: ModalCharacteristicDirection<'_, M, R, P>,
                sampling: &FieldSampling<f64>,
                amplitude: C,
            ) -> Result<DirectionalModeFields, Arc<dyn Error + Send + Sync>>
            where
                ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
                R: ConstitutiveRates,
                P: Into<Channel> + Copy,
            {
                let ($b, $d, $s, $a) = (self, direction, sampling, amplitude);
                $fields
            }

            fn directional_forms<M, R, P>(
                &self,
                modes: [ModalFormOperand<'_, M, R, P>; 2],
                forms: &[ModalForm<'_>],
            ) -> Result<ModalFormIntegrals, Arc<dyn Error + Send + Sync>>
            where
                ValueJet: TensorConstitutiveLift<ComplexPlane, M>,
                R: ConstitutiveRates,
                P: Into<Channel> + Copy,
            {
                let ($fb, $fm, $ff) = (self, modes, forms);
                $forms
            }

            fn checked_characteristic<M>(
                &self,
                stack: CanonicalStack<M, J>,
                coordinates: &CanonicalCoordinates<J>,
                exterior: &ExteriorWavevectors<J>,
                polarisation: impl Into<Channel>,
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

implementation!(Scatter2, scalar);
implementation!(Transfer2, scalar);
implementation!(Scatter4, coupled);
implementation!(Transfer4, coupled);

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
