//! The 4×4 (Berreman) formulation shared by [`Transfer4`](super::Transfer4)
//! and [`Scatter4`](super::Scatter4).
//!
//! Both backends propagate the tangential state `ψ = (E_x, η_y, E_y, −η_x)`,
//! `η = Z₀H`, restricted to the components of the channel: all four for the
//! coupled channel, `(E_x, η_y)` for TM and `(E_y, −η_x)` for TE on separable
//! stacks. They share the region generators, the boundary subspaces of the
//! two sides, the entries (characteristic and Jones matrices), the retained
//! interface states and the field reconstruction; they differ only in how the
//! layers are composed.
//!
//! # Boundary subspaces
//!
//! An open side contributes the span of its outgoing waves, from the spectral
//! projector of the supplied outward wavenumbers applied to the reference
//! vectors `E_y`, `η_y` (see [`split`]). A wall contributes the states it
//! admits — `span{η_y, −η_x}` for a perfect electric wall, `span{E_x, E_y}`
//! for a perfect magnetic one — carried to the outermost interface through
//! the gap and the PML by `exp(∓iℓK)` of the exterior generator, `ℓ` the
//! complex length from the interface to the wall. With a PML the propagator
//! is normalised by the `k∥`-independent `exp(±iℓ n_ref)`, `n_ref² = −2 e₂/N`
//! of the exterior generator at `k∥ = 0` (`e₂` its second elementary
//! symmetric function of the `N` eigenvalues), so thick PMLs stay bounded.
//!
//! # Characteristic
//!
//! A mode is a nonzero state in the left boundary subspace at the first
//! interface that the stack carries into the right boundary subspace at the
//! last. Transfer4 evaluates `det[V_L | M V_R]` with `ψ_left = M ψ_right`;
//! Scatter4 the product of the pivot determinants of the same linear system
//! in bounded eigen-wave amplitudes. Both are holomorphic where the exterior
//! splits are, and differ by a nonvanishing factor.

pub(crate) mod berreman;
pub(crate) mod exponential;
pub(crate) mod integral;
pub(crate) mod matrix;
pub(crate) mod modal;
pub(crate) mod solution;
pub(crate) mod split;

#[cfg(test)]
mod point_tests;

use nalgebra::ComplexField;
use ndarray::Dimension;
use num_traits::One;

use crate::{
    ComplexScalar, FiniteLayerIndex, IncidentSide, JonesAmplitudes, PlaneWaveDeterminant,
    PlaneWavePower,
    algebra::{EntryAlgebra, InlineJet, RealScalarAlgebra, ScalarAlgebra},
    backend::{
        Channel, ChannelError, ExteriorBranch, ExteriorSelectionError, ExteriorWavevectors,
        ModalSolutionSource, ModeReconstructionError, PlaneWaveEntries, PlaneWaveSolution,
        PlaneWaveSolutionSource, RunMode, SolutionWorkspace, TerminationError,
        termination::check_exterior_pattern,
    },
    input::{CanonicalCoordinates, CanonicalStack, canonical::CanonicalTermination},
    material::{ConstitutiveError, ConstitutiveTensor, TensorConstitutiveLift},
    observable::{ProjectJonesAmplitudes, ProjectPlaneWaveModeDeterminant, ProjectPower},
    stack::{StackRegion, StackSide, WallKind},
};

use berreman::{Region, local_positions, reference_positions};
use matrix::Matrix;
use split::{Direction, Monic, Split};

pub use modal::CoupledModeCandidate;
pub use solution::{CoupledSolution, FieldWeights};

/// The errors a 4×4 backend reports while setting up a stack.
pub(crate) trait CoupledError:
    From<TerminationError> + From<ConstitutiveError> + From<ChannelError>
{
    /// The supplied outgoing waves of `side` do not span an invariant
    /// subspace of its exterior.
    fn exterior(side: StackSide) -> Self;
}

/// The boundary subspace of one side at its outermost interface.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Boundary<A> {
    /// The exterior region.
    region: Region<A>,
    /// The outgoing waves of an open side, or the states a closed side
    /// admits.
    pub(crate) outgoing: Matrix<A>,
    /// The incoming waves of an open side.
    incoming: Option<Matrix<A>>,
    /// `σ_min / σ_max` of the outgoing basis of an open side.
    conditioning: Option<f64>,
    /// The distance to the wall of a closed side without a PML.
    gap: Option<A>,
}

impl<A: EntryAlgebra> Boundary<A> {
    /// The boundary of `side` in the inline form `A` of the caller's jets
    /// `J`.
    fn new<J, E, M, Error: CoupledError>(
        side: StackSide,
        coordinates: &CanonicalCoordinates<J>,
        material: &M,
        termination: &CanonicalTermination<J>,
        branch: Option<&ExteriorBranch<J>>,
        channel: Channel,
    ) -> Result<Self, Error>
    where
        J: InlineJet<Inline = A> + TensorConstitutiveLift<E, M>,
        J::Scalar: ComplexScalar,
        J::Dimension: Dimension,
    {
        let name = match side {
            StackSide::Left => StackRegion::LeftExterior,
            StackSide::Right => StackRegion::RightExterior,
        };
        let region = Region::evaluate_inline::<E, M, Error>(name, material, coordinates, channel)?;
        let generator = &region.generator;
        let like = generator.get(0, 0);

        let Some(branch) = branch else {
            let (outgoing, gap) = closed_subspace::<J, E, M, Error>(
                side,
                name,
                coordinates,
                material,
                termination,
                &region,
                channel,
            )?;
            return Ok(Self {
                region,
                outgoing,
                incoming: None,
                conditioning: None,
                gap,
            });
        };

        let outward: Vec<A> = match (branch, channel) {
            (ExteriorBranch::Pair(pair), Channel::Coupled) => {
                pair.iter().map(J::to_inline).collect()
            }
            (ExteriorBranch::Scalar(kappa), _) if channel != Channel::Coupled => {
                vec![kappa.to_inline()]
            }
            _ => return Err(ChannelError::BranchKind { side, channel }.into()),
        };
        let roots: Vec<A> = match side {
            StackSide::Right => outward,
            StackSide::Left => outward.iter().map(A::negate).collect(),
        };
        let split = Split::from_factor(generator, Monic::from_roots(&roots, like))
            .ok_or_else(|| Error::exterior(side))?;
        let reference = Matrix::unit_columns(generator.rows(), &reference_positions(channel), like);
        let outgoing = split.projector.multiply(&reference);

        Ok(Self {
            conditioning: Some(split::conditioning(&outgoing)),
            incoming: Some(reference.subtract(&outgoing)),
            region,
            outgoing,
            gap: None,
        })
    }

    /// The boundary with every jet mapped, for moving between jet forms.
    fn into_map<B>(self, f: &mut impl FnMut(A) -> B) -> Boundary<B> {
        Boundary {
            region: self.region.into_map(f),
            outgoing: self.outgoing.into_map(&mut *f),
            incoming: self.incoming.map(|incoming| incoming.into_map(&mut *f)),
            conditioning: self.conditioning,
            gap: self.gap.map(&mut *f),
        }
    }
}

/// The states a closed side admits at its outermost interface, and the
/// distance to its wall when it has no PML.
#[allow(clippy::type_complexity)]
fn closed_subspace<J, E, M, Error: CoupledError>(
    side: StackSide,
    name: StackRegion,
    coordinates: &CanonicalCoordinates<J>,
    material: &M,
    termination: &CanonicalTermination<J>,
    region: &Region<J::Inline>,
    channel: Channel,
) -> Result<(Matrix<J::Inline>, Option<J::Inline>), Error>
where
    J: InlineJet + TensorConstitutiveLift<E, M>,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    let (Some(kind), Some(gap)) = (termination.wall_kind(), termination.gap_cm()) else {
        return Err(TerminationError::ExteriorBranchMismatch { side, open: true }.into());
    };
    let generator = &region.generator;
    let (n, like) = (generator.rows(), generator.get(0, 0));
    let admitted = match kind {
        WallKind::PerfectElectric => local_positions(channel, &[1, 3]),
        WallKind::PerfectMagnetic => local_positions(channel, &[0, 2]),
    };
    let i = J::Scalar::i();

    // Right: ψ = exp(−iℓ(K − n_ref)) w; left: ψ = exp(iℓ(K + n_ref)) w.
    let mut exponent = match side {
        StackSide::Right => generator.scale(-i),
        StackSide::Left => generator.scale(i),
    };
    let mut length = gap.to_inline();
    let mut wall_gap = Some(length.clone());
    if let Some(pml) = termination.pml() {
        wall_gap = None;
        let k0 = coordinates.vacuum_angular_wavenumber();
        let at_normal = CanonicalCoordinates::new(k0.clone(), ScalarAlgebra::zero_like(k0));
        let normal = Region::evaluate_inline::<E, M, Error>(name, material, &at_normal, channel)?;
        let e2 = &Monic::characteristic(&normal.generator).0[n - 2];
        let reference = e2.scale(nalgebra::convert(-2.0 / n as f64)).sqrt();
        exponent = exponent.shift(&reference.scale(i));
        length = length.add(&pml.optical_thickness_cm.to_inline());
    }
    Ok((
        exponential::exponential(&exponent.scale_by(&length))
            .multiply(&Matrix::unit_columns(n, &admitted, like)),
        wall_gap,
    ))
}

/// The exterior data of a solved stack: the channel and each side's boundary
/// subspace.
#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledContext<A> {
    channel: Channel,
    /// `k₀`.
    vacuum_angular_wavenumber: A,
    left: Boundary<A>,
    right: Boundary<A>,
}

impl<A> CoupledContext<A> {
    /// `σ_min / σ_max` over every sample of the outgoing basis of an open
    /// side, built from the reference vectors `E_y`, `Z₀H_y`; `None` for a
    /// closed side. A small value flags a reference pair that nearly fails
    /// to span the outgoing waves, where the characteristic gauge has a
    /// spurious zero.
    pub fn basis_conditioning(&self, side: StackSide) -> Option<f64> {
        self.boundary(side).conditioning
    }

    pub(crate) fn boundary(&self, side: StackSide) -> &Boundary<A> {
        match side {
            StackSide::Left => &self.left,
            StackSide::Right => &self.right,
        }
    }
}

impl<A: EntryAlgebra> CoupledContext<A> {
    /// The context with every jet mapped, for moving between jet forms.
    fn into_map<B>(self, f: &mut impl FnMut(A) -> B) -> CoupledContext<B> {
        CoupledContext {
            channel: self.channel,
            vacuum_angular_wavenumber: f(self.vacuum_angular_wavenumber),
            left: self.left.into_map(f),
            right: self.right.into_map(f),
        }
    }

    /// The incoming states at the first and at the last interface of the
    /// `2m` excitations: unit amplitude in each incoming wave of the left
    /// side, then of the right side (zero for a closed side).
    pub(crate) fn excitations(&self) -> (Matrix<A>, Matrix<A>) {
        let zeros = Matrix::zeros(
            self.left.outgoing.rows(),
            self.left.outgoing.columns(),
            self.left.outgoing.get(0, 0),
        );
        let incoming =
            |boundary: &Boundary<A>| boundary.incoming.clone().unwrap_or_else(|| zeros.clone());
        (
            incoming(&self.left).hstack(&zeros),
            zeros.hstack(&incoming(&self.right)),
        )
    }
}

/// The principal outward pair of a full exterior for the coupled channel: the
/// roots of each invariant block of a separable exterior (TE then TM), or the
/// two outgoing roots of its generator.
pub(crate) fn principal_pair<E, M, A>(
    side: StackSide,
    name: StackRegion,
    material: &M,
    coordinates: &CanonicalCoordinates<A>,
) -> Result<ExteriorBranch<A>, ExteriorSelectionError>
where
    A: ScalarAlgebra + TensorConstitutiveLift<E, M>,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let (direction, sign) = match side {
        StackSide::Right => (Direction::Forward, A::Scalar::one()),
        StackSide::Left => (Direction::Backward, -A::Scalar::one()),
    };
    let outgoing = |channel| {
        let region =
            Region::evaluate::<E, M, ExteriorSelectionError>(name, material, coordinates, channel)?;
        Split::select(&region.generator, &region.generator, direction)
            .map(|split| split.outgoing.0)
            .ok_or(ExteriorSelectionError::BranchPoint { side })
    };

    let roots = if A::is_separable(material) {
        let root = |channel| outgoing(channel).map(|factor: Vec<A>| factor[0].negate());
        [
            root(Channel::TransverseElectric)?,
            root(Channel::TransverseMagnetic)?,
        ]
    } else {
        // x² − s x + p with s = −f₁, p = f₀.
        let factor = outgoing(Channel::Coupled)?;
        let sum = factor[1].negate();
        let discriminant = sum
            .square()
            .subtract(&factor[0].scale(nalgebra::convert(4.0)))
            .sqrt();
        [sum.add(&discriminant), sum.subtract(&discriminant)]
            .map(|root| root.scale(nalgebra::convert(0.5)))
    };
    Ok(ExteriorBranch::Pair(roots.map(|root| root.scale(sign))))
}

/// The solution of the interface system for the `2m` excitations of
/// [`CoupledContext::excitations`].
pub(crate) struct Solved<A> {
    pub(crate) characteristic: A,
    /// The outgoing amplitudes on the left, `m × 2m`.
    pub(crate) left: Matrix<A>,
    /// The outgoing amplitudes on the right, `m × 2m`.
    pub(crate) right: Matrix<A>,
    /// The state at every interface, `n × 2m`, and the homogeneous system,
    /// when retained.
    pub(crate) states: Option<(Vec<Matrix<A>>, modal::Modal<A>)>,
}

/// A finite layer: its region and thickness.
pub(crate) type Layer<A> = (Region<A>, A);

/// A stack solved in the inline form of the caller's jets `J`, before the
/// results are converted back.
pub(crate) struct Accumulated<J: InlineJet> {
    context: CoupledContext<J::Inline>,
    layers: Vec<Layer<J::Inline>>,
    solved: Solved<J::Inline>,
}

/// Set up `stack` and solve it with `compose`, which composes the layers
/// between the two boundary subspaces and retains the interface states when
/// asked.
///
/// The boundaries, layers and composition run in the inline form of the
/// caller's jets ([`InlineJet`]); [`Accumulated::into_workspace`] converts the
/// results back once, [`Accumulated::into_characteristic`] only the
/// characteristic.
#[allow(clippy::type_complexity)]
pub(crate) fn accumulate<J, E, M, Error: CoupledError>(
    coordinates: &CanonicalCoordinates<J>,
    stack: &CanonicalStack<M, J>,
    exterior: &ExteriorWavevectors<J>,
    channel: Channel,
    request: RunMode,
    compose: impl FnOnce(
        &CoupledContext<J::Inline>,
        &[Layer<J::Inline>],
        bool,
    ) -> Result<Solved<J::Inline>, Error>,
) -> Result<Accumulated<J>, Error>
where
    J: InlineJet + TensorConstitutiveLift<E, M>,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    check_exterior_pattern(stack, exterior)?;
    let boundary = |side, material, termination, branch| {
        Boundary::new::<J, E, M, Error>(side, coordinates, material, termination, branch, channel)
    };
    let context = CoupledContext {
        channel,
        vacuum_angular_wavenumber: coordinates.vacuum_angular_wavenumber().to_inline(),
        left: boundary(
            StackSide::Left,
            stack.left_exterior(),
            stack.left_termination(),
            exterior.left(),
        )?,
        right: boundary(
            StackSide::Right,
            stack.right_exterior(),
            stack.right_termination(),
            exterior.right(),
        )?,
    };
    let layers = layer_regions::<J, E, M, Error>(coordinates, stack, channel)?;
    let solved = compose(&context, &layers, request.is_requested())?;
    Ok(Accumulated {
        context,
        layers,
        solved,
    })
}

impl<J> Accumulated<J>
where
    J: InlineJet,
    J::Scalar: ComplexScalar,
{
    /// The characteristic in the caller's jets.
    pub(crate) fn into_characteristic(self) -> J {
        J::from_inline(self.solved.characteristic)
    }

    /// The solution, and the retained layers and states when they were
    /// requested, in the caller's jets.
    pub(crate) fn into_workspace(self) -> CoupledWorkspace<J> {
        let Self {
            context,
            layers,
            solved,
        } = self;
        let back = &mut J::from_inline;
        let context = context.into_map(back);
        let solved = Solved {
            characteristic: back(solved.characteristic),
            left: solved.left.into_map(&mut *back),
            right: solved.right.into_map(&mut *back),
            states: solved.states.map(|(states, modal)| {
                (
                    states
                        .into_iter()
                        .map(|state| state.into_map(&mut *back))
                        .collect::<Vec<_>>(),
                    modal.into_map(&mut *back),
                )
            }),
        };

        let m = context.left.outgoing.columns();
        let response = |incident: &Boundary<J>,
                        other: &Boundary<J>,
                        columns: std::ops::Range<usize>,
                        own: &Matrix<J>,
                        far: &Matrix<J>| {
            let like = own.get(0, 0);
            if incident.incoming.is_none() {
                let nan = Matrix::from_fn(m, m, |_, _| {
                    like.constant(J::Scalar::from_real(nalgebra::convert(f64::NAN)))
                });
                return Response {
                    reflection: nan.clone(),
                    transmission: nan,
                };
            }
            Response {
                reflection: own.columns_in(columns.clone()),
                transmission: if other.incoming.is_some() {
                    far.columns_in(columns)
                } else {
                    Matrix::zeros(m, m, like)
                },
            }
        };
        let entries = CoupledEntries {
            characteristic: solved.characteristic,
            responses: [
                response(
                    &context.left,
                    &context.right,
                    0..m,
                    &solved.left,
                    &solved.right,
                ),
                response(
                    &context.right,
                    &context.left,
                    m..2 * m,
                    &solved.right,
                    &solved.left,
                ),
            ],
        };
        let retained = solved.states.map(|(states, modal)| Retained {
            layers: layers
                .into_iter()
                .map(|(region, thickness)| (region.into_map(&mut *back), back(thickness)))
                .collect(),
            states,
            modal,
        });

        CoupledWorkspace {
            solution: PlaneWaveSolution::new(entries, context),
            retained,
        }
    }
}

/// The finite layers of a stack for a channel, with their thicknesses, in
/// the inline form of the jets.
pub(crate) fn layer_regions<J, E, M, Error>(
    coordinates: &CanonicalCoordinates<J>,
    stack: &CanonicalStack<M, J>,
    channel: Channel,
) -> Result<Vec<Layer<J::Inline>>, Error>
where
    J: InlineJet + TensorConstitutiveLift<E, M>,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
    Error: From<ConstitutiveError> + From<ChannelError>,
{
    stack
        .layers()
        .iter()
        .enumerate()
        .map(|(index, layer)| {
            let name = StackRegion::Layer(FiniteLayerIndex::new(index));
            Region::evaluate_inline::<E, M, Error>(name, layer.material(), coordinates, channel)
                .map(|region| (region, layer.thickness_cm().to_inline()))
        })
        .collect()
}

/// Implement [`Backend`](crate::Backend) for a 4×4 backend with an
/// `accumulate` method.
macro_rules! coupled_backend {
    ($backend:ty, $error:ty) => {
        impl<J, Domain> $crate::backend::Backend<J, Domain> for $backend
        where
            J: $crate::algebra::InlineJet,
            J::Scalar: $crate::ComplexScalar,
            J::Dimension: ndarray::Dimension,
            <J::Scalar as nalgebra::ComplexField>::RealField: Copy,
        {
            type Error = $error;
            type Entries = $crate::backend::coupled::CoupledEntries<J>;
            type Workspace = $crate::backend::coupled::CoupledWorkspace<J>;

            fn solve<M>(
                &self,
                coordinates: &$crate::CanonicalCoordinates<J>,
                stack: &$crate::input::CanonicalStack<M, J>,
                exterior: &$crate::ExteriorWavevectors<J>,
                channel: $crate::Channel,
            ) -> Result<$crate::backend::PlaneWaveSolution<Self::Entries>, Self::Error>
            where
                J: $crate::material::TensorConstitutiveLift<Domain, M>,
            {
                use $crate::backend::SolutionWorkspace;
                Self::accumulate::<J, Domain, M>(
                    self,
                    coordinates,
                    stack,
                    exterior,
                    channel,
                    $crate::backend::RunMode::ResponseOnly,
                )
                .map(|accumulated| accumulated.into_workspace().into_solution())
            }

            /// The characteristic alone: the context and responses stay in
            /// the inline form and are dropped.
            fn solve_determinant<M>(
                &self,
                coordinates: &$crate::CanonicalCoordinates<J>,
                stack: &$crate::input::CanonicalStack<M, J>,
                exterior: &$crate::ExteriorWavevectors<J>,
                channel: $crate::Channel,
            ) -> Result<
                <Self::Entries as $crate::observable::ProjectPlaneWaveModeDeterminant>::Determinant,
                Self::Error,
            >
            where
                J: $crate::material::TensorConstitutiveLift<Domain, M>,
            {
                Self::accumulate::<J, Domain, M>(
                    self,
                    coordinates,
                    stack,
                    exterior,
                    channel,
                    $crate::backend::RunMode::ResponseOnly,
                )
                .map(|accumulated| {
                    $crate::PlaneWaveDeterminant::new(accumulated.into_characteristic())
                })
            }

            fn retain<M>(
                &self,
                coordinates: &$crate::CanonicalCoordinates<J>,
                stack: &$crate::input::CanonicalStack<M, J>,
                exterior: &$crate::ExteriorWavevectors<J>,
                channel: $crate::Channel,
            ) -> Result<Self::Workspace, Self::Error>
            where
                J: $crate::material::TensorConstitutiveLift<Domain, M>,
            {
                Self::accumulate::<J, Domain, M>(
                    self,
                    coordinates,
                    stack,
                    exterior,
                    channel,
                    $crate::backend::RunMode::InternalFields,
                )
                .map($crate::backend::coupled::Accumulated::into_workspace)
            }
        }
    };
}
pub(crate) use coupled_backend;

/// The reflection and transmission matrices for incidence from one side.
#[derive(Clone, Debug, PartialEq)]
struct Response<A> {
    reflection: Matrix<A>,
    transmission: Matrix<A>,
}

/// The solved response of a 4×4 backend: the characteristic and the Jones
/// matrices for incidence from either side.
///
/// Scalar reflection and transmission amplitudes belong to the scalar
/// channels of the 2×2 backends:
///
/// ```
/// # use lamina_core::{
/// #     CoordinateInput, Coordinates, InPlaneCoordinate, IncidentSide, MaterialStack,
/// #     Polarisation, RealAxisEvaluator, SpectralCoordinate, Stack, Transfer2,
/// #     material::Constant, units::{InverseLengthUnit, Length},
/// # };
/// # let stack: MaterialStack<num_complex::Complex64> =
/// #     Stack::from_materials(Constant::new(1.0, 1.0), Constant::new(2.0, 1.0))
/// #         .material_layer(Constant::new(2.25, 1.0), Length::nanometres(500.0f64))
/// #         .finalise();
/// # let coordinates = CoordinateInput::intrinsic(
/// #     Coordinates::new(
/// #         SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre),
/// #         InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre),
/// #     ),
/// #     ndarray::arr0(1.0e4),
/// #     ndarray::arr0(0.0),
/// # )?;
/// let response = RealAxisEvaluator::new(Transfer2::new()).evaluate(
///     coordinates,
///     &stack,
///     Polarisation::TransverseElectric,
/// )?;
/// let amplitudes = response.amplitudes(IncidentSide::Left)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// The 4×4 entries do not provide them, in any channel; their responses are
/// Jones matrices ([`ProjectJonesAmplitudes`]):
///
/// ```compile_fail
/// # use lamina_core::{
/// #     CoordinateInput, Coordinates, InPlaneCoordinate, IncidentSide, MaterialStack,
/// #     Polarisation, RealAxisEvaluator, SpectralCoordinate, Stack, Transfer4,
/// #     material::Constant, units::{InverseLengthUnit, Length},
/// # };
/// # let stack: MaterialStack<num_complex::Complex64> =
/// #     Stack::from_materials(Constant::new(1.0, 1.0), Constant::new(2.0, 1.0))
/// #         .material_layer(Constant::new(2.25, 1.0), Length::nanometres(500.0f64))
/// #         .finalise();
/// # let coordinates = CoordinateInput::intrinsic(
/// #     Coordinates::new(
/// #         SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre),
/// #         InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre),
/// #     ),
/// #     ndarray::arr0(1.0e4),
/// #     ndarray::arr0(0.0),
/// # )?;
/// let response = RealAxisEvaluator::new(Transfer4::new()).evaluate(
///     coordinates,
///     &stack,
///     Polarisation::TransverseElectric,
/// )?;
/// let amplitudes = response.amplitudes(IncidentSide::Left)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledEntries<A> {
    characteristic: A,
    responses: [Response<A>; 2],
}

impl<A> CoupledEntries<A> {
    fn response(&self, side: IncidentSide) -> &Response<A> {
        &self.responses[usize::from(side == IncidentSide::Right)]
    }
}

impl<A> PlaneWaveEntries for CoupledEntries<A> {
    type ExteriorContext = CoupledContext<A>;
    type Algebra = A;
}

impl<A: Clone> ProjectPlaneWaveModeDeterminant for CoupledEntries<A> {
    type Determinant = PlaneWaveDeterminant<A>;

    fn project_determinant(&self, _: &Self::ExteriorContext) -> Self::Determinant {
        PlaneWaveDeterminant::new(self.characteristic.clone())
    }
}

impl<A: Clone> ProjectJonesAmplitudes for CoupledEntries<A> {
    fn project_jones_amplitudes(
        &self,
        _: &Self::ExteriorContext,
        incident_side: IncidentSide,
    ) -> JonesAmplitudes<A> {
        let response = self.response(incident_side);
        JonesAmplitudes::new(
            response.reflection.rows(),
            response.reflection.entries().to_vec(),
            response.transmission.entries().to_vec(),
        )
    }
}

impl<A> ProjectPower for CoupledEntries<A>
where
    A: RealScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
    A::RealJet: ScalarAlgebra,
{
    type Power = PlaneWavePower<Vec<A::RealJet>>;

    /// Reflectance, transmittance and absorptance per incident eigen-wave:
    /// each power is the signed flux `S_z ∝ Re(ψ† Γ ψ)` of the reflected or
    /// transmitted state relative to the incident one.
    fn project_power(
        &self,
        exterior: &Self::ExteriorContext,
        incident_side: IncidentSide,
    ) -> Self::Power {
        let (incident, other) = match incident_side {
            IncidentSide::Left => (&exterior.left, &exterior.right),
            IncidentSide::Right => (&exterior.right, &exterior.left),
        };
        let response = self.response(incident_side);
        let incoming = incident.incoming.as_ref().unwrap_or(&incident.outgoing);
        let reflected = incident.outgoing.multiply(&response.reflection);
        let transmitted = other.outgoing.multiply(&response.transmission);

        let (mut reflectance, mut transmittance, mut absorptance) =
            (Vec::new(), Vec::new(), Vec::new());
        for wave in 0..incoming.columns() {
            let flux = |state: &Matrix<A>| flux(&state.columns_in(wave..wave + 1));
            let incident_flux = flux(incoming);
            let r = flux(&reflected).divide(&incident_flux).negate();
            let t = flux(&transmitted).divide(&incident_flux);
            absorptance.push(r.one_like().subtract(&r).subtract(&t));
            reflectance.push(r);
            transmittance.push(t);
        }
        PlaneWavePower::new(reflectance, transmittance, absorptance)
    }
}

/// `Re(ψ† Γ ψ)` of a state, `Γ` pairing the components `(2k, 2k + 1)`.
fn flux<A>(state: &Matrix<A>) -> A::RealJet
where
    A: RealScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    (1..state.rows() / 2)
        .fold(
            state.get(0, 0).hermitian_product(state.get(1, 0)),
            |sum, pair| {
                sum.add(
                    &state
                        .get(2 * pair, 0)
                        .hermitian_product(state.get(2 * pair + 1, 0)),
                )
            },
        )
        .real()
}

/// A retained 4×4 solution: the response, every finite layer and the states
/// at every interface.
#[doc(hidden)]
#[derive(Clone, Debug)]
pub struct CoupledWorkspace<A> {
    solution: PlaneWaveSolution<CoupledEntries<A>>,
    retained: Option<Retained<A>>,
}

#[derive(Clone, Debug)]
struct Retained<A> {
    layers: Vec<Layer<A>>,
    /// The state at every interface for the `2m` excitations.
    states: Vec<Matrix<A>>,
    /// The homogeneous interface system.
    modal: modal::Modal<A>,
}

impl<A> PlaneWaveSolutionSource for CoupledWorkspace<A> {
    type Entries = CoupledEntries<A>;

    fn solution(&self) -> crate::backend::PlaneWaveSolutionView<'_, Self::Entries> {
        self.solution.as_view()
    }
}

impl<A> ModalSolutionSource for CoupledWorkspace<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    type Algebra = A;
    type Mode = CoupledModeCandidate<A>;

    fn modal_boundary_solution(&self) -> Result<Self::Mode, ModeReconstructionError> {
        self.mode()
    }
}

impl<A> SolutionWorkspace for CoupledWorkspace<A> {
    fn into_solution(self) -> PlaneWaveSolution<CoupledEntries<A>> {
        self.solution
    }
}

impl<A: Clone> crate::backend::RetainedLayers for CoupledWorkspace<A> {
    type Algebra = A;
    type Layer = Region<A>;

    fn retained_layer_count(&self) -> Option<usize> {
        Some(self.retained.as_ref()?.layers.len())
    }

    fn layer(&self, index: usize) -> Option<&Region<A>> {
        Some(&self.retained.as_ref()?.layers.get(index)?.0)
    }

    fn layer_tensors(&self, index: usize) -> Option<[ConstitutiveTensor<A>; 2]> {
        Some(self.layer(index)?.tensors.clone())
    }

    fn layer_thickness(&self, index: usize) -> Option<&A> {
        Some(&self.retained.as_ref()?.layers.get(index)?.1)
    }
}

impl<A> CoupledWorkspace<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// The solution for unit amplitude in incoming eigen-wave `wave` of
    /// `incident_side`; `None` when the workspace did not retain fields, the
    /// incident side is closed, or `wave` does not exist.
    pub fn excitation(
        &self,
        incident_side: IncidentSide,
        wave: usize,
    ) -> Option<CoupledSolution<'_, A>> {
        let retained = self.retained.as_ref()?;
        let context = self.solution.context();
        let m = context.left.outgoing.columns();
        let incident = match incident_side {
            IncidentSide::Left => &context.left,
            IncidentSide::Right => &context.right,
        };
        if wave >= m || incident.incoming.is_none() {
            return None;
        }
        let column = wave + m * usize::from(incident_side == IncidentSide::Right);
        let states = retained
            .states
            .iter()
            .map(|state| state.columns_in(column..column + 1))
            .collect();
        let from_left = incident_side == IncidentSide::Left;
        Some(CoupledSolution::new(self, states, [!from_left, from_left]))
    }

    /// The mode of the homogeneous interface system: its null vector (see
    /// [`modal`]) with its residual. Exteriors of the mode hold only their
    /// outgoing waves.
    pub fn mode(&self) -> Result<CoupledModeCandidate<A>, ModeReconstructionError> {
        self.retained
            .as_ref()
            .ok_or(ModeReconstructionError::ModeDataNotRetained)?
            .modal
            .candidate()
    }

    /// The solution of a mode of this workspace.
    pub fn mode_solution(&self, mode: &CoupledModeCandidate<A>) -> CoupledSolution<'_, A> {
        CoupledSolution::new(self, mode.states.clone(), [true, true])
    }

    /// `[E, Z₀H]`, all six components in the stack frame, at `depth` into
    /// `region` for unit amplitude in incoming eigen-wave `wave` of
    /// `incident_side` (see [`CoupledSolution::fields`]).
    pub fn fields(
        &self,
        incident_side: IncidentSide,
        wave: usize,
        region: StackRegion,
        depth: &A,
    ) -> Option<[[A; 3]; 2]> {
        self.excitation(incident_side, wave)?.fields(region, depth)
    }
}
