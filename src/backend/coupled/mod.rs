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
pub(crate) mod matrix;
pub(crate) mod split;

use nalgebra::ComplexField;
use ndarray::Dimension;
use num_traits::One;

use crate::{
    ComplexScalar, FiniteLayerIndex, IncidentSide, JonesAmplitudes, PlaneWaveDeterminant,
    PlaneWavePower,
    algebra::{RealScalarAlgebra, ScalarAlgebra},
    backend::{
        Channel, ChannelError, ExteriorBranch, ExteriorSelectionError, ExteriorWavevectors,
        PlaneWaveEntries, PlaneWaveSolution, PlaneWaveSolutionSource, RunMode, SolutionWorkspace,
        TerminationError, termination::check_exterior_pattern,
    },
    input::{CanonicalCoordinates, CanonicalStack, canonical::CanonicalTermination},
    material::{ConstitutiveError, TensorConstitutiveLift},
    observable::{ProjectJonesAmplitudes, ProjectPlaneWaveModeDeterminant, ProjectPower},
    stack::{StackRegion, StackSide, WallKind},
};

use berreman::{Region, local_positions, reference_positions};
use matrix::Matrix;
use split::{Direction, Monic, Split};

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
}

impl<A> Boundary<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    fn new<E, M, Error: CoupledError>(
        side: StackSide,
        coordinates: &CanonicalCoordinates<A>,
        material: &M,
        termination: &CanonicalTermination<A>,
        branch: Option<&ExteriorBranch<A>>,
        channel: Channel,
    ) -> Result<Self, Error>
    where
        A: TensorConstitutiveLift<E, M>,
    {
        let name = match side {
            StackSide::Left => StackRegion::LeftExterior,
            StackSide::Right => StackRegion::RightExterior,
        };
        let region = Region::evaluate::<E, M, Error>(name, material, coordinates, channel)?;
        let generator = &region.generator;
        let like = generator.get(0, 0);

        let Some(branch) = branch else {
            let outgoing = closed_subspace::<A, E, M, Error>(
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
            });
        };

        let outward = match (branch, channel) {
            (ExteriorBranch::Pair(pair), Channel::Coupled) => pair.to_vec(),
            (ExteriorBranch::Scalar(kappa), _) if channel != Channel::Coupled => {
                vec![kappa.clone()]
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
        })
    }
}

/// The states a closed side admits at its outermost interface.
fn closed_subspace<A, E, M, Error: CoupledError>(
    side: StackSide,
    name: StackRegion,
    coordinates: &CanonicalCoordinates<A>,
    material: &M,
    termination: &CanonicalTermination<A>,
    region: &Region<A>,
    channel: Channel,
) -> Result<Matrix<A>, Error>
where
    A: ScalarAlgebra + TensorConstitutiveLift<E, M>,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
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
    let i = A::Scalar::i();

    // Right: ψ = exp(−iℓ(K − n_ref)) w; left: ψ = exp(iℓ(K + n_ref)) w.
    let mut exponent = match side {
        StackSide::Right => generator.scale(-i),
        StackSide::Left => generator.scale(i),
    };
    let mut length = gap.clone();
    if let Some(pml) = termination.pml() {
        let at_normal = CanonicalCoordinates::new(
            coordinates.vacuum_angular_wavenumber().clone(),
            like.zero_like(),
        );
        let normal = Region::evaluate::<E, M, Error>(name, material, &at_normal, channel)?;
        let e2 = &Monic::characteristic(&normal.generator).0[n - 2];
        let reference = e2.scale(nalgebra::convert(-2.0 / n as f64)).sqrt();
        exponent = exponent.shift(&reference.scale(i));
        length = length.add(&pml.optical_thickness_cm);
    }
    Ok(exponential::exponential(&exponent.scale_by(&length))
        .multiply(&Matrix::unit_columns(n, &admitted, like)))
}

/// The exterior data of a solved stack: the channel and each side's boundary
/// subspace.
#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledContext<A> {
    channel: Channel,
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

impl<A> CoupledContext<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
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
    /// The state at every interface, `n × 2m`, when retained.
    pub(crate) states: Option<Vec<Matrix<A>>>,
}

/// A finite layer: its region and thickness.
pub(crate) type Layer<A> = (Region<A>, A);

/// Set up `stack` and solve it with `compose`, which composes the layers
/// between the two boundary subspaces and retains the interface states when
/// asked.
pub(crate) fn accumulate<J, E, M, Error: CoupledError>(
    coordinates: &CanonicalCoordinates<J>,
    stack: &CanonicalStack<M, J>,
    exterior: &ExteriorWavevectors<J>,
    channel: Channel,
    request: RunMode,
    compose: impl FnOnce(&CoupledContext<J>, &[Layer<J>], bool) -> Result<Solved<J>, Error>,
) -> Result<CoupledWorkspace<J>, Error>
where
    J: ScalarAlgebra + TensorConstitutiveLift<E, M>,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    check_exterior_pattern(stack, exterior)?;
    let boundary = |side, material, termination, branch| {
        Boundary::new::<E, M, Error>(side, coordinates, material, termination, branch, channel)
    };
    let context = CoupledContext {
        channel,
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
    let retained = solved.states.map(|states| Retained { layers, states });

    Ok(CoupledWorkspace {
        solution: PlaneWaveSolution::new(entries, context),
        retained,
    })
}

/// The finite layers of a stack for a channel, with their thicknesses.
pub(crate) fn layer_regions<A, E, M, Error>(
    coordinates: &CanonicalCoordinates<A>,
    stack: &CanonicalStack<M, A>,
    channel: Channel,
) -> Result<Vec<Layer<A>>, Error>
where
    A: ScalarAlgebra + TensorConstitutiveLift<E, M>,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
    Error: From<ConstitutiveError> + From<ChannelError>,
{
    stack
        .layers()
        .iter()
        .enumerate()
        .map(|(index, layer)| {
            let name = StackRegion::Layer(FiniteLayerIndex::new(index));
            Region::evaluate::<E, M, Error>(name, layer.material(), coordinates, channel)
                .map(|region| (region, layer.thickness_cm().clone()))
        })
        .collect()
}

/// Implement [`Backend`](crate::Backend) for a 4×4 backend with an
/// `accumulate` method.
macro_rules! coupled_backend {
    ($backend:ty, $error:ty) => {
        impl<J, Domain> $crate::backend::Backend<J, Domain> for $backend
        where
            J: $crate::algebra::ScalarAlgebra,
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
                    coordinates,
                    stack,
                    exterior,
                    channel,
                    $crate::backend::RunMode::ResponseOnly,
                )
                .map(SolutionWorkspace::into_solution)
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
                    coordinates,
                    stack,
                    exterior,
                    channel,
                    $crate::backend::RunMode::InternalFields,
                )
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
}

impl<A> PlaneWaveSolutionSource for CoupledWorkspace<A> {
    type Entries = CoupledEntries<A>;

    fn solution(&self) -> crate::backend::PlaneWaveSolutionView<'_, Self::Entries> {
        self.solution.as_view()
    }
}

impl<A> SolutionWorkspace for CoupledWorkspace<A> {
    fn into_solution(self) -> PlaneWaveSolution<CoupledEntries<A>> {
        self.solution
    }
}

impl<A> CoupledWorkspace<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// `[E, Z₀H]`, all six components in the stack frame, at `depth` into
    /// `region` for unit amplitude in incoming eigen-wave `wave` of
    /// `incident_side`.
    ///
    /// The depth is measured from a finite layer's left interface, and
    /// outward from an exterior's interface. `None` when the workspace did
    /// not retain fields, the incident side is closed, or `region` or `wave`
    /// does not exist. Fields are propagated from the interface state by the
    /// entire layer propagator.
    pub fn fields(
        &self,
        incident_side: IncidentSide,
        wave: usize,
        region: StackRegion,
        depth: &A,
    ) -> Option<[[A; 3]; 2]> {
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
        let (region, interface, depth) = match region {
            StackRegion::Layer(layer) => (
                &retained.layers.get(layer.get())?.0,
                layer.get(),
                depth.clone(),
            ),
            StackRegion::LeftExterior => (&context.left.region, 0, depth.negate()),
            StackRegion::RightExterior => {
                (&context.right.region, retained.layers.len(), depth.clone())
            }
        };
        let column = wave + m * usize::from(incident_side == IncidentSide::Right);
        let state =
            exponential::exponential(&region.generator.scale(A::Scalar::i()).scale_by(&depth))
                .multiply(&retained.states[interface].columns_in(column..column + 1));

        // The full tangential state, with zeros outside the channel.
        let components = berreman::channel_components(context.channel);
        let full: [A; 4] = std::array::from_fn(|component| {
            components
                .iter()
                .position(|c| *c == component)
                .map_or_else(|| depth.zero_like(), |row| state.get(row, 0).clone())
        });
        let dot = |row: &berreman::Row<A>| {
            row.iter()
                .zip(&full)
                .map(|(a, b)| a.multiply(b))
                .reduce(|sum, term| sum.add(&term))
                .unwrap_or_else(|| depth.zero_like())
        };
        let [e_x, eta_y, e_y, minus_eta_x] = full.clone();
        Some([
            [e_x, e_y, dot(&region.elimination.electric)],
            [
                minus_eta_x.negate(),
                eta_y,
                dot(&region.elimination.magnetic),
            ],
        ])
    }
}
