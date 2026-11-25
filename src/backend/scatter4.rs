//! 4×4 scattering-matrix backend.
//!
//! [`Scatter4`] writes the field in each finite layer as forward waves,
//! referenced at the layer's left boundary, and backward waves, referenced at
//! its right boundary:
//!
//! ```text
//! ψ(z) = W⁺ exp(i(z − z_L)K⁺) a + W⁻ exp(i(z − z_R)K⁻) b,
//! ```
//!
//! with `W^±` the spectral projections of the reference vectors onto the
//! forward (`Im n > 0`, or real with `S_z > 0`) and backward invariant
//! subspaces, and `K^±` the generator restricted to them. Only the bounded
//! `exp(i d K⁺)` and `exp(−i d K⁻)` occur. Continuity at each interface is
//! eliminated from left to right, a Redheffer recursion for the left
//! reflection, and the characteristic is the product of the pivot
//! determinants of that elimination: the determinant of the whole interface
//! system in bounded amplitudes. A layer whose forward and backward
//! eigenvalues meet (a layer cutoff or a defective generator) has no split;
//! it is refused with [`Scatter4Error::DefectiveLayer`], and
//! [`Transfer4`](crate::backend::Transfer4) is authoritative there.
//!
//! The split is chosen at the real point `(Re k₀, Re k∥)` — by `Im n`, or the
//! sign of `S_z` for a real `n` — and continued to the evaluation point by
//! matching each eigenvalue to the nearest one there. It is therefore
//! continuous over a thin region about the real axis while no forward
//! eigenvalue of a layer meets a backward one.

use nalgebra::ComplexField;
use ndarray::Dimension;
use num_traits::One;

use crate::{
    CanonicalCoordinates, ComplexScalar, FiniteLayerIndex,
    algebra::ScalarAlgebra,
    backend::{
        Channel, ChannelError, ExteriorWavevectors, RunMode, TerminationError,
        coupled::{
            CoupledError, CoupledWorkspace, Solved, accumulate,
            berreman::{Region, reference_positions},
            coupled_backend,
            exponential::exponential,
            layer_regions,
            matrix::Matrix,
            split::{Direction, Split},
        },
    },
    input::CanonicalStack,
    material::{ConstitutiveError, TensorConstitutiveLift},
    stack::StackSide,
};

/// Failure of the 4×4 scattering backend.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum Scatter4Error {
    #[error(transparent)]
    Termination(#[from] TerminationError),

    #[error(transparent)]
    Constitutive(#[from] ConstitutiveError),

    #[error(transparent)]
    Channel(#[from] ChannelError),

    /// The supplied outgoing wavenumbers of a side do not span an invariant
    /// subspace of its exterior.
    #[error("the {side:?} exterior's outgoing waves do not span an invariant subspace")]
    Exterior { side: StackSide },

    /// A layer's forward and backward eigenvalues do not split: a layer
    /// cutoff or a defective generator. The transfer backend is exact there.
    #[error(
        "layer {layer:?} has no forward/backward split of its eigenvalues (a layer \
         cutoff or a defective generator); use the transfer backend there"
    )]
    DefectiveLayer { layer: FiniteLayerIndex },
}

impl CoupledError for Scatter4Error {
    fn exterior(side: StackSide) -> Self {
        Self::Exterior { side }
    }
}

/// 4×4 (Berreman) scattering backend for any constitutive region.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct Scatter4;

/// One finite layer in eigen-wave form.
struct EigenLayer<A> {
    forward: Matrix<A>,
    backward: Matrix<A>,
    /// `exp(i d K⁺)`.
    forward_phase: Matrix<A>,
    /// `exp(−i d K⁻)`.
    backward_phase: Matrix<A>,
}

impl<A> EigenLayer<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// The layer `region` of `thickness`, split as `reference` is at the
    /// real point.
    fn new(
        region: &Region<A>,
        reference: &Region<A>,
        thickness: &A,
        channel: Channel,
    ) -> Option<Self> {
        let generator = &region.generator;
        let split = Split::select(generator, &reference.generator, Direction::Forward)?;
        let reference = reference_positions(channel);
        let basis = Matrix::unit_columns(generator.rows(), &reference, generator.get(0, 0));
        let forward = split.projector.multiply(&basis);
        let backward = basis.subtract(&forward);
        // exp(±i d K_W), with K W = W K_W restricted to the reference rows.
        let all: Vec<usize> = (0..reference.len()).collect();
        let phase = |waves: &Matrix<A>, sign: A::Scalar| {
            let restricted = waves
                .select(&reference, &all)
                .solve(&generator.multiply(waves).select(&reference, &all));
            exponential(&restricted.scale(sign * A::Scalar::i()).scale_by(thickness))
        };

        Some(Self {
            forward_phase: phase(&forward, A::Scalar::one()),
            backward_phase: phase(&backward, -A::Scalar::one()),
            forward,
            backward,
        })
    }
}

/// Eliminate the interface system from left to right, then back-substitute.
///
/// At each interface the state on its left is `U b + Y`, with `b` the
/// backward amplitudes of the region on the left (referenced at this
/// interface) and `Y` the contribution of the left incoming terms, one column
/// per excitation. Eliminating continuity into the next layer gives
/// `[b; a] = Z [c; 1]`, with `a` and `c` the forward and backward amplitudes
/// of the layer at its left boundary.
fn sweep<A>(
    context: &crate::backend::coupled::CoupledContext<A>,
    layers: &[EigenLayer<A>],
    retain: bool,
) -> Solved<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let (left, right) = (
        &context.boundary(StackSide::Left).outgoing,
        &context.boundary(StackSide::Right).outgoing,
    );
    let (left_incoming, right_incoming) = context.excitations();
    let (m, excitations) = (left.columns(), left_incoming.columns());
    let minus_one = -A::Scalar::one();

    let mut characteristic = left.get(0, 0).one_like();
    let mut state = (left.clone(), left_incoming);
    // Each step: Z, exp(−i d K⁻) of its layer, and the state on its left.
    let mut steps = Vec::with_capacity(layers.len());
    for layer in layers {
        let pivot = state.0.hstack(&layer.forward.scale(minus_one));
        characteristic = characteristic.multiply(&pivot.determinant());
        let elimination = pivot.solve(&layer.backward.hstack(&state.1.scale(minus_one)));
        let lower = elimination.rows_in(m..2 * m);
        let carried = layer.forward.multiply(&layer.forward_phase);
        let next = (
            carried
                .multiply(&lower.columns_in(0..m))
                .multiply(&layer.backward_phase)
                .add(&layer.backward),
            carried.multiply(&lower.columns_in(m..m + excitations)),
        );
        steps.push((
            elimination,
            &layer.backward_phase,
            std::mem::replace(&mut state, next),
        ));
    }

    let pivot = state.0.hstack(&right.scale(minus_one));
    characteristic = characteristic.multiply(&pivot.determinant());
    let last = pivot.solve(&right_incoming.subtract(&state.1));
    let mut backward = last.rows_in(0..m);
    let mut states = vec![state.0.multiply(&backward).add(&state.1)];
    let identity = Matrix::identity(excitations, left.get(0, 0));
    for (elimination, backward_phase, (u, y)) in steps.iter().rev() {
        let amplitudes = backward_phase.multiply(&backward).vstack(&identity);
        backward = elimination.multiply(&amplitudes).rows_in(0..m);
        states.push(u.multiply(&backward).add(y));
    }
    states.reverse();

    Solved {
        characteristic,
        left: backward,
        right: last.rows_in(m..2 * m),
        states: retain.then_some(states),
    }
}

impl Scatter4 {
    /// Construct a 4×4 scattering backend.
    pub const fn new() -> Self {
        Self
    }

    fn accumulate<J, E, M>(
        coordinates: &CanonicalCoordinates<J>,
        stack: &CanonicalStack<M, J>,
        exterior: &ExteriorWavevectors<J>,
        channel: Channel,
        request: RunMode,
    ) -> Result<CoupledWorkspace<J>, Scatter4Error>
    where
        J: ScalarAlgebra + TensorConstitutiveLift<E, M>,
        J::Scalar: ComplexScalar,
        J::Dimension: Dimension,
    {
        let real = |jet: &J| {
            J::lift_constant(jet.value().mapv(|value| J::Scalar::from_real(value.real())))
        };
        let reference = CanonicalCoordinates::new(
            real(coordinates.vacuum_angular_wavenumber()),
            real(coordinates.parallel_angular_wavenumber()),
        );
        let references = layer_regions::<J, E, M, Scatter4Error>(&reference, stack, channel)?;

        accumulate(
            coordinates,
            stack,
            exterior,
            channel,
            request,
            |context, layers, retain| {
                let eigen = layers
                    .iter()
                    .zip(&references)
                    .enumerate()
                    .map(|(index, ((region, thickness), (reference, _)))| {
                        EigenLayer::new(region, reference, thickness, channel).ok_or(
                            Scatter4Error::DefectiveLayer {
                                layer: FiniteLayerIndex::new(index),
                            },
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(sweep(context, &eigen, retain))
            },
        )
    }
}

coupled_backend!(Scatter4, Scatter4Error);
