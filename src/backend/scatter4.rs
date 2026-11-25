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
//! By default the split is chosen at the real point `(Re k₀, Re k∥)` — by
//! `Im n`, or the sign of `S_z` for a real `n` — and continued to the
//! evaluation point by matching each eigenvalue to the nearest one there. It
//! is therefore continuous over a thin region about the real axis while no
//! forward eigenvalue of a layer meets a backward one.
//! [`Scatter4::continued_from`] instead chooses it at a fixed reference point
//! and tracks the eigenvalues along the segment to each evaluation point, in
//! steps short against the distance between the forward and the backward
//! eigenvalues: the characteristic is then holomorphic over a region
//! star-shaped about the reference point that holds no layer branch point,
//! and an evaluation whose segment comes near one is refused
//! ([`Scatter4Error::Continuation`]).

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
            matrix::sample_matrices,
            modal::Modal,
            split::{Direction, Split, classify, matched, separation, try_f64},
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

    /// A forward and a backward eigenvalue of a layer approach between the
    /// reference point of a continued split and the evaluation point: a layer
    /// branch point near the region. The transfer backend is exact there.
    #[error(
        "the forward/backward split of layer {layer:?} cannot be continued from the \
         reference point; use the transfer backend there"
    )]
    Continuation { layer: FiniteLayerIndex },
}

impl CoupledError for Scatter4Error {
    fn exterior(side: StackSide) -> Self {
        Self::Exterior { side }
    }
}

/// 4×4 (Berreman) scattering backend for any constitutive region.
///
/// By default each layer's split is classified at the real point of each
/// evaluation. [`Scatter4::continued_from`] fixes it at one reference point
/// instead, for a characteristic holomorphic over a whole chart region.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Scatter4 {
    /// The reference point `(k₀, k∥)` of a continued split.
    reference: Option<[num_complex::Complex64; 2]>,
}

/// The largest step of the continuation, as a fraction of the path.
const LARGEST_STEP: f64 = 0.25;

/// The smallest step before a continuation is refused.
const SMALLEST_STEP: f64 = 1.0 / 4096.0;

/// A step is accepted when no root moves by more than this fraction of the
/// distance between the forward and the backward roots.
const STEP_FRACTION: f64 = 0.25;

/// Relative distance between a forward and a backward root at or below which
/// a continued split is refused.
const SEPARATION_RATIO: f64 = 1.0e-4;

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
    /// The layer `region` of `thickness` with forward waves from `split`.
    fn new(region: &Region<A>, split: &Split<A>, thickness: &A, channel: Channel) -> Self {
        let generator = &region.generator;
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

        Self {
            forward_phase: phase(&forward, A::Scalar::one()),
            backward_phase: phase(&backward, -A::Scalar::one()),
            forward,
            backward,
        }
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
    // The homogeneous solution: backward amplitudes of the last layer `b`
    // carried back with no excitation, as maps from `b`.
    let mut homogeneous = Matrix::identity(m, left.get(0, 0));
    let mut basis = vec![state.0.clone()];
    let identity = Matrix::identity(excitations, left.get(0, 0));
    for (elimination, backward_phase, (u, y)) in steps.iter().rev() {
        let amplitudes = backward_phase.multiply(&backward).vstack(&identity);
        backward = elimination.multiply(&amplitudes).rows_in(0..m);
        states.push(u.multiply(&backward).add(y));
        homogeneous = elimination
            .columns_in(0..m)
            .multiply(backward_phase)
            .multiply(&homogeneous)
            .rows_in(0..m);
        basis.push(u.multiply(&homogeneous));
    }
    states.reverse();
    basis.reverse();

    Solved {
        characteristic,
        left: backward,
        right: last.rows_in(m..2 * m),
        states: retain.then_some((
            states,
            Modal {
                system: pivot,
                basis,
            },
        )),
    }
}

impl Scatter4 {
    /// Construct a 4×4 scattering backend.
    pub const fn new() -> Self {
        Self { reference: None }
    }

    /// A backend that splits every layer at the reference point
    /// `(k₀, k∥)` and continues the split to each evaluation point along
    /// the straight segment between them, tracking the eigenvalues: the
    /// characteristic is then holomorphic on any region star-shaped about the
    /// reference point in which no forward eigenvalue of a layer meets a
    /// backward one. A layer whose split has no classification at the
    /// reference point is refused with [`Scatter4Error::DefectiveLayer`], one
    /// whose forward and backward eigenvalues approach along the segment with
    /// [`Scatter4Error::Continuation`].
    pub const fn continued_from(
        vacuum_angular_wavenumber: num_complex::Complex64,
        parallel_angular_wavenumber: num_complex::Complex64,
    ) -> Self {
        Self {
            reference: Some([vacuum_angular_wavenumber, parallel_angular_wavenumber]),
        }
    }

    fn accumulate<J, E, M>(
        &self,
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
        let defective = |index| Scatter4Error::DefectiveLayer {
            layer: FiniteLayerIndex::new(index),
        };
        let splits: LayerSplits<'_, J> = match self.reference {
            None => {
                let real = |jet: &J| {
                    J::lift_constant(jet.value().mapv(|value| J::Scalar::from_real(value.real())))
                };
                let reference = CanonicalCoordinates::new(
                    real(coordinates.vacuum_angular_wavenumber()),
                    real(coordinates.parallel_angular_wavenumber()),
                );
                let references =
                    layer_regions::<J, E, M, Scatter4Error>(&reference, stack, channel)?;
                Box::new(move |index, region: &Region<J>| {
                    Split::select(
                        &region.generator,
                        &references[index].0.generator,
                        Direction::Forward,
                    )
                    .ok_or(defective(index))
                })
            }
            Some(reference) => {
                let roots = continued_roots::<J, E, M>(reference, coordinates, stack, channel)?;
                Box::new(move |index, region: &Region<J>| {
                    Split::from_ordered(&region.generator, &roots[index]).ok_or(defective(index))
                })
            }
        };

        accumulate(
            coordinates,
            stack,
            exterior,
            channel,
            request,
            |context, layers, retain| {
                let eigen = layers
                    .iter()
                    .enumerate()
                    .map(|(index, (region, thickness))| {
                        let split = splits(index, region)?;
                        Ok(EigenLayer::new(region, &split, thickness, channel))
                    })
                    .collect::<Result<Vec<_>, Scatter4Error>>()?;
                Ok(sweep(context, &eigen, retain))
            },
        )
    }
}

/// The split of finite layer `index`.
type LayerSplits<'a, J> = Box<dyn Fn(usize, &Region<J>) -> Result<Split<J>, Scatter4Error> + 'a>;

/// Each layer's eigenvalues at every sample of `coordinates`, forward first,
/// classified at `reference` and tracked along the segment from it.
#[allow(clippy::type_complexity)]
fn continued_roots<J, E, M>(
    [k0, k_parallel]: [num_complex::Complex64; 2],
    coordinates: &CanonicalCoordinates<J>,
    stack: &CanonicalStack<M, J>,
    channel: Channel,
) -> Result<Vec<Vec<Vec<J::Scalar>>>, Scatter4Error>
where
    J: ScalarAlgebra + TensorConstitutiveLift<E, M>,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    let scalar = |value: num_complex::Complex64| {
        J::Scalar::from_real(nalgebra::convert(value.re))
            + J::Scalar::i() * J::Scalar::from_real(nalgebra::convert(value.im))
    };
    let (k0, k_parallel) = (scalar(k0), scalar(k_parallel));
    // The generators of every layer, sample by sample, at fraction `t` of
    // the segment.
    let generators = |t: f64| {
        let t = J::Scalar::from_real(nalgebra::convert(t));
        let along = |jet: &J, start: J::Scalar| {
            J::lift_constant(jet.value().mapv(|end| start + (end - start) * t))
        };
        let point = CanonicalCoordinates::new(
            along(coordinates.vacuum_angular_wavenumber(), k0),
            along(coordinates.parallel_angular_wavenumber(), k_parallel),
        );
        layer_regions::<J, E, M, Scatter4Error>(&point, stack, channel).map(|layers| {
            layers
                .iter()
                .map(|(region, _)| sample_matrices(&region.generator))
                .collect::<Vec<_>>()
        })
    };
    let stuck = |layer| Scatter4Error::Continuation {
        layer: FiniteLayerIndex::new(layer),
    };

    let mut roots = generators(0.0)?
        .iter()
        .enumerate()
        .map(|(layer, samples)| {
            samples
                .iter()
                .map(|sample| {
                    classify(sample, Direction::Forward).ok_or(Scatter4Error::DefectiveLayer {
                        layer: FiniteLayerIndex::new(layer),
                    })
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let (mut t, mut step) = (0.0, LARGEST_STEP);
    while t < 1.0 {
        let next = (t + step).min(1.0);
        let mut tracked = Vec::with_capacity(roots.len());
        let mut too_long = None;
        'layers: for (layer, (samples, previous)) in
            generators(next)?.iter().zip(&roots).enumerate()
        {
            let mut layer_roots = Vec::with_capacity(samples.len());
            for (sample, previous) in samples.iter().zip(previous) {
                let now = matched(sample, previous).ok_or(stuck(layer))?;
                let moved = now
                    .iter()
                    .zip(previous)
                    .map(|(a, b)| try_f64((*a - *b).modulus()))
                    .fold(0.0, f64::max);
                if moved.is_nan() || moved > STEP_FRACTION * separation(previous) {
                    too_long = Some(layer);
                    break 'layers;
                }
                let scale = now
                    .iter()
                    .map(|root| try_f64(root.modulus()))
                    .fold(f64::MIN_POSITIVE, f64::max);
                let gap = separation(&now);
                if gap.is_nan() || gap <= SEPARATION_RATIO * scale {
                    return Err(stuck(layer));
                }
                layer_roots.push(now);
            }
            tracked.push(layer_roots);
        }
        match too_long {
            None => {
                roots = tracked;
                t = next;
                step = (2.0 * step).min(LARGEST_STEP);
            }
            Some(layer) => {
                step /= 2.0;
                if step < SMALLEST_STEP {
                    return Err(stuck(layer));
                }
            }
        }
    }
    Ok(roots)
}

coupled_backend!(Scatter4, Scatter4Error);
