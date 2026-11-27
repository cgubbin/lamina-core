//! 4×4 transfer-matrix backend.
//!
//! [`Transfer4`] propagates the tangential state through each finite layer
//! by its entire propagator `exp(−i d K)` (scaling and squaring, no
//! eigen-decomposition), so values and derivatives are exact through layer
//! cutoffs, degenerate (isotropic) and defective layer eigenvalues. Like
//! [`Transfer2`](crate::Transfer2) it may overflow for optically thick,
//! strongly evanescent stacks, where [`Scatter4`](crate::backend::Scatter4)
//! is preferable.
//!
//! With `ψ_left = M ψ_right`, `M = L₀ L₁ … L_N`, the characteristic is
//! `det[V_L | M V_R]` of the two boundary subspaces.

use ndarray::Dimension;
use num_traits::One;

use crate::{
    CanonicalCoordinates, ComplexScalar,
    algebra::InlineJet,
    backend::{
        Channel, ChannelError, ExteriorWavevectors, RunMode, TerminationError,
        coupled::{
            CoupledError, CoupledWorkspace, Layer, Solved, accumulate, coupled_backend,
            exponential::exponential, matrix::Matrix, modal::Modal,
        },
    },
    input::CanonicalStack,
    material::{ConstitutiveError, TensorConstitutiveLift},
    stack::StackSide,
};

/// Failure of the 4×4 transfer backend.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum Transfer4Error {
    #[error(transparent)]
    Termination(#[from] TerminationError),

    #[error(transparent)]
    Constitutive(#[from] ConstitutiveError),

    #[error(transparent)]
    Channel(#[from] ChannelError),

    /// The supplied outgoing wavenumbers of a side do not span an invariant
    /// subspace of its exterior: an exterior branch point, or values that are
    /// not roots of the exterior.
    #[error("the {side:?} exterior's outgoing waves do not span an invariant subspace")]
    Exterior { side: StackSide },

    /// The transfer matrix is not finite; use the scattering backend for
    /// optically thick, strongly absorbing or strongly evanescent stacks.
    #[error(
        "the 4×4 transfer matrix is not finite; use the scattering backend for \
         optically thick, strongly absorbing, or strongly evanescent stacks"
    )]
    NonFinite,
}

impl CoupledError for Transfer4Error {
    fn exterior(side: StackSide) -> Self {
        Self::Exterior { side }
    }
}

/// 4×4 (Berreman) transfer backend for any constitutive region.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct Transfer4;

impl Transfer4 {
    /// Construct a 4×4 transfer backend.
    pub const fn new() -> Self {
        Self
    }

    fn accumulate<J, E, M>(
        &self,
        coordinates: &CanonicalCoordinates<J>,
        stack: &CanonicalStack<M, J>,
        exterior: &ExteriorWavevectors<J>,
        channel: Channel,
        request: RunMode,
    ) -> Result<CoupledWorkspace<J>, Transfer4Error>
    where
        J: InlineJet + TensorConstitutiveLift<E, M>,
        J::Scalar: ComplexScalar,
        J::Dimension: Dimension,
    {
        accumulate(
            coordinates,
            stack,
            exterior,
            channel,
            request,
            |context, layers, retain| {
                let (left, right) = (
                    &context.boundary(StackSide::Left).outgoing,
                    &context.boundary(StackSide::Right).outgoing,
                );
                let (n, m) = (left.rows(), left.columns());
                let i = J::Scalar::i();
                let propagator = |sign: J::Scalar, (region, thickness): &Layer<J::Inline>| {
                    exponential(&region.generator.scale(sign * i).scale_by(thickness))
                };
                let transfer = layers
                    .iter()
                    .fold(Matrix::identity(n, left.get(0, 0)), |transfer, layer| {
                        transfer.multiply(&propagator(-J::Scalar::one(), layer))
                    });
                if !transfer.all_finite() {
                    return Err(Transfer4Error::NonFinite);
                }

                // V_L b + incoming_L = M (V_R t + incoming_R).
                let carried = transfer.multiply(right);
                let (from_left, from_right) = context.excitations();
                let solution = left
                    .hstack(&carried.scale(-J::Scalar::one()))
                    .solve(&transfer.multiply(&from_right).subtract(&from_left));
                let outgoing = solution.rows_in(0..m);
                let system = left.hstack(&carried);
                let states = retain.then(|| {
                    // The excitations, then the left boundary subspace.
                    let start = left.multiply(&outgoing).add(&from_left).hstack(left);
                    let carried: Vec<_> = std::iter::once(start.clone())
                        .chain(layers.iter().scan(start, |state, layer| {
                            *state = propagator(J::Scalar::one(), layer).multiply(state);
                            Some(state.clone())
                        }))
                        .collect();
                    let states = carried.iter().map(|state| state.columns_in(0..2 * m));
                    let basis = carried.iter().map(|state| state.columns_in(2 * m..3 * m));
                    (
                        states.collect(),
                        Modal {
                            system: system.clone(),
                            basis: basis.collect(),
                        },
                    )
                });

                Ok(Solved {
                    characteristic: system.determinant(),
                    left: outgoing,
                    right: solution.rows_in(m..2 * m),
                    states,
                })
            },
        )
    }
}

coupled_backend!(Transfer4, Transfer4Error);
