use thiserror::Error;

use crate::{
    algebra::{ScalarAlgebra, ScaleBy},
    observable::BoundaryState,
    waves::{ExteriorBoundaryWaves, LayerBoundaryWaves},
};

#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ModeReconstructionError {
    #[error("workspace does not retain the data required for modal reconstruction")]
    ModeDataNotRetained,

    #[error("the outgoing boundary system has no usable modal null vector")]
    NoUsableNullVector,

    #[error("the outgoing boundary system has a degenerate modal null space")]
    DegenerateNullSpace,

    #[error(
        "reconstructed modal boundary-wave count {wave_count} does not match \
         retained finite-layer count {layer_count}"
    )]
    LayerCountMismatch {
        wave_count: usize,
        layer_count: usize,
    },
}

/// Unnormalised candidate solution of the outgoing homogeneous boundary
/// problem.
///
/// A candidate contains:
///
/// - a boundary state representing one chosen null-space gauge;
/// - the corresponding right-exterior outgoing amplitude;
/// - a projective residual measuring failure to satisfy the homogeneous
///   boundary condition exactly.
///
/// The overall complex scale is arbitrary. Physical modal normalisation is
/// applied only after reconstruction.
#[doc(hidden)]
#[derive(Clone, Debug)]
pub struct PlaneWaveModeCandidate<A> {
    state: BoundaryState<A>,
    right_outgoing: A,
    residual: A,
    right_closed_state: Option<BoundaryState<A>>,
}

impl<A> PlaneWaveModeCandidate<A> {
    pub(crate) fn new(state: BoundaryState<A>, right_outgoing: A, residual: A) -> Self {
        Self {
            state,
            right_outgoing,
            residual,
            right_closed_state: None,
        }
    }

    /// Attach the state at the outermost right interface of a closed right
    /// side, where the exterior holds a standing wave rather than only an
    /// outgoing one.
    pub(crate) fn with_right_closed_state(mut self, state: Option<BoundaryState<A>>) -> Self {
        self.right_closed_state = state;
        self
    }

    /// Return the state at the outermost right interface when the right side
    /// is closed.
    pub(crate) fn right_closed_state(&self) -> Option<&BoundaryState<A>> {
        self.right_closed_state.as_ref()
    }

    pub(crate) fn state(&self) -> &BoundaryState<A> {
        &self.state
    }

    pub(crate) fn right_outgoing(&self) -> &A {
        &self.right_outgoing
    }

    #[cfg(test)]
    pub(crate) fn residual(&self) -> &A {
        &self.residual
    }

    pub(crate) fn into_residual(self) -> A {
        self.residual
    }

    pub(crate) fn scaled(mut self, factor: &A) -> Self
    where
        A: ScalarAlgebra,
    {
        self.right_outgoing = self.right_outgoing.multiply(factor);
        self.residual = self.residual.multiply(factor);
        self.state = self.state.scale_by(factor);
        self.right_closed_state = self.right_closed_state.map(|state| state.scale_by(factor));

        self
    }
}

/// Reconstruct exterior directional waves for an outgoing homogeneous mode.
///
/// The mode has no incident side. The returned exterior waves must satisfy
/// outgoing boundary conditions on both sides:
///
/// - no right-going incident wave in the left exterior;
/// - no left-going incident wave in the right exterior.
///
/// The reconstructed mode retains the arbitrary complex amplitude carried by
/// `seed`. Modal normalization is applied later.
pub trait ReconstructExteriorModeWaves {
    type Algebra;

    fn reconstruct_exterior_mode_waves(
        &self,
        seed: &PlaneWaveModeCandidate<Self::Algebra>,
    ) -> Result<ExteriorBoundaryWaves<Self::Algebra>, ModeReconstructionError>;
}

/// Reconstruct directional waves in every finite layer for an outgoing
/// homogeneous modal solution.
///
/// The returned vector contains one entry per retained finite layer in
/// physical left-to-right order.
///
/// The reconstruction preserves the arbitrary complex scale carried by
/// `seed`; modal normalisation is applied later.
pub trait ReconstructLayerModeWaves {
    type Algebra;

    fn reconstruct_layer_mode_waves(
        &self,
        seed: &PlaneWaveModeCandidate<Self::Algebra>,
    ) -> Result<Vec<LayerBoundaryWaves<Self::Algebra>>, ModeReconstructionError>;
}

/// Source of an unnormalised outgoing homogeneous boundary solution.
///
/// Implemented by retained backend workspaces that contain sufficient data to
/// construct a representative vector in the outgoing boundary-system null
/// space.
pub trait ModalSolutionSource {
    type Algebra;

    /// The unnormalised mode: a [`PlaneWaveModeCandidate`] for the 2×2
    /// backends, the interface states for the 4×4 ones.
    type Mode;

    /// Construct an arbitrary nonzero representative of the outgoing
    /// homogeneous boundary solution.
    fn modal_boundary_solution(&self) -> Result<Self::Mode, ModeReconstructionError>;
}
