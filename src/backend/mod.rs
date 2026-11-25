//! Numerical backend interfaces and implementations.
//!
//! Backends consume canonical coordinates and canonical stacks and propagate
//! the selected jet algebra through transfer- or scattering-matrix
//! calculations.
//!
//!
//!
//! This module contains:
//!
//! - backend-independent solution and workspace contracts;
//! - exterior-wavevector and per-polarisation medium quantities;
//! - 2×2 transfer- and scattering-matrix backends for regions whose ε and μ
//!   are isotropic or diagonal in the stack frame;
//! - mode-reconstruction support shared by those backends.
//!
//! Callers normally select a concrete backend such as [`Scatter2`] or
//! [`Transfer2`] through the evaluator API rather than invoking backend
//! methods directly.

use crate::{
    CanonicalCoordinates, ComplexScalar, algebra::Jet, input::CanonicalStack,
    material::TensorConstitutiveLift,
};

use ndarray::Dimension;

mod channel;
mod exterior;
mod medium;
mod mode;
pub(crate) mod scatter2;
mod solution;
mod termination;
pub(crate) mod transfer2;
mod workspace;

pub use channel::{Channel, ChannelError};
pub(crate) use exterior::ExteriorContext;
#[cfg(test)]
pub(crate) use exterior::evaluate_exterior_wavevectors;
pub(crate) use exterior::evaluate_stack_exterior_wavevectors;
pub use exterior::{ExteriorBranch, ExteriorWavevectors};
pub(crate) use medium::{ChannelEntries, LayerQuantities, PrincipalConstitutive};
pub(crate) use mode::{
    ModalSolutionSource, ModeReconstructionError, PlaneWaveModeCandidate,
    ReconstructExteriorModeWaves, ReconstructLayerModeWaves,
};
pub use scatter2::Scatter2;
pub use scatter2::Scatter2Error;
pub(crate) use solution::{
    ExteriorContextProvider, PlaneWaveEntries, PlaneWaveSolution, PlaneWaveSolutionSource,
    PlaneWaveSolutionView,
};
#[doc(hidden)]
pub use termination::ClosedSide;
pub use termination::TerminationError;
pub(crate) use termination::{right_exterior_mode_waves, undefined_amplitudes};
pub use transfer2::{Transfer2, Transfer2Error};
pub(crate) use workspace::{RetainedIsotropicLayers, SolutionWorkspace};

/// The polarisation of a scalar channel; the 2×2 backends refuse the coupled
/// channel.
fn scalar_polarisation(channel: Channel) -> Result<crate::Polarisation, ChannelError> {
    channel
        .polarisation()
        .ok_or(ChannelError::Unsupported { channel })
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum RunMode {
    ResponseOnly,
    InternalFields,
}

impl RunMode {
    fn is_requested(&self) -> bool {
        *self == RunMode::InternalFields
    }
}

/// Numerical backend for a canonical plane-wave problem.
///
/// `J` determines the value/derivative algebra propagated through the solve,
/// while `Domain` determines how constitutive material data are evaluated
/// (for example on the real spectral axis or in the complex plane).
///
/// Backends receive the canonical problem as two borrowed components:
///
/// - canonical spectral and in-plane coordinates, which may vary between
///   solves;
/// - a canonical stack, which may be compiled once and reused.
///
/// Exterior longitudinal wavevectors are supplied separately so callers such
/// as mode solvers may control their complex-plane branch choices. `channel`
/// selects the field subspace: a backend that does not solve it is refused
/// with a [`ChannelError`].
///
/// This trait is public because it appears in bounds on public evaluator
/// implementations. Backend implementations are supplied by `lamina-core`;
/// downstream implementation is not currently a supported extension point.
#[doc(hidden)]
pub trait Backend<J, Domain>
where
    J: Jet,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    /// Backend-specific entries stored in a completed plane-wave solution.
    type Entries: PlaneWaveEntries;

    /// Retained backend state used for field and mode reconstruction.
    type Workspace: SolutionWorkspace<Entries = Self::Entries>;

    /// Backend failure type.
    type Error: std::error::Error + Send + Sync + 'static + Clone;

    /// Solve the problem and return only the completed response state.
    ///
    /// No intermediate layer state is retained.
    fn solve<M>(
        &self,
        coordinates: &CanonicalCoordinates<J>,
        stack: &CanonicalStack<M, J>,
        exterior: &ExteriorWavevectors<J>,
        channel: Channel,
    ) -> Result<PlaneWaveSolution<Self::Entries>, Self::Error>
    where
        J: TensorConstitutiveLift<Domain, M>;

    /// Solve the problem while retaining intermediate state required for
    /// subsequent field or mode reconstruction.
    fn retain<M>(
        &self,
        coordinates: &CanonicalCoordinates<J>,
        stack: &CanonicalStack<M, J>,
        exterior: &ExteriorWavevectors<J>,
        channel: Channel,
    ) -> Result<Self::Workspace, Self::Error>
    where
        J: TensorConstitutiveLift<Domain, M>;
}
