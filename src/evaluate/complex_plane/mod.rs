//! Canonical complex-plane evaluation.
//!
//! Complex-plane evaluation operates directly on canonical vacuum and
//! in-plane angular wavenumbers. Unlike real-axis evaluation, it performs no
//! caller-coordinate compilation and attaches no physical meaning to jet
//! derivative slots.
//!
//! The evaluator owns a precompiled canonical stack. Callers provide:
//!
//! - canonical complex coordinates;
//! - explicitly branch-selected exterior longitudinal wavevectors;
//! - the polarization;
//! - the jet algebra used to propagate any desired analytic derivatives.
//!
//! This path is intended for complex continuation, outgoing-mode
//! determinants, argument-principle searches, mode refinement, modal
//! reconstruction, and continuation calculations.
//!
//! [`ComplexPlaneEvaluator::checked_determinant`] checks all input and output
//! jet components for finiteness. It does not certify a search region or accept
//! a candidate as a mode. A pole-free zero search needs a holomorphic
//! characteristic on a neighborhood of its closed domain, with no zeros on
//! integration contours. Chart validity and finite samples alone cannot
//! establish this: material singularities, constitutive-factor zeros, and
//! numerical representation failures must be considered separately.
//!
//! The transfer characteristic continues analytically through regular interior
//! cutoffs. Scattering composition continues removable opposite-admittance
//! interfaces and intermediate projective poles; it can still fail at degenerate
//! directional bases or zero transmission numerators.
//! Such failures do not by themselves identify physical determinant poles.
//! Neither backend guarantees finite precision or representability for every
//! optically thick stack. Checked evaluation reports non-finite results without
//! changing the requested point, backend, or determinant gauge.

mod derivatives;
pub(super) mod mode;
mod normalization;
mod state;
mod unscaled;
pub use derivatives::*;
pub use normalization::*;

use std::fmt::Debug;

use nalgebra::ComplexField;
use ndarray::{Dimension, Ix0};
use num_traits::{Float, FromPrimitive};
use thiserror::Error;

use crate::{
    CanonicalCoordinates, ComplexPlane, ComplexScalar, ExteriorWavevectors, Polarisation, Stack,
    algebra::{Jet, ScalarAlgebra},
    backend::{Backend, PlaneWaveSolution},
    input::{
        CanonicalStack, StackCompileError, StackThicknessJet, ValidationConfig,
        compile_canonical_constant_stack,
    },
    material::TensorConstitutiveLift,
    observable::{PlaneWaveDeterminant, ProjectPlaneWaveModeDeterminant},
};

pub use mode::{ComplexPlaneMode, QnmCreationError};
pub use state::ComplexPlaneState;
pub use unscaled::*;

/// Pointwise numerical failure during checked characteristic evaluation.
///
/// This error does not classify material poles or certify search domains.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum DeterminantEvaluationError<E> {
    #[error("a coordinate or exterior-wavevector jet has non-finite components")]
    NonFiniteInput,
    #[error("determinant backend failed")]
    Backend {
        #[source]
        source: E,
    },
    #[error("the determinant jet has non-finite components")]
    NonFiniteOutput,
}

/// Evaluates canonical complex-plane problems using a precompiled stack.
///
/// `J` determines the analytic value/derivative algebra propagated through
/// coordinates, exterior wavevectors, constitutive models, and layer
/// thicknesses.
///
/// [`compile`](Self::compile) constructs a stack with constant geometric jets.
/// [`from_canonical_stack`](Self::from_canonical_stack) accepts an already
/// compiled stack and therefore permits advanced callers to seed geometry
/// derivatives explicitly.
///
/// The evaluator does not perform coordinate conversion, derivative mapping,
/// or differential-response crystallisation.
#[derive(Clone, Debug)]
pub struct ComplexPlaneEvaluator<J, M, B> {
    backend: B,
    stack: CanonicalStack<M, J>,
}

impl<J, M, B> ComplexPlaneEvaluator<J, M, B> {
    pub fn compile(
        stack: &Stack<M, <J::Scalar as ComplexField>::RealField>,
        backend: B,
    ) -> Result<Self, StackCompileError<<J::Scalar as ComplexField>::RealField>>
    where
        J: Jet<Dimension = Ix0> + StackThicknessJet,
        J::Scalar: ComplexScalar,
        <J::Scalar as ComplexField>::RealField: Float + FromPrimitive + Copy + Debug,
        J::Dimension: Dimension + Clone,
        M: Clone,
    {
        let stack = compile_canonical_constant_stack::<M, _>(
            stack,
            Ix0(),
            &ValidationConfig::permissive(),
        )?;

        Ok(Self { stack, backend })
    }

    pub fn from_canonical_stack(stack: CanonicalStack<M, J>, backend: B) -> Self {
        Self { stack, backend }
    }

    // pub(crate) fn backend(&self) -> &B {
    //     &self.backend
    // }

    #[cfg(test)]
    pub(crate) fn stack(&self) -> &CanonicalStack<M, J> {
        &self.stack
    }

    pub fn into_parts(self) -> (CanonicalStack<M, J>, B) {
        (self.stack, self.backend)
    }
}

impl<J, M, B> ComplexPlaneEvaluator<J, M, B> {
    /// Solve a canonical complex-plane problem without retaining internal
    /// finite-layer state.
    #[allow(clippy::type_complexity)]
    pub fn solve(
        &self,
        coordinates: &CanonicalCoordinates<J>,
        exterior: &ExteriorWavevectors<J>,
        polarisation: Polarisation,
    ) -> Result<
        PlaneWaveSolution<<B as Backend<J, ComplexPlane>>::Entries>,
        <B as Backend<J, ComplexPlane>>::Error,
    >
    where
        J: ScalarAlgebra + TensorConstitutiveLift<ComplexPlane, M>,
        J::Scalar: ComplexScalar,
        J::Dimension: Dimension,
        B: Backend<J, ComplexPlane>,
    {
        self.backend
            .solve(coordinates, &self.stack, exterior, polarisation)
    }

    #[allow(clippy::type_complexity)]
    pub fn determinant(
        &self,
        coordinates: &CanonicalCoordinates<J>,
        exterior: &ExteriorWavevectors<J>,
        polarisation: Polarisation,
    ) -> Result<
        RawModeDeterminant<<B as Backend<J, ComplexPlane>>::Entries>,
        <B as Backend<J, ComplexPlane>>::Error,
    >
    where
        J: ScalarAlgebra + TensorConstitutiveLift<ComplexPlane, M>,
        J::Scalar: ComplexScalar,
        J::Dimension: Dimension,
        M: Clone,
        B: Backend<J, ComplexPlane>,
        <B as Backend<J, ComplexPlane>>::Entries: ProjectPlaneWaveModeDeterminant,
    {
        Ok(self
            .solve(coordinates, exterior, polarisation)?
            .determinant())
    }

    /// Evaluate the outgoing characteristic, checking every input and output
    /// jet component for finiteness. A finite zero is a successful result.
    /// Backend errors retain their source; the raw [`Self::determinant`] API
    /// remains available without these additional checks.
    ///
    /// The caller must supply branch-consistent exterior wavevectors satisfying
    /// `kappa² = epsilon * mu * k0² - k_parallel²`, including their derivatives.
    /// This method neither checks that relation nor selects radiation sheets.
    ///
    /// Success is a **pointwise numerical check**, not root acceptance or domain
    /// certification. Pole-free argument-principle searches require a
    /// holomorphic characteristic on a neighborhood of the closed search region
    /// and zero-free integration contours. Finite contour samples cannot prove
    /// those conditions. Materials must be valid at the fixed spectral point,
    /// or throughout the mapped domain when the spectral coordinate varies.
    /// Material poles and unhandled constitutive-factor zeros must be excluded.
    ///
    /// With `xi = -i * admittance` and right-to-left transfer matrix `M`, the
    /// gauge is `D = [xi_left, -1] M [1, -xi_right]^T`. Where transmission
    /// projection is regular this equals `2 * xi_left / t_left`. Scattering
    /// transmission-chart failures and transfer overflow can be numerical
    /// failures at otherwise regular physical points. No backend substitution,
    /// coordinate perturbation, or denominator replacement is performed.
    /// Underflow and finite loss of precision are not detected.
    #[allow(clippy::type_complexity)]
    pub fn checked_determinant(
        &self,
        coordinates: &CanonicalCoordinates<J>,
        exterior: &ExteriorWavevectors<J>,
        polarisation: Polarisation,
    ) -> Result<PlaneWaveDeterminant<J>, DeterminantEvaluationError<B::Error>>
    where
        J: ScalarAlgebra + TensorConstitutiveLift<ComplexPlane, M>,
        J::Scalar: ComplexScalar,
        M: Clone,
        B: Backend<J, ComplexPlane>,
        B::Entries: ProjectPlaneWaveModeDeterminant<Determinant = PlaneWaveDeterminant<J>>,
    {
        if !coordinates.vacuum_angular_wavenumber().all_finite()
            || !coordinates.parallel_angular_wavenumber().all_finite()
            || !exterior.all(|kappa| kappa.all_finite())
        {
            return Err(DeterminantEvaluationError::NonFiniteInput);
        }
        let determinant = self
            .determinant(coordinates, exterior, polarisation)
            .map_err(|source| DeterminantEvaluationError::Backend { source })?;
        if !determinant.value().all_finite() {
            return Err(DeterminantEvaluationError::NonFiniteOutput);
        }
        Ok(determinant)
    }

    pub fn retain(
        &self,
        coordinates: CanonicalCoordinates<J>,
        exterior: ExteriorWavevectors<J>,
        polarisation: Polarisation,
    ) -> Result<ComplexPlaneState<'_, J, M, B::Workspace>, B::Error>
    where
        J: ScalarAlgebra + TensorConstitutiveLift<ComplexPlane, M>,
        J::Scalar: ComplexScalar,
        J::Dimension: Dimension,
        B: Backend<J, ComplexPlane>,
    {
        let workspace = self
            .backend
            .retain(&coordinates, &self.stack, &exterior, polarisation)?;

        Ok(ComplexPlaneState::new(
            coordinates,
            exterior,
            &self.stack,
            workspace,
            polarisation,
        ))
    }
}

pub(crate) type RawModeDeterminant<E> = <E as ProjectPlaneWaveModeDeterminant>::Determinant;

/// Share a backend error behind a trait object.
///
/// A function rather than a closure: the macros that call it are expanded for
/// Scatter2, whose error type is uninhabited, and a closure body calling
/// `Arc::new(e)` on such a value is flagged as unreachable.
pub(super) fn shared_error<E>(error: E) -> std::sync::Arc<dyn std::error::Error + Send + Sync>
where
    E: std::error::Error + Send + Sync + 'static,
{
    std::sync::Arc::new(error)
}
