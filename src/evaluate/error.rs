use thiserror::Error;

use crate::{
    IncidentSide, ProjectionConstraintError,
    backend::{ExteriorContextProvider, PlaneWaveEntries, PlaneWaveSolutionSource},
    parameter::{DerivativeMappingError, ParameterValidationError},
};

/// Invalid real-axis evaluation request.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum SolveRequestError {
    #[error(transparent)]
    ParameterValidation(#[from] ParameterValidationError),

    #[error(transparent)]
    DerivativeMapping(#[from] DerivativeMappingError),

    #[error(transparent)]
    ProjectionConstraint(#[from] ProjectionConstraintError),

    /// A plane wave cannot be incident from a side closed by a wall or a PML.
    #[error("a plane wave cannot be incident from the {side:?} side: that side is closed")]
    ClosedIncidence { side: IncidentSide },
}

/// Refuse incidence from a closed side.
pub(crate) fn check_open_incidence<S>(
    source: &S,
    side: IncidentSide,
) -> Result<(), SolveRequestError>
where
    S: PlaneWaveSolutionSource,
    <S::Entries as PlaneWaveEntries>::ExteriorContext: ExteriorContextProvider,
{
    let solution = source.solution();
    let context = solution.context();

    let closed = match side {
        IncidentSide::Left => context.left_closure().is_some(),
        IncidentSide::Right => context.right_closure().is_some(),
    };

    if closed {
        return Err(SolveRequestError::ClosedIncidence { side });
    }

    Ok(())
}

/// Failure while compiling or solving a real-axis problem.
#[derive(Debug, Error)]
pub enum RealAxisEvaluationError<C, B> {
    #[error(transparent)]
    Request(#[from] SolveRequestError),

    #[error("failed to compile the plane-wave problem")]
    Compile {
        #[source]
        source: C,
    },

    #[error("plane-wave backend failed")]
    Backend {
        #[source]
        source: B,
    },

    #[error(transparent)]
    Constitutive(crate::material::ConstitutiveError),

    /// The principal exterior branches could not be selected.
    #[error(transparent)]
    Exterior(crate::backend::ExteriorSelectionError),
}

impl<C, B> From<crate::backend::ExteriorSelectionError> for RealAxisEvaluationError<C, B> {
    fn from(error: crate::backend::ExteriorSelectionError) -> Self {
        match error {
            crate::backend::ExteriorSelectionError::Constitutive(error) => {
                Self::Constitutive(error)
            }
            error => Self::Exterior(error),
        }
    }
}

impl<C, B> RealAxisEvaluationError<C, B> {
    pub(crate) fn compile(source: C) -> Self {
        Self::Compile { source }
    }
}
