use thiserror::Error;

use crate::algebra::UnsupportedDerivativeSlot;
use crate::input::compile::coordinates::{
    CoordinateVariable,
    in_plane::{InPlaneCanonicalisationError, InPlaneInputError},
    spectral::SpectralInputError,
};

#[derive(Debug, Error)]
pub enum CoordinateCompileError<R> {
    #[error(transparent)]
    Spectral(#[from] SpectralInputError<R>),

    #[error(transparent)]
    InPlane(#[from] InPlaneInputError<R>),

    #[error(transparent)]
    InPlaneCanonicalisation(#[from] InPlaneCanonicalisationError),

    #[error("failed to seed {variable:?}: {source}")]
    Seed {
        variable: CoordinateVariable,

        #[source]
        source: UnsupportedDerivativeSlot,
    },

    #[error("incident side must be specified to compile with an incident angle")]
    MissingIncidentSide,

    #[error("incident-angle coordinates are not supported for complex input")]
    ComplexIncidentAngleUnsupported,

    /// The incident exterior's refractive index is not a scalar, or cannot be
    /// evaluated; an incident angle is defined for isotropic exteriors only.
    #[error(transparent)]
    IncidentMedium(#[from] crate::material::ConstitutiveError),
}
