//! Narrow supported-backend characteristic evaluation for local derivatives.
use super::ComplexPlaneEvaluator;
use crate::{
    AnalyticalMaterialHandle, CanonicalCoordinates, ExteriorWavevectors, ModeJet1,
    PlaneWaveDeterminant, Polarisation,
    backend::{Scatter2, Transfer2},
    input::CanonicalStack,
};
use ndarray::Ix0;
use num_complex::Complex64 as C;
use std::{error::Error, sync::Arc};

type J = ModeJet1<C, Ix0>;

/// Checked characteristic evaluation on a caller-supplied canonical stack.
///
/// Implemented only by the built-in isotropic modal backends. This keeps
/// derivative postprocessing explicit about its numerical backend while
/// allowing a selected finite-layer thickness jet to reuse core compilation.
pub trait ModalDerivativeBackend: Clone {
    fn checked_characteristic(
        &self,
        stack: crate::CanonicalStack<AnalyticalMaterialHandle<C>, J>,
        coordinates: &CanonicalCoordinates<J>,
        exterior: &ExteriorWavevectors<J>,
        polarisation: Polarisation,
    ) -> Result<J, Arc<dyn Error + Send + Sync>>;
}

macro_rules! implementation {
    ($backend:ty) => {
        impl ModalDerivativeBackend for $backend {
            fn checked_characteristic(
                &self,
                stack: CanonicalStack<AnalyticalMaterialHandle<C>, J>,
                coordinates: &CanonicalCoordinates<J>,
                exterior: &ExteriorWavevectors<J>,
                polarisation: Polarisation,
            ) -> Result<J, Arc<dyn Error + Send + Sync>> {
                let evaluator = ComplexPlaneEvaluator::from_canonical_stack(stack, self.clone());
                evaluator
                    .checked_determinant(coordinates, exterior, polarisation)
                    .map(PlaneWaveDeterminant::into_inner)
                    .map_err(|error| Arc::new(error) as Arc<dyn Error + Send + Sync>)
            }
        }
    };
}

implementation!(Scatter2);
implementation!(Transfer2);
