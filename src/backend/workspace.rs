//! Backend-independent retained-workspace capabilities.
//!
//! Numerical backends may retain intermediate state required for field,
//! layer, and modal reconstruction. These traits expose only the common
//! capabilities needed by backend-independent evaluator and observable code.

use crate::{backend::LayerQuantities, material::ConstitutiveTensor};

use super::{PlaneWaveSolution, PlaneWaveSolutionSource};

/// A retained backend workspace that can be consumed into its completed
/// plane-wave solution.
pub trait SolutionWorkspace: PlaneWaveSolutionSource {
    /// Consume the retained workspace, discarding reconstruction state and
    /// returning the completed solution.
    fn into_solution(self) -> PlaneWaveSolution<Self::Entries>;
}

/// Access to retained finite layers: the stack-frame constitutive tensors,
/// the thickness jet and the backend's evaluated layer (`κ` and admittance of
/// one polarisation for the 2×2 backends, the Berreman generator and the
/// normal-field elimination for the 4×4 backends).
///
/// `None` from [`Self::retained_layer_count`] indicates that the workspace was
/// evaluated without layer retention. When retention is available, valid
/// indices are `0..retained_layer_count`.
pub trait RetainedLayers {
    type Algebra;

    /// The backend's evaluated finite layer.
    type Layer;

    /// Return the number of retained finite layers, or `None` if layer data
    /// were not retained.
    fn retained_layer_count(&self) -> Option<usize>;

    /// Return the evaluated finite layer `index`.
    fn layer(&self, index: usize) -> Option<&Self::Layer>;

    /// Return the stack-frame `[ε, μ]` of finite layer `index`.
    fn layer_tensors(&self, index: usize) -> Option<[ConstitutiveTensor<Self::Algebra>; 2]>;

    /// Return the canonical thickness jet for finite layer `index`.
    fn layer_thickness(&self, index: usize) -> Option<&Self::Algebra>;
}

/// Retained layers of one polarisation of an isotropic or diagonal stack.
pub trait RetainedIsotropicLayers:
    RetainedLayers<Layer = LayerQuantities<<Self as RetainedLayers>::Algebra>>
{
    /// Return the evaluated quantities for finite layer `index`.
    fn layer_quantities(&self, index: usize) -> Option<&LayerQuantities<Self::Algebra>> {
        self.layer(index)
    }
}

impl<W> RetainedIsotropicLayers for W where
    W: RetainedLayers<Layer = LayerQuantities<<W as RetainedLayers>::Algebra>>
{
}
