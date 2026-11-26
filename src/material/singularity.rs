//! Declared poles and zeros of constitutive entries in the complex k0 plane.
//!
//! A complex-frequency search assumes a characteristic that is holomorphic
//! and pole-free in its domain. A pole of any entry, or a zero of an entry a
//! polarisation divides by, breaks that assumption. Models that know where
//! these points lie declare them, so that a solver can exclude them.

/// Whether an entry is singular or vanishes at a point.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SingularityKind {
    /// The entry has a pole.
    Pole,
    /// The entry vanishes.
    Zero,
}

/// The constitutive quantity an entry belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ConstitutiveQuantity {
    /// The relative permittivity ε.
    Permittivity,
    /// The relative permeability μ.
    Permeability,
}

/// A pole or zero of one constitutive entry at a complex vacuum angular
/// wavenumber `k0` (cm⁻¹).
///
/// A declaration lists every such point of a model in the finite complex
/// plane, each location once; multiplicity is not recorded.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralSingularity<C> {
    /// Complex vacuum angular wavenumber in cm⁻¹.
    pub k0: C,
    /// Pole or zero.
    pub kind: SingularityKind,
    /// ε or μ.
    pub quantity: ConstitutiveQuantity,
    /// The diagonal entry (0, 1, 2 for xx, yy, zz in the frame of the
    /// declaration); `None` for every entry of a scalar material.
    pub axis: Option<usize>,
}

impl<C> SpectralSingularity<C> {
    /// A pole or zero of every entry of a scalar material.
    pub fn scalar(k0: C, kind: SingularityKind, quantity: ConstitutiveQuantity) -> Self {
        Self {
            k0,
            kind,
            quantity,
            axis: None,
        }
    }

    /// The same point on diagonal entry `axis`.
    pub fn on_axis(self, axis: usize) -> Self {
        Self {
            axis: Some(axis),
            ..self
        }
    }
}
