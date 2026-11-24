//! Canonical terminations of a planar stack.
//!
//! A canonical termination stores every length in centimetres as a jet with
//! the sampled shape of the canonical coordinates. Positions are absolute in
//! the stack frame (`z = 0` at the first interface). The gap between the
//! outermost interface and the wall is derived at compilation:
//!
//! ```text
//! gap_L = −z_L,        gap_R = z_R − Σ d_k.
//! ```
//!
//! The right gap is a jet carrying the finite-layer thickness derivatives,
//! because the walls stay fixed in the stack frame while a layer thickens.
//! Positions are constant jets.

use crate::{Polarisation, stack::WallKind};

/// The scalar boundary condition a wall imposes on the primary field `ψ`
/// (`E_y` for TE, `H_y` for TM).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum WallCondition {
    /// `ψ = 0`.
    Dirichlet,
    /// `∂ψ/∂z = 0`.
    Neumann,
}

impl WallCondition {
    /// Return the condition imposed by a wall of `kind` in `polarisation`.
    pub(crate) const fn new(kind: WallKind, polarisation: Polarisation) -> Self {
        match (kind, polarisation) {
            (WallKind::PerfectElectric, Polarisation::TransverseElectric)
            | (WallKind::PerfectMagnetic, Polarisation::TransverseMagnetic) => Self::Dirichlet,
            (WallKind::PerfectElectric, Polarisation::TransverseMagnetic)
            | (WallKind::PerfectMagnetic, Polarisation::TransverseElectric) => Self::Neumann,
        }
    }
}

/// How one side of a canonical stack is terminated.
#[derive(Clone, Debug, PartialEq, Default)]
pub(crate) enum CanonicalTermination<J> {
    /// The semi-infinite exterior.
    #[default]
    Open,
    /// A wall at a fixed position.
    Wall {
        kind: WallKind,
        position_cm: J,
        gap_cm: J,
    },
}

impl<J> CanonicalTermination<J> {
    pub(crate) const fn is_open(&self) -> bool {
        matches!(self, Self::Open)
    }

    /// Return the gap between the outermost interface and the wall.
    pub(crate) fn gap_cm(&self) -> Option<&J> {
        match self {
            Self::Open => None,
            Self::Wall { gap_cm, .. } => Some(gap_cm),
        }
    }

    /// Return the kind of the closing wall.
    pub(crate) fn wall_kind(&self) -> Option<WallKind> {
        match self {
            Self::Open => None,
            Self::Wall { kind, .. } => Some(*kind),
        }
    }

    pub(crate) fn try_map<K, E>(
        &self,
        mut f: impl FnMut(&J) -> Result<K, E>,
    ) -> Result<CanonicalTermination<K>, E> {
        Ok(match self {
            Self::Open => CanonicalTermination::Open,
            Self::Wall {
                kind,
                position_cm,
                gap_cm,
            } => CanonicalTermination::Wall {
                kind: *kind,
                position_cm: f(position_cm)?,
                gap_cm: f(gap_cm)?,
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wall_conditions_follow_polarisation() {
        use Polarisation::{TransverseElectric as Te, TransverseMagnetic as Tm};
        use WallKind::{PerfectElectric as Pec, PerfectMagnetic as Pmc};
        assert_eq!(WallCondition::new(Pec, Te), WallCondition::Dirichlet);
        assert_eq!(WallCondition::new(Pec, Tm), WallCondition::Neumann);
        assert_eq!(WallCondition::new(Pmc, Te), WallCondition::Neumann);
        assert_eq!(WallCondition::new(Pmc, Tm), WallCondition::Dirichlet);
    }
}
