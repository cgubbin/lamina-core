//! Terminations closing a planar stack on either side.
//!
//! By default a stack is bounded by two semi-infinite exterior media. Either
//! side may instead be closed
//!
//! by a perfect electric (PEC) or perfect magnetic (PMC) [`Wall`] at an
//! absolute position `z` in the stack's own frame.
//!
//! The stack frame places `z = 0` at the first interface, with the layers
//! running towards `+z` and ending at `z = L`, the sum of the finite-layer
//! thicknesses. A left termination must lie at or to the left of `z = 0` and
//! a right termination at or to the right of `z = L`; a termination exactly on
//! the outermost interface of its side is allowed. The exterior material fills
//! the gap between the outermost interface and the wall.
//!
//! # Wall conditions
//!
//! A PEC wall imposes vanishing tangential electric field and a PMC wall
//! vanishing tangential magnetic field. In the scalar TE (`ψ = E_y`) and TM
//! (`ψ = H_y`) channels this is a Dirichlet or Neumann condition on `ψ`:
//!
//! | wall | TE                  | TM                  |
//! | ---- | ------------------- | ------------------- |
//! | PEC  | `E_y = 0` (`ψ = 0`) | `E_x = 0` (`∂ψ = 0`) |
//! | PMC  | `H_x = 0` (`∂ψ = 0`) | `H_y = 0` (`ψ = 0`) |

use lamina_units::Length;

/// Kind of a perfectly conducting wall.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum WallKind {
    /// Perfect electric conductor: vanishing tangential electric field.
    PerfectElectric,

    /// Perfect magnetic conductor: vanishing tangential magnetic field.
    PerfectMagnetic,
}

/// A perfect electric or perfect magnetic wall at an absolute position in the
/// stack frame.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Wall<F> {
    kind: WallKind,
    position: Length<F>,
}

impl<F> Wall<F> {
    /// Construct a wall of `kind` at absolute `position` in the stack frame.
    pub const fn new(kind: WallKind, position: Length<F>) -> Self {
        Self { kind, position }
    }

    /// Construct a perfect electric wall at `position`.
    pub const fn perfect_electric(position: Length<F>) -> Self {
        Self::new(WallKind::PerfectElectric, position)
    }

    /// Construct a perfect magnetic wall at `position`.
    pub const fn perfect_magnetic(position: Length<F>) -> Self {
        Self::new(WallKind::PerfectMagnetic, position)
    }

    /// Return the wall kind.
    pub const fn kind(&self) -> WallKind {
        self.kind
    }

    /// Return the absolute wall position in the stack frame.
    pub fn position(&self) -> Length<F>
    where
        F: Copy,
    {
        self.position
    }
}

/// How one side of a stack is terminated.
///
/// Further ways of closing a side may be added, so matches need a wildcard
/// arm.
#[derive(Copy, Clone, Debug, PartialEq, Default)]
#[non_exhaustive]
pub enum Termination<F> {
    /// The semi-infinite exterior medium.
    #[default]
    Open,

    /// A perfect electric or perfect magnetic wall.
    Wall(Wall<F>),
}

impl<F> Termination<F> {
    /// Return `true` for the semi-infinite exterior.
    pub const fn is_open(&self) -> bool {
        matches!(self, Self::Open)
    }

    /// Return `true` for a closed side.
    pub const fn is_closed(&self) -> bool {
        !self.is_open()
    }
}

impl<F> From<Wall<F>> for Termination<F> {
    fn from(wall: Wall<F>) -> Self {
        Self::Wall(wall)
    }
}

/// One side of a stack.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum StackSide {
    /// The side at and below `z = 0`.
    Left,

    /// The side at and above `z = L`.
    Right,
}
