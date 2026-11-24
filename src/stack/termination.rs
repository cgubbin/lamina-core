//! Terminations closing a planar stack on either side.
//!
//! By default a stack is bounded by two semi-infinite exterior media. Either
//! side may instead be closed:
//!
//! - by a perfect electric (PEC) or perfect magnetic (PMC) [`Wall`] at an
//!   absolute position `z` in the stack's own frame;
//! - by a perfectly matched layer ([`Pml`]) by complex coordinate stretching,
//!   starting at an absolute position, extending outward by a stated
//!   thickness and backed by a PEC or PMC wall.
//!
//! The stack frame places `z = 0` at the first interface, with the layers
//! running towards `+z` and ending at `z = L`, the sum of the finite-layer
//! thicknesses. A left termination must lie at or to the left of `z = 0` and
//! a right termination at or to the right of `z = L`; a termination exactly on
//! the outermost interface of its side is allowed. The exterior material fills
//! the gap between the outermost interface and the wall or the PML start, and
//! also fills the PML itself.
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
//!
//! # Stretching
//!
//! With `ζ ∈ [0, d]` the outward distance into a PML of thickness `d`, the
//! frequency-independent polynomial stretch is
//!
//! ```text
//! s(ζ) = 1 + c (ζ/d)^m,    ζ̃(ζ) = ζ + c d/(m+1) (ζ/d)^(m+1)
//! ```
//!
//! with real order `m ≥ 0` and complex strength `c`, `Re c ≥ 0`, `Im c > 0`.
//! Fields inside the PML are the analytic continuation of the exterior
//! solution to the stretched coordinate `ζ̃`, and for the characteristic the
//! PML is a homogeneous layer of the exterior material with complex optical
//! thickness `τ = d (1 + c/(m+1))`.

use lamina_units::Length;
use num_complex::Complex;

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

/// Complex stretching profile of a PML.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum StretchProfile<F> {
    /// `s(ζ) = 1 + strength (ζ/d)^order`.
    Polynomial {
        /// Real, non-negative polynomial order `m`.
        order: F,
        /// Complex strength `c` with `Re c ≥ 0` and `Im c > 0`.
        strength: Complex<F>,
    },
}

impl<F> StretchProfile<F> {
    /// Construct a polynomial profile of the given order and strength.
    pub const fn polynomial(order: F, strength: Complex<F>) -> Self {
        Self::Polynomial { order, strength }
    }
}

/// A wall-backed perfectly matched layer by complex coordinate stretching.
///
/// The layer starts at absolute position `start` in the stack frame and
/// extends outward, away from the layers, by `thickness`. Its backing wall
/// therefore lies at `start − thickness` on the left and `start + thickness`
/// on the right.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Pml<F> {
    start: Length<F>,
    thickness: Length<F>,
    profile: StretchProfile<F>,
    backing: WallKind,
}

impl<F> Pml<F> {
    /// Construct a PML. Every field is required.
    pub const fn new(
        start: Length<F>,
        thickness: Length<F>,
        profile: StretchProfile<F>,
        backing: WallKind,
    ) -> Self {
        Self {
            start,
            thickness,
            profile,
            backing,
        }
    }

    /// Return the absolute position of the PML's inner face.
    pub fn start(&self) -> Length<F>
    where
        F: Copy,
    {
        self.start
    }

    /// Return the PML thickness, measured outward from its start.
    pub fn thickness(&self) -> Length<F>
    where
        F: Copy,
    {
        self.thickness
    }

    /// Return the stretching profile.
    pub fn profile(&self) -> &StretchProfile<F> {
        &self.profile
    }

    /// Return the kind of the backing wall.
    pub const fn backing(&self) -> WallKind {
        self.backing
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

    /// A wall-backed perfectly matched layer.
    Pml(Pml<F>),
}

impl<F> Termination<F> {
    /// Return `true` for the semi-infinite exterior.
    pub const fn is_open(&self) -> bool {
        matches!(self, Self::Open)
    }
}

impl<F> From<Wall<F>> for Termination<F> {
    fn from(wall: Wall<F>) -> Self {
        Self::Wall(wall)
    }
}

impl<F> From<Pml<F>> for Termination<F> {
    fn from(pml: Pml<F>) -> Self {
        Self::Pml(pml)
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
