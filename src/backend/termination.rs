//! Closed sides of a stack in the backends.
//!
//! A closed side replaces the outgoing radiation condition of an open
//! exterior by a wall condition, applied after propagation through the gap
//! between the outermost interface and the wall. Both backends read the same
//! closed-side data: Transfer2 carries the wall condition through the gap's
//! transfer matrix, entire in `κ²`; Scatter2 applies the wall as a load at the
//! outermost interface.
//!
//! The exterior material's wavevector on a closed side, used for the
//! directional decomposition of fields there, is the root `κ̃` of
//! `κ² = εμk₀² − β²` with `Im(κ̃ ℓ) ≥ 0`, `ℓ` the length from the outermost
//! interface to the wall. Every quantity reported by the backends is even in
//! that choice.

use nalgebra::ComplexField;
use ndarray::Dimension;
use num_traits::{One, Zero};
use thiserror::Error;

use crate::{
    ComplexScalar, Polarisation,
    algebra::ScalarAlgebra,
    backend::{ExteriorWavevectors, PlaneWaveModeCandidate, isotropic::IsotropicMediumQuantities},
    input::{
        CanonicalCoordinates, CanonicalStack,
        canonical::{CanonicalTermination, WallCondition},
    },
    material::{ConstitutiveEvaluator, ConstitutiveLift},
    stack::StackSide,
    waves::BidirectionalWaves,
};

/// Failure caused by the terminations of a stack.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum TerminationError {
    /// The supplied exterior wavevectors do not match the terminations: open
    /// sides need a wavevector, closed sides take none.
    #[error(
        "the {side:?} side is {}",
        if *open {
            "open, but no exterior wavevector was supplied"
        } else {
            "closed, but an exterior wavevector was supplied"
        }
    )]
    ExteriorBranchMismatch { side: StackSide, open: bool },

    /// A plane-wave excitation's internal fields were requested for a stack
    /// with a closed side.
    #[error("internal fields of a plane-wave excitation are not available for a closed stack")]
    DrivenFieldsOfClosedStack,
}

/// Check that the supplied exterior wavevectors match the terminations.
pub(crate) fn check_exterior_pattern<M, J, A>(
    stack: &CanonicalStack<M, J>,
    exterior: &ExteriorWavevectors<A>,
) -> Result<(), TerminationError> {
    for (side, termination, supplied) in [
        (StackSide::Left, stack.left_termination(), exterior.left()),
        (
            StackSide::Right,
            stack.right_termination(),
            exterior.right(),
        ),
    ] {
        let open = termination.is_open();

        if open != supplied.is_some() {
            return Err(TerminationError::ExteriorBranchMismatch { side, open });
        }
    }

    Ok(())
}

/// One closed side as seen by the backends.
#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub struct ClosedSide<A> {
    side: StackSide,
    condition: WallCondition,
    kappa: A,
    kappa_squared: A,
    factor: A,
    gap: A,
    reference_slope: A,
}

impl<A> ClosedSide<A> {
    pub(crate) fn try_map<B, E>(
        &self,
        mut f: impl FnMut(&A) -> Result<B, E>,
    ) -> Result<ClosedSide<B>, E> {
        Ok(ClosedSide {
            side: self.side,
            condition: self.condition,
            kappa: f(&self.kappa)?,
            kappa_squared: f(&self.kappa_squared)?,
            factor: f(&self.factor)?,
            gap: f(&self.gap)?,
            reference_slope: f(&self.reference_slope)?,
        })
    }
}

impl<A> ClosedSide<A>
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// Build the closed-side data for `termination`, or `None` for an open
    /// side.
    fn new(
        side: StackSide,
        termination: &CanonicalTermination<A>,
        medium: &IsotropicMediumQuantities<A>,
        coordinates: &CanonicalCoordinates<A>,
        polarisation: Polarisation,
    ) -> Option<Self> {
        let kind = termination.wall_kind()?;
        let gap = termination.gap_cm()?.clone();

        let factor = match polarisation {
            Polarisation::TransverseElectric => medium.mu().clone(),
            Polarisation::TransverseMagnetic => medium.epsilon().clone(),
        };

        Some(Self {
            side,
            condition: WallCondition::new(kind, polarisation),
            kappa: decaying_root(medium.kappa(), &gap),
            kappa_squared: medium.kappa_squared().clone(),
            factor,
            gap,
            reference_slope: coordinates
                .vacuum_angular_wavenumber()
                .scale(-A::Scalar::i()),
        })
    }

    /// Transfer matrix `[[cos, −f S], [q S / f, cos]]` of a homogeneous layer
    /// of the exterior material with thickness `d`, as `(cos, S)`, where
    /// `S = sin(κd)/κ`; entire in `q = κ²`.
    fn layer_entries(&self, thickness: &A) -> (A, A) {
        let phase_squared = self.kappa_squared.multiply(&thickness.square());

        (
            phase_squared.cos_sqrt(),
            thickness.multiply(&phase_squared.sinc_sqrt()),
        )
    }

    /// Apply a transfer matrix `[[c, −f s], [q s / f, c]]` to a state.
    fn apply(&self, (cos, sin): &(A, A), (field, slope): (A, A)) -> (A, A) {
        (
            cos.multiply(&field)
                .subtract(&self.factor.multiply(sin).multiply(&slope)),
            self.kappa_squared
                .divide(&self.factor)
                .multiply(sin)
                .multiply(&field)
                .add(&cos.multiply(&slope)),
        )
    }

    /// Apply the row vector `[a, b]` to a transfer matrix, returning the row.
    fn apply_row(&self, (cos, sin): &(A, A), (a, b): (A, A)) -> (A, A) {
        (
            a.multiply(cos)
                .add(&b.multiply(&self.kappa_squared.divide(&self.factor).multiply(sin))),
            b.multiply(cos)
                .subtract(&a.multiply(&self.factor.multiply(sin))),
        )
    }

    /// The transfer state at the outermost interface of a closed right side:
    /// `T(gap) r_wall`, with `r_wall = [0, −ξ_ref]` (Dirichlet) or `[1, 0]`
    /// (Neumann).
    pub(crate) fn right_state(&self) -> (A, A) {
        let one = self.gap.one_like();
        let zero = self.gap.zero_like();

        let state = match self.condition {
            WallCondition::Dirichlet => (zero, self.reference_slope.negate()),
            WallCondition::Neumann => (one, zero),
        };

        self.apply(&self.layer_entries(&self.gap), state)
    }

    /// The boundary row at the outermost interface of a closed left side:
    /// `ℓ_wall T(gap)`, with `ℓ_wall = [ξ_ref, 0]` (Dirichlet) or `[0, −1]`
    /// (Neumann).
    pub(crate) fn left_row(&self) -> (A, A) {
        let one = self.gap.one_like();
        let zero = self.gap.zero_like();

        let row = match self.condition {
            WallCondition::Dirichlet => (self.reference_slope.clone(), zero),
            WallCondition::Neumann => (zero, one.negate()),
        };

        self.apply_row(&self.layer_entries(&self.gap), row)
    }

    /// The wall's reflection in the exterior material's directional basis:
    /// `−1` (Dirichlet) or `+1` (Neumann).
    fn reflection(&self) -> A::Scalar {
        match self.condition {
            WallCondition::Dirichlet => -A::Scalar::one(),
            WallCondition::Neumann => A::Scalar::one(),
        }
    }

    /// The load seen from the outermost interface, `ρ exp(2iκ̃ℓ)`.
    pub(crate) fn load(&self) -> A {
        let two = A::Scalar::one() + A::Scalar::one();

        self.kappa
            .multiply(&self.gap)
            .scale(A::Scalar::i() * two)
            .exp()
            .scale(self.reflection())
    }

    /// The exponent of this side's scale factor in the scattering form of the
    /// characteristic: `−iκ̃ℓ`.
    pub(crate) fn scale_exponent(&self) -> A {
        self.kappa.multiply(&self.gap).scale(-A::Scalar::i())
    }

    /// The side coefficient `c` in the scattering form of the characteristic,
    /// in the exterior material's directional basis at the wall.
    pub(crate) fn coefficient(&self) -> A {
        let half = (A::Scalar::one() + A::Scalar::one()).recip();
        let slope = self.kappa.divide(&self.factor).scale(-A::Scalar::i());

        match (self.side, self.condition) {
            (StackSide::Left, WallCondition::Dirichlet) => self.reference_slope.clone(),
            (StackSide::Left, WallCondition::Neumann) => slope,
            (StackSide::Right, WallCondition::Dirichlet) => {
                self.reference_slope.divide(&slope).scale(half)
            }
            (StackSide::Right, WallCondition::Neumann) => self.gap.one_like().scale(half),
        }
    }

    /// The distance from the outermost interface of a point at outward
    /// distance `distance` (centimetres), or `None` beyond the wall.
    pub(crate) fn exterior_distance(
        &self,
        distance: <A::Scalar as ComplexField>::RealField,
    ) -> Option<<A::Scalar as ComplexField>::RealField> {
        let gap = self
            .gap
            .value()
            .iter()
            .next()
            .copied()
            .expect("a sampled value")
            .real();
        let tolerance = <A::Scalar as ComplexField>::RealField::one() + nalgebra::convert(1.0e-12);

        if distance > gap.clone() * tolerance {
            return None;
        }

        Some(if distance > gap { gap } else { distance })
    }
}

/// Return `κ` multiplied by `±1` per sample so that `Im(κ ℓ) ≥ 0`.
fn decaying_root<A>(kappa: &A, length: &A) -> A
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let product = kappa.multiply(length);

    if product
        .value()
        .iter()
        .all(|value| value.imaginary() >= <A::Scalar as ComplexField>::RealField::zero())
    {
        return kappa.clone();
    }

    let signs = product.value().mapv(|value| {
        if value.imaginary() < <A::Scalar as ComplexField>::RealField::zero() {
            -A::Scalar::one()
        } else {
            A::Scalar::one()
        }
    });

    kappa.multiply(&A::lift_constant(signs))
}

/// The exterior data shared by both backends' contexts: open sides take
/// their wavevectors from `exterior`, closed sides from their termination.
pub(crate) struct ExteriorParts<A> {
    pub(crate) left_admittance: A,
    pub(crate) right_admittance: A,
    pub(crate) left_kappa: A,
    pub(crate) right_kappa: A,
    pub(crate) left_epsilon: A,
    pub(crate) right_epsilon: A,
    pub(crate) left_mu: A,
    pub(crate) right_mu: A,
    pub(crate) left_closure: Option<ClosedSide<A>>,
    pub(crate) right_closure: Option<ClosedSide<A>>,
}

/// Evaluate the exterior data of `stack`, closing its terminated sides.
pub(crate) fn exterior_parts<E, M, A>(
    coordinates: &CanonicalCoordinates<A>,
    stack: &CanonicalStack<M, A>,
    exterior: &ExteriorWavevectors<A>,
    polarisation: Polarisation,
) -> Result<ExteriorParts<A>, TerminationError>
where
    A: ScalarAlgebra + ConstitutiveLift<E, M> + Clone,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
    E: ConstitutiveEvaluator<A::Scalar, A::Dimension, M>,
{
    check_exterior_pattern(stack, exterior)?;

    let left_medium =
        IsotropicMediumQuantities::evaluate::<E, M>(stack.left_exterior(), coordinates);
    let right_medium =
        IsotropicMediumQuantities::evaluate::<E, M>(stack.right_exterior(), coordinates);

    let closure = |side, termination, medium| {
        ClosedSide::new(side, termination, medium, coordinates, polarisation)
    };
    let left_closure = closure(StackSide::Left, stack.left_termination(), &left_medium);
    let right_closure = closure(StackSide::Right, stack.right_termination(), &right_medium);

    let kappa = |supplied: Option<&A>, closure: &Option<ClosedSide<A>>| match supplied {
        Some(kappa) => kappa.clone(),
        None => closure
            .as_ref()
            .expect("the pattern check pairs every closed side with a closure")
            .kappa
            .clone(),
    };

    let left_kappa = kappa(exterior.left(), &left_closure);
    let right_kappa = kappa(exterior.right(), &right_closure);

    let left_admittance = left_medium.admittance_with_kappa(&left_kappa, polarisation);
    let right_admittance = right_medium.admittance_with_kappa(&right_kappa, polarisation);

    let (left_epsilon, left_mu, _) = left_medium.into_parts();
    let (right_epsilon, right_mu, _) = right_medium.into_parts();

    Ok(ExteriorParts {
        left_admittance,
        right_admittance,
        left_kappa,
        right_kappa,
        left_epsilon,
        right_epsilon,
        left_mu,
        right_mu,
        left_closure,
        right_closure,
    })
}

/// The right exterior waves of a modal seed at the outermost interface: the
/// outgoing wave alone on an open side, the standing wave on a closed one.
pub(crate) fn right_exterior_mode_waves<A>(
    seed: &PlaneWaveModeCandidate<A>,
    right_admittance: &A,
) -> BidirectionalWaves<A>
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    match seed.right_closed_state() {
        Some(state) => {
            let slope = right_admittance.scale(-A::Scalar::i());
            let difference = state.secondary().divide(&slope);
            let half = (A::Scalar::one() + A::Scalar::one()).recip();

            BidirectionalWaves::new(
                state.field().subtract(&difference).scale(half),
                state.field().add(&difference).scale(half),
            )
        }
        None => BidirectionalWaves::new(
            seed.right_outgoing().clone(),
            seed.right_outgoing().zero_like(),
        ),
    }
}

/// Amplitudes for incidence from a closed side, which has no meaning; the
/// evaluators refuse it before projecting.
pub(crate) fn undefined_amplitudes<A>(
    source: &ndarray::ArrayBase<ndarray::OwnedRepr<A::Scalar>, A::Dimension>,
) -> crate::PlaneWaveAmplitudes<A>
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar,
{
    let nan = A::filled_constant_like(source, A::Scalar::from_real(nalgebra::convert(f64::NAN)));

    crate::PlaneWaveAmplitudes::new(nan.clone(), nan)
}
