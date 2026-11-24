//! Validation and canonical compilation of stack terminations.

use std::fmt::Debug;

use lamina_units::{Length, LengthUnit};
use nalgebra::ComplexField;
use ndarray::{Array, ArrayBase, Dimension, OwnedRepr};
use num_traits::{Float, FromPrimitive, Zero};
use thiserror::Error;

use crate::{
    ComplexScalar, SeedJet,
    input::{
        StackThicknessJet,
        canonical::{CanonicalTermination, CanonicalWall},
    },
    stack::{Stack, StackSide, Termination, WallKind},
};

/// The termination field named by a validation error.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TerminationField {
    /// The absolute position of a wall.
    WallPosition,
}

impl std::fmt::Display for TerminationField {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::WallPosition => "wall position",
        })
    }
}

/// Invalid termination of one side of a stack.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum TerminationValidationError<R> {
    #[error("{side:?} termination: the {field} is not finite")]
    NonFinite {
        side: StackSide,
        field: TerminationField,
    },

    #[error(
        "{side:?} termination: the {field} {position:?} lies inside the layers \
         (the outermost interface of that side is at {interface:?})"
    )]
    InsideLayers {
        side: StackSide,
        field: TerminationField,
        position: Length<R>,
        interface: Length<R>,
    },
}

/// The validated left and right terminations.
pub(crate) type ValidatedTerminations<R> = (ValidatedTermination<R>, ValidatedTermination<R>);

/// A validated termination with its lengths in centimetres.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ValidatedTermination<R> {
    Open,
    Wall {
        kind: WallKind,
        position_cm: R,
        gap_cm: R,
    },
}

impl<R: Float> ValidatedTermination<R> {
    /// The distance from the outermost interface to the closing wall, in
    /// centimetres, or `None` for an open side.
    pub(crate) fn exterior_extent_cm(&self) -> Option<R> {
        match self {
            Self::Open => None,
            Self::Wall { gap_cm, .. } => Some(*gap_cm),
        }
    }
}

/// Validate both terminations of `stack`.
///
/// The last interface `L` is the left-to-right sum of the layer thicknesses.
/// When the layers and the position share one length unit, the sum and the
/// comparison are carried out in that unit, so that a position equal to the
/// sum of the same values is recognised exactly as lying on the interface;
/// otherwise both are compared in centimetres.
pub(crate) fn validate_terminations<M, R>(
    stack: &Stack<M, R>,
) -> Result<ValidatedTerminations<R>, TerminationValidationError<R>>
where
    R: Float + FromPrimitive,
{
    let thicknesses = stack
        .layers_left_to_right()
        .iter()
        .map(|layer| layer.thickness())
        .collect::<Vec<_>>();

    let left = validate_side(stack.left_termination(), StackSide::Left, &thicknesses)?;
    let right = validate_side(stack.right_termination(), StackSide::Right, &thicknesses)?;

    Ok((left, right))
}

fn validate_side<R>(
    termination: &Termination<R>,
    side: StackSide,
    thicknesses: &[Length<R>],
) -> Result<ValidatedTermination<R>, TerminationValidationError<R>>
where
    R: Float + FromPrimitive,
{
    match termination {
        Termination::Open => Ok(ValidatedTermination::Open),

        Termination::Wall(wall) => {
            let gap_cm = gap(
                wall.position(),
                side,
                TerminationField::WallPosition,
                thicknesses,
            )?;

            Ok(ValidatedTermination::Wall {
                kind: wall.kind(),
                position_cm: positive_zero(wall.position().as_centimetres()),
                gap_cm,
            })
        }
    }
}

/// Return the gap between the outermost interface of `side` and `position`,
/// in centimetres, rejecting a position inside the layers.
fn gap<R>(
    position: Length<R>,
    side: StackSide,
    field: TerminationField,
    thicknesses: &[Length<R>],
) -> Result<R, TerminationValidationError<R>>
where
    R: Float + FromPrimitive,
{
    if !position.value().is_finite() {
        return Err(TerminationValidationError::NonFinite { side, field });
    }

    let (value, unit) = position.into_parts();

    match side {
        StackSide::Left => {
            if value > R::zero() {
                return Err(TerminationValidationError::InsideLayers {
                    side,
                    field,
                    position,
                    interface: Length::new(R::zero(), unit),
                });
            }

            Ok(positive_zero(-value * unit.scale_to_centimetres::<R>()))
        }

        StackSide::Right => {
            let shared = thicknesses.iter().all(|thickness| thickness.unit() == unit);

            let (interface, gap_cm) = if shared {
                let total = thicknesses
                    .iter()
                    .fold(R::zero(), |total, thickness| total + thickness.value());

                (
                    Length::new(total, unit),
                    (value - total) * unit.scale_to_centimetres::<R>(),
                )
            } else {
                let total = thicknesses.iter().fold(R::zero(), |total, thickness| {
                    total + thickness.as_centimetres()
                });

                (
                    Length::new(total, LengthUnit::Centimetre),
                    position.as_centimetres() - total,
                )
            };

            if gap_cm < R::zero() {
                return Err(TerminationValidationError::InsideLayers {
                    side,
                    field,
                    position,
                    interface,
                });
            }

            Ok(positive_zero(gap_cm))
        }
    }
}

fn positive_zero<R: Float>(value: R) -> R {
    if value == R::zero() { R::zero() } else { value }
}

/// Build canonical terminations from validated ones.
///
/// `thicknesses_cm` are the canonical finite-layer thickness jets, from which
/// the right gap inherits its derivatives: the right wall stays fixed in
/// the stack frame, so `gap_R = z_R − Σ d_k`. The value of each gap is the
/// validated one, so a termination exactly on the outermost interface has an
/// exactly zero gap.
pub(crate) fn canonical_terminations<J>(
    validated: ValidatedTerminations<<J::Scalar as ComplexField>::RealField>,
    thicknesses_cm: &[&J],
    source: &ArrayBase<OwnedRepr<J::Scalar>, J::Dimension>,
) -> (CanonicalTermination<J>, CanonicalTermination<J>)
where
    J: StackThicknessJet,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
    <J::Scalar as ComplexField>::RealField: Float + Copy,
{
    let (left, right) = validated;

    let constant =
        |value: J::Scalar| <J as SeedJet>::constant(Array::from_elem(source.raw_dim(), value));
    let real =
        |value: <J::Scalar as ComplexField>::RealField| constant(J::Scalar::from_real(value));

    let left_gap = |gap: <J::Scalar as ComplexField>::RealField| real(gap);

    let right_gap = |gap: <J::Scalar as ComplexField>::RealField| {
        // The value is the validated gap; the derivatives are those of −Σ d_k.
        let total = thicknesses_cm
            .iter()
            .fold(constant(J::Scalar::zero()), |total, d| total.add_jet(d));

        let fluctuation = total.subtract_jet(&total.value_constant());

        real(gap).subtract_jet(&fluctuation)
    };

    (
        canonical_side(left, &constant, left_gap),
        canonical_side(right, &constant, right_gap),
    )
}

fn canonical_side<J>(
    validated: ValidatedTermination<<J::Scalar as ComplexField>::RealField>,
    constant: &impl Fn(J::Scalar) -> J,
    gap: impl Fn(<J::Scalar as ComplexField>::RealField) -> J,
) -> CanonicalTermination<J>
where
    J: StackThicknessJet,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
    <J::Scalar as ComplexField>::RealField: Float + Copy,
{
    let real = |value| constant(J::Scalar::from_real(value));

    match validated {
        ValidatedTermination::Open => CanonicalTermination::Open,

        ValidatedTermination::Wall {
            kind,
            position_cm,
            gap_cm,
        } => CanonicalTermination::Wall(CanonicalWall::new(kind, real(position_cm), gap(gap_cm))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stack::Wall;

    fn stack() -> Stack<&'static str, f64> {
        Stack::builder("left", "right")
            .layer("a", Length::micrometres(0.1))
            .layer("b", Length::micrometres(0.2))
            .finalise()
    }

    #[test]
    fn open_sides_validate_as_open() {
        let (left, right) = validate_terminations(&stack()).unwrap();
        assert_eq!(left, ValidatedTermination::Open);
        assert_eq!(right, ValidatedTermination::Open);
    }

    #[test]
    fn right_wall_inside_the_layers_is_rejected() {
        let stack =
            stack().with_right_termination(Wall::perfect_electric(Length::micrometres(0.25)));
        let error = validate_terminations(&stack).unwrap_err();
        assert!(matches!(
            error,
            TerminationValidationError::InsideLayers {
                side: StackSide::Right,
                field: TerminationField::WallPosition,
                ..
            }
        ));
    }

    #[test]
    fn left_wall_inside_the_layers_is_rejected() {
        let stack = stack().with_left_termination(Wall::perfect_magnetic(Length::nanometres(1.0)));
        let error = validate_terminations(&stack).unwrap_err();
        assert!(matches!(
            error,
            TerminationValidationError::InsideLayers {
                side: StackSide::Left,
                field: TerminationField::WallPosition,
                ..
            }
        ));
    }

    #[test]
    fn walls_on_the_outermost_interfaces_have_zero_gaps() {
        // 0.1 + 0.2 is not 0.3 in binary; the wall at the left-to-right sum is
        // recognised exactly.
        let total = 0.1 + 0.2;
        let stack = stack()
            .with_left_termination(Wall::perfect_electric(Length::micrometres(-0.0)))
            .with_right_termination(Wall::perfect_magnetic(Length::micrometres(total)));
        let (left, right) = validate_terminations(&stack).unwrap();
        let ValidatedTermination::Wall { gap_cm, .. } = left else {
            panic!("left wall expected");
        };
        assert_eq!(gap_cm.to_bits(), 0.0_f64.to_bits());
        let ValidatedTermination::Wall { gap_cm, .. } = right else {
            panic!("right wall expected");
        };
        assert_eq!(gap_cm.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn gaps_are_measured_from_the_outermost_interfaces() {
        let stack = stack()
            .with_left_termination(Wall::perfect_electric(Length::micrometres(-1.0)))
            .with_right_termination(Wall::perfect_magnetic(Length::micrometres(2.3)));
        let (left, right) = validate_terminations(&stack).unwrap();
        let ValidatedTermination::Wall { gap_cm, .. } = left else {
            panic!("left wall expected");
        };
        assert!((gap_cm - 1.0e-4).abs() < 1e-18);
        let ValidatedTermination::Wall { gap_cm, kind, .. } = right else {
            panic!("right wall expected");
        };
        assert_eq!(kind, WallKind::PerfectMagnetic);
        assert!((gap_cm - 2.0e-4).abs() < 1e-17);
    }

    #[test]
    fn non_finite_wall_positions_are_rejected() {
        for position in [f64::NAN, f64::INFINITY] {
            let stack = stack()
                .with_left_termination(Wall::perfect_electric(Length::micrometres(position)));
            assert_eq!(
                validate_terminations(&stack).unwrap_err(),
                TerminationValidationError::NonFinite {
                    side: StackSide::Left,
                    field: TerminationField::WallPosition,
                }
            );
        }
    }
}
