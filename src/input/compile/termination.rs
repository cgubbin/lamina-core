//! Validation and canonical compilation of stack terminations.

use lamina_units::{Length, LengthUnit};
use nalgebra::ComplexField;
use ndarray::Array;
use num_traits::{Float, FromPrimitive, Zero};
use thiserror::Error;

use crate::{
    ComplexScalar, FiniteLayerIndex, SeedJet,
    input::{
        StackCompileError, StackThicknessJet,
        canonical::{CanonicalLayer, CanonicalTermination},
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
#[allow(clippy::type_complexity)]
pub(crate) fn validate_terminations<M, R>(
    stack: &Stack<M, R>,
) -> Result<(ValidatedTermination<R>, ValidatedTermination<R>), TerminationValidationError<R>>
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

        Termination::Wall(wall) => Ok(ValidatedTermination::Wall {
            kind: wall.kind(),
            position_cm: positive_zero(wall.position().as_centimetres()),
            gap_cm: gap(wall.position(), side, thicknesses)?,
        }),
    }
}

/// Return the gap between the outermost interface of `side` and `position`,
/// in centimetres, rejecting a position inside the layers.
fn gap<R>(
    position: Length<R>,
    side: StackSide,
    thicknesses: &[Length<R>],
) -> Result<R, TerminationValidationError<R>>
where
    R: Float + FromPrimitive,
{
    let field = TerminationField::WallPosition;

    if !position.value().is_finite() {
        return Err(TerminationValidationError::NonFinite { side, field });
    }

    let (value, unit) = position.into_parts();

    let (interface, gap_cm) = match side {
        StackSide::Left => (
            Length::new(R::zero(), unit),
            -value * unit.scale_to_centimetres::<R>(),
        ),

        StackSide::Right if thicknesses.iter().all(|thickness| thickness.unit() == unit) => {
            let total = thicknesses
                .iter()
                .fold(R::zero(), |total, thickness| total + thickness.value());

            (
                Length::new(total, unit),
                (value - total) * unit.scale_to_centimetres::<R>(),
            )
        }

        StackSide::Right => {
            let total = thicknesses.iter().fold(R::zero(), |total, thickness| {
                total + thickness.as_centimetres()
            });

            (
                Length::new(total, LengthUnit::Centimetre),
                position.as_centimetres() - total,
            )
        }
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

fn positive_zero<R: Float>(value: R) -> R {
    if value == R::zero() { R::zero() } else { value }
}

/// Validate the terminations of `stack` and build their canonical form.
///
/// The right gap inherits the derivatives of the canonical finite-layer
/// thicknesses of `layers`: the right wall stays fixed in the stack frame, so
/// `gap_R = z_R − Σ d_k`. The value of each gap is the validated one, so a
/// wall exactly on the outermost interface has an exactly zero gap. There no
/// derivative that moves the last interface outward exists, and
/// `outward_layer`, the first layer whose seeded thickness does so, is
/// refused.
#[allow(clippy::type_complexity)]
pub(crate) fn compile_terminations<L, J>(
    stack: &Stack<impl Sized, <J::Scalar as ComplexField>::RealField>,
    layers: &[CanonicalLayer<L, J>],
    sampled_shape: J::Dimension,
    outward_layer: Option<usize>,
) -> Result<
    (CanonicalTermination<J>, CanonicalTermination<J>),
    StackCompileError<<J::Scalar as ComplexField>::RealField>,
>
where
    J: StackThicknessJet,
    J::Scalar: ComplexScalar,
    <J::Scalar as ComplexField>::RealField: Float + FromPrimitive,
    J::Dimension: Clone,
{
    let (left, right) = validate_terminations(stack)?;

    if let (Some(layer), Some(gap)) = (outward_layer, right.exterior_extent_cm()) {
        if gap.is_zero() {
            return Err(StackCompileError::OutwardThicknessDerivative {
                layer: FiniteLayerIndex::new(layer),
                side: StackSide::Right,
            });
        }
    }

    let constant = |value| {
        <J as SeedJet>::constant(Array::from_elem(
            sampled_shape.clone(),
            J::Scalar::from_real(value),
        ))
    };

    let right_gap = |gap| {
        layers.iter().fold(constant(gap), |gap: J, layer| {
            gap.subtract_derivatives(layer.thickness_cm())
        })
    };

    let canonical = |validated, gap: &dyn Fn(_) -> J| match validated {
        ValidatedTermination::Open => CanonicalTermination::Open,

        ValidatedTermination::Wall {
            kind,
            position_cm,
            gap_cm,
        } => CanonicalTermination::Wall {
            kind,
            position_cm: constant(position_cm),
            gap_cm: gap(gap_cm),
        },
    };

    Ok((canonical(left, &constant), canonical(right, &right_gap)))
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
