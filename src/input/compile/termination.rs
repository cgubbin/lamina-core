//! Validation and canonical compilation of stack terminations.

use lamina_units::{Length, LengthUnit};
use nalgebra::ComplexField;
use ndarray::Array;
use num_complex::Complex;
use num_traits::{Float, FromPrimitive, Zero};
use thiserror::Error;

use crate::{
    ComplexScalar, FiniteLayerIndex, SeedJet,
    input::{
        StackCompileError, StackThicknessJet,
        canonical::{CanonicalLayer, CanonicalPml, CanonicalTermination},
    },
    stack::{Stack, StackSide, StretchProfile, Termination, WallKind},
};

/// The termination field named by a validation error.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TerminationField {
    /// The absolute position of a wall.
    WallPosition,
    /// The absolute position of a PML's inner face.
    PmlStart,
    /// The thickness of a PML.
    PmlThickness,
    /// The polynomial order of a PML profile.
    ProfileOrder,
    /// The complex strength of a PML profile.
    ProfileStrength,
}

impl std::fmt::Display for TerminationField {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::WallPosition => "wall position",
            Self::PmlStart => "PML start",
            Self::PmlThickness => "PML thickness",
            Self::ProfileOrder => "profile order",
            Self::ProfileStrength => "profile strength",
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

    #[error("{side:?} termination: the PML thickness must be positive, got {thickness:?}")]
    NonPositiveThickness {
        side: StackSide,
        thickness: Length<R>,
    },

    #[error("{side:?} termination: the profile order must be non-negative, got {order:?}")]
    NegativeOrder { side: StackSide, order: R },

    #[error(
        "{side:?} termination: the profile strength must have a non-negative real part \
         and a positive imaginary part, got {strength:?}"
    )]
    InvalidStrength {
        side: StackSide,
        strength: Complex<R>,
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
    Pml {
        start_cm: R,
        thickness_cm: R,
        gap_cm: R,
        order: R,
        strength: Complex<R>,
        backing: WallKind,
    },
}

impl<R: Float> ValidatedTermination<R> {
    /// The distance from the outermost interface to the closing wall, in
    /// centimetres, or `None` for an open side.
    pub(crate) fn exterior_extent_cm(&self) -> Option<R> {
        match self {
            Self::Open => None,
            Self::Wall { gap_cm, .. } => Some(*gap_cm),
            Self::Pml {
                gap_cm,
                thickness_cm,
                ..
            } => Some(*gap_cm + *thickness_cm),
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
            gap_cm: gap(
                wall.position(),
                side,
                TerminationField::WallPosition,
                thicknesses,
            )?,
        }),

        Termination::Pml(pml) => {
            let gap_cm = gap(pml.start(), side, TerminationField::PmlStart, thicknesses)?;
            let thickness = pml.thickness();
            let StretchProfile::Polynomial { order, strength } = *pml.profile();
            let non_finite = |field| Err(TerminationValidationError::NonFinite { side, field });

            if !thickness.value().is_finite() {
                return non_finite(TerminationField::PmlThickness);
            }

            if thickness.value() <= R::zero() {
                return Err(TerminationValidationError::NonPositiveThickness { side, thickness });
            }

            if !order.is_finite() {
                return non_finite(TerminationField::ProfileOrder);
            }

            if order < R::zero() {
                return Err(TerminationValidationError::NegativeOrder { side, order });
            }

            if !strength.re.is_finite() || !strength.im.is_finite() {
                return non_finite(TerminationField::ProfileStrength);
            }

            if strength.re < R::zero() || strength.im <= R::zero() {
                return Err(TerminationValidationError::InvalidStrength { side, strength });
            }

            Ok(ValidatedTermination::Pml {
                start_cm: positive_zero(pml.start().as_centimetres()),
                thickness_cm: thickness.as_centimetres(),
                gap_cm,
                order: positive_zero(order),
                strength: Complex::new(positive_zero(strength.re), strength.im),
                backing: pml.backing(),
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
/// thicknesses of `layers`: the right wall or PML stays fixed in the stack
/// frame, so `gap_R = z_R − Σ d_k`. The value of each gap is the validated one,
/// so a termination exactly on the outermost interface has an exactly zero
/// gap. There no derivative that moves the last interface outward exists, and
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

    if let (
        Some(layer),
        ValidatedTermination::Wall { gap_cm, .. } | ValidatedTermination::Pml { gap_cm, .. },
    ) = (outward_layer, &right)
    {
        if gap_cm.is_zero() {
            return Err(StackCompileError::OutwardThicknessDerivative {
                layer: FiniteLayerIndex::new(layer),
                side: StackSide::Right,
            });
        }
    }

    let complex = |value| <J as SeedJet>::constant(Array::from_elem(sampled_shape.clone(), value));
    let constant = |value| complex(J::Scalar::from_real(value));

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

        ValidatedTermination::Pml {
            start_cm,
            thickness_cm,
            gap_cm,
            order,
            strength,
            backing,
        } => {
            let one = <<J::Scalar as ComplexField>::RealField as num_traits::One>::one();
            let strength = J::Scalar::from_parts(strength.re, strength.im);
            let stretch =
                J::Scalar::from_real(one) + strength * J::Scalar::from_real(one / (order + one));

            CanonicalTermination::Pml(CanonicalPml {
                start_cm: constant(start_cm),
                thickness_cm: constant(thickness_cm),
                gap_cm: gap(gap_cm),
                optical_thickness_cm: complex(J::Scalar::from_real(thickness_cm) * stretch),
                order: constant(order),
                strength: complex(strength),
                backing,
            })
        }
    };

    Ok((canonical(left, &constant), canonical(right, &right_gap)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stack::{Pml, Wall};

    fn stack() -> Stack<&'static str, f64> {
        Stack::builder("left", "right")
            .layer("a", Length::micrometres(0.1))
            .layer("b", Length::micrometres(0.2))
            .finalise()
    }

    fn pml(start: Length<f64>, thickness: Length<f64>, order: f64, c: Complex<f64>) -> Pml<f64> {
        Pml::new(
            start,
            thickness,
            StretchProfile::polynomial(order, c),
            WallKind::PerfectElectric,
        )
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
            .with_right_termination(pml(
                Length::micrometres(2.3),
                Length::micrometres(1.0),
                2.0,
                Complex::new(0.5, 3.0),
            ));
        let (left, right) = validate_terminations(&stack).unwrap();
        let ValidatedTermination::Wall { gap_cm, .. } = left else {
            panic!("left wall expected");
        };
        assert!((gap_cm - 1.0e-4).abs() < 1e-18);
        let ValidatedTermination::Pml {
            gap_cm,
            thickness_cm,
            ..
        } = right
        else {
            panic!("right PML expected");
        };
        assert!((gap_cm - 2.0e-4).abs() < 1e-17);
        assert!((thickness_cm - 1.0e-4).abs() < 1e-18);
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

    #[test]
    fn invalid_pml_fields_are_rejected_by_side_and_field() {
        let start = Length::micrometres(0.5);
        let cases = [
            (
                pml(start, Length::micrometres(0.0), 2.0, Complex::new(0.0, 1.0)),
                "thickness",
            ),
            (
                pml(start, Length::micrometres(1.0), 2.0, Complex::new(1.0, 0.0)),
                "strength",
            ),
            (
                pml(
                    start,
                    Length::micrometres(1.0),
                    2.0,
                    Complex::new(-0.1, 1.0),
                ),
                "strength",
            ),
            (
                pml(
                    start,
                    Length::micrometres(1.0),
                    -1.0,
                    Complex::new(0.0, 1.0),
                ),
                "order",
            ),
            (
                pml(
                    start,
                    Length::micrometres(f64::NAN),
                    2.0,
                    Complex::new(0.0, 1.0),
                ),
                "nonfinite",
            ),
            (
                pml(
                    Length::micrometres(f64::INFINITY),
                    Length::micrometres(1.0),
                    2.0,
                    Complex::new(0.0, 1.0),
                ),
                "nonfinite start",
            ),
        ];
        for (pml, case) in cases {
            let stack = stack().with_right_termination(pml);
            let error = validate_terminations(&stack).unwrap_err();
            let ok = match (case, &error) {
                ("thickness", TerminationValidationError::NonPositiveThickness { side, .. }) => {
                    *side == StackSide::Right
                }
                ("strength", TerminationValidationError::InvalidStrength { side, .. }) => {
                    *side == StackSide::Right
                }
                ("order", TerminationValidationError::NegativeOrder { side, .. }) => {
                    *side == StackSide::Right
                }
                ("nonfinite", TerminationValidationError::NonFinite { side, field }) => {
                    *side == StackSide::Right && *field == TerminationField::PmlThickness
                }
                ("nonfinite start", TerminationValidationError::NonFinite { side, field }) => {
                    *side == StackSide::Right && *field == TerminationField::PmlStart
                }
                _ => false,
            };
            assert!(ok, "{case}: {error:?}");
        }
    }
}
