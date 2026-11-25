//! Orientation of a material frame in the stack frame.

use num_traits::{Float, FloatConst};
use thiserror::Error;

use super::Tensor3;
use crate::ComplexScalar;

/// A proper rotation `R` from a material's own (crystal) frame to the stack
/// frame.
///
/// Column `k` of `R` holds the stack-frame components of crystal axis `k`
/// (a, b, c), so a tensor transforms as `T_stack = R T_crystal Rᵀ`. The stack
/// frame has `z` along the stack normal, `x` along the in-plane wavevector
/// and `y` invariant.
///
/// Every constructor resolves to the matrix, which is the recorded form of
/// the orientation. A matrix is accepted when `‖RᵀR − I‖_max ≤ 1e-12` and
/// `|det R − 1| ≤ 1e-12`; its entries are kept as given.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rotation<R> {
    matrix: [[R; 3]; 3],
}

/// A signed stack axis, the image of one crystal axis under a signed
/// permutation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SignedAxis {
    PlusX,
    MinusX,
    PlusY,
    MinusY,
    PlusZ,
    MinusZ,
}

/// A matrix that is not a proper rotation.
#[derive(Clone, Copy, Debug, PartialEq, Error)]
pub enum RotationError<R> {
    /// `‖RᵀR − I‖_max` exceeds `1e-12`, or an entry is not finite.
    #[error("the matrix is not orthonormal: max |RᵀR − I| = {deviation:?}")]
    NotOrthonormal { deviation: R },

    /// The matrix is orthonormal but its determinant is not `+1`.
    #[error("the matrix is not a proper rotation: det R = {determinant:?}")]
    Improper { determinant: R },
}

impl<R> Rotation<R> {
    /// The recorded crystal-to-stack matrix.
    pub fn matrix(&self) -> &[[R; 3]; 3] {
        &self.matrix
    }

    /// The stack-frame image `R T Rᵀ` of a crystal-frame tensor `T`, such as
    /// a crystal-frame rate of a tensor material.
    pub fn rotate<C>(&self, tensor: &Tensor3<C>) -> Tensor3<C>
    where
        C: ComplexScalar<RealField = R>,
        R: Copy,
    {
        let r = &self.matrix;
        [0, 1, 2].map(|row| {
            [0, 1, 2].map(|column| {
                let mut entry = C::zero();
                for (k, tensor_row) in tensor.iter().enumerate() {
                    for (l, &value) in tensor_row.iter().enumerate() {
                        entry += value * C::from_real(r[row][k]) * C::from_real(r[column][l]);
                    }
                }
                entry
            })
        })
    }
}

impl<R: Float> Rotation<R> {
    /// The identity: crystal axes a, b, c along stack axes x, y, z.
    pub fn identity() -> Self {
        Self::from_columns(SignedAxis::PlusX, SignedAxis::PlusY, SignedAxis::PlusZ)
    }

    /// Validate and record a crystal-to-stack rotation matrix.
    pub fn from_matrix(matrix: [[R; 3]; 3]) -> Result<Self, RotationError<R>> {
        let tolerance = R::from(1e-12).unwrap_or_else(R::epsilon);

        let mut deviation = R::zero();
        for i in 0..3 {
            for j in 0..3 {
                let product = (0..3).fold(R::zero(), |sum, k| sum + matrix[k][i] * matrix[k][j]);
                let identity = if i == j { R::one() } else { R::zero() };
                deviation = deviation.max((product - identity).abs());
            }
        }
        if deviation > tolerance || matrix.iter().flatten().any(|entry| !entry.is_finite()) {
            return Err(RotationError::NotOrthonormal { deviation });
        }

        let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
        let determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
        if (determinant - R::one()).abs() > tolerance {
            return Err(RotationError::Improper { determinant });
        }

        Ok(Self { matrix })
    }

    /// A signed permutation of the stack axes: crystal axes a, b, c map onto
    /// `a`, `b`, `c`. Of the 48 signed permutations the 24 with determinant
    /// `+1` are accepted.
    pub fn signed_permutation(
        a: SignedAxis,
        b: SignedAxis,
        c: SignedAxis,
    ) -> Result<Self, RotationError<R>> {
        Self::from_matrix(Self::from_columns(a, b, c).matrix)
    }

    /// `true` when every entry is exactly `0` or `±1`.
    pub(crate) fn is_signed_permutation(&self) -> bool {
        self.matrix
            .iter()
            .flatten()
            .all(|&entry| entry == R::zero() || entry.abs() == R::one())
    }

    /// `true` when one crystal axis lies exactly on the stack `y` axis and the
    /// other two exactly in the `xz` plane, so that a tensor diagonal in the
    /// crystal frame has no `xy`, `yx`, `yz` or `zy` component.
    pub(crate) fn has_crystal_axis_on_y(&self) -> bool {
        let m = &self.matrix;
        (0..3).any(|k| {
            m[0][k] == R::zero()
                && m[2][k] == R::zero()
                && (0..3).all(|j| j == k || m[1][j] == R::zero())
        })
    }

    fn from_columns(a: SignedAxis, b: SignedAxis, c: SignedAxis) -> Self {
        let mut matrix = [[R::zero(); 3]; 3];
        for (column, axis) in [a, b, c].into_iter().enumerate() {
            let (row, sign) = match axis {
                SignedAxis::PlusX => (0, R::one()),
                SignedAxis::MinusX => (0, -R::one()),
                SignedAxis::PlusY => (1, R::one()),
                SignedAxis::MinusY => (1, -R::one()),
                SignedAxis::PlusZ => (2, R::one()),
                SignedAxis::MinusZ => (2, -R::one()),
            };
            matrix[row][column] = sign;
        }
        Self { matrix }
    }
}

impl<R: Float + FloatConst> Rotation<R> {
    /// Intrinsic z-y′-z″ Euler angles: `R = R_z(α) R_y(β) R_z(γ)`, in radians.
    ///
    /// Like [`from_optic_axis`](Self::from_optic_axis), angles within an ulp
    /// of a multiple of `π/2` take exact sines and cosines, so quarter turns
    /// give exact signed permutations.
    pub fn from_euler(alpha: R, beta: R, gamma: R) -> Self {
        let (sa, ca) = sin_cos(alpha);
        let (sb, cb) = sin_cos(beta);
        let (sg, cg) = sin_cos(gamma);

        Self {
            matrix: [
                [ca * cb * cg - sa * sg, -ca * cb * sg - sa * cg, ca * sb],
                [sa * cb * cg + ca * sg, -sa * cb * sg + ca * cg, sa * sb],
                [-sb * cg, sb * sg, cb],
            ],
        }
    }

    /// Map the crystal c axis onto the direction
    /// `â = (sin θ cos φ, sin θ sin φ, cos θ)`: polar angle `θ` from the stack
    /// normal `z`, azimuth `φ` from the propagation axis `x` towards `y`, in
    /// radians.
    ///
    /// The rotation is the minimal one taking `ẑ` to `â`: a turn by `θ` about
    /// `û = (−sin φ, cos φ, 0)`. This completes the frame, which matters only
    /// for biaxial models. With `φ = 0` it is a rotation about `y` and leaves
    /// `y` exactly invariant; an axis along `x`, `y` or `z` gives a signed
    /// permutation.
    pub fn from_optic_axis(theta: R, phi: R) -> Self {
        let (s, c) = sin_cos(theta);
        let h = R::one() - c;
        let (sp, cp) = sin_cos(phi);
        let (ux, uy) = (-sp, cp);

        // Rodrigues: R = I + sin θ [u]× + (1 − cos θ)(u uᵀ − I), u_z = 0.
        Self {
            matrix: [
                [R::one() + h * (ux * ux - R::one()), h * ux * uy, s * uy],
                [h * uy * ux, R::one() + h * (uy * uy - R::one()), -s * ux],
                [-s * uy, s * ux, c],
            ],
        }
    }
}

/// `(sin, cos)`, exact (`0`, `±1`) within an ulp of a multiple of `π/2`.
fn sin_cos<R: Float + FloatConst>(angle: R) -> (R, R) {
    let quarters = (angle / R::FRAC_PI_2()).round();
    if (angle - quarters * R::FRAC_PI_2()).abs() > R::epsilon() * angle.abs() {
        return angle.sin_cos();
    }
    let two = R::one() + R::one();
    let four = two + two;
    let quadrant = ((quarters % four) + four) % four;
    let (zero, one) = (R::zero(), R::one());
    match quadrant.to_u8() {
        Some(0) => (zero, one),
        Some(1) => (one, zero),
        Some(2) => (zero, -one),
        Some(3) => (-one, zero),
        _ => angle.sin_cos(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const AXES: [SignedAxis; 6] = [
        SignedAxis::PlusX,
        SignedAxis::MinusX,
        SignedAxis::PlusY,
        SignedAxis::MinusY,
        SignedAxis::PlusZ,
        SignedAxis::MinusZ,
    ];

    fn assert_proper(rotation: &Rotation<f64>) {
        let validated = Rotation::from_matrix(*rotation.matrix()).unwrap();
        assert_eq!(&validated, rotation);
    }

    fn column(rotation: &Rotation<f64>, k: usize) -> [f64; 3] {
        let m = rotation.matrix();
        [m[0][k], m[1][k], m[2][k]]
    }

    fn assert_close(actual: [f64; 3], expected: [f64; 3]) {
        for (a, e) in actual.iter().zip(expected) {
            assert!((a - e).abs() <= 1e-15, "{actual:?} != {expected:?}");
        }
    }

    #[test]
    fn exactly_twenty_four_signed_permutations_are_proper() {
        let mut accepted = 0;
        for a in AXES {
            for b in AXES {
                for c in AXES {
                    if let Ok(rotation) = Rotation::<f64>::signed_permutation(a, b, c) {
                        assert!(rotation.is_signed_permutation());
                        accepted += 1;
                    }
                }
            }
        }
        assert_eq!(accepted, 24);
    }

    #[test]
    fn improper_and_non_orthonormal_matrices_are_rejected() {
        let reflection = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]];
        assert_eq!(
            Rotation::from_matrix(reflection),
            Err(RotationError::Improper { determinant: -1.0 })
        );

        let repeated = Rotation::<f64>::signed_permutation(
            SignedAxis::PlusX,
            SignedAxis::PlusX,
            SignedAxis::PlusZ,
        );
        assert!(matches!(
            repeated,
            Err(RotationError::NotOrthonormal { .. })
        ));

        let mut nearly = *Rotation::<f64>::identity().matrix();
        nearly[0][1] = 2e-12;
        assert!(matches!(
            Rotation::from_matrix(nearly),
            Err(RotationError::NotOrthonormal { .. })
        ));
        nearly[0][1] = 5e-13;
        assert!(Rotation::from_matrix(nearly).is_ok());

        let mut non_finite = *Rotation::<f64>::identity().matrix();
        non_finite[1][0] = f64::NAN;
        assert!(matches!(
            Rotation::from_matrix(non_finite),
            Err(RotationError::NotOrthonormal { .. })
        ));
    }

    #[test]
    fn euler_angles_compose_intrinsic_z_y_z() {
        let (alpha, beta, gamma) = (0.3, 1.1, -0.7);
        let rotation = Rotation::from_euler(alpha, beta, gamma);
        assert_proper(&rotation);

        // The c axis is turned by β from z, then azimuthally by α.
        assert_close(
            column(&rotation, 2),
            [
                alpha.cos() * beta.sin(),
                alpha.sin() * beta.sin(),
                beta.cos(),
            ],
        );
        // Pure β is a rotation about y.
        let about_y = Rotation::from_euler(0.0, beta, 0.0);
        assert!(about_y.has_crystal_axis_on_y());
        assert_eq!(column(&about_y, 1), [0.0, 1.0, 0.0]);
    }

    #[test]
    fn optic_axis_constructor_maps_c_onto_the_requested_direction() {
        for (theta, phi) in [(0.0, 0.0), (0.4, 0.0), (1.2, 0.9), (2.5, -2.0), (3.0, 1.0)] {
            let rotation = Rotation::from_optic_axis(theta, phi);
            assert_proper(&rotation);
            let (s, c): (f64, f64) = (f64::sin(theta), f64::cos(theta));
            assert_close(column(&rotation, 2), [s * phi.cos(), s * phi.sin(), c]);
        }

        assert_eq!(
            Rotation::from_optic_axis(0.0, 0.7),
            Rotation::<f64>::identity()
        );
        let tilted = Rotation::from_optic_axis(0.6, 0.0);
        assert!(tilted.has_crystal_axis_on_y());
        assert!(!tilted.is_signed_permutation());
        assert!(!Rotation::from_optic_axis(0.6, 0.2).has_crystal_axis_on_y());

        // θ = π stays a proper rotation: a half turn about y.
        assert_proper(&Rotation::from_optic_axis(std::f64::consts::PI, 0.0));
    }

    #[test]
    fn quarter_turns_are_exact_signed_permutations() {
        use std::f64::consts::{FRAC_PI_2, PI};

        // Optic axis along x, y and z, and back along −z.
        for (theta, phi, c_axis) in [
            (FRAC_PI_2, 0.0, [1.0, 0.0, 0.0]),
            (FRAC_PI_2, FRAC_PI_2, [0.0, 1.0, 0.0]),
            (0.0, 3.0 * FRAC_PI_2, [0.0, 0.0, 1.0]),
            (PI, 0.0, [0.0, 0.0, -1.0]),
        ] {
            let rotation = Rotation::from_optic_axis(theta, phi);
            assert!(rotation.is_signed_permutation(), "{theta}, {phi}");
            assert_eq!(column(&rotation, 2), c_axis);
            assert_proper(&rotation);
        }

        let euler = Rotation::from_euler(-FRAC_PI_2, FRAC_PI_2, 2.0 * PI);
        assert!(euler.is_signed_permutation());
        assert_proper(&euler);

        // Away from a quarter turn the angles are used as given.
        assert!(!Rotation::from_optic_axis(FRAC_PI_2 + 1e-9, 0.0).is_signed_permutation());
    }

    #[test]
    fn a_crystal_axis_on_y_need_not_be_b() {
        // c along y, then a quarter turn of the xz plane: still separable.
        let c_on_y = Rotation::<f64>::signed_permutation(
            SignedAxis::PlusZ,
            SignedAxis::PlusX,
            SignedAxis::PlusY,
        )
        .unwrap();
        assert!(c_on_y.has_crystal_axis_on_y());
        assert!(
            !Rotation::from_euler(0.2, 0.0, 0.0).has_crystal_axis_on_y(),
            "a turn about z moves a and b off y"
        );
    }
}
