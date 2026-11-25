//! The matrix exponential of a jet-valued matrix by scaling and squaring.
//!
//! `exp(A) = (exp(A / 2^s))^(2^s)` with `exp(A / 2^s)` from its Taylor
//! series. The exponent `s` is chosen once from the values of every sample, so
//! the same arithmetic is applied to every derivative component: the result is
//! entire in the entries of `A`, with no eigenvalues, and so holds at
//! degenerate and defective matrices alike.

use ndarray::Dimension;

use crate::{ComplexScalar, algebra::ScalarAlgebra};

use super::matrix::Matrix;

/// Taylor degree of the scaled kernel; with `‖A / 2^s‖ ≤ 1/2` the truncation
/// error is below `2⁻¹⁷ / 17! ≈ 2·10⁻²⁰`.
const TAYLOR_DEGREE: u32 = 16;

/// Largest squaring count; beyond it every representable result has
/// overflowed.
const MAX_SQUARINGS: u32 = 1100;

/// `exp(A)` of a square matrix.
pub(crate) fn exponential<A>(matrix: &Matrix<A>) -> Matrix<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let norm = matrix.value_norm();
    let squarings = if norm.is_finite() && norm > 0.5 {
        ((norm / 0.5).log2().ceil() as u32).min(MAX_SQUARINGS)
    } else {
        0
    };
    let like = matrix.get(0, 0);
    let half: A::Scalar = nalgebra::convert(0.5);
    let scaled = (0..squarings).fold(matrix.clone(), |scaled, _| scaled.scale(half));

    // Horner form of Σ_k X^k / k!: I + X (I + X/2 (I + X/3 (…))).
    let identity = Matrix::identity(matrix.rows(), like);
    let mut result = identity.clone();
    for k in (1..=TAYLOR_DEGREE).rev() {
        let reciprocal: A::Scalar = nalgebra::convert(1.0 / f64::from(k));
        result = identity.add(&scaled.multiply(&result).scale(reciprocal));
    }

    (0..squarings).fold(result, |result, _| result.multiply(&result))
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use ndarray::{Ix0, arr0};
    use num_complex::Complex64 as C;

    use super::*;
    use crate::algebra::{ArrayJet0, ArrayJet1, RealParameter, SeedJet};

    type J0 = ArrayJet0<C, Ix0, RealParameter>;
    type J1 = ArrayJet1<C, Ix0, RealParameter>;

    fn constant(rows: &[&[C]]) -> Matrix<J0> {
        Matrix::from_fn(rows.len(), rows[0].len(), |row, column| {
            J0::lift_constant(arr0(rows[row][column]))
        })
    }

    fn value(matrix: &Matrix<J0>, row: usize, column: usize) -> C {
        matrix.get(row, column).value()[()]
    }

    #[test]
    fn diagonalisable_matrix_matches_its_eigen_decomposition() {
        // A = V diag(λ) V⁻¹ with V = [[1, 1], [1, −1]] / √2 symmetric.
        let (a, b) = (C::new(0.3, 1.7), C::new(-2.1, 0.4));
        let half = C::new(0.5, 0.0);
        let matrix = constant(&[
            &[(a + b) * half, (a - b) * half],
            &[(a - b) * half, (a + b) * half],
        ]);
        let result = exponential(&matrix.scale(C::new(0.0, 3.0)));
        let (ea, eb) = ((a * C::new(0.0, 3.0)).exp(), (b * C::new(0.0, 3.0)).exp());

        for (row, column, expected) in [
            (0, 0, (ea + eb) * half),
            (0, 1, (ea - eb) * half),
            (1, 1, (ea + eb) * half),
        ] {
            let actual = value(&result, row, column);
            assert_relative_eq!(
                actual.re,
                expected.re,
                epsilon = 1e-12,
                max_relative = 1e-12
            );
            assert_relative_eq!(
                actual.im,
                expected.im,
                epsilon = 1e-12,
                max_relative = 1e-12
            );
        }
    }

    #[test]
    fn defective_matrix_gives_the_jordan_exponential() {
        // exp([[λ, 1], [0, λ]] t) = e^{λt} [[1, t], [0, 1]], with t large.
        let lambda = C::new(-0.2, 1.1);
        let t = C::new(40.0, 0.0);
        let matrix = constant(&[&[lambda * t, t], &[C::new(0.0, 0.0), lambda * t]]);
        let result = exponential(&matrix);
        let e = (lambda * t).exp();

        for (row, column, expected) in [(0, 0, e), (0, 1, e * t), (1, 0, C::new(0.0, 0.0))] {
            let actual = value(&result, row, column);
            assert!((actual - expected).norm() <= 1e-12 * (1.0 + expected.norm()));
        }
    }

    #[test]
    fn derivatives_are_continuous_through_a_degenerate_matrix() {
        // exp(t A(p)) with A(p) = [[0, 1], [p, 0]]: eigenvalues ±√p meet at p = 0.
        // d/dp exp(tA) at p has the closed form of the entire series in p.
        let t = 2.5;
        let derivative_at = |p: f64| {
            let parameter = <J1 as SeedJet>::variable(arr0(C::new(p, 0.0)), 0).unwrap();
            let one = parameter.one_like();
            let zero = parameter.zero_like();
            let matrix = Matrix::from_fn(2, 2, |row, column| match (row, column) {
                (0, 1) => one.clone(),
                (1, 0) => parameter.clone(),
                _ => zero.clone(),
            })
            .scale(C::new(t, 0.0));
            exponential(&matrix).get(0, 0).clone()
        };
        // cosh(t√p) = Σ (t²p)^k / (2k)!; d/dp at 0 = t²/2.
        let at_zero = derivative_at(0.0);
        let expected = t * t / 2.0;
        assert!((at_zero.first()[()] - C::new(expected, 0.0)).norm() < 1e-12);

        for p in [1e-6, -1e-6] {
            let nearby = derivative_at(p);
            assert!((nearby.first()[()] - at_zero.first()[()]).norm() < 1e-4);
        }
    }
}
