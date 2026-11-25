//! Closed-form integrals of quadratic forms of the tangential state.
//!
//! In a homogeneous region `ψ(t) = exp(tA) ψ₀`, with `A = iK` along `+z` and
//! `A = −iK` outward on the left. A Hermitian form `ψ_a† Q ψ_b` or a bilinear
//! form `ψ_aᵀ Q ψ_b` integrates over a finite layer of thickness `d` to
//! `ψ_a,0^{†,ᵀ} I ψ_b,0` with
//!
//! ```text
//! I = ∫₀^d e^{Xt} Q e^{Yt} dt = e^{Xd} G,   exp([[−X, Q], [0, Y]] d) = [[e^{−Xd}, G], [0, e^{Yd}]],
//! ```
//!
//! `X = A_a†` (Hermitian) or `X = A_aᵀ` (bilinear) and `Y = A_b` (Van Loan).
//! Both exponentials are the entire scaling-and-squaring kernel, so `I` is
//! entire in the generators and the thickness: exact at degenerate and
//! defective layer eigenvalues, including a layer cutoff. Over a half-line on
//! which every wave decays, `I = ∫₀^∞ e^{Xt} Q e^{Yt} dt` is the solution of
//! the Sylvester equation `X I + I Y = −Q`.

use ndarray::Dimension;

use crate::{ComplexScalar, algebra::ScalarAlgebra};

use super::{exponential::exponential, matrix::Matrix};

/// `∫₀^d e^{Xt} Q e^{Yt} dt`.
pub(crate) fn layer<A>(x: &Matrix<A>, q: &Matrix<A>, y: &Matrix<A>, thickness: &A) -> Matrix<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let (rows, columns) = (x.rows(), y.rows());
    let minus_one = -<A::Scalar as num_traits::One>::one();
    let block = x
        .scale(minus_one)
        .hstack(q)
        .vstack(&Matrix::zeros(columns, rows, thickness).hstack(y))
        .scale_by(thickness);
    let exponential_block = exponential(&block);
    let all: Vec<usize> = (0..rows).collect();
    let right: Vec<usize> = (rows..rows + columns).collect();
    exponential(&x.scale_by(thickness)).multiply(&exponential_block.select(&all, &right))
}

/// `∫₀^∞ e^{Xt} Q e^{Yt} dt` of decaying `e^{Xt}`, `e^{Yt}`: the solution of
/// `X I + I Y = −Q`, by its Kronecker form.
pub(crate) fn half_line<A>(x: &Matrix<A>, q: &Matrix<A>, y: &Matrix<A>) -> Matrix<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let (rows, columns) = (x.rows(), y.rows());
    let like = q.get(0, 0);
    // vec(X I + I Y) = (1 ⊗ X + Yᵀ ⊗ 1) vec(I), column-major.
    let system = Matrix::from_fn(rows * columns, rows * columns, |row, column| {
        let (i, j) = (row % rows, row / rows);
        let (k, l) = (column % rows, column / rows);
        let mut entry = like.zero_like();
        if j == l {
            entry = entry.add(x.get(i, k));
        }
        if i == k {
            entry = entry.add(y.get(l, j));
        }
        entry
    });
    let rhs = Matrix::from_fn(rows * columns, 1, |row, _| {
        q.get(row % rows, row / rows).negate()
    });
    let solution = system.solve(&rhs);
    Matrix::from_fn(rows, columns, |row, column| {
        solution.get(column * rows + row, 0).clone()
    })
}

#[cfg(test)]
mod tests {
    use ndarray::{Ix0, arr0};
    use num_complex::Complex64 as C;

    use super::*;
    use crate::algebra::{ArrayJet0, ArrayJet1, RealParameter, SeedJet};

    type J0 = ArrayJet0<C, Ix0, RealParameter>;
    type J1 = ArrayJet1<C, Ix0, RealParameter>;

    fn matrix(entries: &[&[C]]) -> Matrix<J0> {
        Matrix::from_fn(entries.len(), entries[0].len(), |row, column| {
            J0::lift_constant(arr0(entries[row][column]))
        })
    }

    fn values(matrix: &Matrix<J0>) -> nalgebra::DMatrix<C> {
        nalgebra::DMatrix::from_fn(matrix.rows(), matrix.columns(), |row, column| {
            matrix.get(row, column).value()[()]
        })
    }

    fn expm(matrix: &nalgebra::DMatrix<C>, t: f64) -> nalgebra::DMatrix<C> {
        values(&exponential(&Matrix::from_fn(
            matrix.nrows(),
            matrix.ncols(),
            |row, column| J0::lift_constant(arr0(matrix[(row, column)] * t)),
        )))
    }

    /// Composite Gauss–Legendre quadrature of `e^{Xt} Q e^{Yt}` on `[0, d]`.
    fn quadrature(
        x: &nalgebra::DMatrix<C>,
        q: &nalgebra::DMatrix<C>,
        y: &nalgebra::DMatrix<C>,
        d: f64,
    ) -> nalgebra::DMatrix<C> {
        let nodes = [
            (-0.906_179_845_938_664, 0.236_926_885_056_189),
            (-0.538_469_310_105_683, 0.478_628_670_499_366),
            (0.0, 0.568_888_888_888_889),
            (0.538_469_310_105_683, 0.478_628_670_499_366),
            (0.906_179_845_938_664, 0.236_926_885_056_189),
        ];
        let panels = 400;
        let h = d / panels as f64;
        let mut sum = nalgebra::DMatrix::zeros(q.nrows(), q.ncols());
        for panel in 0..panels {
            for (node, weight) in nodes {
                let t = h * (panel as f64 + 0.5 + 0.5 * node);
                sum += (expm(x, t) * q * expm(y, t)) * C::from(0.5 * h * weight);
            }
        }
        sum
    }

    fn assert_close(
        actual: &nalgebra::DMatrix<C>,
        expected: &nalgebra::DMatrix<C>,
        tolerance: f64,
    ) {
        let error = (actual - expected).norm();
        assert!(
            error <= tolerance * expected.norm(),
            "{actual} vs {expected}: {error:e}"
        );
    }

    fn c(re: f64, im: f64) -> C {
        C::new(re, im)
    }

    #[test]
    fn layer_integral_matches_quadrature_and_eigen_wave_pairs() {
        let a = matrix(&[
            &[c(0.1, 0.3), c(1.2, -0.2), c(0.0, 0.4)],
            &[c(-0.7, 0.1), c(0.2, -0.5), c(0.3, 0.0)],
            &[c(0.2, 0.2), c(-0.1, 0.6), c(-0.3, 0.9)],
        ]);
        let q = matrix(&[
            &[c(1.0, 0.0), c(0.2, 0.1), c(0.0, -0.3)],
            &[c(0.2, -0.1), c(0.7, 0.0), c(0.4, 0.0)],
            &[c(0.0, 0.3), c(0.4, 0.0), c(1.5, 0.0)],
        ]);
        let b = matrix(&[
            &[c(0.0, 1.1), c(0.5, 0.0), c(0.1, 0.1)],
            &[c(0.3, 0.0), c(-0.2, -0.4), c(0.0, 0.2)],
            &[c(0.0, -0.6), c(0.2, 0.3), c(0.4, 0.5)],
        ]);
        let d = 1.7;
        let integral = layer(&a, &q, &b, &J0::lift_constant(arr0(c(d, 0.0))));
        let expected = quadrature(&values(&a), &values(&q), &values(&b), d);
        assert_close(&values(&integral), &expected, 1e-12);

        // Eigen-wave pairs: with X = V diag(λ) V⁻¹, Y = W diag(ν) W⁻¹,
        // I = V [(V⁻¹ Q W)_{jk} (e^{(λ_j + ν_k) d} − 1)/(λ_j + ν_k)] W⁻¹.
        let (x, y) = (values(&a), values(&b));
        let eigen = |m: &nalgebra::DMatrix<C>| {
            let roots = m.clone().schur().eigenvalues().unwrap();
            let vectors = nalgebra::DMatrix::from_fn(3, 3, |row, column| {
                let shifted = m - nalgebra::DMatrix::identity(3, 3) * (roots[column] + 1e-9);
                let v = shifted
                    .lu()
                    .solve(&nalgebra::DVector::from_element(3, c(1.0, 0.0)))
                    .unwrap();
                v[row] / v.norm()
            });
            (roots, vectors)
        };
        let ((lambda, v), (nu, w)) = (eigen(&x), eigen(&y));
        let core = v.clone().try_inverse().unwrap() * values(&q) * &w;
        let pairs = nalgebra::DMatrix::from_fn(3, 3, |j, k| {
            let s = lambda[j] + nu[k];
            core[(j, k)] * ((s * d).exp() - 1.0) / s
        });
        assert_close(&(v * pairs * w.try_inverse().unwrap()), &expected, 1e-7);
    }

    #[test]
    fn defective_generator_integrates_exactly() {
        // X = Y = [[λ, 1], [0, λ]]: e^{Xt} = e^{λt} [[1, t], [0, 1]].
        let lambda = c(0.0, 2.0);
        let x = matrix(&[&[lambda, c(1.0, 0.0)], &[c(0.0, 0.0), lambda]]);
        let q = matrix(&[&[c(1.0, 0.0), c(0.0, 0.0)], &[c(0.0, 0.0), c(0.0, 0.0)]]);
        let d = 0.9;
        let integral = values(&layer(&x, &q, &x, &J0::lift_constant(arr0(c(d, 0.0)))));
        // e^{Xt} Q e^{Xt} = e^{2λt} [[1, t], [0, 0]].
        let s = lambda * 2.0;
        let e = (s * d).exp();
        let first = (e - 1.0) / s;
        let second = (e * (s * d - 1.0) + 1.0) / (s * s);
        assert_close(
            &integral,
            &nalgebra::DMatrix::from_row_slice(2, 2, &[first, second, c(0.0, 0.0), c(0.0, 0.0)]),
            1e-13,
        );
    }

    #[test]
    fn thickness_derivative_is_the_integrand_at_the_far_end() {
        let x = [[c(0.2, 0.5), c(1.0, 0.0)], [c(-0.4, 0.0), c(0.0, -0.3)]];
        let q = [[c(1.0, 0.0), c(0.3, 0.2)], [c(0.3, -0.2), c(2.0, 0.0)]];
        let jet = |value: C| J1::lift_constant(arr0(value));
        let lift = |m: [[C; 2]; 2]| Matrix::from_fn(2, 2, |row, column| jet(m[row][column]));
        let d = 0.8;
        let thickness = <J1 as SeedJet>::variable(arr0(c(d, 0.0)), 0).unwrap();
        let integral = layer(&lift(x), &lift(q), &lift(x), &thickness);
        let x = nalgebra::DMatrix::from_fn(2, 2, |row, column| x[row][column]);
        let q = nalgebra::DMatrix::from_fn(2, 2, |row, column| q[row][column]);
        let expected = expm(&x, d) * q * expm(&x, d);
        let first =
            nalgebra::DMatrix::from_fn(2, 2, |row, column| integral.get(row, column).first()[()]);
        assert_close(&first, &expected, 1e-12);
    }

    #[test]
    fn half_line_integral_matches_a_long_layer() {
        let x = matrix(&[&[c(-0.8, 0.3), c(0.2, 0.0)], &[c(0.1, 0.1), c(-0.5, -1.0)]]);
        let y = matrix(&[&[c(-0.4, 2.0), c(0.0, 0.3)], &[c(0.5, 0.0), c(-1.1, 0.2)]]);
        let q = matrix(&[&[c(1.0, 0.2), c(0.3, 0.0)], &[c(0.0, -0.7), c(0.6, 0.1)]]);
        let tail = values(&half_line(&x, &q, &y));
        let truncated = quadrature(&values(&x), &values(&q), &values(&y), 40.0);
        assert_close(&tail, &truncated, 1e-12);
    }
}
