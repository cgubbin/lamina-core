//! Small dense matrices of jets.
//!
//! The 4×4 backends work with matrices of at most 8×8 jet-valued entries,
//! so determinants and inverses are formed by cofactor expansion: no pivot
//! has to be chosen from one sample and applied to every sample and
//! derivative component.

use nalgebra::ComplexField;
use ndarray::Dimension;

use crate::{ComplexScalar, algebra::ScalarAlgebra};

/// A dense row-major matrix of jets.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Matrix<A> {
    rows: usize,
    columns: usize,
    entries: Vec<A>,
}

impl<A> Matrix<A> {
    pub(crate) fn from_fn(rows: usize, columns: usize, f: impl FnMut(usize, usize) -> A) -> Self {
        let mut f = f;
        let entries = (0..rows * columns)
            .map(|index| f(index / columns, index % columns))
            .collect();
        Self {
            rows,
            columns,
            entries,
        }
    }

    pub(crate) const fn rows(&self) -> usize {
        self.rows
    }

    pub(crate) const fn columns(&self) -> usize {
        self.columns
    }

    pub(crate) fn get(&self, row: usize, column: usize) -> &A {
        &self.entries[row * self.columns + column]
    }

    pub(crate) fn entries(&self) -> &[A] {
        &self.entries
    }

    pub(crate) fn map<B>(&self, f: impl FnMut(&A) -> B) -> Matrix<B> {
        Matrix {
            rows: self.rows,
            columns: self.columns,
            entries: self.entries.iter().map(f).collect(),
        }
    }
}

impl<A: Clone> Matrix<A> {
    /// The submatrix of the given rows and columns.
    pub(crate) fn select(&self, rows: &[usize], columns: &[usize]) -> Self {
        Self::from_fn(rows.len(), columns.len(), |row, column| {
            self.get(rows[row], columns[column]).clone()
        })
    }

    /// The given contiguous columns.
    pub(crate) fn columns_in(&self, columns: std::ops::Range<usize>) -> Self {
        let rows: Vec<usize> = (0..self.rows).collect();
        self.select(&rows, &columns.collect::<Vec<_>>())
    }

    /// The given contiguous rows.
    pub(crate) fn rows_in(&self, rows: std::ops::Range<usize>) -> Self {
        let columns: Vec<usize> = (0..self.columns).collect();
        self.select(&rows.collect::<Vec<_>>(), &columns)
    }

    /// `selfᵀ`.
    pub(crate) fn transpose(&self) -> Self {
        Self::from_fn(self.columns, self.rows, |row, column| {
            self.get(column, row).clone()
        })
    }

    /// `[self | rhs]`.
    pub(crate) fn hstack(&self, rhs: &Self) -> Self {
        Self::from_fn(self.rows, self.columns + rhs.columns, |row, column| {
            if column < self.columns {
                self.get(row, column).clone()
            } else {
                rhs.get(row, column - self.columns).clone()
            }
        })
    }

    /// `[[self], [rhs]]`.
    pub(crate) fn vstack(&self, rhs: &Self) -> Self {
        Self::from_fn(self.rows + rhs.rows, self.columns, |row, column| {
            if row < self.rows {
                self.get(row, column).clone()
            } else {
                rhs.get(row - self.rows, column).clone()
            }
        })
    }
}

impl<A> Matrix<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// The `rows × columns` zero matrix with the sampled shape of `like`.
    pub(crate) fn zeros(rows: usize, columns: usize, like: &A) -> Self {
        Self::from_fn(rows, columns, |_, _| like.zero_like())
    }

    /// The `n × n` identity with the sampled shape of `like`.
    pub(crate) fn identity(n: usize, like: &A) -> Self {
        Self::from_fn(n, n, |row, column| {
            if row == column {
                like.one_like()
            } else {
                like.zero_like()
            }
        })
    }

    /// The unit columns `e_k` of the `n × n` identity for each `k` in
    /// `columns`.
    pub(crate) fn unit_columns(n: usize, columns: &[usize], like: &A) -> Self {
        Self::from_fn(n, columns.len(), |row, column| {
            if row == columns[column] {
                like.one_like()
            } else {
                like.zero_like()
            }
        })
    }

    pub(crate) fn add(&self, rhs: &Self) -> Self {
        self.zip(rhs, A::add)
    }

    pub(crate) fn subtract(&self, rhs: &Self) -> Self {
        self.zip(rhs, A::subtract)
    }

    fn zip(&self, rhs: &Self, f: impl Fn(&A, &A) -> A) -> Self {
        Self {
            rows: self.rows,
            columns: self.columns,
            entries: self
                .entries
                .iter()
                .zip(&rhs.entries)
                .map(|(a, b)| f(a, b))
                .collect(),
        }
    }

    pub(crate) fn scale(&self, coefficient: A::Scalar) -> Self {
        self.map(|entry| entry.scale(coefficient))
    }

    pub(crate) fn scale_by(&self, factor: &A) -> Self {
        self.map(|entry| entry.multiply(factor))
    }

    /// `self + factor I`.
    pub(crate) fn shift(&self, factor: &A) -> Self {
        Self::from_fn(self.rows, self.columns, |row, column| {
            if row == column {
                self.get(row, column).add(factor)
            } else {
                self.get(row, column).clone()
            }
        })
    }

    pub(crate) fn multiply(&self, rhs: &Self) -> Self {
        Self::from_fn(self.rows, rhs.columns, |row, column| {
            (1..self.columns).fold(
                self.get(row, 0).multiply(rhs.get(0, column)),
                |sum, inner| sum.add(&self.get(row, inner).multiply(rhs.get(inner, column))),
            )
        })
    }

    /// The trace of a square matrix.
    pub(crate) fn trace(&self) -> A {
        (1..self.rows).fold(self.get(0, 0).clone(), |sum, index| {
            sum.add(self.get(index, index))
        })
    }

    /// The determinant of a square matrix by cofactor expansion along the
    /// first row.
    pub(crate) fn determinant(&self) -> A {
        match self.rows {
            1 => self.get(0, 0).clone(),
            2 => self
                .get(0, 0)
                .multiply(self.get(1, 1))
                .subtract(&self.get(0, 1).multiply(self.get(1, 0))),
            n => (0..n)
                .map(|column| self.get(0, column).multiply(&self.cofactor(0, column)))
                .reduce(|sum, term| sum.add(&term))
                .unwrap_or_else(|| self.get(0, 0).one_like()),
        }
    }

    /// The signed minor `(−1)^{row+column} det(self without row, column)`.
    fn cofactor(&self, row: usize, column: usize) -> A {
        let keep = |skip: usize| (0..self.rows).filter(move |&index| index != skip);
        let rows: Vec<usize> = keep(row).collect();
        let columns: Vec<usize> = keep(column).collect();
        let minor = self.select(&rows, &columns).determinant();

        if (row + column) % 2 == 0 {
            minor
        } else {
            minor.negate()
        }
    }

    /// The adjugate of a square matrix: `adj(A) A = det(A) I`.
    pub(crate) fn adjugate(&self) -> Self {
        if self.rows == 1 {
            return Self::identity(1, self.get(0, 0));
        }
        Self::from_fn(self.rows, self.columns, |row, column| {
            self.cofactor(column, row)
        })
    }

    /// `self⁻¹ rhs`, by the adjugate.
    pub(crate) fn solve(&self, rhs: &Self) -> Self {
        let reciprocal = self.determinant().reciprocal();
        self.adjugate().multiply(rhs).scale_by(&reciprocal)
    }

    /// The largest absolute row sum of the values over every sample; NaN
    /// when a value is not finite.
    pub(crate) fn value_norm(&self) -> f64 {
        sample_matrices(self)
            .iter()
            .map(|matrix| {
                matrix
                    .row_iter()
                    .map(|row| {
                        row.iter()
                            .map(|value| {
                                nalgebra::try_convert::<_, f64>(value.abs()).unwrap_or(f64::NAN)
                            })
                            .sum::<f64>()
                    })
                    .fold(0.0, max_nan)
            })
            .fold(0.0, max_nan)
    }

    pub(crate) fn all_finite(&self) -> bool {
        self.entries.iter().all(A::all_finite)
    }
}

impl<A: crate::algebra::RealScalarAlgebra> Matrix<A> {
    /// `self†`, the conjugate transpose.
    pub(crate) fn adjoint(&self) -> Self {
        self.transpose().map(A::conjugated)
    }
}

/// `max` that propagates NaN.
fn max_nan(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

/// The values of a matrix of jets, one matrix per sample in the logical order
/// of the sampled arrays.
pub(crate) fn sample_matrices<A>(matrix: &Matrix<A>) -> Vec<nalgebra::DMatrix<A::Scalar>>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
{
    let values: Vec<Vec<A::Scalar>> = matrix
        .entries()
        .iter()
        .map(|entry| entry.value().iter().copied().collect())
        .collect();
    let samples = values.first().map_or(0, Vec::len);

    (0..samples)
        .map(|sample| {
            nalgebra::DMatrix::from_fn(matrix.rows(), matrix.columns(), |row, column| {
                values[row * matrix.columns() + column][sample]
            })
        })
        .collect()
}
