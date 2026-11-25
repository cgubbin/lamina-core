//! The modal null vector of the 4×4 interface system.
//!
//! A mode is a nonzero solution of the homogeneous interface system: for
//! Transfer4 `[V_L | M V_R] (a, b) = 0`, for Scatter4 the last pivot of its
//! elimination, `[U | −V_R] (b, t) = 0`, back-substituted through the earlier
//! pivots. Either backend retains the square system `N` and, at every
//! interface, the `n × m` map from the first `m` components of a null vector
//! to the state there.
//!
//! The null vector is the right singular vector `v` of the smallest singular
//! value of `N`'s values, with `u` its left singular vector, carried through
//! the jets by the bordered system
//!
//! ```text
//! [[N, u], [v†, 0]] (x, s) = (0, 1),
//! ```
//!
//! nonsingular while the null space is one-dimensional: `x` is analytic and
//! `N x = −s u`, `s` vanishing at the mode. The residual reported is
//! `σ_min / σ_max`; a second singular value below [`DEGENERATE_RATIO`]
//! `σ_max` is a two-dimensional null space and is refused.

use ndarray::{Array, Dimension};

use crate::{ComplexScalar, algebra::ScalarAlgebra, backend::ModeReconstructionError};

use super::matrix::{Matrix, sample_matrices};

/// Relative size of the second-smallest singular value at or below which the
/// null space is two-dimensional.
pub(crate) const DEGENERATE_RATIO: f64 = 1.0e-8;

/// The homogeneous interface system of a retained 4×4 solve.
#[derive(Clone, Debug)]
pub(crate) struct Modal<A> {
    /// The square system `N`.
    pub(crate) system: Matrix<A>,
    /// At every interface, the map from the first `m` components of a null
    /// vector to the state.
    pub(crate) basis: Vec<Matrix<A>>,
}

/// An unnormalised mode of a retained 4×4 workspace: the state at every
/// interface, with the residual `σ_min / σ_max` of the homogeneous system.
#[doc(hidden)]
#[derive(Clone, Debug)]
pub struct CoupledModeCandidate<A> {
    pub(crate) states: Vec<Matrix<A>>,
    residual: f64,
}

impl<A> CoupledModeCandidate<A> {
    /// `σ_min / σ_max` of the homogeneous interface system over every
    /// sample: zero at a mode.
    pub fn residual(&self) -> f64 {
        self.residual
    }
}

impl<A> Modal<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// The mode of the system.
    pub(crate) fn candidate(&self) -> Result<CoupledModeCandidate<A>, ModeReconstructionError> {
        let (vector, residual) = null_vector(&self.system)?;
        let m = self.basis.first().map_or(0, Matrix::columns);
        let leading = vector.rows_in(0..m);
        Ok(CoupledModeCandidate {
            states: self
                .basis
                .iter()
                .map(|basis| basis.multiply(&leading))
                .collect(),
            residual,
        })
    }
}

/// The null vector of a square system and `σ_min / σ_max`.
fn null_vector<A>(system: &Matrix<A>) -> Result<(Matrix<A>, f64), ModeReconstructionError>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let n = system.rows();
    let like = system.get(0, 0);
    // left[k][sample], right[k][sample]: u and v† of each sample.
    let mut left = vec![Vec::new(); n];
    let mut right = vec![Vec::new(); n];
    let mut residual: f64 = 0.0;
    for sample in sample_matrices(system) {
        let svd =
            nalgebra::SVD::try_new(sample, true, true, nalgebra::convert(f64::EPSILON), 10_000)
                .ok_or(ModeReconstructionError::NoUsableNullVector)?;
        let values: Vec<f64> = svd
            .singular_values
            .iter()
            .map(|value| nalgebra::try_convert::<_, f64>(value.clone()).unwrap_or(f64::NAN))
            .collect();
        let largest = values.iter().copied().fold(0.0, f64::max);
        if !(largest.is_finite() && largest > 0.0) || values.iter().any(|v| !v.is_finite()) {
            return Err(ModeReconstructionError::NoUsableNullVector);
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
        let (smallest, second) = (order[0], order.get(1).copied());
        if second.is_some_and(|second| values[second] <= DEGENERATE_RATIO * largest) {
            return Err(ModeReconstructionError::DegenerateNullSpace);
        }
        residual = residual.max(values[smallest] / largest);
        let (Some(u), Some(v_t)) = (svd.u, svd.v_t) else {
            return Err(ModeReconstructionError::NoUsableNullVector);
        };
        for k in 0..n {
            left[k].push(u[(k, smallest)]);
            right[k].push(v_t[(smallest, k)]);
        }
    }
    let lift = |values: Vec<A::Scalar>| {
        Array::from_shape_vec(like.value().raw_dim(), values)
            .map(A::lift_constant)
            .map_err(|_| ModeReconstructionError::NoUsableNullVector)
    };
    let left = left.into_iter().map(lift).collect::<Result<Vec<_>, _>>()?;
    let right = right.into_iter().map(lift).collect::<Result<Vec<_>, _>>()?;

    let bordered = Matrix::from_fn(n + 1, n + 1, |row, column| match (row < n, column < n) {
        (true, true) => system.get(row, column).clone(),
        (true, false) => left[row].clone(),
        (false, true) => right[column].clone(),
        (false, false) => like.zero_like(),
    });
    let unit = Matrix::from_fn(n + 1, 1, |row, _| {
        if row == n {
            like.one_like()
        } else {
            like.zero_like()
        }
    });
    let solution = bordered.solve(&unit);
    if !solution.all_finite() {
        return Err(ModeReconstructionError::NoUsableNullVector);
    }
    Ok((solution.rows_in(0..n), residual))
}
