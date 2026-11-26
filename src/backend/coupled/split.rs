//! Spectral splitting of a Berreman generator.
//!
//! The `2m` eigenvalues of a generator `K` split into an *outgoing* (or
//! forward) set and its complement. Both sets enter only through monic
//! polynomials `f` (outgoing roots) and `g` (the others), with `f g` the
//! characteristic polynomial `p` of `K`. The spectral projector onto the
//! outgoing invariant subspace is a polynomial in `K`,
//!
//! ```text
//! P = a(K) g(K),   a g + b f = 1,   deg a, deg b < m,
//! ```
//!
//! from the Sylvester system of `(f, g)`, whose determinant is the resultant
//! of `f` and `g`. `P` is therefore analytic in the coefficients of `f` while
//! no outgoing root meets a remaining one, including where two outgoing roots
//! coincide (an isotropic region), and singular exactly where an outgoing and
//! a remaining root meet (a branch point or a layer cutoff).
//!
//! Callers that choose the split (real-axis exteriors, scattering layers)
//! classify the roots of each sample's values and refine the factors
//! through the jet algebra by Newton's method on `f g = p`, whose Jacobian is
//! the same Sylvester matrix.

use nalgebra::{DMatrix, DVector};
use ndarray::{Array, Dimension};

use crate::{ComplexScalar, algebra::ScalarAlgebra};

use super::matrix::{Matrix, sample_matrices};

/// Relative size of the imaginary part below which a root is classified by
/// the sign of its flux instead.
const REAL_ROOT_TOLERANCE: f64 = 1.0e-9;

/// Relative idempotency residual above which a projector is refused.
const PROJECTOR_TOLERANCE: f64 = 1.0e-6;

/// Newton steps refining factor jets from exact values: each step doubles
/// the derivative order that is exact, so three cover second-order jets.
const NEWTON_STEPS: usize = 3;

/// A monic polynomial `xᵐ + c_{m−1} xᵐ⁻¹ + … + c₀`, stored as `[c₀, …, c_{m−1}]`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Monic<A>(pub(crate) Vec<A>);

/// The direction a split selects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Direction {
    /// Roots with `Im n > 0`, or real roots carrying positive flux `S_z`.
    Forward,
    /// Roots with `Im n < 0`, or real roots carrying negative flux.
    Backward,
}

impl<A> Monic<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    pub(crate) fn degree(&self) -> usize {
        self.0.len()
    }

    /// `Π (x − root)`.
    pub(crate) fn from_roots(roots: &[A], like: &A) -> Self {
        roots.iter().fold(Self(Vec::new()), |product, root| {
            product.multiply(&Self(vec![root.negate()]), like)
        })
    }

    /// The characteristic polynomial `det(x I − K)` by Faddeev–LeVerrier.
    pub(crate) fn characteristic(generator: &Matrix<A>) -> Self {
        let n = generator.rows();
        let like = generator.get(0, 0);
        let mut coefficients = vec![like.zero_like(); n];
        let mut previous = Matrix::zeros(n, n, like);
        let mut leading = like.one_like();
        for k in 1..=n {
            let current = generator.multiply(&previous).shift(&leading);
            let reciprocal: A::Scalar = nalgebra::convert(-1.0 / k as f64);
            leading = generator.trace_of_product(&current).scale(reciprocal);
            coefficients[n - k] = leading.clone();
            previous = current;
        }
        Self(coefficients)
    }

    /// The full coefficient list `[c₀, …, c_{m−1}, 1]`.
    fn full(&self, like: &A) -> Vec<A> {
        let mut full = self.0.clone();
        full.push(like.one_like());
        full
    }

    fn multiply(&self, rhs: &Self, like: &A) -> Self {
        let (lhs, rhs) = (self.full(like), rhs.full(like));
        let mut product = vec![like.zero_like(); lhs.len() + rhs.len() - 1];
        for (i, a) in lhs.iter().enumerate() {
            for (j, b) in rhs.iter().enumerate() {
                product[i + j] = product[i + j].add(&a.multiply(b));
            }
        }
        product.pop();
        Self(product)
    }

    /// The monic quotient of `self` by `divisor`; the remainder is dropped.
    pub(crate) fn quotient(&self, divisor: &Self, like: &A) -> Self {
        let (n, m) = (self.degree(), divisor.degree());
        let mut remainder = self.full(like);
        let divisor = divisor.full(like);
        let mut quotient = vec![like.zero_like(); n - m + 1];
        for k in (0..=n - m).rev() {
            let coefficient = remainder[k + m].clone();
            for (j, d) in divisor.iter().enumerate() {
                remainder[k + j] = remainder[k + j].subtract(&coefficient.multiply(d));
            }
            quotient[k] = coefficient;
        }
        quotient.pop();
        Self(quotient)
    }

    /// `f(K)` of a square matrix.
    pub(crate) fn of(&self, generator: &Matrix<A>) -> Matrix<A> {
        polynomial(&self.full(generator.get(0, 0)), generator)
    }
}

/// `Σ cᵢ Kⁱ` by Horner's rule.
fn polynomial<A>(coefficients: &[A], generator: &Matrix<A>) -> Matrix<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let n = generator.rows();
    coefficients.iter().rev().fold(
        Matrix::zeros(n, n, generator.get(0, 0)),
        |sum, coefficient| generator.multiply(&sum).shift(coefficient),
    )
}

/// An outgoing (or forward) factor `f`, its complement `g`, and the
/// projector onto the outgoing invariant subspace.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Split<A> {
    pub(crate) outgoing: Monic<A>,
    pub(crate) projector: Matrix<A>,
}

impl<A> Split<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// The split of `generator` whose outgoing roots are the roots of
    /// `outgoing`; `None` when its projector is not idempotent (a branch
    /// point, or roots that are not eigenvalues).
    pub(crate) fn from_factor(generator: &Matrix<A>, outgoing: Monic<A>) -> Option<Self> {
        let like = generator.get(0, 0);
        let complement = Monic::characteristic(generator).quotient(&outgoing, like);
        Self::new(generator, outgoing, &complement)
    }

    /// The split of `generator` into roots of the given direction, chosen
    /// sample by sample from the values of `reference` and continued to
    /// `generator` (see [`continued`]), then refined through the jets;
    /// `None` when a sample does not split evenly or has no split.
    pub(crate) fn select(
        generator: &Matrix<A>,
        reference: &Matrix<A>,
        direction: Direction,
    ) -> Option<Self> {
        let ordered = sample_matrices(generator)
            .iter()
            .zip(&sample_matrices(reference))
            .map(|(sample, reference)| continued(sample, reference, direction))
            .collect::<Option<Vec<_>>>()?;
        Self::from_ordered(generator, &ordered)
    }

    /// The split of `generator` whose chosen roots are the first half of
    /// each sample's `ordered` roots, refined through the jets.
    pub(crate) fn from_ordered(generator: &Matrix<A>, ordered: &[Vec<A::Scalar>]) -> Option<Self> {
        let like = generator.get(0, 0);
        let characteristic = Monic::characteristic(generator);
        let (mut outgoing, mut complement) = factors(generator, ordered)?;

        for _ in 0..NEWTON_STEPS {
            let residual = outgoing.multiply(&complement, like);
            let residual: Vec<A> = residual
                .0
                .iter()
                .zip(&characteristic.0)
                .map(|(product, target)| product.subtract(target))
                .collect();
            let step = sylvester(&outgoing, &complement, like).solve(&Matrix::from_fn(
                residual.len(),
                1,
                |row, _| residual[row].clone(),
            ));
            let m = outgoing.degree();
            for index in 0..m {
                outgoing.0[index] = outgoing.0[index].subtract(step.get(index, 0));
                complement.0[index] = complement.0[index].subtract(step.get(m + index, 0));
            }
        }

        Self::new(generator, outgoing, &complement)
    }

    fn new(generator: &Matrix<A>, outgoing: Monic<A>, complement: &Monic<A>) -> Option<Self> {
        let like = generator.get(0, 0);
        let m = outgoing.degree();
        let unit = Matrix::from_fn(2 * m, 1, |row, _| {
            if row == 0 {
                like.one_like()
            } else {
                like.zero_like()
            }
        });
        let bezout = sylvester(&outgoing, complement, like).solve(&unit);
        let interpolant: Vec<A> = (0..m).map(|index| bezout.get(index, 0).clone()).collect();
        let projector = polynomial(&interpolant, generator).multiply(&complement.of(generator));

        // Idempotency is a property of the values: test it sample by sample
        // without carrying the jets through the product.
        let residual = sample_matrices(&projector)
            .iter()
            .map(|p| row_sum_norm(&(p * p - p)))
            .fold(0.0, |a: f64, b| {
                if a.is_nan() || b.is_nan() {
                    f64::NAN
                } else {
                    a.max(b)
                }
            });
        let scale = projector.value_norm().max(1.0);
        (residual <= PROJECTOR_TOLERANCE * scale).then_some(Self {
            outgoing,
            projector,
        })
    }
}

/// The Sylvester matrix `S` with `S [a; b] = a g + b f` for `deg a, deg b < m`.
fn sylvester<A>(f: &Monic<A>, g: &Monic<A>, like: &A) -> Matrix<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let m = f.degree();
    let (f, g) = (f.full(like), g.full(like));
    Matrix::from_fn(2 * m, 2 * m, |row, column| {
        let (polynomial, shift) = if column < m {
            (&g, column)
        } else {
            (&f, column - m)
        };
        row.checked_sub(shift)
            .and_then(|degree| polynomial.get(degree))
            .map_or_else(|| like.zero_like(), Clone::clone)
    })
}

/// The factors of the first and the second half of each sample's ordered
/// roots, with the roots lifted as constants.
fn factors<A>(generator: &Matrix<A>, ordered: &[Vec<A::Scalar>]) -> Option<(Monic<A>, Monic<A>)>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let like = generator.get(0, 0);
    let lifted = (0..generator.rows())
        .map(|k| {
            let values = ordered
                .iter()
                .map(|roots| roots.get(k).copied())
                .collect::<Option<Vec<_>>>()?;
            Array::from_shape_vec(like.value().raw_dim(), values)
                .ok()
                .map(A::lift_constant)
        })
        .collect::<Option<Vec<A>>>()?;
    let (chosen, rest) = lifted.split_at(generator.rows() / 2);
    Some((
        Monic::from_roots(chosen, like),
        Monic::from_roots(rest, like),
    ))
}

/// The roots of one sample's generator, those in `direction` first, or
/// `None` when they do not split evenly.
pub(crate) fn classify<C: ComplexScalar>(
    generator: &DMatrix<C>,
    direction: Direction,
) -> Option<Vec<C>> {
    let roots = eigenvalues(generator)?;
    let scale = roots
        .iter()
        .map(|root| try_f64(root.modulus()))
        .fold(f64::MIN_POSITIVE, f64::max);
    let mut forward = Vec::new();
    for root in roots.iter().copied() {
        let imaginary = try_f64(root.imaginary());
        forward.push(if imaginary.abs() > REAL_ROOT_TOLERANCE * scale {
            imaginary > 0.0
        } else {
            flux(generator, root, scale)? > 0.0
        });
    }
    let wanted = direction == Direction::Forward;
    let (chosen, rest): (Vec<_>, Vec<_>) = roots
        .iter()
        .zip(forward)
        .partition(|(_, forward)| *forward == wanted);
    (chosen.len() == roots.len() / 2).then(|| {
        chosen
            .into_iter()
            .chain(rest)
            .map(|(root, _)| *root)
            .collect()
    })
}

/// The roots of `generator` continued from the split of `reference` (the
/// same matrix on the real axis): each root takes the class of the reference
/// root it is matched to by the assignment of least total distance, so the
/// split is continuous while no root of one class meets one of the other.
fn continued<C: ComplexScalar>(
    generator: &DMatrix<C>,
    reference: &DMatrix<C>,
    direction: Direction,
) -> Option<Vec<C>> {
    let targets = classify(reference, direction)?;
    matched(generator, &targets)
}

/// The roots of `generator` in the order of `targets`, by the assignment of
/// least total distance.
pub(crate) fn matched<C: ComplexScalar>(generator: &DMatrix<C>, targets: &[C]) -> Option<Vec<C>> {
    let roots = eigenvalues(generator)?;
    let cost = |order: &Vec<usize>| -> f64 {
        order
            .iter()
            .zip(targets)
            .map(|(&k, target)| try_f64((roots[k] - *target).modulus()))
            .sum()
    };
    let best = permutations(roots.len())
        .into_iter()
        .min_by(|a, b| cost(a).total_cmp(&cost(b)))?;
    Some(best.into_iter().map(|k| roots[k]).collect())
}

/// The least distance between a root of the first half of `roots` and one
/// of the second half.
pub(crate) fn separation<C: ComplexScalar>(roots: &[C]) -> f64 {
    let (chosen, rest) = roots.split_at(roots.len() / 2);
    chosen
        .iter()
        .flat_map(|a| rest.iter().map(move |b| try_f64((*a - *b).modulus())))
        .fold(f64::INFINITY, f64::min)
}

/// Every ordering of `0..n`.
fn permutations(n: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![Vec::new()];
    }
    permutations(n - 1)
        .into_iter()
        .flat_map(|shorter| {
            (0..n).map(move |position| {
                let mut order = shorter.clone();
                order.insert(position, n - 1);
                order
            })
        })
        .collect()
}

fn eigenvalues<C: ComplexScalar>(generator: &DMatrix<C>) -> Option<DVector<C>> {
    nalgebra::Schur::try_new(generator.clone(), nalgebra::convert(f64::EPSILON), 10_000)?
        .eigenvalues()
}

/// `Re(x† Γ x)` of an eigenvector `x` of the real root `root`, by one step
/// of inverse iteration; `Γ` pairs the components `(2k, 2k + 1)`.
fn flux<C: ComplexScalar>(generator: &DMatrix<C>, root: C, scale: f64) -> Option<f64> {
    let n = generator.nrows();
    let shift = root + nalgebra::convert::<f64, C>(1.0e-7 * scale);
    let shifted = generator - DMatrix::identity(n, n) * shift;
    let start = DVector::from_element(n, C::one());
    let mut vector = shifted.clone().lu().solve(&start)?;
    vector = shifted.lu().solve(&vector)?;
    let flux = (0..n / 2).fold(C::zero(), |sum, pair| {
        sum + vector[2 * pair].conjugate() * vector[2 * pair + 1]
    });
    let flux = try_f64(flux.real());
    flux.is_finite().then_some(flux)
}

/// The largest absolute row sum; NaN when a value is not finite.
fn row_sum_norm<C: ComplexScalar>(matrix: &DMatrix<C>) -> f64 {
    matrix
        .row_iter()
        .map(|row| row.iter().map(|value| try_f64(value.abs())).sum::<f64>())
        .fold(0.0, |a, b| {
            if a.is_nan() || b.is_nan() {
                f64::NAN
            } else {
                a.max(b)
            }
        })
}

pub(crate) fn try_f64<R: nalgebra::RealField>(value: R) -> f64 {
    nalgebra::try_convert::<_, f64>(value).unwrap_or(f64::NAN)
}

/// The smallest ratio of the smallest to the largest singular value of a
/// basis over every sample: the conditioning of a reference-built basis.
pub(crate) fn conditioning<A>(basis: &Matrix<A>) -> f64
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
{
    sample_matrices(basis)
        .into_iter()
        .map(|sample| {
            nalgebra::SVD::try_new(
                sample,
                false,
                false,
                nalgebra::convert(f64::EPSILON),
                10_000,
            )
            .map_or(0.0, |svd| {
                let values = svd.singular_values.map(try_f64);
                values.min() / values.max()
            })
        })
        .fold(f64::INFINITY, f64::min)
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64 as C;

    use super::*;
    use crate::{
        backend::coupled::berreman::tests::{J, full, jet, random, random_tensor, region},
        material::ConstitutiveTensor,
    };

    fn isotropic(epsilon: C) -> [ConstitutiveTensor<J>; 2] {
        [
            ConstitutiveTensor::Isotropic(jet(epsilon)),
            ConstitutiveTensor::Isotropic(jet(C::new(1.0, 0.0))),
        ]
    }

    #[test]
    fn characteristic_polynomial_vanishes_at_the_eigenvalues() {
        let mut seed = 3;
        let tensors = [
            full(std::array::from_fn(|_| {
                std::array::from_fn(|_| random(&mut seed) + C::new(2.0, 0.0))
            })),
            full(std::array::from_fn(|_| {
                std::array::from_fn(|_| random(&mut seed) + C::new(1.0, 0.0))
            })),
        ];
        let generator = region(&tensors, 0.3).generator;
        let characteristic = Monic::characteristic(&generator);
        let sample = &sample_matrices(&generator)[0];
        for root in sample.clone().schur().eigenvalues().unwrap().iter() {
            let value = characteristic
                .0
                .iter()
                .rev()
                .fold(C::new(1.0, 0.0), |sum, c| sum * root + c.value()[()]);
            assert!(value.norm() < 1e-12, "{value}");
        }
    }

    #[test]
    fn degenerate_isotropic_exterior_projects_onto_te_and_tm_waves() {
        // An isotropic exterior has the double outgoing root κ: the
        // projector is analytic there and takes E_y and η_y to the TE and
        // TM outgoing waves.
        let epsilon = C::new(2.25, 0.0);
        let xi = 0.4;
        let generator = region(&isotropic(epsilon), xi).generator;
        let kappa = (epsilon - xi * xi).sqrt();
        let outgoing = Monic::from_roots(&[jet(kappa), jet(kappa)], &jet(kappa));
        let split = Split::from_factor(&generator, outgoing).unwrap();
        let p = |row, column| split.projector.get(row, column).value()[()];
        let half = C::new(0.5, 0.0);

        // TE: (E_y, −η_x) = (1, κ)/2; TM: (E_x, η_y) = (κ/ε, 1)/2.
        assert!((p(2, 2) - half).norm() < 1e-14);
        assert!((p(3, 2) - kappa * half).norm() < 1e-14);
        assert!((p(1, 1) - half).norm() < 1e-14);
        assert!((p(0, 1) - kappa / epsilon * half).norm() < 1e-14);
        for (row, column) in [(0, 2), (1, 2), (2, 1), (3, 1)] {
            assert!(p(row, column).norm() < 1e-14);
        }
    }

    #[test]
    fn inconsistent_branch_is_refused() {
        let generator = region(&isotropic(C::new(2.25, 0.0)), 0.4).generator;
        let wrong = Monic::from_roots(
            &[jet(C::new(1.0, 0.0)), jet(C::new(1.3, 0.0))],
            generator.get(0, 0),
        );
        assert!(Split::from_factor(&generator, wrong).is_none());
    }

    #[test]
    fn selected_split_matches_the_supplied_split() {
        let mut seed = 5;
        let tensors = [
            random_tensor(&mut seed, true),
            random_tensor(&mut seed, true),
        ];
        let generator = region(&tensors.map(full), 0.3).generator;
        let selected = Split::select(&generator, &generator, Direction::Forward).unwrap();
        let supplied = Split::from_factor(&generator, selected.outgoing.clone()).unwrap();
        let difference = selected
            .projector
            .subtract(&supplied.projector)
            .value_norm();
        assert!(difference < 1e-10, "{difference}");
    }

    #[test]
    fn layer_cutoff_is_degenerate() {
        // At ξ² = ε the forward and backward roots meet at zero.
        let generator = region(&isotropic(C::new(2.25, 0.0)), 1.5).generator;
        assert!(Split::select(&generator, &generator, Direction::Forward).is_none());
    }
}
