//! The Berreman system of a homogeneous region.
//!
//! With `η = Z₀H`, `∂_x = i k∥` and `∂_y = 0`, Maxwell's equations
//! `∇×E = i k₀ μ η` and `∇×η = −i k₀ ε E` split into two algebraic equations,
//! which eliminate the normal components,
//!
//! ```text
//! E_z = −(ξ η_y + ε_zx E_x + ε_zy E_y) / ε_zz
//! η_z =  (ξ E_y − μ_zx η_x − μ_zy η_y) / μ_zz,        ξ = k∥ / k₀,
//! ```
//!
//! and four first-order equations for the tangential state
//! `ψ = (E_x, η_y, E_y, −η_x)`, continuous across interfaces:
//!
//! ```text
//! dψ/dz = i K ψ,   K = k₀ Δ,
//! Δ ψ = ((μη)_y + ξ E_z,  (εE)_x,  −(μη)_x,  (εE)_y − ξ η_z).
//! ```
//!
//! The eigenvalues of `K` are the normal wavenumbers of the region's four
//! plane waves. With μ = I the rows are those of the 4×4 matrix in the
//! design; a diagonal region gives the TM block `(E_x, η_y)` and the TE block
//! `(E_y, −η_x)`, with eigenvalues `±κ_TM` and `±κ_TE`. A region without
//! `xy`, `yz` components keeps the two blocks invariant.

use ndarray::Dimension;

use crate::{
    ComplexScalar,
    algebra::ScalarAlgebra,
    backend::{Channel, ChannelError},
    input::CanonicalCoordinates,
    material::{ConstitutiveError, ConstitutiveTensor, TensorConstitutiveLift},
    stack::StackRegion,
};

use super::matrix::Matrix;

/// Index of `E_x`, `η_y`, `E_y` and `−η_x` in the tangential state.
const EX: usize = 0;
const HY: usize = 1;
const EY: usize = 2;
const HX: usize = 3;

/// The components of the tangential state a channel carries.
pub(crate) const fn channel_components(channel: Channel) -> &'static [usize] {
    match channel {
        Channel::TransverseMagnetic => &[EX, HY],
        Channel::TransverseElectric => &[EY, HX],
        Channel::Coupled => &[EX, HY, EY, HX],
    }
}

/// Positions within a channel's state of the given global components, in the
/// given order, skipping those the channel does not carry.
pub(crate) fn local_positions(channel: Channel, components: &[usize]) -> Vec<usize> {
    let carried = channel_components(channel);
    components
        .iter()
        .filter_map(|component| carried.iter().position(|c| c == component))
        .collect()
}

/// The reference vectors of the exterior and layer bases: `E_y` then `η_y`,
/// which in an isotropic region project onto the TE and the TM wave.
pub(crate) fn reference_positions(channel: Channel) -> Vec<usize> {
    local_positions(channel, &[EY, HY])
}

/// A row vector over the full tangential state `ψ`.
pub(crate) type Row<A> = [A; 4];

/// The linear maps from `ψ` to the eliminated normal components.
#[derive(Clone, Debug, PartialEq)]
pub struct Elimination<A> {
    pub(crate) electric: Row<A>,
    pub(crate) magnetic: Row<A>,
}

/// One region's Berreman generator `K` restricted to a channel, with the
/// elimination of the normal fields.
#[derive(Clone, Debug, PartialEq)]
pub struct Region<A> {
    pub(crate) generator: Matrix<A>,
    pub(crate) elimination: Elimination<A>,
    /// The stack-frame `[ε, μ]`.
    pub(crate) tensors: [ConstitutiveTensor<A>; 2],
}

impl<A> Region<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// Evaluate `material`, named `region`, for `channel`. A scalar channel
    /// needs a separable region.
    pub(crate) fn evaluate<E, M, Error>(
        region: StackRegion,
        material: &M,
        coordinates: &CanonicalCoordinates<A>,
        channel: Channel,
    ) -> Result<Self, Error>
    where
        A: TensorConstitutiveLift<E, M>,
        Error: From<ConstitutiveError> + From<ChannelError>,
    {
        if channel != Channel::Coupled && !A::is_separable(material) {
            return Err(ChannelError::NotInvariant { region, channel }.into());
        }
        let tensors =
            A::constitutive_tensors(region, material, coordinates.vacuum_angular_wavenumber())?;

        Ok(Self::new(
            &tensors,
            coordinates.vacuum_angular_wavenumber(),
            coordinates.parallel_angular_wavenumber(),
            channel,
        ))
    }

    /// The region of stack-frame `[ε, μ]` at `k₀` and `k∥`.
    pub(crate) fn new(
        tensors: &[ConstitutiveTensor<A>; 2],
        k0: &A,
        k_parallel: &A,
        channel: Channel,
    ) -> Self {
        let [epsilon, mu] = tensors;
        let (generator, elimination) = generator(&epsilon.full(), &mu.full(), k0, k_parallel);
        let components = channel_components(channel);

        Self {
            generator: generator.select(components, components),
            elimination,
            tensors: tensors.clone(),
        }
    }
}

/// `K = k₀ Δ` of full `ε`, `μ`, and the elimination rows.
fn generator<A>(
    epsilon: &[[A; 3]; 3],
    mu: &[[A; 3]; 3],
    k0: &A,
    k_parallel: &A,
) -> (Matrix<A>, Elimination<A>)
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let xi = k_parallel.divide(k0);
    let unit = |component: usize, sign: bool| -> Row<A> {
        std::array::from_fn(|index| match (index == component, sign) {
            (false, _) => xi.zero_like(),
            (true, true) => xi.one_like(),
            (true, false) => xi.one_like().negate(),
        })
    };
    let combine = |terms: &[(&A, &Row<A>)]| -> Row<A> {
        std::array::from_fn(|index| {
            terms
                .iter()
                .map(|(weight, row)| weight.multiply(&row[index]))
                .reduce(|sum, term| sum.add(&term))
                .unwrap_or_else(|| xi.zero_like())
        })
    };
    let scale = |row: Row<A>, factor: &A| row.map(|entry| entry.multiply(factor));

    let e_x = unit(EX, true);
    let e_y = unit(EY, true);
    let eta_x = unit(HX, false);
    let eta_y = unit(HY, true);
    let minus_one = xi.one_like().negate();
    let minus_xi = xi.negate();

    let e_z = scale(
        combine(&[
            (&xi, &eta_y),
            (&epsilon[2][0], &e_x),
            (&epsilon[2][1], &e_y),
        ]),
        &epsilon[2][2].reciprocal().negate(),
    );
    let eta_z = scale(
        combine(&[
            (&xi, &e_y),
            (&mu[2][0].negate(), &eta_x),
            (&mu[2][1].negate(), &eta_y),
        ]),
        &mu[2][2].reciprocal(),
    );

    let electric = [&e_x, &e_y, &e_z];
    let magnetic = [&eta_x, &eta_y, &eta_z];
    let product = |tensor: &[[A; 3]; 3], field: [&Row<A>; 3], row: usize| {
        combine(&[
            (&tensor[row][0], field[0]),
            (&tensor[row][1], field[1]),
            (&tensor[row][2], field[2]),
        ])
    };
    let mu_eta_x = product(mu, magnetic, 0);
    let mu_eta_y = product(mu, magnetic, 1);
    let epsilon_e_x = product(epsilon, electric, 0);
    let epsilon_e_y = product(epsilon, electric, 1);

    let one = xi.one_like();
    let rows = [
        combine(&[(&one, &mu_eta_y), (&xi, &e_z)]),
        epsilon_e_x,
        combine(&[(&minus_one, &mu_eta_x)]),
        combine(&[(&one, &epsilon_e_y), (&minus_xi, &eta_z)]),
    ];
    let generator = Matrix::from_fn(4, 4, |row, column| rows[row][column].multiply(k0));

    (
        generator,
        Elimination {
            electric: e_z,
            magnetic: eta_z,
        },
    )
}

#[cfg(test)]
pub(crate) mod tests {
    use ndarray::{Ix0, arr0};
    use num_complex::Complex64 as C;

    use super::*;
    use crate::algebra::{ArrayJet0, RealParameter};

    pub(crate) type J = ArrayJet0<C, Ix0, RealParameter>;

    pub(crate) fn jet(value: C) -> J {
        J::lift_constant(arr0(value))
    }

    pub(crate) fn full(entries: [[C; 3]; 3]) -> ConstitutiveTensor<J> {
        ConstitutiveTensor::Full(entries.map(|row| row.map(jet)))
    }

    /// A deterministic pseudo-random complex number in the unit square.
    pub(crate) fn random(seed: &mut u64) -> C {
        let mut next = || {
            *seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((*seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
        };
        C::new(next(), next())
    }

    /// The coupled region of `[ε, μ]` at `k₀ = 1`, `k∥ = ξ`.
    pub(crate) fn region(tensors: &[ConstitutiveTensor<J>; 2], xi: f64) -> Region<J> {
        Region::new(
            tensors,
            &jet(C::new(1.0, 0.0)),
            &jet(C::new(xi, 0.0)),
            Channel::Coupled,
        )
    }

    fn delta(tensors: &[ConstitutiveTensor<J>; 2], xi: f64) -> nalgebra::Matrix4<C> {
        let generator = region(tensors, xi).generator;
        nalgebra::Matrix4::from_fn(|row, column| generator.get(row, column).value()[()])
    }

    /// `[[0, 1], [sign, 0]]` on `(E_x, η_y)` and on `(E_y, −η_x)`: `Γ` for
    /// `sign = 1`, `J` for `−1`.
    fn pairing(sign: f64) -> nalgebra::Matrix4<C> {
        nalgebra::Matrix4::from_fn(|row, column| match (row % 2, column) {
            (0, column) if column == row + 1 => C::new(1.0, 0.0),
            (1, column) if column + 1 == row => C::new(sign, 0.0),
            _ => C::new(0.0, 0.0),
        })
    }

    /// A random tensor with `T_ji = conjugate(T_ij)`, or `T_ij` without
    /// `conjugate`, and diagonal shifted by 2.
    pub(crate) fn random_tensor(seed: &mut u64, conjugate: bool) -> [[C; 3]; 3] {
        let mut entries = [[C::new(0.0, 0.0); 3]; 3];
        for row in 0..3 {
            for column in row..3 {
                let value = random(seed);
                entries[row][column] = value;
                entries[column][row] = if conjugate { value.conj() } else { value };
            }
            entries[row][row] += C::new(2.0, 0.0);
            if conjugate {
                entries[row][row].im = 0.0;
            }
        }
        entries
    }

    #[test]
    fn isotropic_and_diagonal_regions_are_block_anti_diagonal() {
        let diagonal = |values: [f64; 3]| {
            ConstitutiveTensor::Diagonal(values.map(|value| jet(C::new(value, 0.1))))
        };
        let tensors = [diagonal([2.0, 3.0, 4.0]), diagonal([1.1, 1.2, 1.3])];
        let xi = 0.7;
        let delta = delta(&tensors, xi);
        let c = |value: f64| C::new(value, 0.1);
        let xi2 = C::new(xi * xi, 0.0);

        let expected = [
            (0, 1, c(1.2) - xi2 / c(4.0)),
            (1, 0, c(2.0)),
            (2, 3, c(1.1)),
            (3, 2, c(3.0) - xi2 / c(1.3)),
        ];
        for row in 0..4 {
            for column in 0..4 {
                let value = expected
                    .iter()
                    .find(|(r, c, _)| *r == row && *c == column)
                    .map_or(C::new(0.0, 0.0), |(_, _, value)| *value);
                assert!(
                    (delta[(row, column)] - value).norm() < 1e-14,
                    "Δ[{row}, {column}]"
                );
            }
        }
    }

    #[test]
    fn mu_identity_matches_the_closed_form() {
        let mut seed = 7;
        let epsilon: [[C; 3]; 3] =
            std::array::from_fn(|_| std::array::from_fn(|_| random(&mut seed) + C::new(2.0, 0.0)));
        let identity = ConstitutiveTensor::Isotropic(jet(C::new(1.0, 0.0)));
        let xi = 0.45;
        let delta = delta(&[full(epsilon), identity], xi);
        let e = |row: usize, column: usize| epsilon[row][column];
        let (x, y, z) = (0, 1, 2);
        let xi = C::new(xi, 0.0);
        let zero = C::new(0.0, 0.0);
        let one = C::new(1.0, 0.0);
        let expected = [
            [
                -xi * e(z, x) / e(z, z),
                one - xi * xi / e(z, z),
                -xi * e(z, y) / e(z, z),
                zero,
            ],
            [
                e(x, x) - e(x, z) * e(z, x) / e(z, z),
                -xi * e(x, z) / e(z, z),
                e(x, y) - e(x, z) * e(z, y) / e(z, z),
                zero,
            ],
            [zero, zero, zero, one],
            [
                e(y, x) - e(y, z) * e(z, x) / e(z, z),
                -xi * e(y, z) / e(z, z),
                e(y, y) - xi * xi - e(y, z) * e(z, y) / e(z, z),
                zero,
            ],
        ];
        for row in 0..4 {
            for column in 0..4 {
                assert!((delta[(row, column)] - expected[row][column]).norm() < 1e-13);
            }
        }
    }

    #[test]
    fn hermitian_tensors_conserve_flux() {
        // Γ Δ = Δ† Γ for Hermitian ε, μ and real ξ.
        let mut seed = 11;
        for _ in 0..5 {
            let tensors = [
                full(random_tensor(&mut seed, true)),
                full(random_tensor(&mut seed, true)),
            ];
            let delta = delta(&tensors, 0.83);
            let residual = pairing(1.0) * delta - delta.adjoint() * pairing(1.0);
            assert!(residual.norm() < 1e-13 * delta.norm(), "{residual}");
        }
    }

    #[test]
    fn symmetric_tensors_are_reciprocal() {
        // Δ(ξ)ᵀ J + J Δ(−ξ) = 0 for symmetric ε, μ.
        let mut seed = 13;
        for _ in 0..5 {
            let tensors = [
                full(random_tensor(&mut seed, false)),
                full(random_tensor(&mut seed, false)),
            ];
            let forward = delta(&tensors, 0.61);
            let backward = delta(&tensors, -0.61);
            let residual = forward.transpose() * pairing(-1.0) + pairing(-1.0) * backward;
            assert!(residual.norm() < 1e-13 * forward.norm(), "{residual}");
        }
    }
}
