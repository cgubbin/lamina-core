//! Per-polarisation medium quantities shared by the 2×2 backends.
//!
//! A region whose ε and μ are diagonal in the stack frame keeps TE
//! `(E_y, H_x, H_z)` and TM `(H_y, E_x, E_z)` decoupled, and each is a scalar
//! two-wave problem with a `±κ` pair of normal wavenumbers:
//!
//! ```text
//! TE: κ² = μ_xx (ε_yy k₀² − k∥²/μ_zz),   Y = κ / μ_xx
//! TM: κ² = ε_xx (μ_yy k₀² − k∥²/ε_zz),   Y = κ / ε_xx
//! ```
//!
//! An isotropic region gives `κ² = ε μ k₀² − k∥²` and `Y = κ/μ` (TE) or `κ/ε`
//! (TM), evaluated by exactly these operations.
//!
//! Each polarisation sees three entries: the *canonical* one weighting the
//! canonical field (`ε_yy` for TE, `μ_yy` for TM), the *transverse* one
//! weighting the in-plane reconstructed component and dividing the
//! admittance (`μ_xx`, `ε_xx`), and the *normal* one dividing the normal
//! reconstructed component (`μ_zz`, `ε_zz`). An isotropic region has one
//! transverse-and-normal entry.
//!
//! # Normal-wavenumber branch
//!
//! `κ` is the principal complex square root of `κ²`, with no further
//! pointwise sign correction, so derivatives are local derivatives on that
//! branch. For real passive problems `κ ≥ 0` for propagating and `Im κ ≥ 0`
//! for evanescent waves. Finite layers always use this convention; exterior
//! values reach a backend through
//! [`ExteriorWavevectors`](crate::backend::ExteriorWavevectors) and may carry
//! a caller's branch choice.

use ndarray::Dimension;

use crate::{
    ComplexScalar,
    algebra::ScalarAlgebra,
    input::{CanonicalCoordinates, Polarisation},
    material::{ConstitutiveClass, ConstitutiveError, ConstitutiveTensor, TensorConstitutiveLift},
    stack::StackRegion,
};

/// Stack-frame ε and μ of a region whose principal axes are the stack axes:
/// scalar, or the diagonals `(xx, yy, zz)`.
#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub enum PrincipalConstitutive<A> {
    Isotropic { epsilon: A, mu: A },
    Diagonal { epsilon: [A; 3], mu: [A; 3] },
}

/// The constitutive entries one polarisation sees; `normal` is `None` when
/// it equals `transverse` (isotropic regions).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ChannelEntries<T> {
    pub(crate) canonical: T,
    pub(crate) transverse: T,
    pub(crate) normal: Option<T>,
}

impl<T> ChannelEntries<T> {
    pub(crate) fn map<U>(self, mut f: impl FnMut(T) -> U) -> ChannelEntries<U> {
        ChannelEntries {
            canonical: f(self.canonical),
            transverse: f(self.transverse),
            normal: self.normal.map(f),
        }
    }

    pub(crate) fn zip<U, V>(
        self,
        other: ChannelEntries<U>,
        mut f: impl FnMut(T, U) -> V,
    ) -> ChannelEntries<V> {
        ChannelEntries {
            canonical: f(self.canonical, other.canonical),
            transverse: f(self.transverse, other.transverse),
            normal: self.normal.zip(other.normal).map(|(a, b)| f(a, b)),
        }
    }
}

impl<A> PrincipalConstitutive<A> {
    /// Keep an isotropic or diagonal `[ε, μ]`; refuse a full one.
    pub(crate) fn new(
        region: StackRegion,
        [epsilon, mu]: [ConstitutiveTensor<A>; 2],
    ) -> Result<Self, ConstitutiveError>
    where
        A: Clone,
    {
        let diagonal = |tensor: ConstitutiveTensor<A>| {
            tensor
                .into_diagonal()
                .ok_or(ConstitutiveError::Unsupported {
                    region,
                    class: ConstitutiveClass::Full,
                })
        };
        match (epsilon, mu) {
            (ConstitutiveTensor::Isotropic(epsilon), ConstitutiveTensor::Isotropic(mu)) => {
                Ok(Self::Isotropic { epsilon, mu })
            }
            (epsilon, mu) => Ok(Self::Diagonal {
                epsilon: diagonal(epsilon)?,
                mu: diagonal(mu)?,
            }),
        }
    }

    /// The entries `polarisation` sees.
    pub(crate) fn channel(&self, polarisation: Polarisation) -> ChannelEntries<&A> {
        match (self, polarisation) {
            (Self::Isotropic { epsilon, mu }, Polarisation::TransverseElectric) => ChannelEntries {
                canonical: epsilon,
                transverse: mu,
                normal: None,
            },
            (Self::Isotropic { epsilon, mu }, Polarisation::TransverseMagnetic) => ChannelEntries {
                canonical: mu,
                transverse: epsilon,
                normal: None,
            },
            (Self::Diagonal { epsilon, mu }, Polarisation::TransverseElectric) => ChannelEntries {
                canonical: &epsilon[1],
                transverse: &mu[0],
                normal: Some(&mu[2]),
            },
            (Self::Diagonal { epsilon, mu }, Polarisation::TransverseMagnetic) => ChannelEntries {
                canonical: &mu[1],
                transverse: &epsilon[0],
                normal: Some(&epsilon[2]),
            },
        }
    }

    /// `[ε, μ]` as diagonals `(xx, yy, zz)`.
    pub(crate) fn diagonal(&self) -> [[&A; 3]; 2] {
        match self {
            Self::Isotropic { epsilon, mu } => [[epsilon; 3], [mu; 3]],
            Self::Diagonal { epsilon, mu } => [epsilon.each_ref(), mu.each_ref()],
        }
    }

    /// `κ²/k₀²` at `k∥ = 0`: `ε μ`, or the canonical times the transverse
    /// entry.
    pub(crate) fn index_squared(&self, polarisation: Polarisation) -> A
    where
        A: ScalarAlgebra,
    {
        match self {
            Self::Isotropic { epsilon, mu } => epsilon.multiply(mu),
            Self::Diagonal { .. } => {
                let channel = self.channel(polarisation);
                channel.transverse.multiply(channel.canonical)
            }
        }
    }

    pub(crate) fn try_map<B, E>(
        &self,
        mut f: impl FnMut(&A) -> Result<B, E>,
    ) -> Result<PrincipalConstitutive<B>, E> {
        Ok(match self {
            Self::Isotropic { epsilon, mu } => PrincipalConstitutive::Isotropic {
                epsilon: f(epsilon)?,
                mu: f(mu)?,
            },
            Self::Diagonal { epsilon, mu } => {
                let [ex, ey, ez] = epsilon;
                let [ux, uy, uz] = mu;
                PrincipalConstitutive::Diagonal {
                    epsilon: [f(ex)?, f(ey)?, f(ez)?],
                    mu: [f(ux)?, f(uy)?, f(uz)?],
                }
            }
        })
    }
}

/// Medium quantities of one region for one polarisation: its constitutive
/// entries, normal wavenumber `κ` (with `κ²` retained before the square
/// root) and characteristic admittance `Y = κ / transverse`.
#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub struct LayerQuantities<A> {
    constitutive: PrincipalConstitutive<A>,
    polarisation: Polarisation,
    kappa: A,
    // Preserve q before sqrt: squaring a sqrt jet at q=0 cannot recover its
    // finite derivatives after the intermediate square-root singularity.
    kappa_squared: A,
    admittance: A,
}

impl<A> LayerQuantities<A>
where
    A: ScalarAlgebra,
{
    /// Evaluate the quantities of `material`, named `region` in an error.
    pub(crate) fn evaluate<E, M>(
        region: StackRegion,
        material: &M,
        coordinates: &CanonicalCoordinates<A>,
        polarisation: Polarisation,
    ) -> Result<Self, ConstitutiveError>
    where
        A: TensorConstitutiveLift<E, M> + Clone,
        A::Scalar: ComplexScalar,
        A::Dimension: Dimension,
    {
        let tensors =
            A::constitutive_tensors(region, material, coordinates.vacuum_angular_wavenumber())?;

        Ok(Self::new(
            PrincipalConstitutive::new(region, tensors)?,
            coordinates,
            polarisation,
        ))
    }

    pub(crate) fn new(
        constitutive: PrincipalConstitutive<A>,
        coordinates: &CanonicalCoordinates<A>,
        polarisation: Polarisation,
    ) -> Self {
        let k0 = coordinates.vacuum_angular_wavenumber();
        let k0_squared = k0.multiply(k0);
        let beta = coordinates.parallel_angular_wavenumber();
        let k_parallel_squared = beta.multiply(beta);

        // κ² = n² k₀² − (transverse/normal) k∥², with n² = κ²/k₀² at k∥ = 0.
        let channel = constitutive.channel(polarisation);
        let k_parallel_squared = match channel.normal {
            None => k_parallel_squared,
            Some(normal) => k_parallel_squared.multiply(&channel.transverse.divide(normal)),
        };
        let kappa_squared = constitutive
            .index_squared(polarisation)
            .multiply(&k0_squared)
            .subtract(&k_parallel_squared);
        let kappa = kappa_squared.sqrt();
        let admittance = kappa.divide(constitutive.channel(polarisation).transverse);

        Self {
            constitutive,
            polarisation,
            kappa,
            kappa_squared,
            admittance,
        }
    }

    /// Replace `κ` by a branch-selected value, with its admittance.
    pub(crate) fn with_kappa(mut self, kappa: A) -> Self {
        self.admittance = kappa.divide(self.factor());
        self.kappa = kappa;
        self
    }
}

impl<A> LayerQuantities<A> {
    /// Return the stack-frame constitutive entries.
    pub(crate) fn constitutive(&self) -> &PrincipalConstitutive<A> {
        &self.constitutive
    }

    /// Return the entries the selected polarisation sees.
    pub(crate) fn channel(&self) -> ChannelEntries<&A> {
        self.constitutive.channel(self.polarisation)
    }

    /// Squared longitudinal wavevector, including regular cutoff derivatives.
    pub(crate) fn kappa_squared(&self) -> &A {
        &self.kappa_squared
    }

    /// Return the selected normal angular wavenumber `κ`.
    pub(crate) fn kappa(&self) -> &A {
        &self.kappa
    }

    /// Consume the quantities and return `κ`.
    pub(crate) fn into_kappa(self) -> A {
        self.kappa
    }

    /// Return the selected polarization.
    pub(crate) fn polarisation(&self) -> Polarisation {
        self.polarisation
    }

    /// Return the transverse entry: `μ_xx` for TE and `ε_xx` for TM.
    pub(crate) fn factor(&self) -> &A {
        self.channel().transverse
    }

    /// Return the normal entry: `μ_zz` for TE and `ε_zz` for TM.
    pub(crate) fn normal_factor(&self) -> &A {
        let channel = self.channel();
        channel.normal.unwrap_or(channel.transverse)
    }

    /// Return the cached characteristic admittance.
    pub(crate) fn admittance(&self) -> &A {
        &self.admittance
    }

    /// Consume the quantities and return the cached characteristic admittance.
    pub(crate) fn into_admittance(self) -> A {
        self.admittance
    }

    pub(crate) fn try_map<B, E>(
        &self,
        mut f: impl FnMut(&A) -> Result<B, E>,
    ) -> Result<LayerQuantities<B>, E> {
        Ok(LayerQuantities {
            constitutive: self.constitutive.try_map(&mut f)?,
            polarisation: self.polarisation,
            kappa: f(&self.kappa)?,
            kappa_squared: f(&self.kappa_squared)?,
            admittance: f(&self.admittance)?,
        })
    }
}

#[cfg(test)]
impl<A> LayerQuantities<A>
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// Isotropic quantities with a given `κ`.
    pub(crate) fn test_fixture(kappa: A, epsilon: A, mu: A, polarisation: Polarisation) -> Self {
        let admittance = match polarisation {
            Polarisation::TransverseElectric => kappa.divide(&mu),
            Polarisation::TransverseMagnetic => kappa.divide(&epsilon),
        };
        Self {
            constitutive: PrincipalConstitutive::Isotropic { epsilon, mu },
            polarisation,
            kappa_squared: kappa.square(),
            kappa,
            admittance,
        }
    }

    pub(crate) fn with_admittance(mut self, admittance: A) -> Self {
        self.admittance = admittance;
        self
    }

    /// The permittivity of isotropic quantities.
    pub(crate) fn epsilon(&self) -> &A {
        match &self.constitutive {
            PrincipalConstitutive::Isotropic { epsilon, .. } => epsilon,
            PrincipalConstitutive::Diagonal { .. } => panic!("expected isotropic quantities"),
        }
    }

    /// The permeability of isotropic quantities.
    pub(crate) fn mu(&self) -> &A {
        match &self.constitutive {
            PrincipalConstitutive::Isotropic { mu, .. } => mu,
            PrincipalConstitutive::Diagonal { .. } => panic!("expected isotropic quantities"),
        }
    }

    /// Evaluate on the real spectral axis.
    pub(crate) fn real_axis<M>(
        material: &M,
        coordinates: &CanonicalCoordinates<A>,
        polarisation: Polarisation,
    ) -> Self
    where
        A: TensorConstitutiveLift<crate::RealAxis, M>,
    {
        Self::evaluate::<crate::RealAxis, M>(
            StackRegion::LeftExterior,
            material,
            coordinates,
            polarisation,
        )
        .unwrap()
    }

    /// Evaluate in the complex spectral plane.
    pub(crate) fn complex_plane<M>(
        material: &M,
        coordinates: &CanonicalCoordinates<A>,
        polarisation: Polarisation,
    ) -> Self
    where
        A: TensorConstitutiveLift<crate::ComplexPlane, M>,
    {
        Self::evaluate::<crate::ComplexPlane, M>(
            StackRegion::LeftExterior,
            material,
            coordinates,
            polarisation,
        )
        .unwrap()
    }
}
#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use ndarray::{Array0, Array1, arr0, array};

    use crate::{
        algebra::Jet0,
        backend::medium::LayerQuantities,
        input::{CanonicalCoordinates, Polarisation},
        test_support::{
            C,
            assertions::{assert_complex_close, assert_dispersion_relation},
            c,
            expected::{
                linear_admittance, linear_epsilon, linear_kappa, linear_mu, quadratic_admittance,
                quadratic_epsilon, quadratic_kappa, quadratic_mu,
            },
            materials::{constant, linear, quadratic, vacuum},
        },
    };

    fn scalar_coordinates(
        vacuum_angular_wavenumber: f64,
        parallel_angular_wavenumber: f64,
    ) -> CanonicalCoordinates<Jet0<Array0<C>>> {
        CanonicalCoordinates::new(
            Jet0::new(arr0(c(vacuum_angular_wavenumber))),
            Jet0::new(arr0(c(parallel_angular_wavenumber))),
        )
    }

    fn sampled_coordinates(
        vacuum_angular_wavenumber: Array1<C>,
        parallel_angular_wavenumber: Array1<C>,
    ) -> CanonicalCoordinates<Jet0<Array1<C>>> {
        CanonicalCoordinates::new(
            Jet0::new(vacuum_angular_wavenumber),
            Jet0::new(parallel_angular_wavenumber),
        )
    }

    // ---------------------------------------------------------------------
    // Polarization-independent medium quantities
    // ---------------------------------------------------------------------

    #[test]
    fn squared_wavevector_retains_cutoff_derivatives_before_sqrt() {
        use crate::{ScalarAlgebra, test_support::jet::HoloJ2};
        let coordinates = CanonicalCoordinates::new(
            HoloJ2::variable(arr0(c(1.0))),
            HoloJ2::constant(arr0(c(2.0))),
        );
        for pol in [
            Polarisation::TransverseElectric,
            Polarisation::TransverseMagnetic,
        ] {
            let quantities = LayerQuantities::complex_plane(&constant(4.0, 1.0), &coordinates, pol);
            let q = quantities.kappa_squared();
            assert!(q.all_finite());
            assert_eq!(q.value()[()], c(0.0));
            assert_eq!(q.first()[()], c(8.0));
            assert_eq!(q.second()[()], c(8.0));
            assert!(!quantities.kappa().all_finite());
        }
    }

    #[test]
    fn medium_quantities_store_evaluated_material_properties() {
        let material = constant(4.0, 2.0);
        let coordinates = scalar_coordinates(3.0, 1.0);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_complex_close(quantities.epsilon()[()], c(4.0), 1e-12);

        assert_complex_close(quantities.mu()[()], c(2.0), 1e-12);
    }

    #[test]
    fn computes_normal_wavenumber() {
        let material = constant(4.0, 2.0);
        let coordinates = scalar_coordinates(3.0, 1.0);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        let expected = c((4.0_f64 * 2.0 * 3.0_f64.powi(2) - 1.0_f64.powi(2)).sqrt());

        assert_complex_close(quantities.kappa()[()], expected, 1e-12);
    }

    #[test]
    fn normal_wavenumber_satisfies_dispersion_relation() {
        let material = constant(3.5, 1.4);

        let k0 = 2.3;
        let k_parallel = 0.7;

        let coordinates = scalar_coordinates(k0, k_parallel);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_dispersion_relation(
            quantities.epsilon()[()],
            quantities.mu()[()],
            quantities.kappa()[()],
            k0,
            k_parallel,
            1e-12,
        );
    }

    #[test]
    fn grazing_incidence_has_zero_normal_wavenumber() {
        let material = vacuum();
        let coordinates = scalar_coordinates(2.0, 2.0);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_complex_close(quantities.kappa()[()], C::new(0.0, 0.0), 1e-12);
    }

    #[test]
    fn evanescent_wave_uses_positive_imaginary_branch() {
        let material = vacuum();
        let coordinates = scalar_coordinates(1.0, 2.0);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_relative_eq!(quantities.kappa()[()].re, 0.0, epsilon = 1e-12,);

        assert!(quantities.kappa()[()].im > 0.0);

        assert_complex_close(quantities.kappa()[()], C::new(0.0, 3.0_f64.sqrt()), 1e-12);
    }

    #[test]
    fn sampled_coordinates_are_evaluated_pointwise() {
        let material = constant(4.0, 1.0);

        let coordinates = sampled_coordinates(
            array![c(1.0), c(2.0), c(3.0)],
            array![c(0.0), c(1.0), c(2.0)],
        );

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        let expected = [
            (4.0_f64 * 1.0_f64.powi(2) - 0.0_f64.powi(2)).sqrt(),
            (4.0_f64 * 2.0_f64.powi(2) - 1.0_f64.powi(2)).sqrt(),
            (4.0_f64 * 3.0_f64.powi(2) - 2.0_f64.powi(2)).sqrt(),
        ];

        for (actual, expected) in quantities.kappa().iter().zip(expected) {
            assert_complex_close(*actual, c(expected), 1e-12);
        }
    }

    // ---------------------------------------------------------------------
    // Dispersive material evaluation
    // ---------------------------------------------------------------------

    #[test]
    fn linear_dispersion_is_evaluated_at_vacuum_wavenumber() {
        let material = linear(2.0, 0.5, 3.0, -0.2);

        let k0 = 4.0;

        let coordinates = scalar_coordinates(k0, 0.7);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_complex_close(
            quantities.epsilon()[()],
            linear_epsilon(&material, k0),
            1e-12,
        );

        assert_complex_close(quantities.mu()[()], linear_mu(&material, k0), 1e-12);
    }

    #[test]
    fn linear_dispersion_changes_with_vacuum_wavenumber() {
        let material = linear(2.0, 0.5, 3.0, -0.2);

        let first_coordinates = scalar_coordinates(1.0, 0.25);

        let second_coordinates = scalar_coordinates(3.0, 0.25);

        let first = LayerQuantities::real_axis(
            &material,
            &first_coordinates,
            Polarisation::TransverseElectric,
        );

        let second = LayerQuantities::real_axis(
            &material,
            &second_coordinates,
            Polarisation::TransverseElectric,
        );

        assert_ne!(first.epsilon()[()], second.epsilon()[()],);

        assert_ne!(first.mu()[()], second.mu()[()],);
    }

    #[test]
    fn linear_dispersion_is_used_when_computing_kappa() {
        let material = linear(2.0, 0.5, 1.5, -0.1);

        let k0 = 2.0;
        let k_parallel = 0.7;

        let coordinates = scalar_coordinates(k0, k_parallel);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_complex_close(
            quantities.kappa()[()],
            linear_kappa(&material, k0, k_parallel),
            1e-12,
        );

        assert_dispersion_relation(
            quantities.epsilon()[()],
            quantities.mu()[()],
            quantities.kappa()[()],
            k0,
            k_parallel,
            1e-12,
        );
    }

    #[test]
    fn quadratic_dispersion_is_evaluated_at_vacuum_wavenumber() {
        let material = quadratic(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);

        let k0 = 2.0;

        let coordinates = scalar_coordinates(k0, 0.5);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_complex_close(
            quantities.epsilon()[()],
            quadratic_epsilon(&material, k0),
            1e-12,
        );

        assert_complex_close(quantities.mu()[()], quadratic_mu(&material, k0), 1e-12);

        assert_complex_close(quantities.epsilon()[()], c(17.0), 1e-12);

        assert_complex_close(quantities.mu()[()], c(38.0), 1e-12);
    }

    #[test]
    fn quadratic_dispersion_is_used_when_computing_kappa() {
        let material = quadratic(2.0, 0.4, 0.1, 1.5, -0.2, 0.05);

        let k0 = 2.5;
        let k_parallel = 0.8;

        let coordinates = scalar_coordinates(k0, k_parallel);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_complex_close(
            quantities.kappa()[()],
            quadratic_kappa(&material, k0, k_parallel),
            1e-12,
        );

        assert_dispersion_relation(
            quantities.epsilon()[()],
            quantities.mu()[()],
            quantities.kappa()[()],
            k0,
            k_parallel,
            1e-12,
        );
    }

    #[test]
    fn medium_real_axis_and_complex_plane_agree_on_real_coordinates() {
        let material = quadratic(2.0, 0.4, 0.1, 1.5, -0.2, 0.05);

        let coordinates = scalar_coordinates(2.5, 0.8);

        let real_axis =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        let complex_plane = LayerQuantities::complex_plane(
            &material,
            &coordinates,
            Polarisation::TransverseElectric,
        );

        assert_complex_close(real_axis.epsilon()[()], complex_plane.epsilon()[()], 1e-12);

        assert_complex_close(real_axis.mu()[()], complex_plane.mu()[()], 1e-12);

        assert_complex_close(real_axis.kappa()[()], complex_plane.kappa()[()], 1e-12);
    }

    // ---------------------------------------------------------------------
    // Polarization-specialized layer quantities
    // ---------------------------------------------------------------------

    #[test]
    fn layer_quantities_preserve_selected_polarisation() {
        let material = vacuum();
        let coordinates = scalar_coordinates(2.0, 0.5);

        let te =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        let tm =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseMagnetic);

        assert_eq!(te.polarisation(), Polarisation::TransverseElectric,);

        assert_eq!(tm.polarisation(), Polarisation::TransverseMagnetic,);
    }

    #[test]
    fn te_factor_is_permeability() {
        let material = constant(4.0, 2.0);
        let coordinates = scalar_coordinates(3.0, 1.0);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_complex_close(quantities.factor()[()], quantities.mu()[()], 1e-12);
    }

    #[test]
    fn tm_factor_is_permittivity() {
        let material = constant(4.0, 2.0);
        let coordinates = scalar_coordinates(3.0, 1.0);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseMagnetic);

        assert_complex_close(quantities.factor()[()], quantities.epsilon()[()], 1e-12);
    }

    #[test]
    fn te_admittance_is_cached_kappa_over_mu() {
        let material = constant(4.0, 2.0);
        let coordinates = scalar_coordinates(3.0, 1.0);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        let expected = quantities.kappa()[()] / quantities.mu()[()];

        assert_complex_close(quantities.admittance()[()], expected, 1e-12);
    }

    #[test]
    fn tm_admittance_is_cached_kappa_over_epsilon() {
        let material = constant(4.0, 2.0);
        let coordinates = scalar_coordinates(3.0, 1.0);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseMagnetic);

        let expected = quantities.kappa()[()] / quantities.epsilon()[()];

        assert_complex_close(quantities.admittance()[()], expected, 1e-12);
    }

    #[test]
    fn cached_admittance_matches_expected_dispersive_te_value() {
        let material = linear(2.0, 0.5, 1.5, -0.1);

        let k0 = 2.0;
        let k_parallel = 0.7;

        let coordinates = scalar_coordinates(k0, k_parallel);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseElectric);

        assert_complex_close(
            quantities.admittance()[()],
            linear_admittance(&material, k0, k_parallel, Polarisation::TransverseElectric),
            1e-12,
        );
    }

    #[test]
    fn cached_admittance_matches_expected_dispersive_tm_value() {
        let material = quadratic(2.0, 0.4, 0.1, 1.5, -0.2, 0.05);

        let k0 = 2.5;
        let k_parallel = 0.8;

        let coordinates = scalar_coordinates(k0, k_parallel);

        let quantities =
            LayerQuantities::real_axis(&material, &coordinates, Polarisation::TransverseMagnetic);

        assert_complex_close(
            quantities.admittance()[()],
            quadratic_admittance(&material, k0, k_parallel, Polarisation::TransverseMagnetic),
            1e-12,
        );
    }

    #[test]
    fn layer_real_axis_and_complex_plane_agree_on_real_coordinates() {
        let material = quadratic(2.0, 0.4, 0.1, 1.5, -0.2, 0.05);

        let coordinates = scalar_coordinates(2.5, 0.8);

        for polarisation in [
            Polarisation::TransverseElectric,
            Polarisation::TransverseMagnetic,
        ] {
            let real_axis = LayerQuantities::real_axis(&material, &coordinates, polarisation);

            let complex_plane =
                LayerQuantities::complex_plane(&material, &coordinates, polarisation);

            assert_complex_close(real_axis.epsilon()[()], complex_plane.epsilon()[()], 1e-12);

            assert_complex_close(real_axis.mu()[()], complex_plane.mu()[()], 1e-12);

            assert_complex_close(real_axis.kappa()[()], complex_plane.kappa()[()], 1e-12);

            assert_complex_close(
                real_axis.admittance()[()],
                complex_plane.admittance()[()],
                1e-12,
            );
        }
    }

    #[test]
    fn diagonal_quantities_follow_the_closed_forms() {
        use super::PrincipalConstitutive;

        let (k0, beta) = (2.3, 1.7);
        let epsilon = [c(2.1), C::new(3.4, 0.2), c(4.6)];
        let mu = [c(1.3), c(0.8), C::new(1.1, 0.05)];
        let coordinates = scalar_coordinates(k0, beta);
        let constitutive = PrincipalConstitutive::Diagonal {
            epsilon: epsilon.map(|entry| Jet0::new(arr0(entry))),
            mu: mu.map(|entry| Jet0::new(arr0(entry))),
        };
        for (polarisation, transverse, canonical, normal) in [
            (Polarisation::TransverseElectric, mu[0], epsilon[1], mu[2]),
            (
                Polarisation::TransverseMagnetic,
                epsilon[0],
                mu[1],
                epsilon[2],
            ),
        ] {
            let quantities = LayerQuantities::new(constitutive.clone(), &coordinates, polarisation);
            let kappa_squared = transverse * (canonical * k0 * k0 - beta * beta / normal);
            assert_complex_close(quantities.kappa_squared()[()], kappa_squared, 1e-12);
            assert_complex_close(quantities.kappa()[()], kappa_squared.sqrt(), 1e-12);
            assert_complex_close(
                quantities.admittance()[()],
                kappa_squared.sqrt() / transverse,
                1e-12,
            );
            assert_eq!(quantities.factor()[()], transverse);
            assert_eq!(quantities.normal_factor()[()], normal);
        }
    }
}
