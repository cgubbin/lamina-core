use num_traits::{Float, One};

use crate::{
    ElectromagneticDissipation, ElectromagneticEnergy,
    algebra::{
        CartesianVectorAlgebra, ComplexJet, Jet, RealCartesianVectorAlgebra, RealScalarAlgebra,
        ScalarAlgebra,
    },
    observable::field::constitutive::{
        ElectromagneticDissipationCoefficients, ElectromagneticEnergyCoefficients,
    },
};

use super::{ConstitutiveFields, IsotropicConstitutiveParameters, material::SampledConstitutive};

#[derive(Clone, Debug, PartialEq)]
pub struct ElectromagneticIntensities<V> {
    electric: V,
    magnetic: V,
}

impl<V> ElectromagneticIntensities<V> {
    pub(crate) fn new(electric: V, magnetic: V) -> Self {
        Self { electric, magnetic }
    }

    pub fn electric(&self) -> &V {
        &self.electric
    }

    pub fn magnetic(&self) -> &V {
        &self.magnetic
    }

    pub fn into_parts(self) -> (V, V) {
        (self.electric, self.magnetic)
    }

    pub fn map_vectors<U>(self, f: impl Fn(V) -> U) -> ElectromagneticIntensities<U> {
        ElectromagneticIntensities {
            electric: f(self.electric),
            magnetic: f(self.magnetic),
        }
    }

    pub(crate) fn into_energy(
        self,
        coefficients: &ElectromagneticEnergyCoefficients<V>,
    ) -> ElectromagneticEnergy<V>
    where
        V: ScalarAlgebra,
        V::Scalar: Float,
    {
        let (electric_norm, magnetic_norm) = self.into_parts();

        let electric = electric_norm
            .multiply(coefficients.electric())
            .scale(quarter());

        let magnetic = magnetic_norm
            .multiply(coefficients.magnetic())
            .scale(quarter());

        let total = electric.add(&magnetic);

        ElectromagneticEnergy::new(electric, magnetic, total)
    }

    pub(crate) fn into_dissipation(
        self,
        coefficients: &ElectromagneticDissipationCoefficients<V>,
    ) -> ElectromagneticDissipation<V>
    where
        V: ScalarAlgebra,
    {
        let (electric_norm, magnetic_norm) = self.into_parts();

        let electric = electric_norm.multiply(coefficients.electric());

        let magnetic = magnetic_norm.multiply(coefficients.magnetic());

        let total = electric.add(&magnetic);

        ElectromagneticDissipation::new(electric, magnetic, total)
    }
}

/// Pointwise Cartesian electric and magnetic phasor fields.
///
/// The field uses the electromagnetic normalization chosen by the producing
/// backend. The electric and magnetic vectors share the same ndarray sampling
/// shape.
///
/// The complex Poynting vector uses:
///
/// ```text
/// S = 1/2 E × H*
/// ```
///
/// and the time-averaged Poynting vector is its real part.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectromagneticFields<V> {
    electric: V,
    magnetic: V,
}

impl<V> ElectromagneticFields<V> {
    pub(crate) fn new(electric: V, magnetic: V) -> Self {
        Self { electric, magnetic }
    }

    pub fn electric(&self) -> &V {
        &self.electric
    }

    pub fn magnetic(&self) -> &V {
        &self.magnetic
    }

    pub fn into_parts(self) -> (V, V) {
        (self.electric, self.magnetic)
    }

    /// `D = ε E` and `B = μ H`, component by component for diagonal media.
    pub(crate) fn into_constitutive_fields(
        self,
        constitutive: &SampledConstitutive<IsotropicConstitutiveParameters<V::ScalarAlgebra>, V>,
    ) -> ConstitutiveFields<V>
    where
        V: CartesianVectorAlgebra,
    {
        let (electric, magnetic) = self.into_parts();

        match constitutive {
            SampledConstitutive::Isotropic(constitutive) => ConstitutiveFields::new(
                electric.multiply_by_scalar(constitutive.epsilon()),
                magnetic.multiply_by_scalar(constitutive.mu()),
            ),
            SampledConstitutive::Diagonal {
                electric: epsilon,
                magnetic: mu,
            } => ConstitutiveFields::new(
                electric.multiply_components(epsilon),
                magnetic.multiply_components(mu),
            ),
        }
    }

    /// `Σ_i w_ii |E_i|²` and `Σ_i v_ii |H_i|²` for per-component weights `w`
    /// and `v`.
    fn weighted_intensities(&self, electric: &V, magnetic: &V) -> [V::ScalarAlgebra; 2]
    where
        V: RealCartesianVectorAlgebra,
    {
        [
            self.electric
                .hermitian_dot(&self.electric.multiply_components(electric)),
            self.magnetic
                .hermitian_dot(&self.magnetic.multiply_components(magnetic)),
        ]
    }

    /// Energy density `¼ Re Σ_i w_ii |F_i|²` for per-component dispersive
    /// weights `w = ∂(k₀ ε)/∂k₀` (electric) and `∂(k₀ μ)/∂k₀` (magnetic).
    pub(crate) fn into_diagonal_energy(
        self,
        electric: &V,
        magnetic: &V,
    ) -> ElectromagneticEnergy<V::RealScalarAlgebra>
    where
        V: RealCartesianVectorAlgebra,
        V::RealScalarAlgebra: ScalarAlgebra,
        <V::RealScalarAlgebra as Jet>::Scalar: Float,
    {
        let [electric, magnetic] = self
            .weighted_intensities(electric, magnetic)
            .map(|weighted| V::scalar_real(weighted).scale(quarter()));
        let total = electric.add(&magnetic);

        ElectromagneticEnergy::new(electric, magnetic, total)
    }

    /// Dissipation density `½ Re(k₀) Im Σ_i ε_ii |E_i|²` (electric) and
    /// likewise with `μ` (magnetic).
    pub(crate) fn into_diagonal_dissipation(
        self,
        epsilon: &V,
        mu: &V,
        vacuum_angular_wavenumber: &V::ScalarAlgebra,
    ) -> ElectromagneticDissipation<<V::ScalarAlgebra as ComplexJet>::RealJet>
    where
        V: RealCartesianVectorAlgebra,
        V::ScalarAlgebra: RealScalarAlgebra,
        <V::ScalarAlgebra as ComplexJet>::RealJet: ScalarAlgebra,
        <<V::ScalarAlgebra as ComplexJet>::RealJet as Jet>::Scalar: Float,
    {
        let common = vacuum_angular_wavenumber.real().scale(half());
        let [electric, magnetic] = self
            .weighted_intensities(epsilon, mu)
            .map(|weighted| weighted.imaginary().multiply(&common));
        let total = electric.add(&magnetic);

        ElectromagneticDissipation::new(electric, magnetic, total)
    }

    pub fn map_vectors<U>(self, f: impl Fn(V) -> U) -> ElectromagneticFields<U> {
        ElectromagneticFields {
            electric: f(self.electric),
            magnetic: f(self.magnetic),
        }
    }
}

impl<V> ElectromagneticFields<V> {
    pub fn into_magnitude_squared(self) -> ElectromagneticIntensities<V::RealScalarAlgebra>
    where
        V: RealCartesianVectorAlgebra,
    {
        ElectromagneticIntensities {
            electric: self.electric.magnitude_squared(),
            magnetic: self.magnetic.magnitude_squared(),
        }
    }

    /// Return the pointwise squared electric-field magnitude.
    pub fn electric_magnitude_squared(&self) -> V::RealScalarAlgebra
    where
        V: RealCartesianVectorAlgebra,
    {
        self.electric.magnitude_squared()
    }

    /// Return the pointwise squared magnetic-field magnitude.
    pub fn magnetic_magnitude_squared(&self) -> V::RealScalarAlgebra
    where
        V: RealCartesianVectorAlgebra,
    {
        self.magnetic.magnitude_squared()
    }

    /// Return the pointwise complex Poynting vector.
    ///
    /// This evaluates `1/2 E × H*`.
    ///
    /// Fields use the exp(-iωt) phasor convention.
    pub fn complex_poynting_vector(&self) -> V
    where
        V: RealCartesianVectorAlgebra,
        V::Coefficient: One
            + Copy
            + std::ops::Add<Output = V::Coefficient>
            + std::ops::Div<Output = V::Coefficient>,
    {
        self.electric
            .cross(&self.magnetic.conjugated())
            .scale_by_constant(half())
    }

    /// Return the time average poynting vector
    ///
    /// This evaluates `1/2 Re(E × H*)`.
    ///
    /// Fields use the exp(-iωt) phasor convention.
    pub fn time_averaged_poynting_vector(&self) -> V::RealVector
    where
        V: RealCartesianVectorAlgebra,
        V::Coefficient: One
            + Copy
            + std::ops::Add<Output = V::Coefficient>
            + std::ops::Div<Output = V::Coefficient>,
    {
        self.complex_poynting_vector().real()
    }
}

fn quarter<C>() -> C
where
    C: One + Copy + std::ops::Add<Output = C> + std::ops::Div<Output = C>,
{
    let one = C::one();
    let two = one + one;

    one / (two + two)
}

fn half<C>() -> C
where
    C: One + Copy + std::ops::Add<Output = C> + std::ops::Div<Output = C>,
{
    let one = C::one();

    one / (one + one)
}

#[cfg(test)]
mod tests {
    use ndarray::{Array1, Ix1, arr1};
    use num_complex::Complex64;

    use crate::{
        algebra::{Jet0, RealParameter},
        field::VectorField,
    };

    use super::*;

    type C = Complex64;
    type D = Ix1;
    type Vector = Jet0<VectorField<C, D>, RealParameter>;
    type Field = ElectromagneticFields<Vector>;

    const TOLERANCE: f64 = 1.0e-12;

    fn c(real: f64, imaginary: f64) -> C {
        C::new(real, imaginary)
    }

    fn scalar_vector(value: f64) -> Vector {
        Jet0::new(VectorField::new_unchecked(
            arr1(&[c(value, 0.0)]),
            arr1(&[c(value + 0.1, 0.0)]),
            arr1(&[c(value + 0.2, 0.0)]),
        ))
    }

    fn vector(x: &[C], y: &[C], z: &[C]) -> Vector {
        Jet0::new(VectorField::new_unchecked(
            Array1::from_vec(x.to_vec()),
            Array1::from_vec(y.to_vec()),
            Array1::from_vec(z.to_vec()),
        ))
    }

    fn assert_real_close(actual: f64, expected: f64) {
        let error = (actual - expected).abs();

        assert!(
            error <= TOLERANCE,
            "expected {expected:e}, \
             got {actual:e}; \
             absolute error = {error:e}",
        );
    }

    fn assert_complex_close(actual: C, expected: C) {
        let error = (actual - expected).norm();

        assert!(
            error <= TOLERANCE,
            "expected {expected:?}, \
             got {actual:?}; \
             absolute error = {error:e}",
        );
    }

    fn assert_real_array_close(actual: &Array1<f64>, expected: &Array1<f64>) {
        assert_eq!(actual.raw_dim(), expected.raw_dim(),);

        for (&actual, &expected) in actual.iter().zip(expected.iter()) {
            assert_real_close(actual, expected);
        }
    }

    fn assert_complex_array_close(actual: &Array1<C>, expected: &Array1<C>) {
        assert_eq!(actual.raw_dim(), expected.raw_dim(),);

        for (&actual, &expected) in actual.iter().zip(expected.iter()) {
            assert_complex_close(actual, expected);
        }
    }

    #[test]
    fn construction_preserves_electric_and_magnetic_fields() {
        let electric = scalar_vector(1.0);
        let magnetic = scalar_vector(2.0);

        let field = Field::new(electric.clone(), magnetic.clone());

        assert_eq!(field.electric(), &electric,);

        assert_eq!(field.magnetic(), &magnetic,);

        assert_eq!(field.into_parts(), (electric, magnetic),);
    }

    #[test]
    fn magnitude_methods_use_corresponding_vector() {
        let electric = vector(
            &[c(3.0, 4.0), c(1.0, 0.0)],
            &[c(0.0, 2.0), c(0.0, 2.0)],
            &[c(1.0, 0.0), c(2.0, 0.0)],
        );

        let magnetic = vector(
            &[c(1.0, 0.0), c(2.0, 0.0)],
            &[c(2.0, 0.0), c(3.0, 0.0)],
            &[c(3.0, 0.0), c(4.0, 0.0)],
        );

        let field = Field::new(electric, magnetic);

        assert_real_array_close(&field.electric_magnitude_squared(), &arr1(&[30.0, 9.0]));

        assert_real_array_close(&field.magnetic_magnitude_squared(), &arr1(&[14.0, 29.0]));
    }

    #[test]
    fn complex_poynting_is_half_e_cross_conjugate_h() {
        let zero = c(0.0, 0.0);

        let electric = vector(&[c(1.0, 1.0)], &[zero], &[zero]);

        let magnetic = vector(&[zero], &[c(2.0, 1.0)], &[zero]);

        let field = Field::new(electric, magnetic);

        let result = field.complex_poynting_vector();

        let expected_z = c(0.5, 0.0) * c(1.0, 1.0) * c(2.0, 1.0).conj();

        assert_complex_array_close(result.x(), &arr1(&[zero]));

        assert_complex_array_close(result.y(), &arr1(&[zero]));

        assert_complex_array_close(result.z(), &arr1(&[expected_z]));
    }

    #[test]
    fn time_averaged_poynting_is_real_part_of_complex_poynting() {
        let zero = c(0.0, 0.0);

        let field = Field::new(
            vector(&[c(1.0, 1.0)], &[zero], &[zero]),
            vector(&[zero], &[c(2.0, 1.0)], &[zero]),
        );

        let complex = field.complex_poynting_vector();

        let averaged = field.time_averaged_poynting_vector();

        assert_real_array_close(averaged.x(), &complex.x().mapv(|value| value.re));

        assert_real_array_close(averaged.y(), &complex.y().mapv(|value| value.re));

        assert_real_array_close(averaged.z(), &complex.z().mapv(|value| value.re));
    }
}
