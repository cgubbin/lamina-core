//! Shape-preserving 2×2 transfer matrices.
//!
//! This module provides an entry representation:
//!
//! [`Transfer2Entries`] is an internal algebraic representation used while
//!   evaluating value-only, first-order, and second-order transfer matrices.
//!
//! [`Transfer2Entries`] is generic over its scalar-entry representation. Its
//! entries may therefore be sampled arrays, first-order jets, or second-order
//! jets.
//!
//! Transfer matrices compose by ordinary matrix multiplication. Each layer
//! matrix maps the state at its right boundary to the state at its left
//! boundary:
//!
//! ```text
//! ψ_left = L ψ_right.
//! ```
//!
//! Layers are appended in physical left-to-right order. If the accumulated
//! matrix contains the layers already encountered and `L` is the next layer to
//! the right, accumulation is:
//!
//! ```text
//! M_total <- M_total L.
//! ```
//!
//! Consequently, for layers `L₀, L₁, …, Lₙ` stored from left to right,
//!
//! ```text
//! ψ_left = L₀ L₁ … Lₙ ψ_right.
//! ```
//!
//! Layer-specific construction is implemented in the private [`layer`] module.

use nalgebra::ComplexField;
use ndarray::{ArrayBase, Dimension, IntoDimension, OwnedRepr};
use num_traits::{One, Zero};

use crate::{
    ComplexScalar, IncidentSide, PlaneWaveAmplitudes, PlaneWaveDeterminant, PlaneWavePower,
    Polarisation,
    algebra::{Jet, RealScalarAlgebra, ScalarAlgebra},
    backend::{
        ClosedSide, ExteriorContextProvider, ExteriorWavevectors, PlaneWaveEntries,
        PlaneWaveModeCandidate, TerminationError, exterior_parts,
        isotropic::IsotropicLayerQuantities,
        transfer2::{error::Transfer2Entry, state::right_outgoing_transfer_state},
        undefined_amplitudes,
    },
    input::{CanonicalCoordinates, CanonicalStack},
    material::{ConstitutiveEvaluator, ConstitutiveLift},
    observable::{BoundaryState, ProjectAmplitudes, ProjectPlaneWaveModeDeterminant, ProjectPower},
};

use super::{
    projection::{right_incoming_column, right_outgoing_column},
    state::{TransferState, transfer_state_slope},
};

/// Internal algebraic representation of a 2×2 transfer matrix.
///
/// `A` determines the sampled shape, derivative order, and parameter policy
/// propagated through each matrix entry. All four entries share the same
/// algebraic representation.
#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub struct Transfer2Entries<A> {
    pub(super) m11: A,
    pub(super) m12: A,
    pub(super) m21: A,
    pub(super) m22: A,
}

impl<A> Transfer2Entries<A> {
    /// Construct an internal matrix from its four entries.
    pub(crate) fn new(m11: A, m12: A, m21: A, m22: A) -> Self {
        Self { m11, m12, m21, m22 }
    }

    /// Return entry `(1, 1)`.
    pub(crate) fn m11(&self) -> &A {
        &self.m11
    }

    /// Return entry `(1, 2)`.
    pub(crate) fn m12(&self) -> &A {
        &self.m12
    }

    /// Return entry `(2, 1)`.
    pub(crate) fn m21(&self) -> &A {
        &self.m21
    }

    /// Return entry `(2, 2)`.
    pub(crate) fn m22(&self) -> &A {
        &self.m22
    }

    // /// Consume the matrix and return its entries in row-major order.
    // pub(crate) fn into_parts(self) -> (A, A, A, A) {
    //     (self.m11, self.m12, self.m21, self.m22)
    // }

    pub(crate) fn first_non_finite(&self) -> Option<(Transfer2Entry, Vec<usize>)>
    where
        A: ScalarAlgebra,
        A::Scalar: ComplexField,
        A::Dimension: Dimension,
    {
        [
            (Transfer2Entry::M11, self.m11.value()),
            (Transfer2Entry::M12, self.m12.value()),
            (Transfer2Entry::M21, self.m21.value()),
            (Transfer2Entry::M22, self.m22.value()),
        ]
        .into_iter()
        .find_map(|(entry, values)| {
            values.indexed_iter().find_map(|(index, value)| {
                (!value.is_finite()).then(|| {
                    (
                        entry,
                        index.into_dimension().as_array_view().to_owned().to_vec(),
                    )
                })
            })
        })
    }

    // pub(super) fn sample_source(&self) -> &ArrayBase<OwnedRepr<A::Scalar>, A::Dimension>
    // where
    //     A: ScalarAlgebra,
    // {
    //     self.m11().value()
    // }
}

impl<A> Transfer2Entries<A> {
    pub(crate) fn identity_like(source: &ArrayBase<OwnedRepr<A::Scalar>, A::Dimension>) -> Self
    where
        A: ScalarAlgebra + Clone,
        A::Scalar: ComplexField + One + Zero,
    {
        let zero = A::filled_constant_like(source, <A::Scalar as Zero>::zero());
        let one = A::filled_constant_like(source, <A::Scalar as One>::one());

        Self::new(one.clone(), zero.clone(), zero, one)
    }

    /// Multiply two matrices using the supplied scalar algebra.
    ///
    /// This computes `self * rhs`:
    ///
    /// ```text
    /// [a b] [e f]   [ae + bg  af + bh]
    /// [c d] [g h] = [ce + dg  cf + dh]
    /// ```
    ///
    /// When `A` is a jet representation, all matrix-product derivatives are
    /// generated by the scalar product and sum rules.
    pub(crate) fn multiply(&self, rhs: &Self) -> Self
    where
        A: ScalarAlgebra,
        A::Scalar: ComplexField,
        A::Dimension: Dimension,
    {
        let m11 = self
            .m11()
            .multiply(rhs.m11())
            .add(&self.m12().multiply(rhs.m21()));

        let m12 = self
            .m11()
            .multiply(rhs.m12())
            .add(&self.m12().multiply(rhs.m22()));

        let m21 = self
            .m21()
            .multiply(rhs.m11())
            .add(&self.m22().multiply(rhs.m21()));

        let m22 = self
            .m21()
            .multiply(rhs.m12())
            .add(&self.m22().multiply(rhs.m22()));

        Self::new(m11, m12, m21, m22)
    }

    pub(crate) fn from_layer(quantities: &IsotropicLayerQuantities<A>, thickness: &A) -> Self
    where
        A: ScalarAlgebra,
        A::Scalar: ComplexField + One,
        A::Dimension: Dimension,
    {
        // Entire functions of q avoid the removable interior cutoff and the
        // singular derivatives of sqrt(q). Do not use the cached kappa jet.
        let q = quantities.kappa_squared();
        let phase_squared = q.multiply(&thickness.square());
        let cos = phase_squared.cos_sqrt();
        let integrated_sinc = thickness.multiply(&phase_squared.sinc_sqrt());

        Self::new(
            cos.clone(),
            integrated_sinc.multiply(quantities.factor()).negate(),
            integrated_sinc.multiply(&q.divide(quantities.factor())),
            cos,
        )
    }
}

#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub struct Transfer2ExteriorContext<A> {
    left_admittance: A,
    right_admittance: A,
    left_kappa: A,
    right_kappa: A,
    left_epsilon: A,
    right_epsilon: A,
    left_mu: A,
    right_mu: A,
    vacuum_angular_wavenumber: A,
    parallel_angular_wavenumber: A,
    polarisation: Polarisation,

    left_closure: Option<ClosedSide<A>>,
    right_closure: Option<ClosedSide<A>>,
}

impl<A> ExteriorContextProvider for Transfer2ExteriorContext<A> {
    type Algebra = A;

    fn left_admittance(&self) -> &A {
        &self.left_admittance
    }

    fn right_admittance(&self) -> &A {
        &self.right_admittance
    }

    fn left_kappa(&self) -> &A {
        &self.left_kappa
    }

    fn right_kappa(&self) -> &A {
        &self.right_kappa
    }

    fn left_epsilon(&self) -> &Self::Algebra {
        &self.left_epsilon
    }

    fn right_epsilon(&self) -> &Self::Algebra {
        &self.right_epsilon
    }

    fn left_mu(&self) -> &Self::Algebra {
        &self.left_mu
    }

    fn right_mu(&self) -> &Self::Algebra {
        &self.right_mu
    }

    fn vacuum_angular_wavenumber(&self) -> &Self::Algebra {
        &self.vacuum_angular_wavenumber
    }

    fn parallel_angular_wavenumber(&self) -> &Self::Algebra {
        &self.parallel_angular_wavenumber
    }

    fn polarisation(&self) -> Polarisation {
        self.polarisation
    }

    fn left_closure(&self) -> Option<&ClosedSide<A>> {
        self.left_closure.as_ref()
    }

    fn right_closure(&self) -> Option<&ClosedSide<A>> {
        self.right_closure.as_ref()
    }
}

impl<J> Transfer2ExteriorContext<J> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_parts(
        left_admittance: J,
        right_admittance: J,
        left_kappa: J,
        right_kappa: J,
        left_epsilon: J,
        right_epsilon: J,
        left_mu: J,
        right_mu: J,
        vacuum_angular_wavenumber: J,
        parallel_angular_wavenumber: J,
        polarisation: Polarisation,
        left_closure: Option<ClosedSide<J>>,
        right_closure: Option<ClosedSide<J>>,
    ) -> Self {
        Self {
            left_admittance,
            right_admittance,
            left_kappa,
            right_kappa,
            left_epsilon,
            right_epsilon,
            left_mu,
            right_mu,
            vacuum_angular_wavenumber,
            parallel_angular_wavenumber,
            polarisation,
            left_closure,
            right_closure,
        }
    }

    #[cfg(test)]
    pub(super) fn new<E, M>(
        coordinates: &CanonicalCoordinates<J>,
        left_exterior: &M,
        right_exterior: &M,
        exterior: &ExteriorWavevectors<J>,
        polarisation: Polarisation,
    ) -> Self
    where
        J: ScalarAlgebra + ConstitutiveLift<E, M> + Clone,
        J::Scalar: ComplexScalar,
        J::Dimension: Dimension,
        E: ConstitutiveEvaluator<J::Scalar, J::Dimension, M>,
        M: Clone,
    {
        let stack = CanonicalStack::new(left_exterior.clone(), right_exterior.clone(), Vec::new());

        Self::for_stack::<E, M>(coordinates, &stack, exterior, polarisation)
            .expect("an open stack takes both exterior wavevectors")
    }

    /// Construct the context for `stack`, whose open sides take their
    /// wavevectors from `exterior` and whose closed sides take theirs from
    /// the termination.
    pub(crate) fn for_stack<E, M>(
        coordinates: &CanonicalCoordinates<J>,
        stack: &CanonicalStack<M, J>,
        exterior: &ExteriorWavevectors<J>,
        polarisation: Polarisation,
    ) -> Result<Self, TerminationError>
    where
        J: ScalarAlgebra + ConstitutiveLift<E, M> + Clone,
        J::Scalar: ComplexScalar,
        J::Dimension: Dimension,
        E: ConstitutiveEvaluator<J::Scalar, J::Dimension, M>,
    {
        let parts = exterior_parts::<E, M, J>(coordinates, stack, exterior, polarisation)?;

        Ok(Self::from_parts(
            parts.left_admittance,
            parts.right_admittance,
            parts.left_kappa,
            parts.right_kappa,
            parts.left_epsilon,
            parts.right_epsilon,
            parts.left_mu,
            parts.right_mu,
            coordinates.vacuum_angular_wavenumber().clone(),
            coordinates.parallel_angular_wavenumber().clone(),
            polarisation,
            parts.left_closure,
            parts.right_closure,
        ))
    }
}

impl<A> PlaneWaveEntries for Transfer2Entries<A> {
    type ExteriorContext = Transfer2ExteriorContext<A>;
    type Algebra = A;
}

impl<J> ProjectAmplitudes for Transfer2Entries<J>
where
    J: ScalarAlgebra + Clone,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    type Amplitudes = PlaneWaveAmplitudes<J>;

    fn project_amplitudes(
        &self,
        exterior: &Self::ExteriorContext,
        incident_side: IncidentSide,
    ) -> Self::Amplitudes {
        if exterior.has_closure() {
            return closed_amplitudes(self, exterior, incident_side);
        }

        let left_slope = transfer_state_slope(exterior.left_admittance());

        let right_slope = transfer_state_slope(exterior.right_admittance());

        let (right_outgoing_field, right_outgoing_slope) =
            right_outgoing_column(self, &right_slope);

        let denominator = left_slope
            .multiply(&right_outgoing_field)
            .subtract(&right_outgoing_slope);

        let two = J::filled_constant_like(
            denominator.value(),
            <J::Scalar as One>::one() + <J::Scalar as One>::one(),
        );

        match incident_side {
            IncidentSide::Left => {
                /*
                 * At the left boundary:
                 *
                 *   field = 1 + r
                 *   slope = ξL (r - 1)
                 *
                 * At the right boundary only the transmitted right-going wave
                 * remains. This gives:
                 *
                 *   t = 2 ξL / D
                 *   r = t p - 1.
                 */
                let transmission = two.multiply(&left_slope).divide(&denominator);

                let reflection = transmission
                    .multiply(&right_outgoing_field)
                    .subtract(&transmission.constant(<J::Scalar as One>::one()));

                PlaneWaveAmplitudes::new(reflection, transmission)
            }

            IncidentSide::Right => {
                /*
                 * At the right boundary the incident left-going basis state is
                 * [1, +ξR], while the reflected right-going basis state is
                 * [1, -ξR].
                 */
                let (right_incoming_field, right_incoming_slope) =
                    right_incoming_column(self, &right_slope);

                let reflection = right_incoming_slope
                    .subtract(&left_slope.multiply(&right_incoming_field))
                    .divide(&denominator);

                /*
                 * Compute transmission from the propagated field rather than
                 * assuming det(M) = 1:
                 *
                 *   t = a + r p.
                 *
                 * This remains correct for any transfer representation using
                 * the documented state convention.
                 */
                let transmission =
                    right_incoming_field.add(&reflection.multiply(&right_outgoing_field));

                PlaneWaveAmplitudes::new(reflection, transmission)
            }
        }
    }
}

/// Plane-wave amplitudes of a stack with a closed side.
///
/// Incidence comes from the open side; nothing is transmitted, and the
/// reflection follows from the closed side's boundary condition. Incidence
/// from a closed side has no meaning and gives non-finite amplitudes; the
/// evaluators refuse it before projecting.
fn closed_amplitudes<J>(
    entries: &Transfer2Entries<J>,
    exterior: &Transfer2ExteriorContext<J>,
    incident_side: IncidentSide,
) -> PlaneWaveAmplitudes<J>
where
    J: ScalarAlgebra + Clone,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    let source = entries.m11().value();
    let zero = J::filled_constant_like(source, <J::Scalar as Zero>::zero());
    let one = J::filled_constant_like(source, <J::Scalar as One>::one());

    match (
        incident_side,
        exterior.left_closure(),
        exterior.right_closure(),
    ) {
        (IncidentSide::Left, None, Some(right)) => {
            /*
             * The right boundary state is fixed by the closure up to an
             * amplitude a: at the left boundary
             *
             *   a p = 1 + r,    a q = ξL (r − 1).
             */
            let left_slope = transfer_state_slope(exterior.left_admittance());
            let (field, slope) = right.right_state();
            let left = entries.apply_state(&TransferState::new(field, slope));
            let two = one.add(&one);
            let amplitude = two
                .multiply(&left_slope)
                .divide(&left_slope.multiply(left.field()).subtract(left.slope()));
            let reflection = amplitude.multiply(left.field()).subtract(&one);

            PlaneWaveAmplitudes::new(reflection, zero)
        }

        (IncidentSide::Right, Some(left), None) => {
            let right_slope = transfer_state_slope(exterior.right_admittance());
            let (incoming_field, incoming_slope) = right_incoming_column(entries, &right_slope);
            let (outgoing_field, outgoing_slope) = right_outgoing_column(entries, &right_slope);
            let (field_weight, slope_weight) = left.left_row();
            let incoming = field_weight
                .multiply(&incoming_field)
                .add(&slope_weight.multiply(&incoming_slope));
            let outgoing = field_weight
                .multiply(&outgoing_field)
                .add(&slope_weight.multiply(&outgoing_slope));

            PlaneWaveAmplitudes::new(incoming.divide(&outgoing).negate(), zero)
        }

        _ => undefined_amplitudes(source),
    }
}

impl<J> ProjectPower for Transfer2Entries<J>
where
    J: Clone + RealScalarAlgebra,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
    J::RealJet: ScalarAlgebra,
    <J::RealJet as Jet>::Scalar: One,
{
    type Power = PlaneWavePower<J::RealJet>;

    fn project_power(
        &self,
        exterior: &Self::ExteriorContext,
        incident_side: IncidentSide,
    ) -> Self::Power {
        let (reflection, transmission) = self
            .project_amplitudes(exterior, incident_side)
            .into_parts();

        let (incident_admittance, transmitted_admittance) = match incident_side {
            IncidentSide::Left => (exterior.left_admittance(), exterior.right_admittance()),

            IncidentSide::Right => (exterior.right_admittance(), exterior.left_admittance()),
        };

        PlaneWavePower::from_amplitudes_and_admittance(
            &reflection,
            &transmission,
            incident_admittance,
            transmitted_admittance,
        )
    }
}

impl<J> ProjectPlaneWaveModeDeterminant for Transfer2Entries<J>
where
    J: ScalarAlgebra + Clone,
    J::Scalar: ComplexScalar + One,
    J::Dimension: Dimension,
{
    type Determinant = PlaneWaveDeterminant<J>;

    fn project_determinant(&self, exterior: &Self::ExteriorContext) -> Self::Determinant {
        let candidate = right_gauged_mode_candidate(self, exterior);

        PlaneWaveDeterminant::new(candidate.into_residual())
    }
}

/// The outgoing-mode candidate gauged by a unit right boundary state.
///
/// With `ξ = −iY`, the characteristic is `D = ℓ_L · M · r_R`: an open right
/// side contributes `r_R = [1, −ξ_R]ᵀ` (a unit outgoing wave) and an open left
/// side `ℓ_L = [ξ_L, −1]`; a closed side contributes the wall condition
/// carried through its gap (see [`ClosedSide`]).
pub(crate) fn right_gauged_mode_candidate<A>(
    entries: &Transfer2Entries<A>,
    context: &Transfer2ExteriorContext<A>,
) -> PlaneWaveModeCandidate<A>
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar + One,
    A::Dimension: Dimension,
{
    let one = A::filled_constant_like(entries.m11().value(), A::Scalar::one());

    let (right_state, right_closed_state) = match context.right_closure() {
        None => (
            right_outgoing_transfer_state(&one, context.right_admittance()),
            None,
        ),
        Some(closure) => {
            let (field, slope) = closure.right_state();
            (
                TransferState::new(field.clone(), slope.clone()),
                Some(BoundaryState::new(field, slope)),
            )
        }
    };

    let left_state = entries.apply_state(&right_state);

    let residual = match context.left_closure() {
        None => {
            let left_slope = transfer_state_slope(context.left_admittance());

            left_slope
                .multiply(left_state.field())
                .subtract(left_state.slope())
        }
        Some(closure) => {
            let (field_weight, slope_weight) = closure.left_row();

            field_weight
                .multiply(left_state.field())
                .add(&slope_weight.multiply(left_state.slope()))
        }
    };

    let (field, secondary) = left_state.into_parts();

    PlaneWaveModeCandidate::new(BoundaryState::new(field, secondary), one, residual)
        .with_right_closed_state(right_closed_state)
}

#[cfg(test)]
mod tests {
    use ndarray::{Ix0, arr0};

    use super::*;
    use crate::{
        algebra::{ArrayJet0, RealParameter},
        input::{CanonicalCoordinates, Polarisation},
        test_support::{
            C, TOLERANCE,
            assertions::assert_complex_close,
            c,
            jet::real_j0,
            materials::{constant, linear},
        },
    };

    type A = ArrayJet0<C, Ix0, RealParameter>;

    fn entries(m11: C, m12: C, m21: C, m22: C) -> Transfer2Entries<A> {
        Transfer2Entries::new(real_j0(m11), real_j0(m12), real_j0(m21), real_j0(m22))
    }

    fn coordinates(k0: f64, k_parallel: f64) -> CanonicalCoordinates<A> {
        CanonicalCoordinates::new(real_j0(c(k0)), real_j0(c(k_parallel)))
    }

    #[test]
    fn identity_has_unit_diagonal_and_zero_off_diagonal() {
        let identity = Transfer2Entries::<A>::identity_like(&arr0(c(0.0)));

        assert_complex_close(identity.m11()[()], c(1.0), TOLERANCE);
        assert_complex_close(identity.m12()[()], c(0.0), TOLERANCE);
        assert_complex_close(identity.m21()[()], c(0.0), TOLERANCE);
        assert_complex_close(identity.m22()[()], c(1.0), TOLERANCE);
    }

    #[test]
    fn cutoff_layer_has_analytic_spectral_and_thickness_hessian() {
        use crate::{SeedJet, test_support::jet::HoloJB2};
        type J = HoloJB2;
        let d = 0.3_f64;
        let coordinates = CanonicalCoordinates::new(
            <J as SeedJet>::variable(arr0(c(1.0)), 0).unwrap(),
            <J as SeedJet>::constant(arr0(c(2.0))),
        );
        let thickness = <J as SeedJet>::variable(arr0(c(d)), 1).unwrap();
        for (pol, f) in [
            (Polarisation::TransverseElectric, 1.0),
            (Polarisation::TransverseMagnetic, 4.0),
        ] {
            let layer =
                IsotropicLayerQuantities::complex_plane(&constant(4.0, 1.0), &coordinates, pol);
            let m = Transfer2Entries::from_layer(&layer, &thickness);
            // Independent coefficients of cos(d sqrt(q)), sin(d sqrt(q))/sqrt(q),
            // q=4*k0^2-4 at k0=1. Axes are k0 and thickness.
            let expected = [
                (
                    m.m11(),
                    [
                        1.0,
                        -4.0 * d * d,
                        0.0,
                        16.0 / 3.0 * d.powi(4) - 4.0 * d * d,
                        -8.0 * d,
                        0.0,
                    ],
                ),
                (
                    m.m12(),
                    [
                        -f * d,
                        4.0 / 3.0 * f * d.powi(3),
                        -f,
                        f * (4.0 / 3.0 * d.powi(3) - 16.0 / 15.0 * d.powi(5)),
                        4.0 * f * d * d,
                        0.0,
                    ],
                ),
                (
                    m.m21(),
                    [
                        0.0,
                        8.0 * d / f,
                        0.0,
                        (8.0 * d - 64.0 / 3.0 * d.powi(3)) / f,
                        8.0 / f,
                        0.0,
                    ],
                ),
                (
                    m.m22(),
                    [
                        1.0,
                        -4.0 * d * d,
                        0.0,
                        16.0 / 3.0 * d.powi(4) - 4.0 * d * d,
                        -8.0 * d,
                        0.0,
                    ],
                ),
            ];
            for (entry, expected) in expected {
                assert!(entry.all_finite());
                let actual = [
                    entry.value()[()],
                    entry.axis0()[()],
                    entry.axis1()[()],
                    entry.axis0_axis0()[()],
                    entry.axis0_axis1()[()],
                    entry.axis1_axis1()[()],
                ];
                for (actual, expected) in actual.into_iter().zip(expected) {
                    assert_complex_close(actual, c(expected), 1e-12);
                }
            }
        }
    }

    #[test]
    fn identity_is_neutral_on_both_sides() {
        let matrix = entries(c(1.0), c(2.0), c(3.0), c(4.0));

        let identity = Transfer2Entries::<A>::identity_like(&arr0(c(0.0)));

        assert_eq!(identity.multiply(&matrix), matrix,);

        assert_eq!(matrix.multiply(&identity), matrix,);
    }

    #[test]
    fn multiplication_uses_standard_matrix_order() {
        let left = entries(c(1.0), c(2.0), c(3.0), c(4.0));

        let right = entries(c(5.0), c(6.0), c(7.0), c(8.0));

        let product = left.multiply(&right);

        assert_complex_close(product.m11()[()], c(1.0 * 5.0 + 2.0 * 7.0), TOLERANCE);
        assert_complex_close(product.m12()[()], c(1.0 * 6.0 + 2.0 * 8.0), TOLERANCE);
        assert_complex_close(product.m21()[()], c(3.0 * 5.0 + 4.0 * 7.0), TOLERANCE);
        assert_complex_close(product.m22()[()], c(3.0 * 6.0 + 4.0 * 8.0), TOLERANCE);
    }

    #[test]
    fn multiplication_is_not_accidentally_reversed() {
        let left = entries(c(1.0), c(2.0), c(0.0), c(1.0));

        let right = entries(c(1.0), c(0.0), c(3.0), c(1.0));

        assert_ne!(left.multiply(&right), right.multiply(&left),);
    }

    #[test]
    fn zero_thickness_layer_is_identity() {
        let material = constant(4.0, 1.0);
        let coordinates = coordinates(2.0, 0.3);

        let quantities = IsotropicLayerQuantities::real_axis(
            &material,
            &coordinates,
            Polarisation::TransverseElectric,
        );

        let matrix = Transfer2Entries::from_layer(&quantities, &real_j0(c(0.0)));

        let identity = Transfer2Entries::<A>::identity_like(&arr0(c(0.0)));

        assert_eq!(matrix, identity);
    }

    #[test]
    fn homogeneous_layer_matches_analytic_entries() {
        let material = constant(4.0, 1.0);
        let coordinates = coordinates(2.0, 0.0);

        let quantities = IsotropicLayerQuantities::real_axis(
            &material,
            &coordinates,
            Polarisation::TransverseElectric,
        );

        let thickness = 0.2;

        let matrix = Transfer2Entries::from_layer(&quantities, &real_j0(c(thickness)));

        let kappa = 4.0;
        let admittance = 4.0;
        let theta = kappa * thickness;

        assert_complex_close(matrix.m11()[()], c(theta.cos()), TOLERANCE);

        assert_complex_close(matrix.m12()[()], c(-theta.sin() / admittance), TOLERANCE);

        assert_complex_close(matrix.m21()[()], c(admittance * theta.sin()), TOLERANCE);

        assert_complex_close(matrix.m22()[()], c(theta.cos()), TOLERANCE);
    }

    #[test]
    fn te_and_tm_use_their_respective_admittances() {
        let material = constant(4.0, 2.0);
        let coordinates = coordinates(3.0, 1.0);
        let thickness = real_j0(c(0.1));

        let te = Transfer2Entries::from_layer(
            &IsotropicLayerQuantities::real_axis(
                &material,
                &coordinates,
                Polarisation::TransverseElectric,
            ),
            &thickness,
        );

        let tm = Transfer2Entries::from_layer(
            &IsotropicLayerQuantities::real_axis(
                &material,
                &coordinates,
                Polarisation::TransverseMagnetic,
            ),
            &thickness,
        );

        assert_ne!(te.m12(), tm.m12());
        assert_ne!(te.m21(), tm.m21());
    }

    #[test]
    fn dispersive_material_changes_layer_matrix() {
        let material = linear(2.0, 0.4, 1.0, 0.1);

        let first = Transfer2Entries::from_layer(
            &IsotropicLayerQuantities::real_axis(
                &material,
                &coordinates(2.0, 0.2),
                Polarisation::TransverseElectric,
            ),
            &real_j0(c(0.1)),
        );

        let second = Transfer2Entries::from_layer(
            &IsotropicLayerQuantities::real_axis(
                &material,
                &coordinates(3.0, 0.2),
                Polarisation::TransverseElectric,
            ),
            &real_j0(c(0.1)),
        );

        assert_ne!(first, second);
    }

    #[test]
    fn detects_non_finite_entry_and_index() {
        let matrix = Transfer2Entries::new(
            real_j0(c(1.0)),
            real_j0(C::new(f64::INFINITY, 0.0)),
            real_j0(c(0.0)),
            real_j0(c(1.0)),
        );

        assert_eq!(
            matrix.first_non_finite(),
            Some((Transfer2Entry::M12, Vec::new(),)),
        );
    }
}

#[cfg(test)]
mod projection_tests {
    use approx::assert_relative_eq;
    use ndarray::{Ix0, arr0};
    use num_complex::Complex64;

    use super::*;

    use crate::{
        RealAxis,
        algebra::{ArrayJet0, Jet0, RealParameter},
        backend::IsotropicMediumQuantities,
        input::{CanonicalCoordinates, IncidentSide},
        test_support::materials::constant,
    };

    type C = Complex64;
    type J = ArrayJet0<C, Ix0, RealParameter>;

    const TOLERANCE: f64 = 1.0e-12;

    fn c(real: f64) -> C {
        C::new(real, 0.0)
    }

    fn jet(value: C) -> J {
        Jet0::new(arr0(value))
    }

    fn assert_close(actual: C, expected: C) {
        assert_relative_eq!(
            actual.re,
            expected.re,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );

        assert_relative_eq!(
            actual.im,
            expected.im,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    fn coordinates() -> CanonicalCoordinates<J> {
        CanonicalCoordinates::new(jet(c(2.0)), jet(c(0.0)))
    }

    fn exterior_context(
        left_epsilon: f64,
        right_epsilon: f64,
        polarisation: Polarisation,
    ) -> Transfer2ExteriorContext<J> {
        let left_exterior = constant(left_epsilon, 1.0);
        let right_exterior = constant(right_epsilon, 1.0);
        let coordinates = coordinates();

        let exterior = ExteriorWavevectors::new(
            IsotropicLayerQuantities::evaluate::<RealAxis, _>(
                &left_exterior,
                &coordinates,
                polarisation,
            )
            .kappa()
            .clone(),
            IsotropicLayerQuantities::evaluate::<RealAxis, _>(
                &right_exterior,
                &coordinates,
                polarisation,
            )
            .kappa()
            .clone(),
        );

        Transfer2ExteriorContext::new::<crate::domain::RealAxis, _>(
            &coordinates,
            &left_exterior,
            &right_exterior,
            &exterior,
            polarisation,
        )
    }

    fn identity() -> Transfer2Entries<J> {
        Transfer2Entries::identity_like(coordinates().vacuum_angular_wavenumber().value())
    }

    #[test]
    fn equal_exteriors_and_identity_give_unit_transmission() {
        let entries = identity();

        let exterior = exterior_context(1.0, 1.0, Polarisation::TransverseElectric);

        for side in [IncidentSide::Left, IncidentSide::Right] {
            let amplitudes = entries.project_amplitudes(&exterior, side);

            assert_close(amplitudes.reflection()[()], c(0.0));

            assert_close(amplitudes.transmission()[()], c(1.0));
        }
    }

    #[test]
    fn identity_with_different_exteriors_matches_interface_fresnel_amplitudes() {
        let entries = identity();

        let exterior = exterior_context(1.0, 4.0, Polarisation::TransverseElectric);

        let left = entries.project_amplitudes(&exterior, IncidentSide::Left);

        assert_close(left.reflection()[()], c(-1.0 / 3.0));

        assert_close(left.transmission()[()], c(2.0 / 3.0));

        let right = entries.project_amplitudes(&exterior, IncidentSide::Right);

        assert_close(right.reflection()[()], c(1.0 / 3.0));

        assert_close(right.transmission()[()], c(4.0 / 3.0));
    }

    #[test]
    fn identity_interface_conserves_power() {
        let entries = identity();

        let exterior = exterior_context(1.0, 4.0, Polarisation::TransverseElectric);

        for side in [IncidentSide::Left, IncidentSide::Right] {
            let power = entries.project_power(&exterior, side);

            let total =
                power.reflectance()[()] + power.transmittance()[()] + power.absorptance()[()];

            assert_relative_eq!(total, 1.0, epsilon = TOLERANCE, max_relative = TOLERANCE,);

            assert_relative_eq!(
                power.absorptance()[()],
                0.0,
                epsilon = TOLERANCE,
                max_relative = TOLERANCE,
            );
        }
    }

    #[test]
    fn determinant_projection_uses_both_exterior_admittances() {
        let entries = Transfer2Entries::new(jet(c(1.0)), jet(c(0.2)), jet(c(0.3)), jet(c(1.1)));

        let first = entries.project_determinant(&exterior_context(
            1.0,
            2.25,
            Polarisation::TransverseElectric,
        ));

        let second = entries.project_determinant(&exterior_context(
            1.0,
            4.0,
            Polarisation::TransverseElectric,
        ));

        assert_ne!(first.value(), second.value(),);
    }

    #[test]
    fn determinant_projection_uses_supplied_exterior_wavevector_branches() {
        let entries = Transfer2Entries::new(jet(c(1.0)), jet(c(0.2)), jet(c(0.3)), jet(c(1.1)));

        let left_exterior = constant(1.0, 1.0);
        let right_exterior = constant(2.25, 1.0);
        let coordinates = coordinates();

        let ordinary = ExteriorWavevectors::new(
            IsotropicMediumQuantities::real_axis(&left_exterior, &coordinates).into_kappa(),
            IsotropicMediumQuantities::real_axis(&right_exterior, &coordinates).into_kappa(),
        );

        let flipped = ExteriorWavevectors::new(
            ordinary.left().unwrap().negate(),
            ordinary.right().unwrap().clone(),
        );

        let ordinary_context = Transfer2ExteriorContext::new::<RealAxis, _>(
            &coordinates,
            &left_exterior,
            &right_exterior,
            &ordinary,
            Polarisation::TransverseElectric,
        );

        let flipped_context = Transfer2ExteriorContext::new::<RealAxis, _>(
            &coordinates,
            &left_exterior,
            &right_exterior,
            &flipped,
            Polarisation::TransverseElectric,
        );

        let ordinary = entries.project_determinant(&ordinary_context);

        let flipped = entries.project_determinant(&flipped_context);

        assert_ne!(ordinary.value(), flipped.value());
    }
}
