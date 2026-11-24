//! Projective scattering representation and modal projection.
//!
//! The ordinary scattering entries are represented by a common denominator
//! and four numerators. Redheffer composition is performed directly on this
//! homogeneous representation, avoiding division by intermediate scattering
//! denominators.
//!
//! Physical amplitudes and powers are obtained by returning to ordinary
//! scattering entries. Modal quantities are projected directly from the
//! projective representation so that poles of the ordinary scattering matrix
//! need not be materialised.

use ndarray::{ArrayBase, Dimension, OwnedRepr};
use num_traits::{One, Zero};

use crate::{
    ComplexScalar, IncidentSide, PlaneWaveAmplitudes, PlaneWaveDeterminant, PlaneWavePower,
    algebra::{Jet, RealScalarAlgebra, ScalarAlgebra},
    backend::{
        ExteriorContextProvider, PlaneWaveEntries, PlaneWaveModeCandidate,
        scatter2::{Scatter2ExteriorContext, entries::transfer_state_slope},
        undefined_amplitudes,
    },
    observable::{BoundaryState, ProjectAmplitudes, ProjectPlaneWaveModeDeterminant, ProjectPower},
};

use super::Scatter2Entries;

/// Projective representation of a two-port scattering matrix.
///
/// The represented scattering entries satisfy:
///
/// ```text
/// sᵢⱼ = nᵢⱼ / d.
/// ```
///
/// Multiplying `d` and every numerator by the same nonzero scalar leaves the
/// represented scattering matrix unchanged. The complementary determinant h
/// satisfies n12*n21-n11*n22=d*h and must scale with every other entry. Carry
/// h independently: reconstructing it by division at d=0 loses continuation.
#[derive(Clone, Debug, PartialEq)]
pub struct Scatter2ProjectiveEntries<A> {
    denominator: A,
    n11: A,
    n12: A,
    n21: A,
    n22: A,
    complement: A,
}

impl<A> PlaneWaveEntries for Scatter2ProjectiveEntries<A> {
    type ExteriorContext = Scatter2ExteriorContext<A>;
    type Algebra = A;
}

impl<A> Scatter2ProjectiveEntries<A> {
    pub(crate) const fn from_parts(
        denominator: A,
        n11: A,
        n12: A,
        n21: A,
        n22: A,
        complement: A,
    ) -> Self {
        Self {
            denominator,
            n11,
            n12,
            n21,
            n22,
            complement,
        }
    }

    pub(crate) fn complement(&self) -> &A {
        &self.complement
    }

    pub(crate) fn denominator(&self) -> &A {
        &self.denominator
    }

    pub(crate) fn n11(&self) -> &A {
        &self.n11
    }

    pub(crate) fn n12(&self) -> &A {
        &self.n12
    }

    pub(crate) fn n21(&self) -> &A {
        &self.n21
    }

    pub(crate) fn n22(&self) -> &A {
        &self.n22
    }

    pub(crate) fn sample_source(&self) -> &ArrayBase<OwnedRepr<A::Scalar>, A::Dimension>
    where
        A: ScalarAlgebra,
    {
        self.n11.value()
    }
}

impl<A> Scatter2ProjectiveEntries<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar + One + Zero,
    A::Dimension: Dimension,
{
    pub(crate) fn identity_like(
        source: &ndarray::ArrayBase<ndarray::OwnedRepr<A::Scalar>, A::Dimension>,
    ) -> Self {
        let zero = A::filled_constant_like(source, A::Scalar::zero());

        let one = A::filled_constant_like(source, A::Scalar::one());

        Self {
            denominator: one.clone(),
            n11: zero.clone(),
            n12: one.clone(),
            n21: one.clone(),
            n22: zero,
            complement: one,
        }
    }
}

impl<A> Scatter2ProjectiveEntries<A>
where
    A: ScalarAlgebra,
{
    pub(crate) fn entries(&self) -> Scatter2Entries<A> {
        Scatter2Entries::from_parts(
            self.n11.divide(&self.denominator),
            self.n12.divide(&self.denominator),
            self.n21.divide(&self.denominator),
            self.n22.divide(&self.denominator),
        )
    }
}

impl<A> Scatter2ProjectiveEntries<A>
where
    A: ScalarAlgebra,
    A::Scalar: One,
{
    pub(crate) fn from_entries(entries: &Scatter2Entries<A>) -> Self {
        let denominator = A::filled_constant_like(entries.sample_source(), A::Scalar::one());

        Self {
            denominator,
            n11: entries.s11().clone(),
            n12: entries.s12().clone(),
            n21: entries.s21().clone(),
            n22: entries.s22().clone(),
            complement: entries
                .s12()
                .multiply(entries.s21())
                .subtract(&entries.s11().multiply(entries.s22())),
        }
    }
}

impl<A> From<Scatter2Entries<A>> for Scatter2ProjectiveEntries<A>
where
    A: ScalarAlgebra,
    A::Scalar: One,
{
    fn from(entries: Scatter2Entries<A>) -> Self {
        Self::from_entries(&entries)
    }
}

impl<A> Scatter2ProjectiveEntries<A>
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar + One,
    A::Dimension: Dimension,
{
    pub(crate) fn right_gauged_mode_candidate(
        &self,
        context: &Scatter2ExteriorContext<A>,
    ) -> PlaneWaveModeCandidate<A> {
        let slope = transfer_state_slope(context.left_admittance());

        if context.has_closure() {
            return self.closed_mode_candidate(context, &slope);
        }

        /*
         * Projectively scaled left waves:
         *
         * forward  = denominator
         * backward = n11
         */
        let field = self.denominator().add(self.n11());

        let secondary = self.n11().subtract(self.denominator()).multiply(&slope);

        let two = A::Scalar::one() + A::Scalar::one();

        let residual = slope.multiply(self.denominator()).scale(two);

        let right_outgoing = self.n21().clone();

        PlaneWaveModeCandidate::new(
            BoundaryState::new(field, secondary),
            right_outgoing,
            residual,
        )
    }

    /// The left port's waves `(forward, backward)`, projectively scaled by
    /// the outgoing right-port amplitude `n21`: `(d, n11)`, or, with a closed
    /// right side that reflects the outgoing wave back with its load `ρ̂_R`,
    /// `(d − ρ̂_R n22, n11 + ρ̂_R h)`.
    pub(crate) fn left_port_waves(&self, right_load: Option<&A>) -> (A, A) {
        match right_load {
            Some(load) => self.loaded_port(load, self.n11(), self.n22()),
            None => (self.denominator().clone(), self.n11().clone()),
        }
    }

    /// `(d − ρ̂ far, near + ρ̂ h)`: the projective waves at one port when the
    /// other port carries the load `ρ̂`, with `near` and `far` the reflection
    /// numerators of the two ports.
    fn loaded_port(&self, load: &A, near: &A, far: &A) -> (A, A) {
        (
            self.denominator().subtract(&load.multiply(far)),
            near.add(&load.multiply(self.complement())),
        )
    }

    /// The left waves and the residual of a stack with a closed side. A
    /// closed left side requires `forward = ρ̂_L backward`.
    fn closed_left_waves(&self, context: &Scatter2ExteriorContext<A>, slope: &A) -> (A, A, A) {
        let right_load = context.right_closure().map(|closure| closure.load());

        let (forward, backward) = self.left_port_waves(right_load.as_ref());

        let residual = match context.left_closure() {
            Some(closure) => closure
                .coefficient()
                .multiply(&forward.subtract(&closure.load().multiply(&backward))),
            None => {
                let two = A::Scalar::one() + A::Scalar::one();
                slope.multiply(&forward).scale(two)
            }
        };

        (forward, backward, residual)
    }

    /// The candidate of a stack with a closed side, gauged by `n21`.
    fn closed_mode_candidate(
        &self,
        context: &Scatter2ExteriorContext<A>,
        slope: &A,
    ) -> PlaneWaveModeCandidate<A> {
        let (forward, backward, residual) = self.closed_left_waves(context, slope);

        let field = forward.add(&backward);

        let secondary = backward.subtract(&forward).multiply(slope);

        let right_outgoing = self.n21().clone();

        let right_closed_state = context.right_closure().map(|closure| {
            let load = closure.load();
            let right_slope = transfer_state_slope(context.right_admittance());
            let one = load.one_like();

            BoundaryState::new(
                right_outgoing.multiply(&one.add(&load)),
                right_slope
                    .multiply(&right_outgoing)
                    .multiply(&load.subtract(&one)),
            )
        });

        PlaneWaveModeCandidate::new(
            BoundaryState::new(field, secondary),
            right_outgoing,
            residual,
        )
        .with_right_closed_state(right_closed_state)
    }

    /// The characteristic of a stack with a closed side, in the transfer
    /// gauge `D = ℓ_L · M · r_R`:
    ///
    /// ```text
    /// D̂ = c_L c_R / n21 · exp(Σ e) · [d − ρ̂_L n11 − ρ̂_R n22 − ρ̂_L ρ̂_R h]
    /// ```
    ///
    /// where an open left side has `c_L = 2ξ_L` and an open right side
    /// `c_R = 1`, without load or exponent, and each closed side contributes its
    /// coefficient, its load `ρ̂ = ρ exp(2iκ̃ℓ)` and its exponent `e = −iκ̃ℓ`.
    /// `c_L` times the bracket is the candidate's residual.
    fn closed_determinant(&self, context: &Scatter2ExteriorContext<A>) -> A {
        let slope = transfer_state_slope(context.left_admittance());

        let (_, _, residual) = self.closed_left_waves(context, &slope);

        let mut value = residual.divide(self.n21());

        if let Some(closure) = context.right_closure() {
            value = value.multiply(&closure.coefficient());
        }

        let exponent = [context.left_closure(), context.right_closure()]
            .into_iter()
            .flatten()
            .map(|closure| closure.scale_exponent())
            .reduce(|left, right| left.add(&right))
            .expect("a closed side");

        value.multiply(&exponent.exp())
    }
}

pub(crate) fn cascade_projection<A>(
    left: &Scatter2ProjectiveEntries<A>,
    right: &Scatter2ProjectiveEntries<A>,
) -> Scatter2ProjectiveEntries<A>
where
    A: ScalarAlgebra,
{
    // Substitute n12*n21-n11*n22=d*h before evaluating: the discarded
    // dl/dr factors are removable at intermediate network poles.
    let denominator = left
        .denominator()
        .multiply(right.denominator())
        .subtract(&left.n22().multiply(right.n11()));
    let n11 = left
        .n11()
        .multiply(right.denominator())
        .add(&left.complement().multiply(right.n11()));
    let n12 = left.n12().multiply(right.n12());
    let n21 = left.n21().multiply(right.n21());
    let n22 = right
        .n22()
        .multiply(left.denominator())
        .add(&right.complement().multiply(left.n22()));
    let complement = left
        .complement()
        .multiply(right.complement())
        .subtract(&left.n11().multiply(right.n22()));
    Scatter2ProjectiveEntries::from_parts(denominator, n11, n12, n21, n22, complement)
}

impl<J> ProjectAmplitudes for Scatter2ProjectiveEntries<J>
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
            return self.closed_amplitudes(exterior, incident_side);
        }

        self.entries().project_amplitudes(exterior, incident_side)
    }
}

impl<J> ProjectPower for Scatter2ProjectiveEntries<J>
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
        if exterior.has_closure() {
            let (reflection, transmission) =
                self.closed_amplitudes(exterior, incident_side).into_parts();

            let admittance = match incident_side {
                IncidentSide::Left => exterior.left_admittance(),
                IncidentSide::Right => exterior.right_admittance(),
            };

            return PlaneWavePower::from_amplitudes_and_admittance(
                &reflection,
                &transmission,
                admittance,
                admittance,
            );
        }

        self.entries().project_power(exterior, incident_side)
    }
}

impl<J> Scatter2ProjectiveEntries<J>
where
    J: ScalarAlgebra + Clone,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    /// Plane-wave amplitudes of a stack with a closed side: the closed side
    /// is a load `ρ̂` at its outermost interface, so nothing is transmitted and
    /// the reflection from the open side is `(n11 + ρ̂_R h)/(d − ρ̂_R n22)`
    /// (left incidence) or `(n22 + ρ̂_L h)/(d − ρ̂_L n11)` (right incidence).
    /// Incidence from a closed side gives non-finite amplitudes; the
    /// evaluators refuse it before projecting.
    fn closed_amplitudes(
        &self,
        exterior: &Scatter2ExteriorContext<J>,
        incident_side: IncidentSide,
    ) -> PlaneWaveAmplitudes<J> {
        let zero = self.denominator().zero_like();

        match (
            incident_side,
            exterior.left_closure(),
            exterior.right_closure(),
        ) {
            (IncidentSide::Left, None, Some(right)) => {
                let (forward, backward) = self.left_port_waves(Some(&right.load()));

                PlaneWaveAmplitudes::new(backward.divide(&forward), zero)
            }
            (IncidentSide::Right, Some(left), None) => {
                let (forward, backward) = self.loaded_port(&left.load(), self.n22(), self.n11());

                PlaneWaveAmplitudes::new(backward.divide(&forward), zero)
            }
            _ => undefined_amplitudes(self.sample_source()),
        }
    }
}

/// The projection uses `n21` as its transmission chart and is therefore
/// singular where `n21 = 0`. Such points are distinct from poles represented
/// by `denominator = 0`.
impl<J> ProjectPlaneWaveModeDeterminant for Scatter2ProjectiveEntries<J>
where
    J: ScalarAlgebra + Clone,
    J::Scalar: ComplexScalar + One,
    J::Dimension: Dimension,
{
    type Determinant = PlaneWaveDeterminant<J>;

    fn project_determinant(&self, exterior: &Self::ExteriorContext) -> Self::Determinant {
        if exterior.has_closure() {
            return PlaneWaveDeterminant::new(self.closed_determinant(exterior));
        }

        let left_slope = transfer_state_slope(exterior.left_admittance());

        let two = J::Scalar::one() + J::Scalar::one();

        let value = left_slope
            .multiply(self.denominator())
            .scale(two)
            .divide(self.n21());

        PlaneWaveDeterminant::new(value)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    use crate::{
        Polarisation, RealAxis,
        algebra::ScalarAlgebra,
        backend::{
            ExteriorWavevectors, IsotropicLayerQuantities, scatter2::Scatter2ExteriorContext,
        },
        input::CanonicalCoordinates,
        material::Constant,
        observable::BoundaryState,
        test_support::{
            C, TOLERANCE,
            assertions::{assert_array_close, assert_complex_close},
            c,
            jet::{RealJ0, real_j0},
        },
    };

    type Projection = Scatter2ProjectiveEntries<RealJ0>;

    fn jet(value: C) -> RealJ0 {
        real_j0(value)
    }

    fn scalar(value: &RealJ0) -> C {
        value.value()[()]
    }

    fn projection(denominator: C, n11: C, n12: C, n21: C, n22: C) -> Projection {
        Scatter2ProjectiveEntries::from_parts(
            jet(denominator),
            jet(n11),
            jet(n12),
            jet(n21),
            jet(n22),
            jet(if denominator == c(0.0) {
                c(0.0)
            } else {
                (n12 * n21 - n11 * n22) / denominator
            }),
        )
    }

    fn make_context() -> Scatter2ExteriorContext<RealJ0> {
        let source = jet(c(1.7));

        let coordinates = CanonicalCoordinates::new(source.clone(), source);

        Scatter2ExteriorContext::new::<RealAxis, _>(
            &coordinates,
            &Constant::vacuum(),
            &Constant::vacuum(),
            &ExteriorWavevectors::new(
                IsotropicLayerQuantities::evaluate::<RealAxis, _>(
                    &Constant::vacuum(),
                    &coordinates,
                    Polarisation::TransverseElectric,
                )
                .kappa()
                .clone(),
                IsotropicLayerQuantities::evaluate::<RealAxis, _>(
                    &Constant::vacuum(),
                    &coordinates,
                    Polarisation::TransverseElectric,
                )
                .kappa()
                .clone(),
            ),
            Polarisation::TransverseElectric,
        )
    }

    fn scale_projection(data: &Projection, scale: C) -> Projection {
        Scatter2ProjectiveEntries::from_parts(
            data.denominator().scale(scale),
            data.n11().scale(scale),
            data.n12().scale(scale),
            data.n21().scale(scale),
            data.n22().scale(scale),
            data.complement().scale(scale),
        )
    }

    fn assert_state_scaled(
        actual: &BoundaryState<RealJ0>,
        expected: &BoundaryState<RealJ0>,
        scale: C,
    ) {
        assert_complex_close(
            scalar(actual.field()),
            scale * scalar(expected.field()),
            TOLERANCE,
        );

        assert_complex_close(
            scalar(actual.secondary()),
            scale * scalar(expected.secondary()),
            TOLERANCE,
        );
    }

    #[test]
    fn entries_divides_every_numerator_by_common_denominator() {
        let data = projection(c(2.0), c(3.0), c(5.0), c(7.0), c(11.0));

        let entries = data.entries();

        assert_complex_close(scalar(entries.s11()), c(1.5), TOLERANCE);

        assert_complex_close(scalar(entries.s12()), c(2.5), TOLERANCE);

        assert_complex_close(scalar(entries.s21()), c(3.5), TOLERANCE);

        assert_complex_close(scalar(entries.s22()), c(5.5), TOLERANCE);
    }

    #[test]
    fn intermediate_pole_has_finite_analytic_continuation() {
        use crate::test_support::jet::HoloJ2;
        use ndarray::arr0;
        type J = HoloJ2;
        let constant = |x: f64| J::constant(arr0(c(x)));
        let right = Scatter2ProjectiveEntries::from_entries(&Scatter2Entries::from_parts(
            constant(0.2),
            constant(0.7),
            constant(0.8),
            constant(0.1),
        ));
        for delta in [0.0, 1e-7, -1e-7] {
            // h=0 is supplied by the rank-one numerator's parameterization,
            // not reconstructed through a 0/0 division at delta=0.
            let left = Scatter2ProjectiveEntries::from_parts(
                J::variable(arr0(c(delta))),
                constant(1.0),
                constant(2.0),
                constant(3.0),
                constant(6.0),
                constant(0.0),
            );
            let combined = cascade_projection(&left, &right);
            let entries = combined.entries();
            for (entry, numerator, offset) in [
                (entries.s11(), 1.0, 0.0),
                (entries.s12(), 1.4, 0.0),
                (entries.s21(), 2.4, 0.0),
                (entries.s22(), 3.36, 0.1),
            ] {
                let denominator = delta - 1.2;
                assert!(entry.all_finite());
                assert_complex_close(
                    entry.value()[()],
                    c(offset + numerator / denominator),
                    1e-12,
                );
                assert_complex_close(
                    entry.first()[()],
                    c(-numerator / denominator.powi(2)),
                    1e-12,
                );
                assert_complex_close(
                    entry.second()[()],
                    c(2.0 * numerator / denominator.powi(3)),
                    1e-12,
                );
            }
            assert!(combined.denominator().value()[()].norm() > 1.0);
        }
    }

    #[test]
    fn common_projective_scaling_does_not_change_into_scattering_entries() {
        let data = projection(c(2.0), c(3.0), c(5.0), c(7.0), c(11.0));

        let scale = c(1.7) + C::i() * c(-0.4);

        let scaled = scale_projection(&data, scale);

        let actual = scaled.entries();
        let expected = data.entries();

        assert_array_close(actual.s11().value(), expected.s11().value(), TOLERANCE);

        assert_array_close(actual.s12().value(), expected.s12().value(), TOLERANCE);

        assert_array_close(actual.s21().value(), expected.s21().value(), TOLERANCE);

        assert_array_close(actual.s22().value(), expected.s22().value(), TOLERANCE);
    }

    #[test]
    fn determinant_is_invariant_under_common_projective_scaling() {
        let data = projection(c(2.0), c(3.0), c(5.0), c(7.0), c(11.0));

        let scale = c(1.7) + C::i() * c(-0.4);

        let scaled = scale_projection(&data, scale);

        let context = make_context();

        let actual = scaled.project_determinant(&context);
        let expected = data.project_determinant(&context);

        assert_array_close(actual.value().value(), expected.value().value(), TOLERANCE);
    }

    #[test]
    fn projective_determinant_matches_inverse_transmission_characteristic() {
        let data = projection(c(2.0), c(3.0), c(5.0), c(7.0), c(11.0));

        let context = make_context();

        let determinant = data.project_determinant(&context);

        let entries = data.entries();

        let slope = transfer_state_slope(context.left_admittance());

        let two = C::one() + C::one();

        let expected = slope.scale(two).divide(entries.s21());

        assert_array_close(determinant.value().value(), expected.value(), TOLERANCE);
    }

    #[test]
    fn right_gauged_candidate_uses_denominator_and_n11_as_left_waves() {
        let data = projection(c(2.0), c(3.0), c(5.0), c(7.0), c(11.0));

        let context = make_context();

        let candidate = data.right_gauged_mode_candidate(&context);

        let slope = transfer_state_slope(context.left_admittance());

        let expected_field = data.denominator().add(data.n11());

        let expected_secondary = data.n11().subtract(data.denominator()).multiply(&slope);

        assert_array_close(
            candidate.state().field().value(),
            expected_field.value(),
            TOLERANCE,
        );

        assert_array_close(
            candidate.state().secondary().value(),
            expected_secondary.value(),
            TOLERANCE,
        );
    }

    #[test]
    fn candidate_residual_is_left_outgoing_boundary_residual() {
        let data = projection(c(2.0), c(3.0), c(5.0), c(7.0), c(11.0));

        let context = make_context();

        let candidate = data.right_gauged_mode_candidate(&context);

        let slope = transfer_state_slope(context.left_admittance());

        let expected = slope
            .multiply(candidate.state().field())
            .subtract(candidate.state().secondary());

        assert_array_close(candidate.residual().value(), expected.value(), TOLERANCE);
    }

    #[test]
    fn chart_normalized_candidate_residual_equals_determinant() {
        let data = projection(c(2.0), c(3.0), c(5.0), c(7.0), c(11.0));

        let context = make_context();

        let candidate = data.right_gauged_mode_candidate(&context);

        let determinant = data.project_determinant(&context);

        let chart_normalized = candidate.residual().divide(data.n21());

        assert_array_close(
            chart_normalized.value(),
            determinant.value().value(),
            TOLERANCE,
        );
    }

    #[test]
    fn projective_candidate_scales_with_representation() {
        let data = projection(c(2.0), c(3.0), c(5.0), c(7.0), c(11.0));

        let scale = c(1.7) + C::i() * c(-0.4);

        let scaled = scale_projection(&data, scale);

        let context = make_context();

        let base = data.right_gauged_mode_candidate(&context);

        let scaled = scaled.right_gauged_mode_candidate(&context);

        assert_state_scaled(scaled.state(), base.state(), scale);

        assert_complex_close(
            scalar(scaled.residual()),
            scale * scalar(base.residual()),
            TOLERANCE,
        );
    }

    #[test]
    fn candidate_remains_finite_at_exact_projective_pole() {
        let data = projection(c(0.0), c(3.0), c(5.0), c(7.0), c(35.0 / 3.0));

        let context = make_context();

        let candidate = data.right_gauged_mode_candidate(&context);

        assert!(scalar(candidate.state().field()).re.is_finite(),);

        assert!(scalar(candidate.state().field()).im.is_finite(),);

        assert!(scalar(candidate.state().secondary()).re.is_finite(),);

        assert!(scalar(candidate.state().secondary()).im.is_finite(),);

        assert_complex_close(scalar(candidate.residual()), c(0.0), TOLERANCE);

        let determinant = data.project_determinant(&context);

        assert_complex_close(scalar(determinant.value()), c(0.0), TOLERANCE);
    }

    #[test]
    fn regularized_candidate_is_n21_times_unit_right_outgoing_candidate() {
        let data = projection(c(2.0), c(3.0), c(5.0), c(7.0), c(11.0));

        let context = make_context();

        let regularized = data.right_gauged_mode_candidate(&context);

        let entries = data.entries();

        let left_forward = entries.s21().reciprocal();

        let left_backward = entries.s11().divide(entries.s21());

        let slope = transfer_state_slope(context.left_admittance());

        let unit_state = BoundaryState::new(
            left_forward.add(&left_backward),
            left_backward.subtract(&left_forward).multiply(&slope),
        );

        assert_state_scaled(regularized.state(), &unit_state, scalar(data.n21()));
    }
}

#[cfg(test)]
mod cascade_tests {
    use super::*;

    use crate::{
        backend::scatter2::entries::cascade,
        test_support::{
            C, TOLERANCE,
            assertions::assert_array_close,
            c,
            jet::{RealJ0, real_j0},
        },
    };

    type Projection = Scatter2ProjectiveEntries<RealJ0>;
    type Entries = Scatter2Entries<RealJ0>;

    fn jet(value: C) -> RealJ0 {
        real_j0(value)
    }

    fn entries(s11: C, s12: C, s21: C, s22: C) -> Entries {
        Scatter2Entries::from_parts(jet(s11), jet(s12), jet(s21), jet(s22))
    }

    #[test]
    fn complementary_determinant_preserves_component_jets() {
        use crate::test_support::jet::HoloJ2;
        use ndarray::arr0;
        type J = HoloJ2;
        let p = J::variable(arr0(c(0.3)));
        let zero = J::constant(arr0(c(0.0)));
        let ordinary = Scatter2Entries::from_parts(zero.clone(), p.clone().exp(), p.exp(), zero);
        let projective = Scatter2ProjectiveEntries::from_entries(&ordinary);
        let h = projective.complement();
        let expected = (0.6_f64).exp();
        for (actual, expected) in [
            (h.value()[()], expected),
            (h.first()[()], 2.0 * expected),
            (h.second()[()], 4.0 * expected),
        ] {
            assert!((actual - c(expected)).norm() < 1e-12);
        }
        let identity = Scatter2ProjectiveEntries::<J>::identity_like(&arr0(c(0.0)));
        for component in [&identity, &projective] {
            let residual = component
                .n12()
                .multiply(component.n21())
                .subtract(&component.n11().multiply(component.n22()))
                .subtract(&component.denominator().multiply(component.complement()));
            assert!(residual.value()[()].norm() < 1e-12);
            assert!(residual.first()[()].norm() < 1e-12);
            assert!(residual.second()[()].norm() < 1e-12);
        }
    }

    #[test]
    fn reduced_cascade_preserves_bivariate_arrays_associativity_and_analytic_scaling() {
        use crate::{
            SeedJet,
            algebra::{ArrayJetBivariate2, HolomorphicParameter},
        };
        use ndarray::{Ix1, array};
        type J = ArrayJetBivariate2<C, Ix1, HolomorphicParameter>;
        let x = <J as SeedJet>::variable(array![c(0.0), c(0.2), C::new(-0.1, 0.2)], 0).unwrap();
        let y = <J as SeedJet>::variable(array![c(0.1), c(0.0), C::new(0.3, -0.1)], 1).unwrap();
        let network = |offset: f64| {
            Scatter2Entries::from_parts(
                x.scale(c(0.1))
                    .add(&y.scale(c(0.03)))
                    .add(&x.constant(c(offset))),
                x.multiply(&y).scale(c(0.02)).add(&x.constant(c(0.7))),
                x.scale(c(0.04)).add(&x.constant(c(0.6))),
                y.square().scale(c(0.02)).subtract(&x.constant(c(0.1))),
            )
        };
        let assert_jet = |a: &J, b: &J| {
            for (a, b) in [
                (a.value(), b.value()),
                (a.axis0(), b.axis0()),
                (a.axis1(), b.axis1()),
                (a.axis0_axis0(), b.axis0_axis0()),
                (a.axis0_axis1(), b.axis0_axis1()),
                (a.axis1_axis1(), b.axis1_axis1()),
            ] {
                for (a, b) in a.iter().zip(b) {
                    assert!((*a - *b).norm() < 1e-12, "{a:?} != {b:?}");
                }
            }
        };
        let assert_entries = |a: &Scatter2Entries<J>, b: &Scatter2Entries<J>| {
            for (a, b) in [
                (a.s11(), b.s11()),
                (a.s12(), b.s12()),
                (a.s21(), b.s21()),
                (a.s22(), b.s22()),
            ] {
                assert_jet(a, b);
            }
        };
        let a = network(0.1);
        let b = network(0.2);
        let c0 = network(0.05);
        let pa = Scatter2ProjectiveEntries::from_entries(&a);
        let pb = Scatter2ProjectiveEntries::from_entries(&b);
        let pc = Scatter2ProjectiveEntries::from_entries(&c0);
        let combined = cascade_projection(&pa, &pb);
        assert_entries(&combined.entries(), &cascade(&a, &b));
        let invariant = combined
            .n12()
            .multiply(combined.n21())
            .subtract(&combined.n11().multiply(combined.n22()));
        assert_jet(
            &invariant,
            &combined.denominator().multiply(combined.complement()),
        );
        assert_entries(
            &cascade_projection(&combined, &pc).entries(),
            &cascade_projection(&pa, &cascade_projection(&pb, &pc)).entries(),
        );
        // A nonconstant scale tests propagation of gauge derivatives, not just values.
        let scale = x.add(&y).exp();
        let scaled = Scatter2ProjectiveEntries::from_parts(
            pa.denominator().multiply(&scale),
            pa.n11().multiply(&scale),
            pa.n12().multiply(&scale),
            pa.n21().multiply(&scale),
            pa.n22().multiply(&scale),
            pa.complement().multiply(&scale),
        );
        assert_entries(
            &cascade_projection(&scaled, &pb).entries(),
            &combined.entries(),
        );
    }

    fn assert_entries_close(actual: &Entries, expected: &Entries) {
        assert_array_close(actual.s11().value(), expected.s11().value(), TOLERANCE);

        assert_array_close(actual.s12().value(), expected.s12().value(), TOLERANCE);

        assert_array_close(actual.s21().value(), expected.s21().value(), TOLERANCE);

        assert_array_close(actual.s22().value(), expected.s22().value(), TOLERANCE);
    }

    #[test]
    fn projective_cascade_reduces_to_ordinary_redheffer_cascade() {
        let left = entries(c(0.10), c(0.80), c(0.70), c(-0.20));

        let right = entries(c(0.30), c(0.60), c(0.50), c(-0.10));

        let expected = cascade(&left, &right);

        let actual = cascade_projection(
            &Projection::from_entries(&left),
            &Projection::from_entries(&right),
        )
        .entries();

        assert_entries_close(&actual, &expected);
    }

    #[test]
    fn projective_cascade_is_unchanged_by_operand_chart_scaling() {
        let left = Projection::from_entries(&entries(c(0.10), c(0.80), c(0.70), c(-0.20)));

        let right = Projection::from_entries(&entries(c(0.30), c(0.60), c(0.50), c(-0.10)));

        let left_scale = c(1.3) + C::i() * c(0.2);

        let right_scale = c(-0.7) + C::i() * c(0.4);

        let scaled_left = Scatter2ProjectiveEntries::from_parts(
            left.denominator().scale(left_scale),
            left.n11().scale(left_scale),
            left.n12().scale(left_scale),
            left.n21().scale(left_scale),
            left.n22().scale(left_scale),
            left.complement().scale(left_scale),
        );

        let scaled_right = Scatter2ProjectiveEntries::from_parts(
            right.denominator().scale(right_scale),
            right.n11().scale(right_scale),
            right.n12().scale(right_scale),
            right.n21().scale(right_scale),
            right.n22().scale(right_scale),
            right.complement().scale(right_scale),
        );

        let expected = cascade_projection(&left, &right).entries();

        let actual = cascade_projection(&scaled_left, &scaled_right).entries();

        assert_entries_close(&actual, &expected);
    }

    #[test]
    fn projective_identity_reduces_to_redheffer_identity() {
        let source = ndarray::arr0(c(0.0));

        let identity: Projection = Projection::identity_like(&source);

        let entries = identity.entries();

        assert_array_close(entries.s11().value(), &ndarray::arr0(c(0.0)), TOLERANCE);

        assert_array_close(entries.s12().value(), &ndarray::arr0(c(1.0)), TOLERANCE);

        assert_array_close(entries.s21().value(), &ndarray::arr0(c(1.0)), TOLERANCE);

        assert_array_close(entries.s22().value(), &ndarray::arr0(c(0.0)), TOLERANCE);
    }
}
