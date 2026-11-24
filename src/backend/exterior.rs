use ndarray::Dimension;

use crate::{
    ComplexScalar,
    algebra::ScalarAlgebra,
    backend::isotropic::IsotropicMediumQuantities,
    input::{CanonicalCoordinates, CanonicalStack},
    material::{ConstitutiveEvaluator, ConstitutiveLift},
};

/// Canonical longitudinal angular wavevectors in the open exterior media.
///
/// Each side holds the branch-selected out-of-plane angular wavevector `κ` of
/// its semi-infinite exterior medium when that side is open, and nothing when
/// the side is closed by a wall or a PML: a closed side's characteristic does
/// not depend on a choice of exterior branch. A pattern that disagrees with
/// the stack's terminations is rejected by the backends with
/// [`TerminationError::ExteriorBranchMismatch`](crate::backend::TerminationError).
///
/// Branch selection is completed before this value reaches a numerical
/// backend. This allows advanced callers, such as complex-frequency mode
/// solvers, to supply branch choices determined by an external conformal map
/// rather than relying on the backend to reconstruct them.
#[derive(Clone, Debug, PartialEq)]
pub struct ExteriorWavevectors<A> {
    left: Option<A>,
    right: Option<A>,
}

impl<A> ExteriorWavevectors<A> {
    /// Construct wavevectors for a stack open on both sides.
    pub fn new(left: A, right: A) -> Self {
        Self {
            left: Some(left),
            right: Some(right),
        }
    }

    /// Construct wavevectors for a stack closed on the left and open on the
    /// right.
    pub fn closed_left(right: A) -> Self {
        Self {
            left: None,
            right: Some(right),
        }
    }

    /// Construct wavevectors for a stack open on the left and closed on the
    /// right.
    pub fn closed_right(left: A) -> Self {
        Self {
            left: Some(left),
            right: None,
        }
    }

    /// Construct the empty pattern of a stack closed on both sides.
    pub fn closed() -> Self {
        Self {
            left: None,
            right: None,
        }
    }

    /// Construct wavevectors from optional per-side values.
    pub fn from_options(left: Option<A>, right: Option<A>) -> Self {
        Self { left, right }
    }

    /// Return the wavevector in the left exterior medium, if the left side is
    /// open.
    pub fn left(&self) -> Option<&A> {
        self.left.as_ref()
    }

    /// Return the wavevector in the right exterior medium, if the right side
    /// is open.
    pub fn right(&self) -> Option<&A> {
        self.right.as_ref()
    }

    /// Return `true` when every supplied wavevector satisfies `predicate`.
    pub(crate) fn all(&self, mut predicate: impl FnMut(&A) -> bool) -> bool {
        self.left
            .iter()
            .chain(self.right.iter())
            .all(&mut predicate)
    }
}

#[cfg(test)]
pub(crate) fn evaluate_exterior_wavevectors<E, M, J>(
    coordinates: &CanonicalCoordinates<J>,
    left_exterior: &M,
    right_exterior: &M,
) -> ExteriorWavevectors<J>
where
    J: ScalarAlgebra + ConstitutiveLift<E, M> + Clone,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
    E: ConstitutiveEvaluator<J::Scalar, J::Dimension, M>,
{
    let left_kappa =
        IsotropicMediumQuantities::evaluate::<E, M>(left_exterior, coordinates).into_kappa();

    let right_kappa =
        IsotropicMediumQuantities::evaluate::<E, M>(right_exterior, coordinates).into_kappa();

    ExteriorWavevectors::new(left_kappa, right_kappa)
}

/// Principal-root exterior wavevectors for the open sides of `stack`.
pub(crate) fn evaluate_stack_exterior_wavevectors<E, M, J>(
    coordinates: &CanonicalCoordinates<J>,
    stack: &CanonicalStack<M, J>,
) -> ExteriorWavevectors<J>
where
    J: ScalarAlgebra + ConstitutiveLift<E, M> + Clone,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
    E: ConstitutiveEvaluator<J::Scalar, J::Dimension, M>,
{
    let kappa =
        |material| IsotropicMediumQuantities::evaluate::<E, M>(material, coordinates).into_kappa();

    ExteriorWavevectors::from_options(
        stack
            .left_termination()
            .is_open()
            .then(|| kappa(stack.left_exterior())),
        stack
            .right_termination()
            .is_open()
            .then(|| kappa(stack.right_exterior())),
    )
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use ndarray::{Ix0, arr0};
    use num_complex::Complex64;

    use super::*;

    use crate::{
        Polarisation,
        algebra::{ArrayJet0, RealParameter},
        backend::IsotropicLayerQuantities,
        domain::RealAxis,
        input::CanonicalCoordinates,
        material::Constant,
    };

    type C = Complex64;
    type J = ArrayJet0<C, Ix0, RealParameter>;

    const TOLERANCE: f64 = 1.0e-12;

    fn jet(value: C) -> J {
        J::new(arr0(value))
    }

    fn coordinates(
        vacuum_angular_wavenumber: C,
        parallel_angular_wavenumber: C,
    ) -> CanonicalCoordinates<J> {
        CanonicalCoordinates::new(
            jet(vacuum_angular_wavenumber),
            jet(parallel_angular_wavenumber),
        )
    }

    fn assert_jet_close(actual: &J, expected: &J) {
        let actual = actual.value()[()];
        let expected = expected.value()[()];

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

    #[test]
    fn exterior_wavevectors_match_layer_quantities() {
        let coordinates = coordinates(C::new(2.5, 0.0), C::new(0.7, 0.0));

        let left = Constant::new(2.0, 1.0);
        let right = Constant::new(4.0, 1.5);

        let exterior = evaluate_exterior_wavevectors::<RealAxis, _, J>(&coordinates, &left, &right);

        let left_quantities = IsotropicLayerQuantities::evaluate::<RealAxis, _>(
            &left,
            &coordinates,
            Polarisation::TransverseElectric,
        );

        let right_quantities = IsotropicLayerQuantities::evaluate::<RealAxis, _>(
            &right,
            &coordinates,
            Polarisation::TransverseElectric,
        );

        assert_jet_close(exterior.left().unwrap(), left_quantities.kappa());

        assert_jet_close(exterior.right().unwrap(), right_quantities.kappa());
    }

    #[test]
    fn exterior_wavevectors_preserve_left_and_right_materials() {
        let coordinates = coordinates(C::new(2.5, 0.0), C::new(0.3, 0.0));

        let left = Constant::new(1.0, 1.0);
        let right = Constant::new(9.0, 1.0);

        let exterior = evaluate_exterior_wavevectors::<RealAxis, _, J>(&coordinates, &left, &right);

        let left_expected = IsotropicLayerQuantities::evaluate::<RealAxis, _>(
            &left,
            &coordinates,
            Polarisation::TransverseElectric,
        );

        let right_expected = IsotropicLayerQuantities::evaluate::<RealAxis, _>(
            &right,
            &coordinates,
            Polarisation::TransverseElectric,
        );

        assert_jet_close(exterior.left().unwrap(), left_expected.kappa());

        assert_jet_close(exterior.right().unwrap(), right_expected.kappa());

        assert_ne!(
            exterior.left().unwrap().value()[()],
            exterior.right().unwrap().value()[()],
        );
    }

    #[test]
    fn exterior_wavevectors_are_independent_of_polarisation() {
        let coordinates = coordinates(C::new(2.5, 0.0), C::new(0.7, 0.0));

        let left = Constant::new(2.0, 3.0);
        let right = Constant::new(4.0, 5.0);

        let te = evaluate_exterior_wavevectors::<RealAxis, _, J>(&coordinates, &left, &right);

        let tm = evaluate_exterior_wavevectors::<RealAxis, _, J>(&coordinates, &left, &right);

        assert_jet_close(te.left().unwrap(), tm.left().unwrap());
        assert_jet_close(te.right().unwrap(), tm.right().unwrap());
    }
}
