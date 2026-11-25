use ndarray::Dimension;

use crate::{
    ComplexScalar, IncidentSide, Polarisation,
    algebra::ScalarAlgebra,
    backend::{
        ClosedSide, ExteriorContextProvider, TerminationError, medium::LayerQuantities,
        termination::check_exterior_pattern,
    },
    input::{CanonicalCoordinates, CanonicalStack},
    material::{ConstitutiveError, TensorConstitutiveLift},
    stack::{StackRegion, StackSide},
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

/// Exterior data of a solved stack: each side's medium quantities, whose
/// `κ` and admittance are the supplied (open side) or terminated (closed
/// side) ones, each side's closure, and the coordinates and polarisation.
#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub struct ExteriorContext<A> {
    left: LayerQuantities<A>,
    right: LayerQuantities<A>,
    vacuum_angular_wavenumber: A,
    parallel_angular_wavenumber: A,
    polarisation: Polarisation,
    left_closure: Option<ClosedSide<A>>,
    right_closure: Option<ClosedSide<A>>,
}

impl<A> ExteriorContextProvider for ExteriorContext<A> {
    type Algebra = A;

    fn left(&self) -> &LayerQuantities<A> {
        &self.left
    }

    fn right(&self) -> &LayerQuantities<A> {
        &self.right
    }

    fn vacuum_angular_wavenumber(&self) -> &A {
        &self.vacuum_angular_wavenumber
    }

    fn parallel_angular_wavenumber(&self) -> &A {
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

impl<A> ExteriorContext<A> {
    /// Evaluate the exterior data of `stack`: open sides take their
    /// wavevectors from `exterior`, closed sides from their termination.
    pub(crate) fn for_stack<E, M, Error>(
        coordinates: &CanonicalCoordinates<A>,
        stack: &CanonicalStack<M, A>,
        exterior: &ExteriorWavevectors<A>,
        polarisation: Polarisation,
    ) -> Result<Self, Error>
    where
        A: ScalarAlgebra + TensorConstitutiveLift<E, M> + Clone,
        A::Scalar: ComplexScalar,
        A::Dimension: Dimension,
        Error: From<TerminationError> + From<ConstitutiveError>,
    {
        check_exterior_pattern(stack, exterior)?;

        let side = |side, region, material, termination, supplied: Option<&A>| {
            let medium =
                LayerQuantities::evaluate::<E, M>(region, material, coordinates, polarisation)?;
            let closure = ClosedSide::new(side, termination, &medium, coordinates)?;
            let kappa = match (supplied, &closure) {
                (Some(kappa), _) => kappa.clone(),
                (None, Some(closed)) => closed.kappa.clone(),
                (None, None) => {
                    return Err(
                        TerminationError::ExteriorBranchMismatch { side, open: true }.into(),
                    );
                }
            };
            Ok::<_, Error>((medium.with_kappa(kappa), closure))
        };

        let (left, left_closure) = side(
            StackSide::Left,
            StackRegion::LeftExterior,
            stack.left_exterior(),
            stack.left_termination(),
            exterior.left(),
        )?;
        let (right, right_closure) = side(
            StackSide::Right,
            StackRegion::RightExterior,
            stack.right_exterior(),
            stack.right_termination(),
            exterior.right(),
        )?;

        Ok(Self {
            left,
            right,
            vacuum_angular_wavenumber: coordinates.vacuum_angular_wavenumber().clone(),
            parallel_angular_wavenumber: coordinates.parallel_angular_wavenumber().clone(),
            polarisation,
            left_closure,
            right_closure,
        })
    }

    pub(crate) fn try_map<B, Error>(
        &self,
        mut f: impl FnMut(&A) -> Result<B, Error>,
    ) -> Result<ExteriorContext<B>, Error> {
        Ok(ExteriorContext {
            left: self.left.try_map(&mut f)?,
            right: self.right.try_map(&mut f)?,
            vacuum_angular_wavenumber: f(&self.vacuum_angular_wavenumber)?,
            parallel_angular_wavenumber: f(&self.parallel_angular_wavenumber)?,
            polarisation: self.polarisation,
            left_closure: self
                .left_closure
                .as_ref()
                .map(|closure| closure.try_map(&mut f))
                .transpose()?,
            right_closure: self
                .right_closure
                .as_ref()
                .map(|closure| closure.try_map(&mut f))
                .transpose()?,
        })
    }

    /// The admittances of the incident and the transmitting side.
    pub(crate) fn incident_and_transmitted_admittances(&self, side: IncidentSide) -> (&A, &A) {
        match side {
            IncidentSide::Left => (self.left.admittance(), self.right.admittance()),
            IncidentSide::Right => (self.right.admittance(), self.left.admittance()),
        }
    }

    /// A context of isotropic exteriors from its parts.
    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_parts(
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
    ) -> Self
    where
        A: ScalarAlgebra,
        A::Scalar: ComplexScalar,
        A::Dimension: Dimension,
    {
        let side = |kappa, epsilon, mu, admittance| {
            LayerQuantities::test_fixture(kappa, epsilon, mu, polarisation)
                .with_admittance(admittance)
        };
        Self {
            left: side(left_kappa, left_epsilon, left_mu, left_admittance),
            right: side(right_kappa, right_epsilon, right_mu, right_admittance),
            vacuum_angular_wavenumber,
            parallel_angular_wavenumber,
            polarisation,
            left_closure,
            right_closure,
        }
    }

    /// The context of an open stack with exteriors `left` and `right`.
    #[cfg(test)]
    pub(crate) fn new<E, M>(
        coordinates: &CanonicalCoordinates<A>,
        left: &M,
        right: &M,
        exterior: &ExteriorWavevectors<A>,
        polarisation: Polarisation,
    ) -> Self
    where
        A: ScalarAlgebra + TensorConstitutiveLift<E, M> + Clone,
        A::Scalar: ComplexScalar,
        A::Dimension: Dimension,
        M: Clone,
    {
        let stack = CanonicalStack::new(left.clone(), right.clone(), Vec::new());

        Self::for_stack::<E, M, crate::backend::Scatter2Error>(
            coordinates,
            &stack,
            exterior,
            polarisation,
        )
        .unwrap()
    }
}

/// Principal-root TE exterior wavevectors of two exterior materials.
#[cfg(test)]
pub(crate) fn evaluate_exterior_wavevectors<E, M, J>(
    coordinates: &CanonicalCoordinates<J>,
    left_exterior: &M,
    right_exterior: &M,
) -> ExteriorWavevectors<J>
where
    J: ScalarAlgebra + TensorConstitutiveLift<E, M> + Clone,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    let kappa = |region, material| {
        LayerQuantities::evaluate::<E, M>(
            region,
            material,
            coordinates,
            Polarisation::TransverseElectric,
        )
        .unwrap()
        .into_kappa()
    };

    ExteriorWavevectors::new(
        kappa(StackRegion::LeftExterior, left_exterior),
        kappa(StackRegion::RightExterior, right_exterior),
    )
}

/// Principal-root exterior wavevectors of `polarisation` for the open sides
/// of `stack`.
pub(crate) fn evaluate_stack_exterior_wavevectors<E, M, J>(
    coordinates: &CanonicalCoordinates<J>,
    stack: &CanonicalStack<M, J>,
    polarisation: Polarisation,
) -> Result<ExteriorWavevectors<J>, ConstitutiveError>
where
    J: ScalarAlgebra + TensorConstitutiveLift<E, M> + Clone,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    let kappa = |region, material| {
        LayerQuantities::evaluate::<E, M>(region, material, coordinates, polarisation)
            .map(LayerQuantities::into_kappa)
    };

    Ok(ExteriorWavevectors::from_options(
        stack
            .left_termination()
            .is_open()
            .then(|| kappa(StackRegion::LeftExterior, stack.left_exterior()))
            .transpose()?,
        stack
            .right_termination()
            .is_open()
            .then(|| kappa(StackRegion::RightExterior, stack.right_exterior()))
            .transpose()?,
    ))
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
        backend::LayerQuantities,
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

        let left_quantities =
            LayerQuantities::real_axis(&left, &coordinates, Polarisation::TransverseElectric);

        let right_quantities =
            LayerQuantities::real_axis(&right, &coordinates, Polarisation::TransverseElectric);

        assert_jet_close(exterior.left().unwrap(), left_quantities.kappa());

        assert_jet_close(exterior.right().unwrap(), right_quantities.kappa());
    }

    #[test]
    fn exterior_wavevectors_preserve_left_and_right_materials() {
        let coordinates = coordinates(C::new(2.5, 0.0), C::new(0.3, 0.0));

        let left = Constant::new(1.0, 1.0);
        let right = Constant::new(9.0, 1.0);

        let exterior = evaluate_exterior_wavevectors::<RealAxis, _, J>(&coordinates, &left, &right);

        let left_expected =
            LayerQuantities::real_axis(&left, &coordinates, Polarisation::TransverseElectric);

        let right_expected =
            LayerQuantities::real_axis(&right, &coordinates, Polarisation::TransverseElectric);

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
