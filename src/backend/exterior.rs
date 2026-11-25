use ndarray::Dimension;

use crate::{
    ComplexScalar, IncidentSide, Polarisation,
    algebra::ScalarAlgebra,
    backend::{
        Channel, ChannelError, ClosedSide, ExteriorContextProvider, TerminationError,
        medium::LayerQuantities, termination::check_exterior_pattern,
    },
    input::{CanonicalCoordinates, CanonicalStack, canonical::CanonicalTermination},
    material::{ConstitutiveError, TensorConstitutiveLift},
    stack::{StackRegion, StackSide},
};

/// The branch-selected outward normal wavenumbers of one open exterior.
///
/// A TE or TM channel has one outgoing wave per open side, a coupled channel
/// two. Values are *outward*: an outgoing wave varies as `exp(iκ |z|)` away
/// from the stack, so the right exterior's outgoing normal wavenumbers are
/// `κ` and the left exterior's are `−κ`. An isotropic exterior gives `κ` for a
/// scalar channel and `(κ, κ)` for the coupled channel on both sides.
#[derive(Clone, Debug, PartialEq)]
pub enum ExteriorBranch<A> {
    /// The outward normal wavenumber of a TE or TM channel.
    Scalar(A),
    /// The outward normal wavenumbers of the two outgoing waves of the
    /// coupled channel, in either order.
    Pair([A; 2]),
}

impl<A> ExteriorBranch<A> {
    /// The wavenumber of a scalar branch, or `None` for a pair.
    pub fn scalar(&self) -> Option<&A> {
        match self {
            Self::Scalar(kappa) => Some(kappa),
            Self::Pair(_) => None,
        }
    }

    /// Return every supplied wavenumber.
    pub(crate) fn values(&self) -> &[A] {
        match self {
            Self::Scalar(kappa) => std::slice::from_ref(kappa),
            Self::Pair(pair) => pair,
        }
    }
}

/// Canonical outward normal angular wavevectors in the open exterior media.
///
/// Each side holds the branch-selected [`ExteriorBranch`] of its
/// semi-infinite exterior medium when that side is open, and nothing when the
/// side is closed by a wall or a PML: a closed side's characteristic does not
/// depend on a choice of exterior branch. A pattern that disagrees with the
/// stack's terminations is rejected by the backends with
/// [`TerminationError::ExteriorBranchMismatch`](crate::backend::TerminationError),
/// and a branch kind that does not match the channel with
/// [`ChannelError::BranchKind`](crate::backend::ChannelError).
///
/// Branch selection is completed before this value reaches a numerical
/// backend. This allows advanced callers, such as complex-frequency mode
/// solvers, to supply branch choices determined by an external conformal map
/// rather than relying on the backend to reconstruct them.
#[derive(Clone, Debug, PartialEq)]
pub struct ExteriorWavevectors<A> {
    left: Option<ExteriorBranch<A>>,
    right: Option<ExteriorBranch<A>>,
}

impl<A> ExteriorWavevectors<A> {
    /// Construct scalar-channel wavevectors for a stack open on both sides.
    pub fn new(left: A, right: A) -> Self {
        Self::from_options(Some(left), Some(right))
    }

    /// Construct scalar-channel wavevectors for a stack closed on the left
    /// and open on the right.
    pub fn closed_left(right: A) -> Self {
        Self::from_options(None, Some(right))
    }

    /// Construct scalar-channel wavevectors for a stack open on the left and
    /// closed on the right.
    pub fn closed_right(left: A) -> Self {
        Self::from_options(Some(left), None)
    }

    /// Construct the empty pattern of a stack closed on both sides.
    pub fn closed() -> Self {
        Self::from_branches(None, None)
    }

    /// Construct scalar-channel wavevectors from optional per-side values.
    pub fn from_options(left: Option<A>, right: Option<A>) -> Self {
        Self::from_branches(
            left.map(ExteriorBranch::Scalar),
            right.map(ExteriorBranch::Scalar),
        )
    }

    /// Construct wavevectors from optional per-side branches.
    pub fn from_branches(
        left: Option<ExteriorBranch<A>>,
        right: Option<ExteriorBranch<A>>,
    ) -> Self {
        Self { left, right }
    }

    /// Return the branch of the left exterior medium, if the left side is
    /// open.
    pub fn left(&self) -> Option<&ExteriorBranch<A>> {
        self.left.as_ref()
    }

    /// Return the branch of the right exterior medium, if the right side is
    /// open.
    pub fn right(&self) -> Option<&ExteriorBranch<A>> {
        self.right.as_ref()
    }

    /// Return `true` when every supplied wavevector satisfies `predicate`.
    pub(crate) fn all(&self, mut predicate: impl FnMut(&A) -> bool) -> bool {
        self.left
            .iter()
            .chain(self.right.iter())
            .flat_map(ExteriorBranch::values)
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
        Error: From<TerminationError> + From<ConstitutiveError> + From<ChannelError>,
    {
        check_exterior_pattern(stack, exterior)?;

        let side = |side, region, material, termination, supplied: Option<&ExteriorBranch<A>>| {
            let supplied = match supplied {
                None => None,
                Some(branch) => Some(branch.scalar().ok_or(ChannelError::BranchKind {
                    side,
                    channel: polarisation.into(),
                })?),
            };
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

/// Principal exterior branches of `channel` for the open sides of `stack`:
/// the principal root of each polarisation's `κ²`, one for TE or TM and
/// `(κ_TE, κ_TM)` for the coupled channel.
pub(crate) fn evaluate_stack_exterior_wavevectors<E, M, J>(
    coordinates: &CanonicalCoordinates<J>,
    stack: &CanonicalStack<M, J>,
    channel: Channel,
) -> Result<ExteriorWavevectors<J>, ConstitutiveError>
where
    J: ScalarAlgebra + TensorConstitutiveLift<E, M> + Clone,
    J::Scalar: ComplexScalar,
    J::Dimension: Dimension,
{
    let kappa = |region, material: &M, polarisation| {
        LayerQuantities::evaluate::<E, M>(region, material, coordinates, polarisation)
            .map(LayerQuantities::into_kappa)
    };
    let branch = |region, material: &M| {
        Ok(match channel.polarisation() {
            Some(polarisation) => ExteriorBranch::Scalar(kappa(region, material, polarisation)?),
            None => ExteriorBranch::Pair([
                kappa(region, material, Polarisation::TransverseElectric)?,
                kappa(region, material, Polarisation::TransverseMagnetic)?,
            ]),
        })
    };
    let open = |termination: &CanonicalTermination<J>, region, material| {
        termination
            .is_open()
            .then(|| branch(region, material))
            .transpose()
    };

    Ok(ExteriorWavevectors::from_branches(
        open(
            stack.left_termination(),
            StackRegion::LeftExterior,
            stack.left_exterior(),
        )?,
        open(
            stack.right_termination(),
            StackRegion::RightExterior,
            stack.right_exterior(),
        )?,
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

        assert_jet_close(
            exterior.left().unwrap().scalar().unwrap(),
            left_quantities.kappa(),
        );

        assert_jet_close(
            exterior.right().unwrap().scalar().unwrap(),
            right_quantities.kappa(),
        );
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

        assert_jet_close(
            exterior.left().unwrap().scalar().unwrap(),
            left_expected.kappa(),
        );

        assert_jet_close(
            exterior.right().unwrap().scalar().unwrap(),
            right_expected.kappa(),
        );

        assert_ne!(
            exterior.left().unwrap().scalar().unwrap().value()[()],
            exterior.right().unwrap().scalar().unwrap().value()[()],
        );
    }

    #[test]
    fn exterior_wavevectors_are_independent_of_polarisation() {
        let coordinates = coordinates(C::new(2.5, 0.0), C::new(0.7, 0.0));

        let left = Constant::new(2.0, 3.0);
        let right = Constant::new(4.0, 5.0);

        let te = evaluate_exterior_wavevectors::<RealAxis, _, J>(&coordinates, &left, &right);

        let tm = evaluate_exterior_wavevectors::<RealAxis, _, J>(&coordinates, &left, &right);

        assert_jet_close(
            te.left().unwrap().scalar().unwrap(),
            tm.left().unwrap().scalar().unwrap(),
        );
        assert_jet_close(
            te.right().unwrap().scalar().unwrap(),
            tm.right().unwrap().scalar().unwrap(),
        );
    }
}
