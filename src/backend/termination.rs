//! Closed sides of a stack: walls and wall-backed PMLs in the backends.
//!
//! A closed side replaces the outgoing radiation condition of an open
//! exterior by a wall condition, applied after propagation through the gap
//! between the outermost interface and the wall or PML, and through the PML.
//! For the characteristic a PML of thickness `d` is a homogeneous layer of
//! the exterior material with complex thickness `τ = d (1 + c/(m+1))`.
//!
//! Each PML side is normalised by the β-independent factor `exp(iκ_ref τ)`,
//! with `κ_ref` the PML material's wavevector at `β = 0` (`k₀ √(εμ)`, or
//! `k₀ √(ε_yy μ_xx)` for TE and `k₀ √(ε_xx μ_yy)` for TM in a diagonal one), so that
//! the characteristic stays entire and even in `β` while thick PMLs do not
//! overflow. The growing exponentials are never formed on their own: every
//! exponent is combined with the normalisation before exponentiating.
//!
//! The exterior material's wavevector on a closed side, used for the
//! directional decomposition of fields there, is the root `κ̃` of the
//! exterior's `κ²` for the polarisation with `Im(κ̃ ℓ) ≥ 0`, `ℓ` the total complex length from the
//! outermost interface to the wall. Every quantity reported by the backends is
//! even in that choice.

use nalgebra::ComplexField;
use ndarray::Dimension;
use num_traits::{One, Zero};
use thiserror::Error;

use crate::{
    ComplexScalar,
    algebra::ScalarAlgebra,
    backend::{ExteriorWavevectors, LayerQuantities, PlaneWaveModeCandidate},
    input::{
        CanonicalCoordinates, CanonicalStack,
        canonical::{CanonicalPml, CanonicalTermination, WallCondition},
    },
    stack::StackSide,
    waves::BidirectionalWaves,
};

/// Largest magnitude of the natural logarithm of a PML normalisation scale
/// accepted at an evaluation point.
///
/// Double precision represents magnitudes between about `e^−708` and `e^709`;
/// the margin leaves room for the rest of the stack.
pub(crate) const PML_LOG_SCALE_LIMIT: f64 = 600.0;

/// Above this `|Im κτ|` the PML matrix is formed from combined exponentials
/// rather than from its entire cosine and sine series.
const PML_EXPONENTIAL_FORM: f64 = 30.0;

/// Failure caused by the terminations of a stack.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum TerminationError {
    /// The supplied exterior wavevectors do not match the terminations: open
    /// sides need a wavevector, closed (wall or PML) sides take none.
    #[error(
        "the {side:?} side is {}",
        if *open {
            "open, but no exterior wavevector was supplied"
        } else {
            "closed, but an exterior wavevector was supplied"
        }
    )]
    ExteriorBranchMismatch { side: StackSide, open: bool },

    /// The normalised PML contribution at the evaluation point is too large
    /// or too small to be represented: `log_magnitude` is the natural
    /// logarithm of its scale, rounded to an integer.
    #[error(
        "the {side:?} PML normalisation scale exp({log_magnitude}) is outside the \
         representable range"
    )]
    PmlScaleOutOfRange { side: StackSide, log_magnitude: i64 },

    /// A plane-wave excitation's internal fields were requested for a stack
    /// with a closed side.
    #[error("internal fields of a plane-wave excitation are not available for a closed stack")]
    DrivenFieldsOfClosedStack,
}

/// Check that the supplied exterior wavevectors match the terminations.
pub(crate) fn check_exterior_pattern<M, J, A>(
    stack: &CanonicalStack<M, J>,
    exterior: &ExteriorWavevectors<A>,
) -> Result<(), TerminationError> {
    for (side, termination, supplied) in [
        (StackSide::Left, stack.left_termination(), exterior.left()),
        (
            StackSide::Right,
            stack.right_termination(),
            exterior.right(),
        ),
    ] {
        let open = termination.is_open();

        if open != supplied.is_some() {
            return Err(TerminationError::ExteriorBranchMismatch { side, open });
        }
    }

    Ok(())
}

/// A PML as seen by the backends: its canonical layer and the PML material's
/// wavevector `κ_ref` at `β = 0`.
#[derive(Clone, Debug, PartialEq)]
struct ClosedPml<A> {
    layer: CanonicalPml<A>,
    reference_kappa: A,
}

/// One closed side as seen by the backends.
#[doc(hidden)]
#[derive(Clone, Debug, PartialEq)]
pub struct ClosedSide<A> {
    side: StackSide,
    condition: WallCondition,
    pub(super) kappa: A,
    kappa_squared: A,
    factor: A,
    gap: A,
    pml: Option<ClosedPml<A>>,
    reference_slope: A,
}

impl<A> ClosedSide<A> {
    pub(crate) fn try_map<B, E>(
        &self,
        mut f: impl FnMut(&A) -> Result<B, E>,
    ) -> Result<ClosedSide<B>, E> {
        Ok(ClosedSide {
            side: self.side,
            condition: self.condition,
            kappa: f(&self.kappa)?,
            kappa_squared: f(&self.kappa_squared)?,
            factor: f(&self.factor)?,
            gap: f(&self.gap)?,
            pml: match &self.pml {
                Some(pml) => Some(ClosedPml {
                    layer: pml.layer.try_map(&mut f)?,
                    reference_kappa: f(&pml.reference_kappa)?,
                }),
                None => None,
            },
            reference_slope: f(&self.reference_slope)?,
        })
    }
}

impl<A> ClosedSide<A>
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    /// Build the closed-side data for `termination`, or `None` for an open
    /// side, checking the PML normalisation scale.
    pub(super) fn new(
        side: StackSide,
        termination: &CanonicalTermination<A>,
        medium: &LayerQuantities<A>,
        coordinates: &CanonicalCoordinates<A>,
    ) -> Result<Option<Self>, TerminationError> {
        let (Some(kind), Some(gap)) = (termination.wall_kind(), termination.gap_cm()) else {
            return Ok(None);
        };

        let k0 = coordinates.vacuum_angular_wavenumber();

        let mut closed = Self {
            side,
            condition: WallCondition::new(kind, medium.polarisation()),
            kappa: medium.kappa().clone(),
            kappa_squared: medium.kappa_squared().clone(),
            factor: medium.factor().clone(),
            gap: gap.clone(),
            pml: termination.pml().map(|layer| ClosedPml {
                layer: layer.clone(),
                reference_kappa: medium
                    .constitutive()
                    .index_squared(medium.polarisation())
                    .sqrt()
                    .multiply(k0),
            }),
            reference_slope: k0.scale(-A::Scalar::i()),
        };

        closed.kappa = decaying_root(medium.kappa(), &closed.length());

        closed.check_scale()?;

        Ok(Some(closed))
    }

    /// Refuse a PML normalisation scale `exp(−Im(κ_ref τ) + |Im(κ̃ τ)|)` that
    /// is not representable at some sample.
    fn check_scale(&self) -> Result<(), TerminationError> {
        let Some(pml) = &self.pml else {
            return Ok(());
        };

        let tau = &pml.layer.optical_thickness_cm;
        let reference = pml.reference_kappa.multiply(tau);
        let growth = self.kappa.multiply(tau);
        let limit: <A::Scalar as ComplexField>::RealField = nalgebra::convert(PML_LOG_SCALE_LIMIT);

        for (reference, growth) in reference.value().iter().zip(growth.value()) {
            let log_scale = -reference.imaginary() + ComplexField::abs(growth.imaginary());

            if !ComplexField::is_finite(&log_scale) || ComplexField::abs(log_scale.clone()) > limit
            {
                let log_magnitude = nalgebra::try_convert::<_, f64>(log_scale)
                    .filter(|value| !value.is_nan())
                    .map_or(i64::MAX, |value| value.round() as i64);

                return Err(TerminationError::PmlScaleOutOfRange {
                    side: self.side,
                    log_magnitude,
                });
            }
        }

        Ok(())
    }

    /// Total complex length from the outermost interface to the wall.
    fn length(&self) -> A {
        match &self.pml {
            Some(pml) => self.gap.add(&pml.layer.optical_thickness_cm),
            None => self.gap.clone(),
        }
    }

    /// Transfer matrix `[[cos, −f S], [q S / f, cos]]` of a homogeneous layer
    /// of the exterior material with thickness `d`, as `(cos, S)`, where
    /// `S = sin(κd)/κ`; entire in `q = κ²`.
    fn layer_entries(&self, thickness: &A) -> (A, A) {
        let phase_squared = self.kappa_squared.multiply(&thickness.square());

        (
            phase_squared.cos_sqrt(),
            thickness.multiply(&phase_squared.sinc_sqrt()),
        )
    }

    /// `exp(iκ_ref τ) cos(κτ)` and `exp(iκ_ref τ) sin(κτ)/κ`, without forming
    /// the growing exponentials on their own.
    fn normalised_pml_entries(&self, pml: &ClosedPml<A>) -> (A, A) {
        let tau = &pml.layer.optical_thickness_cm;
        let i = A::Scalar::i();
        let reference = pml.reference_kappa.multiply(tau).scale(i).exp();

        let threshold: <A::Scalar as ComplexField>::RealField =
            nalgebra::convert(PML_EXPONENTIAL_FORM);

        let growth = self.kappa.multiply(tau);

        let series_mask = growth
            .value()
            .mapv(|value| ComplexField::abs(value.imaginary()) <= threshold);

        if series_mask.iter().all(|series| *series) {
            let (cos, sin) = self.layer_entries(tau);
            return (cos.multiply(&reference), sin.multiply(&reference));
        }

        let one = A::Scalar::one();
        let zero = A::Scalar::zero();
        let series = A::lift_constant(series_mask.mapv(|series| if series { one } else { zero }));
        let exponential =
            A::lift_constant(series_mask.mapv(|series| if series { zero } else { one }));

        // Series form, with q set to zero where the exponential form is used.
        let series_q = self.kappa_squared.multiply(&series);
        let series_phase = series_q.multiply(&tau.square());
        let series_cos = series_phase.cos_sqrt().multiply(&reference);
        let series_sin = tau.multiply(&series_phase.sinc_sqrt()).multiply(&reference);

        // Exponential form with the root κ̂, Im(κ̂ τ) ≥ 0, and q set to one where
        // the series form is used so that the root stays regular.
        let exponential_q = self.kappa_squared.multiply(&exponential).add(&series);
        let kappa = decaying_root(&exponential_q.sqrt(), tau);
        let two = one + one;
        let decreasing = pml
            .reference_kappa
            .subtract(&kappa)
            .multiply(tau)
            .scale(i)
            .exp();
        let round_trip = kappa.multiply(tau).scale(i * two).exp();
        let exponential_cos = decreasing
            .multiply(&round_trip.add(&round_trip.one_like()))
            .scale(one / two);
        let exponential_sin = decreasing
            .multiply(&round_trip.subtract(&round_trip.one_like()))
            .divide(&kappa.scale(i * two));

        (
            series_cos
                .multiply(&series)
                .add(&exponential_cos.multiply(&exponential)),
            series_sin
                .multiply(&series)
                .add(&exponential_sin.multiply(&exponential)),
        )
    }

    /// Apply a transfer matrix `[[c, −f s], [q s / f, c]]` to a state.
    fn apply(&self, (cos, sin): &(A, A), (field, slope): (A, A)) -> (A, A) {
        (
            cos.multiply(&field)
                .subtract(&self.factor.multiply(sin).multiply(&slope)),
            self.kappa_squared
                .divide(&self.factor)
                .multiply(sin)
                .multiply(&field)
                .add(&cos.multiply(&slope)),
        )
    }

    /// Apply the row vector `[a, b]` to a transfer matrix, returning the row.
    fn apply_row(&self, (cos, sin): &(A, A), (a, b): (A, A)) -> (A, A) {
        (
            a.multiply(cos)
                .add(&b.multiply(&self.kappa_squared.divide(&self.factor).multiply(sin))),
            b.multiply(cos)
                .subtract(&a.multiply(&self.factor.multiply(sin))),
        )
    }

    /// The transfer state at the outermost interface of a closed right side:
    /// `T(gap) T̂(τ) r_wall`, with `r_wall = [0, −ξ_ref]` (Dirichlet) or
    /// `[1, 0]` (Neumann) and `T̂` the normalised PML matrix.
    pub(crate) fn right_state(&self) -> (A, A) {
        let one = self.gap.one_like();
        let zero = self.gap.zero_like();

        let mut state = match self.condition {
            WallCondition::Dirichlet => (zero, self.reference_slope.negate()),
            WallCondition::Neumann => (one, zero),
        };

        if let Some(pml) = &self.pml {
            state = self.apply(&self.normalised_pml_entries(pml), state);
        }

        self.apply(&self.layer_entries(&self.gap), state)
    }

    /// The boundary row at the outermost interface of a closed left side:
    /// `ℓ_wall T̂(τ) T(gap)`, with `ℓ_wall = [ξ_ref, 0]` (Dirichlet) or
    /// `[0, −1]` (Neumann).
    pub(crate) fn left_row(&self) -> (A, A) {
        let one = self.gap.one_like();
        let zero = self.gap.zero_like();

        let mut row = match self.condition {
            WallCondition::Dirichlet => (self.reference_slope.clone(), zero),
            WallCondition::Neumann => (zero, one.negate()),
        };

        if let Some(pml) = &self.pml {
            row = self.apply_row(&self.normalised_pml_entries(pml), row);
        }

        self.apply_row(&self.layer_entries(&self.gap), row)
    }

    /// The wall's reflection in the exterior material's directional basis:
    /// `−1` (Dirichlet) or `+1` (Neumann).
    fn reflection(&self) -> A::Scalar {
        match self.condition {
            WallCondition::Dirichlet => -A::Scalar::one(),
            WallCondition::Neumann => A::Scalar::one(),
        }
    }

    /// The load seen from the outermost interface, `ρ exp(2iκ̃ℓ)`.
    pub(crate) fn load(&self) -> A {
        let two = A::Scalar::one() + A::Scalar::one();

        self.kappa
            .multiply(&self.length())
            .scale(A::Scalar::i() * two)
            .exp()
            .scale(self.reflection())
    }

    /// The exponent of this side's scale factor in the scattering form of the
    /// characteristic: `−iκ̃ℓ + iκ_ref τ`.
    pub(crate) fn scale_exponent(&self) -> A {
        let i = A::Scalar::i();
        let exponent = self.kappa.multiply(&self.length()).scale(-i);

        match &self.pml {
            Some(pml) => exponent.add(
                &pml.reference_kappa
                    .multiply(&pml.layer.optical_thickness_cm)
                    .scale(i),
            ),
            None => exponent,
        }
    }

    /// The side coefficient `c` in the scattering form of the characteristic,
    /// in the exterior material's directional basis at the wall.
    pub(crate) fn coefficient(&self) -> A {
        let half = (A::Scalar::one() + A::Scalar::one()).recip();
        let slope = self.kappa.divide(&self.factor).scale(-A::Scalar::i());

        match (self.side, self.condition) {
            (StackSide::Left, WallCondition::Dirichlet) => self.reference_slope.clone(),
            (StackSide::Left, WallCondition::Neumann) => slope,
            (StackSide::Right, WallCondition::Dirichlet) => {
                self.reference_slope.divide(&slope).scale(half)
            }
            (StackSide::Right, WallCondition::Neumann) => self.gap.one_like().scale(half),
        }
    }

    /// The complex distance from the outermost interface of a point at real
    /// outward distance `distance` (centimetres), or `None` beyond the wall.
    ///
    /// In the gap this is the distance itself; inside the PML it is
    /// `gap + ζ̃(ζ)` with `ζ = distance − gap` and
    /// `ζ̃(ζ) = ζ + c d/(m+1) (ζ/d)^(m+1)`.
    pub(crate) fn complex_distance(
        &self,
        distance: <A::Scalar as ComplexField>::RealField,
    ) -> Option<A::Scalar> {
        type R<A> = <<A as crate::algebra::Jet>::Scalar as ComplexField>::RealField;

        let point = |jet: &A| jet.value().iter().next().copied().expect("a sampled value");
        let one = R::<A>::one();
        let tolerance = one.clone() + nalgebra::convert(1.0e-12);
        let gap: R<A> = point(&self.gap).real();

        let Some(pml) = self.pml.as_ref().filter(|_| distance > gap) else {
            // In the gap; a sample on a wall given in other units may lie
            // beyond it by round-off.
            if distance > gap.clone() * tolerance {
                return None;
            }

            let distance = if distance > gap { gap } else { distance };

            return Some(A::Scalar::from_real(distance));
        };

        let thickness: R<A> = point(&pml.layer.thickness_cm).real();

        if distance > (gap.clone() + thickness.clone()) * tolerance {
            return None;
        }

        let mut zeta: R<A> = distance - gap.clone();

        if zeta > thickness {
            zeta = thickness.clone();
        }

        let order: R<A> = point(&pml.layer.order).real();
        let strength = point(&pml.layer.strength);
        let fraction: R<A> = zeta.clone() / thickness.clone();
        let power: R<A> = ComplexField::powf(fraction, order.clone() + one.clone());
        let stretched = A::Scalar::from_real(zeta)
            + strength * A::Scalar::from_real(thickness / (order + one) * power);

        Some(A::Scalar::from_real(gap) + stretched)
    }
}

/// Return `κ` multiplied by `±1` per sample so that `Im(κ ℓ) ≥ 0`.
fn decaying_root<A>(kappa: &A, length: &A) -> A
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let product = kappa.multiply(length);

    if product
        .value()
        .iter()
        .all(|value| value.imaginary() >= <A::Scalar as ComplexField>::RealField::zero())
    {
        return kappa.clone();
    }

    let signs = product.value().mapv(|value| {
        if value.imaginary() < <A::Scalar as ComplexField>::RealField::zero() {
            -A::Scalar::one()
        } else {
            A::Scalar::one()
        }
    });

    kappa.multiply(&A::lift_constant(signs))
}

/// The right exterior waves of a modal seed at the outermost interface: the
/// outgoing wave alone on an open side, the standing wave on a closed one.
pub(crate) fn right_exterior_mode_waves<A>(
    seed: &PlaneWaveModeCandidate<A>,
    right_admittance: &A,
) -> BidirectionalWaves<A>
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    match seed.right_closed_state() {
        Some(state) => {
            let slope = right_admittance.scale(-A::Scalar::i());
            let difference = state.secondary().divide(&slope);
            let half = (A::Scalar::one() + A::Scalar::one()).recip();

            BidirectionalWaves::new(
                state.field().subtract(&difference).scale(half),
                state.field().add(&difference).scale(half),
            )
        }
        None => BidirectionalWaves::new(
            seed.right_outgoing().clone(),
            seed.right_outgoing().zero_like(),
        ),
    }
}

/// Amplitudes for incidence from a closed side, which has no meaning; the
/// evaluators refuse it before projecting.
pub(crate) fn undefined_amplitudes<A>(
    source: &ndarray::ArrayBase<ndarray::OwnedRepr<A::Scalar>, A::Dimension>,
) -> crate::PlaneWaveAmplitudes<A>
where
    A: ScalarAlgebra + Clone,
    A::Scalar: ComplexScalar,
{
    let nan = A::filled_constant_like(source, A::Scalar::from_real(nalgebra::convert(f64::NAN)));

    crate::PlaneWaveAmplitudes::new(nan.clone(), nan)
}
