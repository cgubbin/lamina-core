//! Single solutions of a retained 4×4 workspace, their fields and their
//! integrated observables.
//!
//! A solution is the tangential state `ψ` at every interface. Inside a region
//! `ψ(t) = exp(tA) ψ₀`, with `t` the depth from a finite layer's left
//! interface (`A = iK`) or outward from an exterior's interface (`A = iK` on
//! the right, `−iK` on the left). Quadratic forms of the six components
//! `f = (E, η) = T ψ`, `T` the region's `6 × n` field map, integrate in closed
//! form over finite layers and wall gaps, and over open exteriors whose
//! field holds only outgoing waves (see [`super::integral`]).
//!
//! Forms are weighted by a pair of `3 × 3` tensors `[W_e, W_m]` per region,
//! `E_a^{†,ᵀ} W_e E_b` and `η_a^{†,ᵀ} W_m η_b`. With `W = [ε, μ]` the real
//! part of the Hermitian form is the energy density and its imaginary part the
//! dissipation density, `E† ((ε − ε†)/2i) E` and likewise for `μ`.

use ndarray::Dimension;

use crate::{
    BilinearLayerOverlap, Channel, ComplexScalar, FiniteLayerIndex, HermitianLayerOverlap,
    IncidentSide, LayerDissipation, LayerEnergy, LayerPower, Layers,
    algebra::{RealScalarAlgebra, ScalarAlgebra},
    material::ConstitutiveTensor,
    stack::StackRegion,
};

use super::{
    CoupledContext, CoupledWorkspace,
    berreman::{Region, reference_positions},
    exponential::exponential,
    integral,
    matrix::Matrix,
};

/// One solution of a retained 4×4 workspace: the state at every interface.
#[doc(hidden)]
#[derive(Clone, Debug)]
pub struct CoupledSolution<'a, A> {
    workspace: &'a CoupledWorkspace<A>,
    /// The `n × 1` state at every interface, left to right.
    states: Vec<Matrix<A>>,
    /// Whether the field of each exterior (left, right) holds only its
    /// outgoing waves.
    outgoing: [bool; 2],
}

/// Weights `[W_e, W_m]` of a region given its stack-frame `[ε, μ]`.
pub trait FieldWeights<A>:
    Fn(StackRegion, &[ConstitutiveTensor<A>; 2]) -> [ConstitutiveTensor<A>; 2]
{
}

impl<A, F> FieldWeights<A> for F where
    F: Fn(StackRegion, &[ConstitutiveTensor<A>; 2]) -> [ConstitutiveTensor<A>; 2]
{
}

/// How a form pairs its left operand: conjugate transpose or transpose.
type Pairing<A> = fn(&Matrix<A>) -> Matrix<A>;

/// One region of a solution: its generator data, the exponent `A` per unit
/// depth, the state at its reference interface and its extent.
struct Segment<'s, A> {
    region: &'s Region<A>,
    exponent: Matrix<A>,
    state: Matrix<A>,
    extent: Extent<'s, A>,
}

enum Extent<'s, A> {
    Finite(&'s A),
    /// An open exterior holding only outgoing waves: the outgoing basis and
    /// the coefficients of the state in it.
    Outgoing {
        basis: &'s Matrix<A>,
        coefficients: Matrix<A>,
    },
}

impl<'a, A> CoupledSolution<'a, A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    pub(crate) fn new(
        workspace: &'a CoupledWorkspace<A>,
        states: Vec<Matrix<A>>,
        outgoing: [bool; 2],
    ) -> Self {
        Self {
            workspace,
            states,
            outgoing,
        }
    }

    fn context(&self) -> &'a CoupledContext<A> {
        self.workspace.solution.context()
    }

    /// The solution multiplied by `factor`.
    pub fn scaled(&self, factor: &A) -> Self {
        Self {
            states: self
                .states
                .iter()
                .map(|state| state.scale_by(factor))
                .collect(),
            ..self.clone()
        }
    }

    /// The region, exponent, reference state and extent of `region`; `None`
    /// when it does not exist, or for an exterior whose field cannot be
    /// integrated (a PML, or incoming waves).
    fn segment(&self, region: StackRegion) -> Option<Segment<'a, A>> {
        let retained = self.workspace.retained.as_ref()?;
        let i = A::Scalar::i();
        let (boundary, side, state) = match region {
            StackRegion::Layer(layer) => {
                let (region, thickness) = retained.layers.get(layer.get())?;
                return Some(Segment {
                    region,
                    exponent: region.generator.scale(i),
                    state: self.states.get(layer.get())?.clone(),
                    extent: Extent::Finite(thickness),
                });
            }
            StackRegion::LeftExterior => (&self.context().left, 0, self.states.first()?),
            StackRegion::RightExterior => (&self.context().right, 1, self.states.last()?),
        };
        let sign = if side == 0 { -i } else { i };
        let exponent = boundary.region.generator.scale(sign);
        let extent = match (&boundary.gap, &boundary.incoming) {
            (Some(gap), _) => Extent::Finite(gap),
            (None, Some(incoming)) if self.outgoing[side] => Extent::Outgoing {
                basis: &boundary.outgoing,
                coefficients: boundary
                    .outgoing
                    .hstack(incoming)
                    .solve(state)
                    .rows_in(0..boundary.outgoing.columns()),
            },
            _ => return None,
        };
        Some(Segment {
            region: &boundary.region,
            exponent,
            state: state.clone(),
            extent,
        })
    }

    /// `[E, Z₀H]`, all six components in the stack frame, at `depth` into
    /// `region`: from a finite layer's left interface, or outward from an
    /// exterior's interface. An exterior that holds only outgoing waves keeps
    /// only those.
    pub fn fields(&self, region: StackRegion, depth: &A) -> Option<[[A; 3]; 2]> {
        let channel = self.context().channel;
        let (region, state) = match self.segment(region) {
            // Propagated in the outgoing subspace, so that rounding in the
            // incoming waves does not grow with depth.
            Some(Segment {
                region,
                exponent,
                extent:
                    Extent::Outgoing {
                        basis,
                        coefficients,
                    },
                ..
            }) => {
                let restricted = restricted(&exponent, basis, channel);
                let carried = exponential(&restricted.scale_by(depth)).multiply(&coefficients);
                (region, basis.multiply(&carried))
            }
            Some(segment) => (
                segment.region,
                exponential(&segment.exponent.scale_by(depth)).multiply(&segment.state),
            ),
            None => {
                // A PML side or an incident side: the full state.
                let (boundary, state, sign) = match region {
                    StackRegion::LeftExterior => {
                        (&self.context().left, self.states.first()?, -A::Scalar::i())
                    }
                    StackRegion::RightExterior => {
                        (&self.context().right, self.states.last()?, A::Scalar::i())
                    }
                    StackRegion::Layer(_) => return None,
                };
                let exponent = boundary.region.generator.scale(sign).scale_by(depth);
                (&boundary.region, exponential(&exponent).multiply(state))
            }
        };
        let fields = region.field_map(channel).multiply(&state);
        let entry = |row: usize| fields.get(row, 0).clone();
        Some([
            [entry(0), entry(1), entry(2)],
            [entry(3), entry(4), entry(5)],
        ])
    }

    /// `∫ (T_a ψ_a)^{†,ᵀ} M (T_b ψ_b) dz` over `region`, for each `6 × 6`
    /// weight `M` produced from the two regions.
    fn integrate(
        &self,
        other: &Self,
        region: StackRegion,
        pairing: Pairing<A>,
        weights: &[Matrix<A>],
    ) -> Option<Vec<A>> {
        let (a, b) = (self.segment(region)?, other.segment(region)?);
        let channel = self.context().channel;
        let (map_a, map_b) = (a.region.field_map(channel), b.region.field_map(channel));
        weights
            .iter()
            .map(|weight| {
                let q = pairing(&map_a).multiply(weight).multiply(&map_b);
                let value = match (&a.extent, &b.extent) {
                    (Extent::Finite(thickness), Extent::Finite(_)) => {
                        let kernel =
                            integral::layer(&pairing(&a.exponent), &q, &b.exponent, thickness);
                        pairing(&a.state).multiply(&kernel).multiply(&b.state)
                    }
                    (
                        Extent::Outgoing {
                            basis: basis_a,
                            coefficients: c_a,
                        },
                        Extent::Outgoing {
                            basis: basis_b,
                            coefficients: c_b,
                        },
                    ) => {
                        let (x, y) = (
                            restricted(&a.exponent, basis_a, channel),
                            restricted(&b.exponent, basis_b, channel),
                        );
                        let q = pairing(basis_a).multiply(&q).multiply(basis_b);
                        let kernel = integral::half_line(&pairing(&x), &q, &y);
                        pairing(c_a).multiply(&kernel).multiply(c_b)
                    }
                    _ => return None,
                };
                Some(value.get(0, 0).clone())
            })
            .collect()
    }

    /// The bilinear overlap `∫ E_aᵀ W_e E_b dz`, `∫ η_aᵀ W_m η_b dz` of this
    /// solution and `other` over `region`; holomorphic in the coordinates.
    pub fn bilinear_overlap(
        &self,
        other: &Self,
        region: StackRegion,
        weights: impl FieldWeights<A>,
    ) -> Option<BilinearLayerOverlap<A>> {
        let [electric, magnetic] = self.overlap(other, region, Matrix::transpose, &weights)?;
        let total = electric.add(&magnetic);
        Some(BilinearLayerOverlap::from_parts(electric, magnetic, total))
    }

    fn overlap(
        &self,
        other: &Self,
        region: StackRegion,
        pairing: Pairing<A>,
        weights: &impl FieldWeights<A>,
    ) -> Option<[A; 2]> {
        let tensors = &self.segment(region)?.region.tensors;
        let [electric, magnetic] = weights(region, tensors);
        let like = self.states.first()?.get(0, 0);
        let block = |tensor: &ConstitutiveTensor<A>, offset: usize| {
            let full = tensor.full();
            Matrix::from_fn(6, 6, |row, column| {
                match (row.checked_sub(offset), column.checked_sub(offset)) {
                    (Some(row), Some(column)) if row < 3 && column < 3 => full[row][column].clone(),
                    _ => like.zero_like(),
                }
            })
        };
        let values = self.integrate(
            other,
            region,
            pairing,
            &[block(&electric, 0), block(&magnetic, 3)],
        )?;
        let [electric, magnetic]: [A; 2] = values.try_into().ok()?;
        Some([electric, magnetic])
    }

    /// The finite-layer regions, left to right.
    pub fn layers(&self) -> impl Iterator<Item = StackRegion> + use<A> {
        (0..self.states.len().saturating_sub(1))
            .map(|index| StackRegion::Layer(FiniteLayerIndex::new(index)))
    }
}

impl<A> CoupledSolution<'_, A>
where
    A: RealScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
    A::RealJet: ScalarAlgebra,
{
    /// The Hermitian overlap `∫ E_a† W_e E_b dz`, `∫ η_a† W_m η_b dz` of this
    /// solution and `other` over `region`. With `other` this solution and
    /// unit weights, the field norms `∫ |E|² dz` and `∫ |η|² dz`.
    pub fn hermitian_overlap(
        &self,
        other: &Self,
        region: StackRegion,
        weights: impl FieldWeights<A>,
    ) -> Option<HermitianLayerOverlap<A>> {
        let [electric, magnetic] = self.overlap(other, region, Matrix::adjoint, &weights)?;
        let total = electric.add(&magnetic);
        Some(HermitianLayerOverlap::new(electric, magnetic, total))
    }

    /// `½ Re ∫ (E_y η_z* − E_z η_y*) dz` over `region`: `Z₀` times the power
    /// per unit width carried along `x` (`η = Z₀H`).
    pub fn longitudinal_power(&self, region: StackRegion) -> Option<A::RealJet> {
        let like = self.states.first()?.get(0, 0);
        // f† M f = η_z* E_y − η_y* E_z, f = (E, η).
        let density = Matrix::from_fn(6, 6, |row, column| match (row, column) {
            (5, 1) => like.one_like(),
            (4, 2) => like.one_like().negate(),
            _ => like.zero_like(),
        });
        let [value] = self
            .integrate(self, region, Matrix::adjoint, &[density])?
            .try_into()
            .ok()?;
        Some(value.real().scale(nalgebra::convert(0.5)))
    }

    /// `½ Re(ψ† Γ ψ)` at interface `index`: `Z₀ S_z`, the power per unit area
    /// crossing it along `+z`.
    pub fn flux(&self, index: usize) -> Option<A::RealJet> {
        Some(super::flux(self.states.get(index)?).scale(nalgebra::convert(0.5)))
    }
}

/// The generator of an invariant subspace with basis `basis`: `R` with
/// `exponent · basis = basis · R`, from the reference rows.
fn restricted<A>(exponent: &Matrix<A>, basis: &Matrix<A>, channel: Channel) -> Matrix<A>
where
    A: ScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
{
    let rows = reference_positions(channel);
    let all: Vec<usize> = (0..basis.columns()).collect();
    basis
        .select(&rows, &all)
        .solve(&exponent.multiply(basis).select(&rows, &all))
}

/// Layer observables of one excitation, per unit incident flux `|Z₀ S_z|` of
/// the incoming wave, as the 2×2 backends report them.
impl<A> CoupledWorkspace<A>
where
    A: RealScalarAlgebra,
    A::Scalar: ComplexScalar,
    A::Dimension: Dimension,
    A::RealJet: ScalarAlgebra,
{
    /// The excitation and the reciprocal of its incident flux.
    fn normalised(
        &self,
        incident_side: IncidentSide,
        wave: usize,
    ) -> Option<(CoupledSolution<'_, A>, A::RealJet)> {
        let solution = self.excitation(incident_side, wave)?;
        let context = self.solution.context();
        let incoming = match incident_side {
            IncidentSide::Left => &context.left,
            IncidentSide::Right => &context.right,
        }
        .incoming
        .as_ref()?;
        let flux = super::flux(&incoming.columns_in(wave..wave + 1)).scale(nalgebra::convert(0.5));
        // An incoming wave carries flux away from its side.
        let magnitude = match incident_side {
            IncidentSide::Left => flux,
            IncidentSide::Right => flux.negate(),
        };
        Some((solution, magnitude.reciprocal()))
    }

    /// Net flux into each finite layer at its two interfaces and the power it
    /// absorbs, signed along `+z`.
    pub fn layer_power(
        &self,
        incident_side: IncidentSide,
        wave: usize,
    ) -> Option<Layers<LayerPower<A::RealJet>>> {
        let (solution, scale) = self.normalised(incident_side, wave)?;
        let fluxes = (0..solution.states.len())
            .map(|index| Some(solution.flux(index)?.multiply(&scale)))
            .collect::<Option<Vec<_>>>()?;
        Some(Layers::new(
            fluxes
                .windows(2)
                .map(|pair| {
                    LayerPower::new(pair[0].clone(), pair[1].clone(), pair[0].subtract(&pair[1]))
                })
                .collect(),
        ))
    }

    /// `¼ ∫ (E† W_e E + η† W_m η) dz` of each finite layer, with the Hermitian
    /// parts of the weights: `[ε, μ]` for the nondispersive energy,
    /// `[∂(k₀ε)/∂k₀, ∂(k₀μ)/∂k₀]` for the dispersive one.
    pub fn layer_energy(
        &self,
        incident_side: IncidentSide,
        wave: usize,
        weights: impl FieldWeights<A>,
    ) -> Option<Layers<LayerEnergy<A::RealJet>>> {
        let (solution, scale) = self.normalised(incident_side, wave)?;
        let scale = scale.scale(nalgebra::convert(0.25));
        solution
            .layers()
            .map(|region| {
                let overlap = solution.hermitian_overlap(&solution, region, &weights)?;
                let (electric, magnetic, total) = overlap.into_parts();
                let part = |value: A| value.real().multiply(&scale);
                Some(LayerEnergy::new(
                    part(electric),
                    part(magnetic),
                    part(total),
                ))
            })
            .collect::<Option<_>>()
            .map(Layers::new)
    }

    /// `(k₀/2) ∫ (E† ε'' E + η† μ'' η) dz` of each finite layer, with
    /// `ε'' = (ε − ε†)/2i` and likewise for `μ`.
    pub fn layer_dissipation(
        &self,
        incident_side: IncidentSide,
        wave: usize,
    ) -> Option<Layers<LayerDissipation<A::RealJet>>> {
        let (solution, scale) = self.normalised(incident_side, wave)?;
        let scale = self
            .solution
            .context()
            .vacuum_angular_wavenumber
            .real()
            .multiply(&scale)
            .scale(nalgebra::convert(0.5));
        solution
            .layers()
            .map(|region| {
                let overlap =
                    solution.hermitian_overlap(&solution, region, |_, tensors| tensors.clone())?;
                let (electric, magnetic, total) = overlap.into_parts();
                let part = |value: A| value.imaginary().multiply(&scale);
                Some(LayerDissipation::new(
                    part(electric),
                    part(magnetic),
                    part(total),
                ))
            })
            .collect::<Option<_>>()
            .map(Layers::new)
    }
}
