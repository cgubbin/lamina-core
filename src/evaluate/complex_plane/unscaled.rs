//! Owned primal homogeneous fields, independent of QNM normalization.
use super::{
    ComplexPlaneEvaluator,
    mode::{ExteriorWaves, reconstruct_mode_fields},
};
use crate::{
    AnalyticalMaterialHandle, CanonicalCoordinates, ComplexPlane, ElectromagneticFields,
    ExteriorWavevectors, FieldSampling, Polarisation, SpatialResponse, Stack, StackRegion,
    TensorConstitutiveLift, VectorField,
    algebra::{ArrayJet0, HolomorphicParameter},
    backend::{LayerQuantities, ModalSolutionSource, Scatter2, Transfer2},
};
use ndarray::{Ix0, Ix1, arr0};
use num_complex::Complex64 as C;
type J = ArrayJet0<C, Ix0, HolomorphicParameter>;
/// Owned E/H components in core Cartesian and impedance conventions.
pub type OwnedModeFields = SpatialResponse<ElectromagneticFields<VectorField<C, Ix1>>, f64>;
#[derive(Clone, Debug)]
pub struct UnscaledFieldResult {
    pub fields: OwnedModeFields,
    /// Backend seed boundary residual divided by a finite nonzero seed scale.
    /// Scale = |secondary| + |xi_left * primary|; dimensionless ratio.
    pub relative_boundary_defect: f64,
    pub seed_residual: C,
    /// Fresh characteristic in the selected backend gauge, before amplitude scaling.
    pub characteristic: C,
    pub seed_scale: f64,
}
#[derive(Debug, thiserror::Error)]
pub enum UnscaledFieldError {
    #[error("invalid or inconsistent physical coordinates, exterior roots or gauge")]
    InvalidInput,
    #[error("stack compilation failed")]
    Stack(#[source] crate::StackCompileError<f64>),
    #[error("retained backend evaluation failed")]
    Backend(#[source] std::sync::Arc<dyn std::error::Error + Send + Sync>),
    #[error("homogeneous seed reconstruction failed")]
    Seed(#[source] crate::backend::ModeReconstructionError),
    #[error("field sampling or reconstruction failed")]
    Fields(#[source] crate::observable::FieldReconstructionError<f64>),
    #[error("degenerate or nonfinite homogeneous seed scale")]
    DegenerateSeed,
    #[error("sampled fields are not representable")]
    NonFiniteOutput,
}
/// High-level recovery capability implemented by the built-in 2×2 backends,
/// for stacks of scalar materials or regions (`M`).
/// Uses the backend's arbitrary seed gauge, then multiplies both E and H by
/// `amplitude`. No energy integral or normalization factor is evaluated.
/// Exterior roots are supplied in physical left/right order, in cm^-1, for
/// open sides only: a side closed by a wall or a PML takes `None`, and its
/// fields are the standing wave between the outermost interface and the wall.
pub trait ModalFieldBackend<M = AnalyticalMaterialHandle<C>> {
    #[allow(clippy::too_many_arguments)] // explicit physical data; no parallel device representation
    fn unscaled_fields(
        &self,
        stack: &Stack<M, f64>,
        k0: C,
        beta: C,
        kappa: [Option<C>; 2],
        polarisation: Polarisation,
        sampling: &FieldSampling<f64>,
        amplitude: C,
    ) -> Result<UnscaledFieldResult, UnscaledFieldError>;
}
macro_rules! implementation {
    ($backend:ty) => {
        impl<M> ModalFieldBackend<M> for $backend
        where
            M: Clone,
            J: TensorConstitutiveLift<ComplexPlane, M>,
        {
            fn unscaled_fields(
                &self,
                stack: &Stack<M, f64>,
                k0: C,
                beta: C,
                kappa: [Option<C>; 2],
                polarisation: Polarisation,
                sampling: &FieldSampling<f64>,
                amplitude: C,
            ) -> Result<UnscaledFieldResult, UnscaledFieldError> {
                let finite = |v: C| v.re.is_finite() && v.im.is_finite();
                if ![k0, beta, amplitude].into_iter().all(finite)
                    || !kappa.iter().flatten().copied().all(finite)
                    || amplitude.norm() == 0.
                {
                    return Err(UnscaledFieldError::InvalidInput);
                }
                let coordinates =
                    CanonicalCoordinates::new(J::constant(arr0(k0)), J::constant(arr0(beta)));
                let mut left_factor = C::new(0., 0.);
                let mut left_kappa = C::new(0., 0.);
                for (i, (region, material, open)) in [
                    (
                        StackRegion::LeftExterior,
                        stack.left_exterior(),
                        stack.left_termination().is_open(),
                    ),
                    (
                        StackRegion::RightExterior,
                        stack.right_exterior(),
                        stack.right_termination().is_open(),
                    ),
                ]
                .into_iter()
                .enumerate()
                {
                    if open != kappa[i].is_some() {
                        return Err(UnscaledFieldError::InvalidInput);
                    }
                    let medium = LayerQuantities::evaluate::<ComplexPlane, _>(
                        region,
                        material,
                        &coordinates,
                        polarisation,
                    )
                    .map_err(super::shared_error)
                    .map_err(UnscaledFieldError::Backend)?;
                    let index_squared =
                        medium.constitutive().index_squared(polarisation).value()[()];
                    let q = medium.kappa_squared().value()[()];
                    // A closed side's wavevector enters only through q.
                    let root = kappa[i].unwrap_or_else(|| q.sqrt());
                    let scale =
                        (index_squared * k0 * k0).norm() + (beta * beta).norm() + root.norm_sqr();
                    if !finite(index_squared)
                        || !finite(q)
                        || !scale.is_finite()
                        || scale == 0.
                        || (root * root - q).norm() > 1e-8 * scale
                    {
                        return Err(UnscaledFieldError::InvalidInput);
                    }
                    if i == 0 {
                        left_kappa = root;
                        left_factor = medium.factor().value()[()];
                    }
                }
                let evaluator = ComplexPlaneEvaluator::<J, _, _>::compile(stack, self.clone())
                    .map_err(UnscaledFieldError::Stack)?;
                let state = evaluator
                    .retain(
                        CanonicalCoordinates::new(J::constant(arr0(k0)), J::constant(arr0(beta))),
                        ExteriorWavevectors::from_options(
                            kappa[0].map(|kappa| J::constant(arr0(kappa))),
                            kappa[1].map(|kappa| J::constant(arr0(kappa))),
                        ),
                        polarisation,
                    )
                    .map_err(super::shared_error)
                    .map_err(UnscaledFieldError::Backend)?;
                let characteristic = state.determinant().value().value()[()];
                if !finite(characteristic) {
                    return Err(UnscaledFieldError::NonFiniteOutput);
                }
                let seed = state
                    .workspace()
                    .modal_boundary_solution()
                    .map_err(UnscaledFieldError::Seed)?;
                let residual = seed.clone().into_residual().value()[()];
                let seed_scale = seed.state().secondary().value()[()].norm()
                    + (C::new(0., -1.) * left_kappa / left_factor
                        * seed.state().field().value()[()])
                        .norm();
                let relative_boundary_defect = residual.norm() / seed_scale;
                if seed_scale == 0.
                    || !seed_scale.is_finite()
                    || !relative_boundary_defect.is_finite()
                {
                    return Err(UnscaledFieldError::DegenerateSeed);
                }
                let resolved = sampling
                    .resolve(stack)
                    .map_err(|e| UnscaledFieldError::Fields(e.into()))?;
                let scaled = seed.scaled(&J::constant(arr0(amplitude)));
                let fields =
                    reconstruct_mode_fields(&state, &scaled, &resolved, ExteriorWaves::Outgoing)
                        .map_err(UnscaledFieldError::Fields)?;
                let (e, h) = fields.into_parts();
                let e = e.into_inner();
                let h = h.into_inner();
                if [e.x(), e.y(), e.z(), h.x(), h.y(), h.z()]
                    .into_iter()
                    .flat_map(|a| a.iter())
                    .any(|v| !finite(*v))
                {
                    return Err(UnscaledFieldError::NonFiniteOutput);
                }
                Ok(UnscaledFieldResult {
                    fields: SpatialResponse::new(ElectromagneticFields::new(e, h), resolved),
                    relative_boundary_defect,
                    seed_residual: residual,
                    characteristic,
                    seed_scale,
                })
            }
        }
    };
}
implementation!(Scatter2);
implementation!(Transfer2);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ScalarAlgebra;
    use crate::{Constant, ExteriorSampling, LayerSampling, Stack};
    use lamina_units::Length;
    fn c(v: f64) -> C {
        C::new(v, 0.)
    }
    #[test]
    fn interface_without_normalisation_and_common_gauge() {
        let stack = Stack::from_analytical_materials(Constant::vacuum(), Constant::dielectric(-2.))
            .finalise();
        let k = [C::new(0., 1.), C::new(0., 2.)];
        let sampling = FieldSampling::new()
            .left_exterior(ExteriorSampling::distances([
                Length::zero(),
                Length::centimetres(0.3),
            ]))
            .right_exterior(ExteriorSampling::distances([
                Length::zero(),
                Length::centimetres(0.3),
            ]));
        let a = Scatter2::new()
            .unscaled_fields(
                &stack,
                c(1.),
                c(2_f64.sqrt()),
                k.map(Some),
                Polarisation::TransverseMagnetic,
                &sampling,
                c(1.),
            )
            .unwrap();
        let b = Transfer2::new()
            .unscaled_fields(
                &stack,
                c(1.),
                c(2_f64.sqrt()),
                k.map(Some),
                Polarisation::TransverseMagnetic,
                &sampling,
                C::new(2., 1.),
            )
            .unwrap();
        let ah = a.fields.quantity().magnetic().y();
        let bh = b.fields.quantity().magnetic().y();
        let scale = ah[0] / bh[0];
        for (av, bv) in [
            a.fields.quantity().electric().x(),
            a.fields.quantity().electric().z(),
            ah,
        ]
        .into_iter()
        .zip([
            b.fields.quantity().electric().x(),
            b.fields.quantity().electric().z(),
            bh,
        ]) {
            for (x, y) in av.iter().zip(bv.iter()) {
                assert!((*x - scale * y).norm() < 1e-12);
            }
        }
        assert!((ah[1] / ah[0] - c((-0.3_f64).exp())).norm() < 1e-12);
        assert!((ah[3] / ah[2] - c((-0.6_f64).exp())).norm() < 1e-12);
        assert_eq!(a.fields.sampling().len(), 4);
        assert!(a.relative_boundary_defect < 1e-12);
        let evaluator =
            ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Transfer2::new()).unwrap();
        let state = evaluator
            .retain(
                CanonicalCoordinates::new(
                    J::constant(arr0(c(1.))),
                    J::constant(arr0(c(2_f64.sqrt()))),
                ),
                ExteriorWavevectors::new(J::constant(arr0(k[0])), J::constant(arr0(k[1]))),
                Polarisation::TransverseMagnetic,
            )
            .unwrap();
        assert!(
            state.mode().is_err(),
            "zero-layer energy normalization is unavailable while unscaled fields work"
        );
    }
    #[test]
    fn legacy_normalized_fields_are_one_common_scale() {
        let stack = Stack::from_analytical_materials(Constant::vacuum(), Constant::vacuum())
            .analytical_layer(Constant::dielectric(4.), Length::centimetres(0.3))
            .finalise();
        let k0 = C::new(std::f64::consts::PI / 0.6, -3_f64.ln() / 0.6);
        let sampling = FieldSampling::new()
            .left_exterior(ExteriorSampling::point(Length::zero()))
            .layer(0, LayerSampling::uniform(3))
            .right_exterior(ExteriorSampling::point(Length::zero()));
        let evaluator =
            ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Transfer2::new()).unwrap();
        let state = evaluator
            .retain(
                CanonicalCoordinates::new(J::constant(arr0(k0)), J::constant(arr0(c(0.)))),
                ExteriorWavevectors::new(J::constant(arr0(k0)), J::constant(arr0(k0))),
                Polarisation::TransverseElectric,
            )
            .unwrap();
        let normalized = state.mode().unwrap();
        let expected = normalized.fields(&sampling).unwrap();
        let result = Transfer2::new()
            .unscaled_fields(
                &stack,
                k0,
                c(0.),
                [Some(k0); 2],
                Polarisation::TransverseElectric,
                &sampling,
                c(1.),
            )
            .unwrap();
        let ey = expected.quantity().electric().value().y();
        let hy = expected.quantity().magnetic().value();
        let raw = result.fields.quantity();
        let scale = ey[0] / raw.electric().y()[0];
        for (a, b) in [
            (ey, raw.electric().y()),
            (hy.x(), raw.magnetic().x()),
            (hy.z(), raw.magnetic().z()),
        ] {
            for (x, y) in a.iter().zip(b.iter()) {
                assert!((*x - scale * y).norm() < 1e-10);
            }
        }
        assert!(normalized.seed_normalisation().total().all_finite());
    }
}
