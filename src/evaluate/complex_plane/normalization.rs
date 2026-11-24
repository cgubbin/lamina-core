//! Complete scalar pole norm. No infinite exterior integral is evaluated.
use super::{
    ComplexPlaneEvaluator, UnscaledFieldError,
    mode::{raw_layer_integration_inputs_unchecked, raw_qnm_normalisation_unchecked},
};
use crate::{
    AnalyticalMaterialHandle, AnalyticalMaterialStack, CanonicalCoordinates, DerivativeOrder,
    EvaluateDifferentiableMeromorphicMaterial, EvaluateMeromorphicMaterial, ExteriorSampling,
    ExteriorWavevectors, FieldSampling, ModalFieldBackend, Polarisation, Scalar,
    algebra::{ArrayJet0, HolomorphicParameter},
    backend::{ModalSolutionSource, Scatter2, Transfer2},
};
use lamina_units::Length;
use ndarray::{Ix0, arr0};
use num_complex::Complex64 as C;
type J = ArrayJet0<C, Ix0, HolomorphicParameter>;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PoleVariable {
    K0AtFixedBeta,
    BetaAtFixedK0,
}
#[derive(Clone, Copy, Debug)]
pub struct PoleNormalizationInput {
    pub k0: C,
    pub beta: C,
    pub kappa: [C; 2],
    pub polarisation: Polarisation,
    pub variable: PoleVariable,
}
#[derive(Clone, Copy, Debug)]
pub struct PoleNormalizationSettings {
    pub amplitude: C,
    /// Outward movement of the left and right finite integration planes, cm.
    pub exterior_offsets_cm: [f64; 2],
    /// Factor units: squared primary field in the core seed gauge.
    pub factor_atol: f64,
    pub factor_rtol: f64,
    pub factor_reference_scale: f64,
    pub max_condition: f64,
    /// Exterior kappa exclusion in cm^-1; constitutive exclusion is dimensionless.
    pub kappa_cutoff: f64,
    pub constitutive_cutoff: f64,
}
impl Default for PoleNormalizationSettings {
    fn default() -> Self {
        Self {
            amplitude: C::new(1., 0.),
            exterior_offsets_cm: [0.; 2],
            factor_atol: 1e-12,
            factor_rtol: 1e-12,
            factor_reference_scale: 1.,
            max_condition: 1e10,
            kappa_cutoff: 1e-10,
            constitutive_cutoff: 1e-12,
        }
    }
}
#[derive(Clone, Debug)]
pub struct PoleNormalizationData {
    pub input: PoleNormalizationInput,
    pub settings: PoleNormalizationSettings,
    pub convention: &'static str,
    pub factor_units: &'static str,
    pub layers: Vec<C>,
    pub exterior_segments: [C; 2],
    pub outgoing_boundaries: [C; 2],
    pub total: C,
    pub contribution_scale: f64,
    pub condition: f64,
    pub exclusion_threshold: f64,
    pub scale: C,
    /// Compatibility evidence only; not added to the scalar factor a second time.
    pub legacy_finite_maxwell: Option<C>,
    pub characteristic: C,
    pub relative_boundary_defect: f64,
}
impl PoleNormalizationData {
    /// Res_omega = c * Res_k0; factor transforms inversely. c uses cm/s.
    pub fn angular_frequency_factor(&self) -> Option<C> {
        if self.input.variable == PoleVariable::K0AtFixedBeta {
            Some(self.total / 2.99792458e10)
        } else {
            None
        }
    }
}
#[derive(Debug, thiserror::Error)]
pub enum PoleNormalizationError {
    #[error("invalid normalization settings or physical input")]
    InvalidInput,
    #[error("unscaled homogeneous recovery failed")]
    Fields(#[source] UnscaledFieldError),
    #[error("core normalization evaluation failed")]
    Evaluation(#[source] std::sync::Arc<dyn std::error::Error + Send + Sync>),
    #[error("singular constitutive factor or exterior cutoff in {region}")]
    Singular { region: String },
    #[error("nonfinite normalization contribution or scale")]
    NonFinite,
    #[error("small or ill-conditioned normalization factor")]
    Unavailable {
        evidence: Box<PoleNormalizationData>,
    },
    #[error("pole normalization with outgoing exterior terms needs a stack open on both sides")]
    ClosedStack,
}
pub trait ModalNormalizationBackend: ModalFieldBackend {
    fn pole_normalization(
        &self,
        stack: &AnalyticalMaterialStack<C>,
        input: PoleNormalizationInput,
        settings: PoleNormalizationSettings,
    ) -> Result<PoleNormalizationData, PoleNormalizationError>;
}
fn finite(v: C) -> bool {
    v.re.is_finite() && v.im.is_finite()
}
/// Exact local material values and spectral derivatives, not finite differences.
fn constitutive(m: &AnalyticalMaterialHandle<C>, k: C, spectral: bool) -> [C; 4] {
    let k = Scalar(k);
    let e = m.evaluate_relative_permittivity_complex(k);
    let u = m.evaluate_relative_permeability_complex(k);
    if spectral {
        [
            e,
            u,
            m.evaluate_relative_permittivity_complex_derivative(k, DerivativeOrder::First),
            m.evaluate_relative_permeability_complex_derivative(k, DerivativeOrder::First),
        ]
    } else {
        [e, u, C::new(0., 0.), C::new(0., 0.)]
    }
}
fn coefficients(v: [C; 4], input: PoleNormalizationInput) -> (C, C, C, C) {
    let [e, u, ed, ud] = v;
    let (f, g, fd, gd) = if input.polarisation == Polarisation::TransverseElectric {
        (u, e, ud, ed)
    } else {
        (e, u, ed, ud)
    };
    let qd = match input.variable {
        PoleVariable::K0AtFixedBeta => {
            2. * input.k0 * g + input.k0 * input.k0 * gd + input.beta * input.beta * fd / (f * f)
        }
        PoleVariable::BetaAtFixedK0 => -2. * input.beta / f,
    };
    (f, fd, -fd / (f * f), qd)
}
macro_rules! implement {
    ($backend:ty)=>{
        impl ModalNormalizationBackend for $backend {
            fn pole_normalization(
                &self,
                stack: &AnalyticalMaterialStack<C>,
                input: PoleNormalizationInput,
                settings: PoleNormalizationSettings,
            ) -> Result<PoleNormalizationData, PoleNormalizationError> {
                if !stack.is_open() {
                    return Err(PoleNormalizationError::ClosedStack);
                }
                let cfg = settings;
                if !finite(cfg.amplitude)
                    || cfg.amplitude.norm() == 0.
                    || [
                        cfg.factor_atol,
                        cfg.factor_rtol,
                        cfg.factor_reference_scale,
                        cfg.kappa_cutoff,
                        cfg.constitutive_cutoff,
                        cfg.max_condition,
                        cfg.exterior_offsets_cm[0],
                        cfg.exterior_offsets_cm[1],
                    ]
                    .into_iter()
                    .any(|v| !v.is_finite() || v < 0.)
                    || cfg.factor_reference_scale == 0.
                    || cfg.kappa_cutoff == 0.
                    || cfg.constitutive_cutoff == 0.
                    || cfg.max_condition < 1.
                {
                    return Err(PoleNormalizationError::InvalidInput);
                }
                let sampling = FieldSampling::new()
                    .left_exterior(ExteriorSampling::point(Length::zero()))
                    .right_exterior(ExteriorSampling::point(Length::zero()));
                let recovered = self
                    .unscaled_fields(
                        stack,
                        input.k0,
                        input.beta,
                        input.kappa.map(Some),
                        input.polarisation,
                        &sampling,
                        cfg.amplitude,
                    )
                    .map_err(PoleNormalizationError::Fields)?;
                let primary = if input.polarisation == Polarisation::TransverseElectric {
                    recovered.fields.quantity().electric().y()
                } else {
                    recovered.fields.quantity().magnetic().y()
                };
                let evaluator = ComplexPlaneEvaluator::<J, _, _>::compile(stack, self.clone())
                    .map_err(super::shared_error)
                    .map_err(PoleNormalizationError::Evaluation)?;
                let state = evaluator
                    .retain(
                        CanonicalCoordinates::new(J::constant(arr0(input.k0)), J::constant(arr0(input.beta))),
                        ExteriorWavevectors::new(
                            J::constant(arr0(input.kappa[0])),
                            J::constant(arr0(input.kappa[1])),
                        ),
                        input.polarisation,
                    )
                    .map_err(super::shared_error)
                    .map_err(PoleNormalizationError::Evaluation)?;
                let seed = state
                    .workspace()
                    .modal_boundary_solution()
                    .map_err(super::shared_error)
                    .map_err(PoleNormalizationError::Evaluation)?
                    .scaled(&J::constant(arr0(cfg.amplitude)));
                let legacy_finite_maxwell = raw_qnm_normalisation_unchecked(&seed, &state)
                    .ok()
                    .map(|n| n.total().value()[()]);
                let integrated = raw_layer_integration_inputs_unchecked(&seed, state.workspace())
                    .map_err(super::shared_error)
                    .map_err(PoleNormalizationError::Evaluation)?
                    .integrate_bilinear()
                    .into_inner();
                let spectral = input.variable == PoleVariable::K0AtFixedBeta;
                let mut layers = Vec::new();
                for (i, (layer, material)) in integrated
                    .into_iter()
                    .zip(stack.layers_left_to_right())
                    .enumerate()
                {
                    let v = constitutive(material.material(), input.k0, spectral);
                    let (f, fd, _, qd) = coefficients(v, input);
                    if !v.into_iter().all(finite) || f.norm() <= cfg.constitutive_cutoff {
                        return Err(PoleNormalizationError::Singular {
                            region: format!("layer {i}"),
                        });
                    }
                    let (products, _) = layer.into_parts();
                    // secondary=psi'/f, hence -p_lambda*(psi')² = f_lambda*secondary².
                    layers.push(
                        fd * products.secondary_secondary().value()[()]
                            + qd * products.field_field().value()[()],
                    );
                }
                let mut outgoing_boundaries = [C::new(0., 0.); 2];
                let mut exterior_segments = outgoing_boundaries;
                for (i, m) in [stack.left_exterior(), stack.right_exterior()]
                    .into_iter()
                    .enumerate()
                {
                    let v = constitutive(m, input.k0, spectral);
                    let (f, fd, pd, qd) = coefficients(v, input);
                    let k = input.kappa[i];
                    if !v.into_iter().all(finite)
                        || f.norm() <= cfg.constitutive_cutoff
                        || k.norm() <= cfg.kappa_cutoff
                    {
                        return Err(PoleNormalizationError::Singular {
                            region: format!("exterior {i}"),
                        });
                    }
                    let qd_wave = if spectral {
                        (v[2] * v[1] + v[0] * v[3]) * input.k0 * input.k0 + 2. * v[0] * v[1] * input.k0
                    } else {
                        -2. * input.beta
                    };
                    let kd = qd_wave / (2. * k);
                    // Opposite Robin signs cancel the opposite endpoint signs.
                    let boundary = C::new(0., 1.) * (kd / f - k * fd / (f * f)) * primary[i] * primary[i];
                    let offset = cfg.exterior_offsets_cm[i];
                    let argument = C::new(0., 1.) * 2. * k * offset;
                    let integral = J::constant(arr0(argument)).exprel().value()[()] * offset;
                    exterior_segments[i] = (pd * k * k + qd) * primary[i] * primary[i] * integral;
                    outgoing_boundaries[i] = boundary * argument.exp();
                }
                let total = layers
                    .iter()
                    .copied()
                    .chain(exterior_segments)
                    .chain(outgoing_boundaries)
                    .sum::<C>();
                let contribution_scale = layers
                    .iter()
                    .chain(exterior_segments.iter())
                    .chain(outgoing_boundaries.iter())
                    .map(|v| v.norm())
                    .sum::<f64>();
                if !finite(total)
                    || !contribution_scale.is_finite()
                    || layers
                        .iter()
                        .chain(exterior_segments.iter())
                        .chain(outgoing_boundaries.iter())
                        .any(|v| !finite(*v))
                {
                    return Err(PoleNormalizationError::NonFinite);
                }
                let condition = contribution_scale / total.norm();
                let exclusion_threshold =
                    cfg.factor_atol * cfg.factor_reference_scale + cfg.factor_rtol * contribution_scale;
                if !exclusion_threshold.is_finite() {
                    return Err(PoleNormalizationError::InvalidInput);
                }
                let scale = C::new(1., 0.) / total.sqrt();
                let data = PoleNormalizationData {
                    input,
                    settings: cfg,
                    convention: "scalar TE/TM outgoing operator; positive unit-delta jump; one planar Fourier component",
                    factor_units: "squared primary seed field (physical pole variable in cm^-1)",
                    layers,
                    exterior_segments,
                    outgoing_boundaries,
                    total,
                    contribution_scale,
                    condition,
                    exclusion_threshold,
                    scale,
                    legacy_finite_maxwell,
                    characteristic: recovered.characteristic,
                    relative_boundary_defect: recovered.relative_boundary_defect,
                };
                if total.norm() <= exclusion_threshold
                    || !condition.is_finite()
                    || condition > cfg.max_condition
                    || !finite(scale)
                {
                    return Err(PoleNormalizationError::Unavailable {
                        evidence: Box::new(data),
                    });
                }
                Ok(data)
            }
        }
    }
}
implement!(Scatter2);
implement!(Transfer2);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Constant, Stack};
    fn c(v: f64) -> C {
        C::new(v, 0.)
    }
    fn cavity(pol: Polarisation) -> (AnalyticalMaterialStack<C>, PoleNormalizationInput) {
        let stack = Stack::from_analytical_materials(Constant::vacuum(), Constant::vacuum())
            .analytical_layer(Constant::dielectric(4.), Length::centimetres(0.3))
            .finalise();
        let k = C::new(std::f64::consts::PI / 0.6, -3_f64.ln() / 0.6);
        (
            stack,
            PoleNormalizationInput {
                k0: k,
                beta: c(0.),
                kappa: [k; 2],
                polarisation: pol,
                variable: PoleVariable::K0AtFixedBeta,
            },
        )
    }
    fn close(a: C, b: C) {
        assert!(
            (a - b).norm() < 1e-9 * (1. + a.norm() + b.norm()),
            "{a} != {b}"
        );
    }
    #[test]
    fn normal_incidence_maxwell_cancellation_is_not_double_counted() {
        for pol in [
            Polarisation::TransverseElectric,
            Polarisation::TransverseMagnetic,
        ] {
            let (stack, input) = cavity(pol);
            let n = Transfer2::new()
                .pole_normalization(&stack, input, Default::default())
                .unwrap();
            let sign = if pol == Polarisation::TransverseElectric {
                1.
            } else {
                -1.
            };
            close(n.legacy_finite_maxwell.unwrap(), c(0.9));
            close(n.total, sign * input.k0 * 0.9);
            close(n.total, sign * input.k0 * n.legacy_finite_maxwell.unwrap());
            // Boundary term is needed in the scalar expression, but adding it
            // once more to the equivalent Maxwell factor is incorrect.
            assert!(
                (n.total + n.outgoing_boundaries.iter().sum::<C>() - sign * input.k0 * 0.9).norm()
                    > 0.1
            );
            close(
                n.angular_frequency_factor().unwrap() * 2.99792458e10,
                n.total,
            );
            let shifted = Transfer2::new()
                .pole_normalization(
                    &stack,
                    input,
                    PoleNormalizationSettings {
                        exterior_offsets_cm: [0.11, 0.23],
                        ..Default::default()
                    },
                )
                .unwrap();
            close(n.total, shifted.total);
            let a = C::new(0.7, -0.8);
            let rescaled = Transfer2::new()
                .pole_normalization(
                    &stack,
                    input,
                    PoleNormalizationSettings {
                        amplitude: a,
                        ..Default::default()
                    },
                )
                .unwrap();
            close(rescaled.total, a * a * n.total);
        }
    }
    #[test]
    fn analytic_integration_includes_electric_and_magnetic_dispersion() {
        use crate::test_support::materials::LinearDispersion;
        for pol in [
            Polarisation::TransverseElectric,
            Polarisation::TransverseMagnetic,
        ] {
            let material = LinearDispersion::new(3., 0.2, 1.2, 0.07);
            let stack = Stack::from_analytical_materials(Constant::vacuum(), Constant::vacuum())
                .analytical_layer(material, Length::centimetres(0.3))
                .finalise();
            let k = C::new(2., -0.1);
            let beta = c(0.4);
            let exterior = (k * k - beta * beta).sqrt();
            let input = PoleNormalizationInput {
                k0: k,
                beta,
                kappa: [exterior; 2],
                polarisation: pol,
                variable: PoleVariable::K0AtFixedBeta,
            };
            let n = Transfer2::new()
                .pole_normalization(&stack, input, Default::default())
                .unwrap();
            let e = c(3.) + 0.2 * k;
            let u = c(1.2) + 0.07 * k;
            let (f, g, fd, gd) = if pol == Polarisation::TransverseElectric {
                (u, e, c(0.07), c(0.2))
            } else {
                (e, u, c(0.2), c(0.07))
            };
            let wave = (e * u * k * k - beta * beta).sqrt();
            let qd = 2. * k * g + k * k * gd + beta * beta * fd / (f * f);
            let integrand = |z: f64| {
                let x = wave * (z - 0.3);
                let psi = x.cos() + C::i() * exterior * f / wave * x.sin();
                let dp = -wave * x.sin() + C::i() * exterior * f * x.cos();
                fd / (f * f) * dp * dp + qd * psi * psi
            };
            let count = 2000;
            let h = 0.3 / count as f64;
            let sum = (0..=count)
                .map(|i| {
                    integrand(i as f64 * h)
                        * if i == 0 || i == count {
                            1.
                        } else if i % 2 == 0 {
                            2.
                        } else {
                            4.
                        }
                })
                .sum::<C>()
                * h
                / 3.;
            close(n.layers[0], sum);
        }
    }
    #[test]
    fn diagnostics_keep_small_factor_evidence() {
        let (stack, input) = cavity(Polarisation::TransverseElectric);
        let error = Transfer2::new()
            .pole_normalization(
                &stack,
                input,
                PoleNormalizationSettings {
                    factor_atol: 1e9,
                    ..Default::default()
                },
            )
            .unwrap_err();
        assert!(
            matches!(error,PoleNormalizationError::Unavailable{evidence} if evidence.layers.len()==1 && evidence.total.norm()>0.)
        );
        assert!(matches!(
            Transfer2::new().pole_normalization(
                &stack,
                input,
                PoleNormalizationSettings {
                    amplitude: c(0.),
                    ..Default::default()
                }
            ),
            Err(PoleNormalizationError::InvalidInput)
        ));
        assert!(
            Transfer2::new()
                .pole_normalization(
                    &stack,
                    input,
                    PoleNormalizationSettings {
                        exterior_offsets_cm: [1e9, 0.],
                        ..Default::default()
                    }
                )
                .is_err()
        );
    }
}
