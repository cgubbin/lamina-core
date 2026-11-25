//! Modes of the 4×4 backends: the null vector of the interface system and
//! its six-component fields, against Maxwell's equations, against each other
//! and against the 2×2 backends on separable stacks.

use lamina_units::{Length, LengthUnit};
use ndarray::arr0;

use crate::{
    AnalyticalTensorStack, Backend, CanonicalCoordinates, Channel, ComplexPlane,
    ComplexPlaneEvaluator, ExteriorBranch, ExteriorSampling, ExteriorWavevectors, FieldSampling,
    FiniteLayerIndex, LayerSampling, ModalFieldBackend, Polarisation, Rotation, Scatter2, Scatter4,
    SeedJet, Stack, StackRegion, Transfer2, Transfer4,
    backend::{CoupledSolution, CoupledWorkspace, ModeReconstructionError, RetainedLayers},
    material::TensorConstitutiveLift,
    test_support::{jet::RealJ0, materials::ComplexQuadratic},
};

use super::coupled::{C, assert_close, c, isotropic, lossless, oriented};

const K0: f64 = 4.0;
const CLADDING: f64 = 1.0;
const SUBSTRATE: f64 = 2.0736;
const CORE: f64 = 1.0;

/// A film of `core` on an isotropic substrate (right) under an isotropic
/// cladding (left).
fn film(core: crate::AnalyticalMaterialRegion<C>) -> AnalyticalTensorStack<C> {
    Stack::builder(
        isotropic(lossless(CLADDING)),
        isotropic(lossless(SUBSTRATE)),
    )
    .layer(core, Length::centimetres(CORE))
    .finalise()
}

fn uniaxial(rotation: Rotation<f64>) -> crate::AnalyticalMaterialRegion<C> {
    oriented([lossless(4.0), lossless(4.4), lossless(3.6)], rotation)
}

fn jet(value: C) -> RealJ0 {
    <RealJ0 as SeedJet>::constant(arr0(value))
}

/// The decaying normal wavenumber `i√(β² − ε k₀²)` of an isotropic exterior.
fn decaying(epsilon: f64, beta: C) -> C {
    c(0.0, 1.0) * (beta * beta - epsilon * K0 * K0).sqrt()
}

fn inputs(
    beta: C,
    channel: Channel,
) -> (CanonicalCoordinates<RealJ0>, ExteriorWavevectors<RealJ0>) {
    let branch = |epsilon| {
        let kappa = jet(decaying(epsilon, beta));
        Some(match channel {
            Channel::Coupled => ExteriorBranch::Pair([kappa.clone(), kappa]),
            _ => ExteriorBranch::Scalar(kappa),
        })
    };
    (
        CanonicalCoordinates::new(jet(c(K0, 0.0)), jet(beta)),
        ExteriorWavevectors::from_branches(branch(CLADDING), branch(SUBSTRATE)),
    )
}

/// The characteristic of `backend` at complex `β`.
macro_rules! characteristic {
    ($stack:expr, $backend:expr, $channel:expr) => {{
        let evaluator = ComplexPlaneEvaluator::<RealJ0, _, _>::compile($stack, $backend).unwrap();
        let channel: Channel = $channel;
        move |beta: C| {
            let (coordinates, exterior) = inputs(beta, channel);
            evaluator
                .checked_determinant(&coordinates, &exterior, channel)
                .unwrap()
                .value()[()]
        }
    }};
}

/// The guided roots of `f` on `(lo, hi)`: local minima of `|f|` on a grid,
/// refined by the secant method.
fn roots(f: impl Fn(C) -> C, lo: f64, hi: f64) -> Vec<f64> {
    let steps = 600;
    let at = |k: usize| lo + (hi - lo) * k as f64 / steps as f64;
    let values: Vec<f64> = (0..=steps).map(|k| f(c(at(k), 0.0)).norm()).collect();
    (1..steps)
        .filter(|&k| values[k] < values[k - 1] && values[k] < values[k + 1])
        .filter_map(|k| {
            let (mut a, mut b) = (c(at(k - 1), 0.0), c(at(k + 1), 0.0));
            for _ in 0..80 {
                let (fa, fb) = (f(a), f(b));
                if fb == fa {
                    break;
                }
                let next = b - fb * (b - a) / (fb - fa);
                a = b;
                b = next;
                if (b - a).norm() < 1e-15 * b.norm() {
                    break;
                }
            }
            (b.im.abs() < 1e-9 && b.re > lo && b.re < hi).then_some(b.re)
        })
        .collect()
}

fn retain<B>(
    stack: &AnalyticalTensorStack<C>,
    backend: B,
    beta: f64,
    channel: Channel,
) -> CoupledWorkspace<RealJ0>
where
    B: Backend<RealJ0, ComplexPlane, Workspace = CoupledWorkspace<RealJ0>> + Clone,
{
    let evaluator = ComplexPlaneEvaluator::<RealJ0, _, _>::compile(stack, backend).unwrap();
    let (coordinates, exterior) = inputs(c(beta, 0.0), channel);
    evaluator
        .retain(coordinates, exterior, channel)
        .unwrap()
        .workspace()
        .clone()
}

type Fields = [[C; 3]; 2];

fn sample(solution: &CoupledSolution<'_, RealJ0>, region: StackRegion, depth: f64) -> Fields {
    solution
        .fields(region, &jet(c(depth, 0.0)))
        .unwrap()
        .map(|field| field.map(|entry| entry.value()[()]))
}

const LAYER: StackRegion = StackRegion::Layer(FiniteLayerIndex::new(0));

/// The sample points of the whole line: exterior depths and layer depths.
fn points() -> Vec<(StackRegion, f64)> {
    let mut points: Vec<_> = [0.0, 0.2, 0.5, 1.0]
        .into_iter()
        .map(|depth| (StackRegion::LeftExterior, depth))
        .collect();
    points.extend((0..=8).map(|k| (LAYER, CORE * k as f64 / 8.0)));
    points.extend([0.0, 0.2, 0.5, 1.0].map(|depth| (StackRegion::RightExterior, depth)));
    points
}

/// The whole-line longitudinal power `½ Re ∫ (E_y η_z* − E_z η_y*) dz`.
fn power(solution: &CoupledSolution<'_, RealJ0>) -> f64 {
    [StackRegion::LeftExterior, LAYER, StackRegion::RightExterior]
        .into_iter()
        .map(|region| solution.longitudinal_power(region).unwrap().value()[()])
        .sum()
}

fn tangential(fields: &Fields) -> f64 {
    [fields[0][0], fields[0][1], fields[1][0], fields[1][1]]
        .iter()
        .map(|v| v.norm_sqr())
        .sum::<f64>()
        .sqrt()
}

/// Maxwell's equations at `depth` in the core, by central differences.
fn assert_maxwell(
    solution: &CoupledSolution<'_, RealJ0>,
    stack: &AnalyticalTensorStack<C>,
    beta: f64,
    label: &str,
) {
    let [epsilon, mu] = <RealJ0 as TensorConstitutiveLift<ComplexPlane, _>>::constitutive_tensors(
        LAYER,
        stack.layers_left_to_right()[0].material(),
        &jet(c(K0, 0.0)),
    )
    .unwrap()
    .map(|tensor| tensor.full().map(|row| row.map(|entry| entry.value()[()])));
    let apply = |tensor: [[C; 3]; 3], field: [C; 3]| {
        std::array::from_fn::<C, 3, _>(|row| (0..3).map(|k| tensor[row][k] * field[k]).sum())
    };
    for depth in [0.2, 0.55, 0.9] {
        let h = 1e-4;
        let [e, eta] = sample(solution, LAYER, depth);
        let [e_plus, eta_plus] = sample(solution, LAYER, depth + h);
        let [e_minus, eta_minus] = sample(solution, LAYER, depth - h);
        let dz = |plus: C, minus: C| (plus - minus) / (2.0 * h);
        let (ikx, ik0) = (c(0.0, beta), c(0.0, K0));
        let (mu_eta, epsilon_e) = (apply(mu, eta), apply(epsilon, e));
        let residuals = [
            -dz(e_plus[1], e_minus[1]) - ik0 * mu_eta[0],
            dz(e_plus[0], e_minus[0]) - ikx * e[2] - ik0 * mu_eta[1],
            ikx * e[1] - ik0 * mu_eta[2],
            -dz(eta_plus[1], eta_minus[1]) + ik0 * epsilon_e[0],
            dz(eta_plus[0], eta_minus[0]) - ikx * eta[2] + ik0 * epsilon_e[1],
            ikx * eta[1] + ik0 * epsilon_e[2],
        ];
        let scale = e.iter().chain(&eta).map(|v| v.norm()).fold(0.0, f64::max) * K0;
        for (index, residual) in residuals.iter().enumerate() {
            assert!(
                residual.norm() < 1e-6 * scale,
                "{label} depth {depth} equation {index}: {residual}"
            );
        }
    }
}

#[test]
fn hybrid_mode_fields_satisfy_maxwell_carry_unit_power_and_agree_between_backends() {
    let stack = film(uniaxial(Rotation::from_euler(0.5, 0.6, 0.3)));
    let (lo, hi) = (K0 * SUBSTRATE.sqrt() + 0.05, K0 * 4.4f64.sqrt() - 0.05);
    let found = roots(
        characteristic!(&stack, Transfer4::new(), Channel::Coupled),
        lo,
        hi,
    );
    assert!(found.len() >= 4, "{found:?}");

    for beta in found {
        let workspaces = [
            retain(&stack, Transfer4::new(), beta, Channel::Coupled),
            retain(&stack, Scatter4::new(), beta, Channel::Coupled),
        ];
        let mut normalised = Vec::new();
        for (backend, workspace) in workspaces.iter().enumerate() {
            let label = format!("β {beta} backend {backend}");
            let mode = workspace.mode().unwrap();
            assert!(
                mode.residual() < 1e-10,
                "{label}: residual {}",
                mode.residual()
            );
            let solution = workspace.mode_solution(&mode);
            assert_maxwell(&solution, &stack, beta, &label);

            // Tangential E and η are continuous at both interfaces.
            for (left, right) in [
                (
                    sample(&solution, StackRegion::LeftExterior, 0.0),
                    sample(&solution, LAYER, 0.0),
                ),
                (
                    sample(&solution, LAYER, CORE),
                    sample(&solution, StackRegion::RightExterior, 0.0),
                ),
            ] {
                for (a, b) in [
                    (left[0][0], right[0][0]),
                    (left[0][1], right[0][1]),
                    (left[1][0], right[1][0]),
                    (left[1][1], right[1][1]),
                ] {
                    assert!(
                        (a - b).norm() < 1e-8 * tangential(&left),
                        "{label} continuity"
                    );
                }
            }
            // |ψ| decays into both exteriors, and the mode is hybrid.
            for region in [StackRegion::LeftExterior, StackRegion::RightExterior] {
                let magnitudes: Vec<f64> = [0.0, 0.2, 0.5, 1.0, 2.0]
                    .map(|depth| tangential(&sample(&solution, region, depth)))
                    .to_vec();
                assert!(
                    magnitudes.windows(2).all(|pair| pair[1] < pair[0]),
                    "{label} {region:?} {magnitudes:?}"
                );
            }
            let [e, eta] = sample(&solution, LAYER, 0.37);
            assert!(
                e.iter()
                    .chain(&eta)
                    .all(|v| v.norm() > 1e-3 * tangential(&[e, eta])),
                "{label} hybrid"
            );

            // Unit power under the stated normalisation.
            let p = power(&solution);
            assert!(p > 0.0, "{label} forward power");
            let scaled = solution.scaled(&jet(c(p.sqrt().recip(), 0.0)));
            assert_close(c(power(&scaled), 0.0), c(1.0, 0.0), 1e-12, &label);
            normalised.push(
                points()
                    .into_iter()
                    .map(|(region, depth)| sample(&scaled, region, depth))
                    .collect::<Vec<_>>(),
            );
        }
        // Transfer4 and Scatter4 agree up to a phase.
        let (transfer, scatter) = (&normalised[0], &normalised[1]);
        let reference = transfer[5][0][1];
        let phase = scatter[5][0][1] / reference;
        assert_close(c(phase.norm(), 0.0), c(1.0, 0.0), 1e-9, "phase");
        for (t, s) in transfer.iter().zip(scatter) {
            for part in 0..2 {
                for k in 0..3 {
                    assert!(
                        (s[part][k] - phase * t[part][k]).norm() < 1e-8 * tangential(t).max(1e-3),
                        "β {beta} fields"
                    );
                }
            }
        }
    }
}

#[test]
fn separable_modes_equal_the_2x2_modes_without_cross_terms() {
    let stack = film(uniaxial(Rotation::identity()));
    let (lo, hi) = (K0 * SUBSTRATE.sqrt() + 0.05, K0 * 4.4f64.sqrt() - 0.05);
    let sampling = FieldSampling::new()
        .left_exterior(ExteriorSampling::distances_in(
            [0.0, 0.3],
            LengthUnit::Centimetre,
        ))
        .layer(
            0,
            LayerSampling::offsets_in([0.0, 0.25, 0.6, 1.0], LengthUnit::Centimetre),
        )
        .right_exterior(ExteriorSampling::distances_in(
            [0.0, 0.4],
            LengthUnit::Centimetre,
        ));
    let points = [
        (StackRegion::LeftExterior, 0.0),
        (StackRegion::LeftExterior, 0.3),
        (LAYER, 0.0),
        (LAYER, 0.25),
        (LAYER, 0.6),
        (LAYER, 1.0),
        (StackRegion::RightExterior, 0.0),
        (StackRegion::RightExterior, 0.4),
    ];
    for (polarisation, cross) in [
        // TE: E_x, E_z, η_y vanish; TM: E_y, η_x, η_z vanish.
        (Polarisation::TransverseElectric, [(0, 0), (0, 2), (1, 1)]),
        (Polarisation::TransverseMagnetic, [(0, 1), (1, 0), (1, 2)]),
    ] {
        let scalar = roots(
            characteristic!(&stack, Transfer2::new(), polarisation.into()),
            lo,
            hi,
        );
        assert!(scalar.len() >= 2, "{polarisation:?} {scalar:?}");
        for beta in scalar {
            let reference = Scatter2::new()
                .unscaled_fields(
                    &stack,
                    c(K0, 0.0),
                    c(beta, 0.0),
                    [
                        Some(decaying(CLADDING, c(beta, 0.0))),
                        Some(decaying(SUBSTRATE, c(beta, 0.0))),
                    ],
                    polarisation,
                    &sampling,
                    c(1.0, 0.0),
                )
                .unwrap()
                .fields;
            let reference = reference.quantity();
            let two: Vec<Fields> = (0..points.len())
                .map(|index| {
                    let at = |field: &crate::VectorField<C, ndarray::Ix1>| {
                        [field.x()[index], field.y()[index], field.z()[index]]
                    };
                    [at(reference.electric()), at(reference.magnetic())]
                })
                .collect();
            for (channel, backend) in [
                (Channel::from(polarisation), 0),
                (Channel::Coupled, 0),
                (Channel::from(polarisation), 1),
                (Channel::Coupled, 1),
            ] {
                let label = format!("{polarisation:?} β {beta} {channel:?} backend {backend}");
                let workspace = if backend == 0 {
                    retain(&stack, Scatter4::new(), beta, channel)
                } else {
                    retain(&stack, Transfer4::new(), beta, channel)
                };
                let mode = workspace.mode().unwrap();
                assert!(mode.residual() < 1e-10, "{label}: {}", mode.residual());
                let solution = workspace.mode_solution(&mode);
                let four: Vec<Fields> = points
                    .iter()
                    .map(|&(region, depth)| sample(&solution, region, depth))
                    .collect();
                // The scale fixed by the largest 2×2 component at the centre.
                let (part, k) = if polarisation == Polarisation::TransverseElectric {
                    (0, 1)
                } else {
                    (1, 1)
                };
                let scale = two[4][part][k] / four[4][part][k];
                for (index, (f, t)) in four.iter().zip(&two).enumerate() {
                    let size = tangential(t);
                    for part in 0..2 {
                        for k in 0..3 {
                            assert!(
                                (f[part][k] * scale - t[part][k]).norm() < 1e-8 * size,
                                "{label} point {index} ({part}, {k}): {} vs {}",
                                f[part][k] * scale,
                                t[part][k]
                            );
                        }
                    }
                    for (part, k) in cross {
                        assert!(
                            f[part][k].norm() < 1e-12 * tangential(f),
                            "{label} cross term"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn a_two_dimensional_null_space_is_refused() {
    // ε = μ in every region: TE and TM are degenerate at every mode.
    let dual = |value: f64| {
        isotropic(ComplexQuadratic {
            epsilon: [c(value, 0.0), C::default(), C::default()],
            mu: [c(value, 0.0), C::default(), C::default()],
        })
    };
    let stack = Stack::builder(dual(1.0), dual(1.0))
        .layer(dual(2.0), Length::centimetres(CORE))
        .finalise();
    let evaluator =
        ComplexPlaneEvaluator::<RealJ0, _, _>::compile(&stack, Transfer2::new()).unwrap();
    let te = |beta: C| {
        let kappa = jet(c(0.0, 1.0) * (beta * beta - K0 * K0).sqrt());
        let exterior = ExteriorWavevectors::from_branches(
            Some(ExteriorBranch::Scalar(kappa.clone())),
            Some(ExteriorBranch::Scalar(kappa)),
        );
        let coordinates = CanonicalCoordinates::new(jet(c(K0, 0.0)), jet(beta));
        evaluator
            .checked_determinant(&coordinates, &exterior, Polarisation::TransverseElectric)
            .unwrap()
            .value()[()]
    };
    let beta = roots(te, K0 + 0.05, 2.0 * K0 - 0.05)[0];
    let kappa = jet(c(0.0, 1.0) * (beta * beta - K0 * K0).sqrt());
    let exterior = ExteriorWavevectors::from_branches(
        Some(ExteriorBranch::Pair([kappa.clone(), kappa.clone()])),
        Some(ExteriorBranch::Pair([kappa.clone(), kappa])),
    );
    let coordinates = CanonicalCoordinates::new(jet(c(K0, 0.0)), jet(c(beta, 0.0)));
    for workspace in [
        ComplexPlaneEvaluator::<RealJ0, _, _>::compile(&stack, Transfer4::new())
            .unwrap()
            .retain(coordinates.clone(), exterior.clone(), Channel::Coupled)
            .unwrap()
            .workspace()
            .clone(),
        ComplexPlaneEvaluator::<RealJ0, _, _>::compile(&stack, Scatter4::new())
            .unwrap()
            .retain(coordinates.clone(), exterior.clone(), Channel::Coupled)
            .unwrap()
            .workspace()
            .clone(),
    ] {
        assert_eq!(
            workspace.mode().unwrap_err(),
            ModeReconstructionError::DegenerateNullSpace
        );
        assert_eq!(workspace.retained_layer_count(), Some(1));
    }
}
