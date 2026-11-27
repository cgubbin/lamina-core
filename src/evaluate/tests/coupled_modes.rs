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

/// `R_x(θ) R₀`, with `R_x` the rotation about the stack's `x` axis.
fn rotated(theta: f64) -> Rotation<f64> {
    let base = *Rotation::from_euler(0.5, 0.6, 0.3).matrix();
    let (s, co) = theta.sin_cos();
    let about_x = [[1.0, 0.0, 0.0], [0.0, co, -s], [0.0, s, co]];
    Rotation::from_matrix(std::array::from_fn(|row| {
        std::array::from_fn(|column| (0..3).map(|k| about_x[row][k] * base[k][column]).sum())
    }))
    .unwrap()
}

#[test]
fn rotation_angle_rates_match_central_differences_of_rotated_stacks() {
    use crate::{
        ConstitutiveTensor, ModalCharacteristicDirection, ModalDerivativeBackend, ModeJet1,
    };
    type J1 = ModeJet1<C, ndarray::Ix0>;
    let theta = 0.4;
    let beta = 7.1;
    let constant = |value: C| J1::from_parts(arr0(value), arr0(C::default()));
    let coordinates = CanonicalCoordinates::new(constant(c(K0, 0.0)), constant(c(beta, 0.0)));
    let pair = |epsilon| {
        let kappa = constant(decaying(epsilon, c(beta, 0.0)));
        Some(ExteriorBranch::Pair([kappa.clone(), kappa]))
    };
    let exterior = ExteriorWavevectors::from_branches(pair(CLADDING), pair(SUBSTRATE));
    // dε/dθ = G ε − ε G, G the generator of rotations about x.
    let epsilon = {
        let [epsilon, _] =
            <RealJ0 as TensorConstitutiveLift<ComplexPlane, _>>::constitutive_tensors(
                LAYER,
                film(uniaxial(rotated(theta))).layers_left_to_right()[0].material(),
                &jet(c(K0, 0.0)),
            )
            .unwrap();
        epsilon.full().map(|row| row.map(|entry| entry.value()[()]))
    };
    let generator = [[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]];
    let rate: [[C; 3]; 3] = std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            (0..3)
                .map(|k| {
                    generator[row][k] * epsilon[k][column] - epsilon[row][k] * generator[k][column]
                })
                .sum()
        })
    });
    let zero = || {
        [
            ConstitutiveTensor::Isotropic(C::default()),
            ConstitutiveTensor::Isotropic(C::default()),
        ]
    };
    let rates = [
        zero(),
        [
            ConstitutiveTensor::Full(rate),
            ConstitutiveTensor::Isotropic(C::default()),
        ],
        zero(),
    ];
    let characteristic = |transfer: bool, theta: f64| {
        let stack = film(uniaxial(rotated(theta)));
        let direction = ModalCharacteristicDirection {
            stack: &stack,
            coordinates: &coordinates,
            exterior: &exterior,
            polarisation: Channel::Coupled,
            thickness_rates: &[0.0],
            constitutive_rates: &rates,
        };
        if transfer {
            Transfer4::new()
                .directional_characteristic(direction)
                .unwrap()
        } else {
            Scatter4::new()
                .directional_characteristic(direction)
                .unwrap()
        }
    };
    for transfer in [true, false] {
        let h = 1e-5;
        let difference = (characteristic(transfer, theta + h).value()[()]
            - characteristic(transfer, theta - h).value()[()])
            / (2.0 * h);
        let derivative = characteristic(transfer, theta).first()[()];
        let name = format!("transfer {transfer}");
        assert!(derivative.norm() > 1e-3, "{name} {derivative}");
        assert_close(derivative, difference, 1e-7, &name);
    }

    // The 2×2 backends refuse a full rate; the 4×4 fields along the
    // direction are the mode's fields.
    let stack = film(uniaxial(Rotation::identity()));
    let direction = ModalCharacteristicDirection {
        stack: &stack,
        coordinates: &CanonicalCoordinates::new(constant(c(K0, 0.0)), constant(c(beta, 0.0))),
        exterior: &ExteriorWavevectors::from_branches(
            Some(ExteriorBranch::Scalar(constant(decaying(
                CLADDING,
                c(beta, 0.0),
            )))),
            Some(ExteriorBranch::Scalar(constant(decaying(
                SUBSTRATE,
                c(beta, 0.0),
            )))),
        ),
        polarisation: Polarisation::TransverseElectric,
        thickness_rates: &[0.0],
        constitutive_rates: &rates,
    };
    assert!(
        Scatter2::new()
            .directional_characteristic(direction)
            .is_err()
    );
    let stack = film(uniaxial(rotated(theta)));
    let direction = ModalCharacteristicDirection {
        stack: &stack,
        coordinates: &coordinates,
        exterior: &exterior,
        polarisation: Channel::Coupled,
        thickness_rates: &[0.0],
        constitutive_rates: &rates,
    };
    let sampling = FieldSampling::new().layer(
        0,
        LayerSampling::offsets_in([0.0, 0.5], LengthUnit::Centimetre),
    );
    for fields in [
        Transfer4::new()
            .directional_fields(direction, &sampling, c(1.0, 0.0))
            .unwrap(),
        Scatter4::new()
            .directional_fields(direction, &sampling, c(1.0, 0.0))
            .unwrap(),
    ] {
        let tangent = fields.tangent.quantity();
        assert!(
            tangent
                .electric()
                .x()
                .iter()
                .chain(tangent.magnetic().z())
                .all(|v| v.norm() > 0.0)
        );
    }
}

#[test]
fn bivariate_derivatives_agree_between_backends_and_with_central_differences() {
    use super::coupled::full_stack_with;
    use crate::{IncidentSide, Parameter, RealAxisEvaluator, observable::ProjectJonesAmplitudes};
    let (k0, beta, thickness) = (2.6, 1.9, 0.13);
    let entry = |k: f64, b: f64, d: f64| {
        let solution = RealAxisEvaluator::new(Transfer4::new())
            .evaluate(
                crate::test_support::planar::scalar_real_input(k, b),
                &full_stack_with(d),
                Channel::Coupled,
            )
            .unwrap()
            .into_parts()
            .0;
        let jones = solution
            .entries()
            .project_jones_amplitudes(solution.context(), IncidentSide::Left);
        jones.reflection(0, 1).unwrap().value()[()]
    };
    let layer = Parameter::LayerThickness(FiniteLayerIndex::new(1));
    for (axis0, axis1) in [
        (Parameter::Spectral, Parameter::InPlane),
        (Parameter::InPlane, layer),
        (Parameter::Spectral, layer),
    ] {
        macro_rules! mixed {
            ($backend:expr) => {{
                let solution = RealAxisEvaluator::new($backend)
                    .evaluate_bivariate_second(
                        crate::test_support::planar::scalar_real_input(k0, beta),
                        &full_stack_with(thickness),
                        Channel::Coupled,
                        axis0,
                        axis1,
                    )
                    .unwrap()
                    .into_parts()
                    .0;
                let jones = solution
                    .entries()
                    .project_jones_amplitudes(solution.context(), IncidentSide::Left);
                let r = jones.reflection(0, 1).unwrap().clone();
                [
                    r.value()[()],
                    r.axis0()[()],
                    r.axis1()[()],
                    r.axis0_axis1()[()],
                ]
            }};
        }
        let (transfer, scatter) = (mixed!(Transfer4::new()), mixed!(Scatter4::new()));
        let label = format!("{axis0:?} {axis1:?}");
        for (t, s) in transfer.iter().zip(&scatter) {
            assert_close(*s, *t, 1e-7, &label);
        }
        let h = 1e-4;
        let shifted = |a: f64, b: f64| {
            let mut point = [k0, beta, thickness];
            for (axis, step) in [(axis0, a), (axis1, b)] {
                let slot = match axis {
                    Parameter::Spectral => 0,
                    Parameter::InPlane => 1,
                    _ => 2,
                };
                point[slot] += step;
            }
            entry(point[0], point[1], point[2])
        };
        let difference =
            (shifted(h, h) - shifted(h, -h) - shifted(-h, h) + shifted(-h, -h)) / (4.0 * h * h);
        assert_close(transfer[3], difference, 1e-5, &format!("{label} mixed"));
    }
}

mod forms {
    //! Closed-form modal forms against quadrature of sampled fields, central
    //! differences and the other backends.

    use std::sync::Arc;

    use ndarray::Ix0;

    use super::*;
    use crate::{
        ConstitutiveTensor, ModalCharacteristicDirection, ModalDerivativeBackend, ModalForm,
        ModalFormIntegrals, ModalFormOperand, ModeJet1, Pml, ScalarAlgebra, StretchProfile, Wall,
        WallKind,
    };

    type J = ModeJet1<C, Ix0>;
    type Weight = [[C; 6]; 6];

    fn seed(value: C, rate: C) -> J {
        <J as SeedJet>::constant(arr0(value)).add(
            &<J as SeedJet>::variable(arr0(c(0.0, 0.0)), 0)
                .unwrap()
                .scale(rate),
        )
    }

    /// `η_z* E_y − η_y* E_z`: twice the longitudinal power density.
    fn power() -> Weight {
        let mut w = [[c(0.0, 0.0); 6]; 6];
        w[5][1] = c(1.0, 0.0);
        w[4][2] = c(-1.0, 0.0);
        w
    }

    /// `|E|² + |η|²`, Hermitian; `E·E + η·η`, bilinear.
    fn identity() -> Weight {
        std::array::from_fn(|i| std::array::from_fn(|j| c(f64::from(u8::from(i == j)), 0.0)))
    }

    /// `E_y η_z − E_z η_y`, bilinear.
    fn pairing() -> Weight {
        let mut w = [[c(0.0, 0.0); 6]; 6];
        w[1][5] = c(1.0, 0.0);
        w[2][4] = c(-1.0, 0.0);
        w
    }

    const WEIGHTS: [(bool, fn() -> Weight); 4] = [
        (true, power),
        (true, identity),
        (false, pairing),
        (false, identity),
    ];

    /// A mode's direction inputs: `β` and its rate along the direction, and
    /// each open side's decaying wavenumber `κ = i√(β² − ε k₀²)` with its
    /// rate `−β β'/κ`.
    struct Inputs {
        coordinates: CanonicalCoordinates<J>,
        exterior: ExteriorWavevectors<J>,
        thickness: Vec<f64>,
        rates: Vec<[ConstitutiveTensor<C>; 2]>,
    }

    fn form_inputs(
        beta: f64,
        beta_rate: f64,
        thickness_rate: f64,
        channel: Channel,
        open: [bool; 2],
    ) -> Inputs {
        let branch = |epsilon: f64, open: bool| {
            open.then(|| {
                let kappa = decaying(epsilon, c(beta, 0.0));
                let kappa = seed(kappa, -c(beta * beta_rate, 0.0) / kappa);
                match channel {
                    Channel::Coupled => ExteriorBranch::Pair([kappa.clone(), kappa]),
                    _ => ExteriorBranch::Scalar(kappa),
                }
            })
        };
        let zero = ConstitutiveTensor::Isotropic(c(0.0, 0.0));
        Inputs {
            coordinates: CanonicalCoordinates::new(
                seed(c(K0, 0.0), c(0.0, 0.0)),
                seed(c(beta, 0.0), c(beta_rate, 0.0)),
            ),
            exterior: ExteriorWavevectors::from_branches(
                branch(CLADDING, open[0]),
                branch(SUBSTRATE, open[1]),
            ),
            thickness: vec![thickness_rate],
            rates: vec![[zero.clone(), zero]; 3],
        }
    }

    impl Inputs {
        fn direction<'a>(
            &'a self,
            stack: &'a AnalyticalTensorStack<C>,
            channel: Channel,
        ) -> ModalCharacteristicDirection<
            'a,
            crate::AnalyticalMaterialRegion<C>,
            [ConstitutiveTensor<C>; 2],
            Channel,
        > {
            ModalCharacteristicDirection {
                stack,
                coordinates: &self.coordinates,
                exterior: &self.exterior,
                polarisation: channel,
                thickness_rates: &self.thickness,
                constitutive_rates: &self.rates,
            }
        }
    }

    fn evaluate<B: ModalDerivativeBackend>(
        backend: &B,
        stack: &AnalyticalTensorStack<C>,
        inputs: &Inputs,
        channel: Channel,
        weights: &[Weight],
    ) -> Result<ModalFormIntegrals, Arc<dyn std::error::Error + Send + Sync>> {
        let operand = ModalFormOperand {
            direction: inputs.direction(stack, channel),
            amplitude: c(1.0, 0.0),
        };
        let forms: Vec<ModalForm<'_>> = WEIGHTS
            .iter()
            .zip(weights)
            .map(|((conjugate, _), weight)| ModalForm {
                conjugate: *conjugate,
                weights: std::slice::from_ref(weight),
            })
            .collect();
        backend.directional_forms([operand, operand], &forms)
    }

    /// Five-point Gauss–Legendre on `panels` uniform panels of `[0, length]`:
    /// nodes and weights, symmetric under reversal.
    fn gauss(length: f64, panels: usize) -> (Vec<f64>, Vec<f64>) {
        let rule = [
            (-0.906_179_845_938_664, 0.236_926_885_056_189),
            (-0.538_469_310_105_683, 0.478_628_670_499_366),
            (0.0, 0.568_888_888_888_889),
            (0.538_469_310_105_683, 0.478_628_670_499_366),
            (0.906_179_845_938_664, 0.236_926_885_056_189),
        ];
        let h = length / panels as f64;
        (0..panels)
            .flat_map(|panel| {
                rule.map(|(node, weight)| (h * (panel as f64 + 0.5 + 0.5 * node), 0.5 * h * weight))
            })
            .unzip()
    }

    /// The forms of [`WEIGHTS`] by quadrature of `directional_fields`, region
    /// by region; exteriors to `reach` cm or the wall.
    fn quadrature<B: ModalDerivativeBackend>(
        backend: &B,
        stack: &AnalyticalTensorStack<C>,
        inputs: &Inputs,
        channel: Channel,
        reach: [f64; 2],
    ) -> Vec<Vec<C>> {
        let direction = inputs.direction(stack, channel);
        let regions = [(reach[0], 4000usize), (CORE, 400), (reach[1], 4000)];
        let mut values = vec![vec![]; WEIGHTS.len()];
        for (index, (length, panels)) in regions.into_iter().enumerate() {
            let (nodes, weights) = gauss(length, panels);
            let sampling = match index {
                0 => FieldSampling::new().left_exterior(ExteriorSampling::distances_in(
                    nodes,
                    LengthUnit::Centimetre,
                )),
                1 => FieldSampling::new()
                    .layer(0, LayerSampling::offsets_in(nodes, LengthUnit::Centimetre)),
                _ => FieldSampling::new().right_exterior(ExteriorSampling::distances_in(
                    nodes,
                    LengthUnit::Centimetre,
                )),
            };
            let fields = backend
                .directional_fields(direction, &sampling, c(1.0, 0.0))
                .unwrap()
                .fields;
            let (e, h) = (fields.quantity().electric(), fields.quantity().magnetic());
            let f = |i: usize| [e.x()[i], e.y()[i], e.z()[i], h.x()[i], h.y()[i], h.z()[i]];
            for (form, (conjugate, weight)) in WEIGHTS.iter().enumerate() {
                let weight = weight();
                let mut sum = c(0.0, 0.0);
                for (i, w) in weights.iter().enumerate() {
                    let f = f(i);
                    for (row, entries) in weight.iter().enumerate() {
                        let left = if *conjugate { f[row].conj() } else { f[row] };
                        for (column, entry) in entries.iter().enumerate() {
                            sum += left * entry * f[column] * *w;
                        }
                    }
                }
                values[form].push(sum);
            }
        }
        values
    }

    fn assert_forms_match(closed: &ModalFormIntegrals, sampled: &[Vec<C>], label: &str) {
        for (form, (closed, sampled)) in closed.iter().zip(sampled).enumerate() {
            let scale = sampled.iter().map(|v| v.norm()).fold(0.0, f64::max);
            for (region, (closed, sampled)) in closed.iter().zip(sampled).enumerate() {
                let closed = closed.value()[()];
                assert!(
                    (closed - sampled).norm() <= 1e-10 * scale,
                    "{label} form {form} region {region}: {closed} vs {sampled}"
                );
            }
        }
    }

    fn weights() -> Vec<Weight> {
        WEIGHTS.iter().map(|(_, weight)| weight()).collect()
    }

    #[test]
    fn hybrid_forms_match_quadrature_of_sampled_fields() {
        let stack = film(uniaxial(Rotation::from_euler(0.5, 0.6, 0.3)));
        let (lo, hi) = (K0 * SUBSTRATE.sqrt() + 0.05, K0 * 4.4f64.sqrt() - 0.05);
        let found = roots(
            characteristic!(&stack, Transfer4::new(), Channel::Coupled),
            lo,
            hi,
        );
        assert!(found.len() >= 4, "{found:?}");
        for beta in [found[0], found[found.len() - 1]] {
            let mut ratios = vec![];
            let inputs = form_inputs(beta, 0.0, 0.0, Channel::Coupled, [true, true]);
            // Out to where |f|² has decayed by e^{-50}.
            let reach =
                [CLADDING, SUBSTRATE].map(|epsilon| 25.0 / decaying(epsilon, c(beta, 0.0)).im);
            let sampled = quadrature(&Scatter4::new(), &stack, &inputs, Channel::Coupled, reach);
            for (name, closed) in [
                (
                    "Scatter4",
                    evaluate(
                        &Scatter4::new(),
                        &stack,
                        &inputs,
                        Channel::Coupled,
                        &weights(),
                    ),
                ),
                (
                    "Transfer4",
                    evaluate(
                        &Transfer4::new(),
                        &stack,
                        &inputs,
                        Channel::Coupled,
                        &weights(),
                    ),
                ),
            ] {
                let closed = closed.unwrap();
                if name == "Scatter4" {
                    assert_forms_match(&closed, &sampled, &format!("β {beta} {name}"));
                }
                // The Hermitian norm is real; power over it is free of the
                // gauge and agrees between the backends.
                let total = |form: usize| -> C { closed[form].iter().map(|v| v.value()[()]).sum() };
                assert!(
                    total(1).im.abs() < 1e-12 * total(1).re,
                    "{name} {}",
                    total(1)
                );
                ratios.push(total(0).re / total(1).re);
            }
            assert!(
                (ratios[0] - ratios[1]).abs() < 1e-10 * ratios[0].abs(),
                "{ratios:?}"
            );
        }
    }

    #[test]
    fn scalar_forms_keep_the_2x2_gauge_and_match_quadrature() {
        let stack = film(uniaxial(Rotation::identity()));
        let (lo, hi) = (K0 * SUBSTRATE.sqrt() + 0.05, K0 * 4.4f64.sqrt() - 0.05);
        for polarisation in [
            Polarisation::TransverseElectric,
            Polarisation::TransverseMagnetic,
        ] {
            let channel = Channel::from(polarisation);
            let found = roots(characteristic!(&stack, Transfer2::new(), channel), lo, hi);
            assert!(!found.is_empty());
            let beta = found[0];
            let inputs = form_inputs(beta, 0.0, 0.0, channel, [true, true]);
            let reach =
                [CLADDING, SUBSTRATE].map(|epsilon| 25.0 / decaying(epsilon, c(beta, 0.0)).im);
            for (name, closed, sampled) in [
                (
                    "Scatter2",
                    evaluate(&Scatter2::new(), &stack, &inputs, channel, &weights()).unwrap(),
                    quadrature(&Scatter2::new(), &stack, &inputs, channel, reach),
                ),
                (
                    "Transfer2",
                    evaluate(&Transfer2::new(), &stack, &inputs, channel, &weights()).unwrap(),
                    quadrature(&Transfer2::new(), &stack, &inputs, channel, reach),
                ),
            ] {
                assert_forms_match(&closed, &sampled, &format!("{polarisation:?} {name}"));
            }
            // The same mode through the coupled channel of Scatter4: the
            // gauge-free power-to-energy ratio agrees.
            let coupled = inputs_coupled(beta);
            let ratio = |closed: &ModalFormIntegrals| {
                let total = |form: usize| -> C { closed[form].iter().map(|v| v.value()[()]).sum() };
                total(0) / total(1)
            };
            let two =
                ratio(&evaluate(&Scatter2::new(), &stack, &inputs, channel, &weights()).unwrap());
            let four = ratio(
                &evaluate(
                    &Scatter4::new(),
                    &stack,
                    &coupled,
                    Channel::Coupled,
                    &weights(),
                )
                .unwrap(),
            );
            assert!(
                (two - four).norm() < 1e-10 * two.norm(),
                "{polarisation:?} {two} {four}"
            );
        }
    }

    fn inputs_coupled(beta: f64) -> Inputs {
        form_inputs(beta, 0.0, 0.0, Channel::Coupled, [true, true])
    }

    #[test]
    fn forms_over_a_wall_gap_match_quadrature() {
        // A perfect electric wall 0.6 cm beyond the film on the right.
        let walled: AnalyticalTensorStack<C> = Stack::builder(
            isotropic(lossless(CLADDING)),
            isotropic(lossless(SUBSTRATE)),
        )
        .layer(uniaxial(Rotation::identity()), Length::centimetres(CORE))
        .right_termination(Wall::new(
            WallKind::PerfectElectric,
            Length::centimetres(CORE + 0.6),
        ))
        .finalise();
        let channel = Channel::TransverseElectric;
        let characteristic = {
            let evaluator =
                ComplexPlaneEvaluator::<RealJ0, _, _>::compile(&walled, Transfer2::new()).unwrap();
            move |beta: C| {
                let coordinates = CanonicalCoordinates::new(jet(c(K0, 0.0)), jet(beta));
                let exterior = ExteriorWavevectors::from_branches(
                    Some(ExteriorBranch::Scalar(jet(decaying(CLADDING, beta)))),
                    None,
                );
                evaluator
                    .checked_determinant(&coordinates, &exterior, channel)
                    .unwrap()
                    .value()[()]
            }
        };
        let (lo, hi) = (K0 * SUBSTRATE.sqrt() + 0.05, K0 * 4.4f64.sqrt() - 0.05);
        let found = roots(characteristic, lo, hi);
        assert!(!found.is_empty());
        let inputs = form_inputs(found[0], 0.0, 0.0, channel, [true, false]);
        let reach = [25.0 / decaying(CLADDING, c(found[0], 0.0)).im, 0.6];
        for backend in [0, 1] {
            let (closed, sampled) = match backend {
                0 => (
                    evaluate(&Scatter2::new(), &walled, &inputs, channel, &weights()).unwrap(),
                    quadrature(&Scatter2::new(), &walled, &inputs, channel, reach),
                ),
                _ => (
                    evaluate(&Transfer2::new(), &walled, &inputs, channel, &weights()).unwrap(),
                    quadrature(&Transfer2::new(), &walled, &inputs, channel, reach),
                ),
            };
            assert_forms_match(&closed, &sampled, &format!("walled backend {backend}"));
        }
    }

    #[test]
    fn form_rates_along_a_thickness_match_central_differences() {
        // Power over energy is free of the mode's gauge; its rate along the
        // film thickness, with β moving along the branch, against central
        // differences of re-solved modes.
        let rotation = Rotation::from_euler(0.5, 0.6, 0.3);
        let stack_at = |thickness: f64| -> AnalyticalTensorStack<C> {
            Stack::builder(
                isotropic(lossless(CLADDING)),
                isotropic(lossless(SUBSTRATE)),
            )
            .layer(uniaxial(rotation), Length::centimetres(thickness))
            .finalise()
        };
        let (lo, hi) = (K0 * SUBSTRATE.sqrt() + 0.05, K0 * 4.4f64.sqrt() - 0.05);
        let modes = |thickness: f64| {
            roots(
                characteristic!(&stack_at(thickness), Transfer4::new(), Channel::Coupled),
                lo,
                hi,
            )
        };
        let ratio = |closed: &ModalFormIntegrals| {
            let total = |form: usize| -> J {
                closed[form]
                    .iter()
                    .skip(1)
                    .fold(closed[form][0].clone(), |sum, v| sum.add(v))
            };
            total(0).divide(&total(1))
        };
        let base = modes(CORE);
        let h = 1e-5;
        let (plus, minus) = (modes(CORE + h), modes(CORE - h));
        assert_eq!(base.len(), plus.len());
        let stack = stack_at(CORE);
        for (index, &beta) in base.iter().enumerate() {
            // β' = −F_d / F_β from the characteristic.
            let characteristic = |inputs: &Inputs| {
                Scatter4::new()
                    .directional_characteristic(inputs.direction(&stack, Channel::Coupled))
                    .unwrap()
                    .first()[()]
            };
            let along_beta =
                characteristic(&form_inputs(beta, 1.0, 0.0, Channel::Coupled, [true, true]));
            let along_d =
                characteristic(&form_inputs(beta, 0.0, 1.0, Channel::Coupled, [true, true]));
            let beta_rate = (-along_d / along_beta).re;
            let seeded = form_inputs(beta, beta_rate, 1.0, Channel::Coupled, [true, true]);
            let rate = ratio(
                &evaluate(
                    &Scatter4::new(),
                    &stack,
                    &seeded,
                    Channel::Coupled,
                    &weights(),
                )
                .unwrap(),
            )
            .first()[()];
            let at = |thickness: f64, beta: f64| {
                let stack = stack_at(thickness);
                ratio(
                    &evaluate(
                        &Scatter4::new(),
                        &stack,
                        &form_inputs(beta, 0.0, 0.0, Channel::Coupled, [true, true]),
                        Channel::Coupled,
                        &weights(),
                    )
                    .unwrap(),
                )
                .value()[()]
            };
            let difference = (at(CORE + h, plus[index]) - at(CORE - h, minus[index])) / (2.0 * h);
            assert!(
                (rate - difference).norm() < 1e-6 * difference.norm().max(1e-3),
                "mode {index}: {rate} vs {difference}"
            );
        }
    }

    #[test]
    fn a_pml_side_is_refused() {
        let stack: AnalyticalTensorStack<C> = Stack::builder(
            isotropic(lossless(CLADDING)),
            isotropic(lossless(SUBSTRATE)),
        )
        .layer(
            uniaxial(Rotation::from_euler(0.5, 0.6, 0.3)),
            Length::centimetres(CORE),
        )
        .right_termination(Pml::new(
            Length::centimetres(CORE + 0.3),
            Length::centimetres(0.5),
            StretchProfile::polynomial(2.0, c(0.5, 3.0)),
            WallKind::PerfectElectric,
        ))
        .finalise();
        let inputs = form_inputs(K0 * 1.6, 0.0, 0.0, Channel::Coupled, [true, false]);
        assert!(
            evaluate(
                &Scatter4::new(),
                &stack,
                &inputs,
                Channel::Coupled,
                &weights()
            )
            .is_err()
        );
    }
}

/// The characteristic-only solve behind `determinant` gives the bits of the
/// full solve's characteristic, for value and first-order jets.
#[test]
fn determinant_only_solve_matches_the_full_solve_bit_for_bit() {
    use crate::test_support::jet::HoloJ1;

    let bits = |value: C| [value.re.to_bits(), value.im.to_bits()];
    let stack = film(uniaxial(Rotation::from_euler(0.4, -1.3, 2.0)));
    for beta in [c(4.3, 0.0), c(5.1, 0.02), c(6.7, -0.1)] {
        let (coordinates, exterior) = inputs(beta, Channel::Coupled);
        macro_rules! check {
            ($backend:expr) => {{
                let evaluator =
                    ComplexPlaneEvaluator::<RealJ0, _, _>::compile(&stack, $backend).unwrap();
                let only = evaluator
                    .determinant(&coordinates, &exterior, Channel::Coupled)
                    .unwrap();
                let full = evaluator
                    .solve(&coordinates, &exterior, Channel::Coupled)
                    .unwrap()
                    .determinant();
                assert_eq!(
                    bits(only.value().value()[()]),
                    bits(full.value().value()[()])
                );
            }};
        }
        check!(Scatter4::new());
        check!(Transfer4::new());

        // β seeded: value and first derivative.
        let seeded = |value: C| <HoloJ1 as SeedJet>::constant(arr0(value));
        let coordinates = CanonicalCoordinates::new(
            seeded(c(K0, 0.0)),
            <HoloJ1 as SeedJet>::variable(arr0(beta), 0).unwrap(),
        );
        let kappa = |epsilon| seeded(decaying(epsilon, beta));
        let exterior = ExteriorWavevectors::from_branches(
            Some(ExteriorBranch::Pair([kappa(CLADDING), kappa(CLADDING)])),
            Some(ExteriorBranch::Pair([kappa(SUBSTRATE), kappa(SUBSTRATE)])),
        );
        macro_rules! check_first {
            ($backend:expr) => {{
                let evaluator =
                    ComplexPlaneEvaluator::<HoloJ1, _, _>::compile(&stack, $backend).unwrap();
                let only = evaluator
                    .determinant(&coordinates, &exterior, Channel::Coupled)
                    .unwrap();
                let full = evaluator
                    .solve(&coordinates, &exterior, Channel::Coupled)
                    .unwrap()
                    .determinant();
                let (only, full) = (only.value(), full.value());
                assert_eq!(bits(only.value()[()]), bits(full.value()[()]));
                assert_eq!(bits(only.first()[()]), bits(full.first()[()]));
            }};
        }
        check_first!(Scatter4::new());
        check_first!(Transfer4::new());
    }
}
