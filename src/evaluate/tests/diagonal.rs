//! Stacks whose permittivity and permeability are diagonal in the stack frame,
//! through both 2×2 backends.

use std::f64::consts::FRAC_PI_2;

use lamina_units::Length;
use ndarray::{Array1, arr0};

use crate::{
    AnalyticalMaterialHandle, AnalyticalMaterialRegion, AnalyticalTensorStack,
    CanonicalCoordinates, Channel, ComplexPlane, ComplexPlaneEvaluator, ConstitutiveClass,
    ConstitutiveError, ExteriorBranch, ExteriorSampling, ExteriorWavevectors, FieldSampling,
    FiniteLayerIndex, IncidentSide, LayerSampling, Parameter, Pml, Polarisation,
    RealAxisEvaluationError, RealAxisEvaluator, Rotation, Scatter2, SeedJet, Stack, StackRegion,
    StretchProfile, Transfer2, Wall, WallKind,
    backend::evaluate_stack_exterior_wavevectors,
    test_support::{
        jet::{HoloJ0, HoloJ1, HoloJB2},
        materials::{ComplexQuadratic, principal_region},
        planar::{principal_exterior_wavevectors, sampled_real_input, scalar_real_input},
    },
};

type C = num_complex::Complex64;
type Region = AnalyticalMaterialRegion<C>;

const POLARISATIONS: [Polarisation; 2] = [
    Polarisation::TransverseElectric,
    Polarisation::TransverseMagnetic,
];

fn c(re: f64, im: f64) -> C {
    C::new(re, im)
}

/// `|actual − expected| ≤ tolerance · max(1, |expected|)`.
fn assert_close(actual: C, expected: C, tolerance: f64, label: &str) {
    let error = (actual - expected).norm();
    assert!(
        error <= tolerance * expected.norm().max(1.0),
        "{label}: expected {expected}, got {actual}, error {error:e}"
    );
}

fn assert_arrays_close(actual: &Array1<C>, expected: &Array1<C>, tolerance: f64, label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}");
    for (index, (&a, &e)) in actual.iter().zip(expected).enumerate() {
        assert_close(a, e, tolerance, &format!("{label}[{index}]"));
    }
}

/// `ε = a₀ + a₁ k₀ + a₂ k₀²`, `μ = b₀ + b₁ k₀`: lossy, dispersive and magnetic.
fn dispersive(epsilon: C, mu: C) -> ComplexQuadratic {
    ComplexQuadratic {
        epsilon: [epsilon, c(0.04, 0.01), c(-0.003, 0.002)],
        mu: [mu, c(-0.02, 0.005), c(0.0, 0.0)],
    }
}

fn constant(epsilon: C, mu: C) -> ComplexQuadratic {
    ComplexQuadratic {
        epsilon: [epsilon, C::default(), C::default()],
        mu: [mu, C::default(), C::default()],
    }
}

fn isotropic(material: ComplexQuadratic) -> Region {
    AnalyticalMaterialHandle::new(material).into()
}

/// Principal axes `materials` along stack x, y and z.
fn diagonal(materials: [ComplexQuadratic; 3]) -> Region {
    principal_region(materials, Rotation::identity())
}

/// A uniaxial region with its optic axis along stack axis `axis` (x, y, z),
/// reached through the optic-axis rotation of a crystal with `c` along z.
fn uniaxial(ordinary: ComplexQuadratic, extraordinary: ComplexQuadratic, axis: usize) -> Region {
    let (theta, phi) = [(FRAC_PI_2, 0.0), (FRAC_PI_2, FRAC_PI_2), (0.0, 0.0)][axis];
    principal_region(
        [ordinary, ordinary, extraordinary],
        Rotation::from_optic_axis(theta, phi),
    )
}

/// Stack-frame diagonals `(ε, μ)` of the constant regions built above.
fn uniaxial_entries(ordinary: f64, extraordinary: f64, axis: usize) -> [C; 3] {
    let mut entries = [c(ordinary, 0.0); 3];
    entries[axis] = c(extraordinary, 0.0);
    entries
}

/// `(κ, Y)` of one polarisation in a medium with diagonal `ε` and `μ`, from
/// the closed forms.
fn closed_form(epsilon: [C; 3], mu: [C; 3], k0: C, beta: C, polarisation: Polarisation) -> (C, C) {
    let (transverse, canonical, normal) = match polarisation {
        Polarisation::TransverseElectric => (mu[0], epsilon[1], mu[2]),
        Polarisation::TransverseMagnetic => (epsilon[0], mu[1], epsilon[2]),
    };
    let kappa = (transverse * (canonical * k0 * k0 - beta * beta / normal)).sqrt();
    (kappa, kappa / transverse)
}

/// A lossy stack of a biaxial dispersive layer, a uniaxial layer with its
/// optic axis along x and an isotropic layer, between lossless exteriors.
fn biaxial_stack() -> AnalyticalTensorStack<C> {
    Stack::builder(
        isotropic(constant(c(1.0, 0.0), c(1.0, 0.0))),
        diagonal([
            constant(c(1.7, 0.0), c(1.0, 0.0)),
            constant(c(2.0, 0.0), c(1.0, 0.0)),
            constant(c(1.4, 0.0), c(1.0, 0.0)),
        ]),
    )
    .layer(
        diagonal([
            dispersive(c(2.2, 0.15), c(1.1, 0.03)),
            dispersive(c(3.1, 0.05), c(0.9, 0.02)),
            dispersive(c(4.0, 0.3), c(1.3, 0.01)),
        ]),
        Length::centimetres(0.21),
    )
    .layer(
        uniaxial(
            dispersive(c(2.5, 0.1), c(1.0, 0.0)),
            dispersive(c(5.0, 0.2), c(1.0, 0.0)),
            0,
        ),
        Length::centimetres(0.13),
    )
    .layer(
        isotropic(dispersive(c(1.9, 0.08), c(1.05, 0.02))),
        Length::centimetres(0.17),
    )
    .finalise()
}

macro_rules! for_each_backend {
    ($evaluator:ident, $body:block) => {{
        {
            let $evaluator = RealAxisEvaluator::new(Scatter2::new());
            $body
        }
        {
            let $evaluator = RealAxisEvaluator::new(Transfer2::new());
            $body
        }
    }};
}

#[test]
fn uniaxial_interface_reflection_matches_admittance_fresnel_coefficients() {
    let (k0, beta) = (2.0, 0.9);
    let mut cases: Vec<(Region, [C; 3], [C; 3])> = (0..3)
        .map(|axis| {
            (
                uniaxial(
                    constant(c(2.25, 0.0), c(1.0, 0.0)),
                    constant(c(4.0, 0.0), c(1.0, 0.0)),
                    axis,
                ),
                uniaxial_entries(2.25, 4.0, axis),
                [c(1.0, 0.0); 3],
            )
        })
        .collect();
    let (epsilon, mu) = (
        [c(2.25, 0.0), c(3.0, 0.0), c(4.0, 0.0)],
        [c(1.2, 0.0), c(0.9, 0.0), c(1.1, 0.0)],
    );
    cases.push((
        diagonal([0, 1, 2].map(|i| constant(epsilon[i], mu[i]))),
        epsilon,
        mu,
    ));

    for (region, epsilon, mu) in cases {
        let stack =
            Stack::builder(isotropic(constant(c(1.0, 0.0), c(1.0, 0.0))), region).finalise();
        for polarisation in POLARISATIONS {
            let (_, incident) = closed_form(
                [c(1.0, 0.0); 3],
                [c(1.0, 0.0); 3],
                c(k0, 0.0),
                c(beta, 0.0),
                polarisation,
            );
            let (_, transmitted) = closed_form(epsilon, mu, c(k0, 0.0), c(beta, 0.0), polarisation);
            let reflection = (incident - transmitted) / (incident + transmitted);
            let transmission = 2.0 * incident / (incident + transmitted);

            for_each_backend!(evaluator, {
                let result = evaluator
                    .evaluate(scalar_real_input(k0, beta), &stack, polarisation)
                    .unwrap();
                let amplitudes = result.amplitudes(IncidentSide::Left).unwrap();
                assert_close(amplitudes.reflection()[()], reflection, 1e-13, "reflection");
                assert_close(
                    amplitudes.transmission()[()],
                    transmission,
                    1e-13,
                    "transmission",
                );
                let power = result.power(IncidentSide::Left).unwrap();
                assert!((power.reflectance()[()] + power.transmittance()[()] - 1.0).abs() < 1e-13);
            });
        }
    }
}

#[test]
fn birefringent_cladding_has_polarisation_dependent_light_lines() {
    // Optic axis along z: TE vanishes at k∥² = ε_yy μ_zz k0², TM at ε_zz μ_yy k0².
    let k0 = 1.7;
    let stack = Stack::builder(
        uniaxial(
            constant(c(2.0, 0.0), c(1.0, 0.0)),
            constant(c(3.0, 0.0), c(1.0, 0.0)),
            2,
        ),
        isotropic(constant(c(1.0, 0.0), c(1.0, 0.0))),
    )
    .finalise();
    let canonical = crate::input::compile_canonical_constant_stack::<_, HoloJ0>(
        &stack,
        ndarray::Ix0(),
        &crate::ValidationConfig::permissive(),
    )
    .unwrap();
    let kappa = |beta: f64, polarisation| {
        let coordinates = CanonicalCoordinates::new(
            HoloJ0::constant(arr0(c(k0, 0.0))),
            HoloJ0::constant(arr0(c(beta, 0.0))),
        );
        evaluate_stack_exterior_wavevectors::<ComplexPlane, _, _>(
            &coordinates,
            &canonical,
            Channel::from(polarisation),
        )
        .unwrap()
        .left()
        .and_then(ExteriorBranch::scalar)
        .unwrap()
        .value()[()]
    };
    let te_line = k0 * 2.0_f64.sqrt();
    let tm_line = k0 * 3.0_f64.sqrt();

    assert!(kappa(te_line, Polarisation::TransverseElectric).norm() < 1e-7);
    assert!(kappa(te_line, Polarisation::TransverseMagnetic).norm() > 0.5);
    assert!(kappa(tm_line, Polarisation::TransverseMagnetic).norm() < 1e-7);
    // Below its line a wave propagates with κ > 0; above, Im κ > 0.
    for (polarisation, line) in [
        (Polarisation::TransverseElectric, te_line),
        (Polarisation::TransverseMagnetic, tm_line),
    ] {
        let below = kappa(0.9 * line, polarisation);
        let above = kappa(1.1 * line, polarisation);
        assert!(
            below.re > 0.0 && below.im.abs() < 1e-14,
            "{polarisation:?} {below}"
        );
        assert!(
            above.im > 0.0 && above.re.abs() < 1e-14,
            "{polarisation:?} {above}"
        );
    }
}

#[test]
fn isotropic_tensor_regions_reproduce_scalar_regions() {
    let materials = [
        dispersive(c(1.0, 0.0), c(1.0, 0.0)),
        dispersive(c(2.2, 0.15), c(1.1, 0.03)),
        dispersive(c(3.5, 0.02), c(0.95, 0.01)),
        dispersive(c(1.8, 0.0), c(1.2, 0.0)),
    ];
    let build = |region: &dyn Fn(ComplexQuadratic) -> Region| {
        Stack::builder(region(materials[0]), region(materials[3]))
            .layer(region(materials[1]), Length::centimetres(0.23))
            .layer(region(materials[2]), Length::centimetres(0.11))
            .finalise()
    };
    let scalar = build(&isotropic);
    let tensor = build(&|material| diagonal([material; 3]));
    let sampling = FieldSampling::new()
        .left_exterior(ExteriorSampling::uniform(2, Length::centimetres(0.1)))
        .layer(0, LayerSampling::uniform_with_boundaries(5, true, true))
        .layer(1, LayerSampling::uniform(4))
        .right_exterior(ExteriorSampling::uniform(2, Length::centimetres(0.1)));
    let tolerance = 1e-12;

    for polarisation in POLARISATIONS {
        for_each_backend!(evaluator, {
            let retain = |stack| {
                evaluator
                    .retain_first(
                        scalar_real_input(2.3, 1.1),
                        stack,
                        polarisation,
                        Parameter::Spectral,
                    )
                    .unwrap()
            };
            let (scalar_state, tensor_state) = (retain(&scalar), retain(&tensor));
            let (scalar_point, tensor_point) = (
                scalar_state.project_point(&()).unwrap(),
                tensor_state.project_point(&()).unwrap(),
            );
            for side in [IncidentSide::Left, IncidentSide::Right] {
                let (s, t) = (
                    scalar_point.excitation(side).unwrap(),
                    tensor_point.excitation(side).unwrap(),
                );
                let complex = |a: C, b: C, label: &str| assert_close(a, b, tolerance, label);
                let real = |a: f64, b: f64, label: &str| {
                    assert_close(c(a, 0.0), c(b, 0.0), tolerance, label)
                };
                let (sa, ta) = (s.amplitudes(), t.amplitudes());
                complex(
                    ta.value().reflection()[()],
                    sa.value().reflection()[()],
                    "r",
                );
                complex(
                    ta.first().reflection()[()],
                    sa.first().reflection()[()],
                    "dr",
                );
                complex(
                    ta.value().transmission()[()],
                    sa.value().transmission()[()],
                    "t",
                );

                let (sf, tf) = (
                    s.evaluate_fields(&sampling).unwrap(),
                    t.evaluate_fields(&sampling).unwrap(),
                );
                let (sd, td) = (
                    s.evaluate_constitutive_fields(&sampling).unwrap(),
                    t.evaluate_constitutive_fields(&sampling).unwrap(),
                );
                for (label, a, b) in [
                    (
                        "E",
                        sf.quantity().value().electric(),
                        tf.quantity().value().electric(),
                    ),
                    (
                        "H",
                        sf.quantity().value().magnetic(),
                        tf.quantity().value().magnetic(),
                    ),
                    (
                        "D",
                        sd.quantity().value().electric_displacement(),
                        td.quantity().value().electric_displacement(),
                    ),
                    (
                        "B",
                        sd.quantity().value().magnetic_induction(),
                        td.quantity().value().magnetic_induction(),
                    ),
                ] {
                    assert_arrays_close(b.x(), a.x(), tolerance, label);
                    assert_arrays_close(b.y(), a.y(), tolerance, label);
                    assert_arrays_close(b.z(), a.z(), tolerance, label);
                }

                macro_rules! densities {
                    ($e:expr) => {
                        (
                            $e.evaluate_energy_density(&sampling)
                                .unwrap()
                                .quantity()
                                .value()
                                .total()
                                .clone(),
                            $e.evaluate_dissipation_density(&sampling)
                                .unwrap()
                                .quantity()
                                .value()
                                .total()
                                .clone(),
                        )
                    };
                }
                let ((se, sq), (te, tq)) = (densities!(s), densities!(t));
                for (a, b) in se.iter().zip(&te).chain(sq.iter().zip(&tq)) {
                    real(*b, *a, "density");
                }

                let (sl, tl) = (
                    s.layer_energy_dispersive().unwrap(),
                    t.layer_energy_dispersive().unwrap(),
                );
                let (sx, tx) = (
                    s.layer_dissipation().unwrap(),
                    t.layer_dissipation().unwrap(),
                );
                for index in 0..2 {
                    let index = FiniteLayerIndex::new(index);
                    let (a, b) = (
                        sl.value().get(index).unwrap(),
                        tl.value().get(index).unwrap(),
                    );
                    real(b.electric()[()], a.electric()[()], "electric energy");
                    real(b.magnetic()[()], a.magnetic()[()], "magnetic energy");
                    let (a, b) = (
                        sl.first().get(index).unwrap(),
                        tl.first().get(index).unwrap(),
                    );
                    real(b.total()[()], a.total()[()], "energy gradient");
                    let (a, b) = (
                        sx.value().get(index).unwrap(),
                        tx.value().get(index).unwrap(),
                    );
                    real(b.electric()[()], a.electric()[()], "electric dissipation");
                    real(b.magnetic()[()], a.magnetic()[()], "magnetic dissipation");
                }
            }
        });
    }

    // The characteristic and its k0 derivative in the complex plane.
    for polarisation in POLARISATIONS {
        let coordinates = CanonicalCoordinates::new(
            HoloJ1::variable(arr0(c(2.3, -0.1))),
            HoloJ1::constant(arr0(c(1.4, 0.2))),
        );
        let determinant = |stack: &AnalyticalTensorStack<C>| {
            let evaluator =
                ComplexPlaneEvaluator::<HoloJ1, _, _>::compile(stack, Scatter2::new()).unwrap();
            let exterior = principal_exterior_wavevectors(evaluator.stack(), &coordinates);
            evaluator
                .checked_determinant(&coordinates, &exterior, polarisation)
                .unwrap()
                .into_inner()
        };
        let (a, b) = (determinant(&scalar), determinant(&tensor));
        assert_close(b.value()[()], a.value()[()], tolerance, "characteristic");
        assert_close(
            b.first()[()],
            a.first()[()],
            tolerance,
            "characteristic k0 derivative",
        );
    }
}

#[test]
fn uniaxial_slab_guided_modes_match_closed_form_conditions() {
    // Symmetric slab of half-thickness a in an isotropic cladding. Even modes
    // satisfy q tan(qa)/p = α/p_c, odd modes −q cot(qa)/p = α/p_c, with
    // q² = p (w k0² − β²/r) and α² = β² − ε_c k0²; (p, w, r) = (μ_xx, ε_yy,
    // μ_zz) for TE and (ε_xx, μ_yy, ε_zz) for TM, p_c = μ_c or ε_c.
    let (k0, cladding, ordinary, extraordinary, half) = (3.0, 1.5, 2.25, 4.0, 0.8);
    let mut worst: f64 = 0.0;
    let mut modes = 0;
    for axis in 0..3 {
        let epsilon = uniaxial_entries(ordinary, extraordinary, axis).map(|entry| entry.re);
        let clad = || isotropic(constant(c(cladding, 0.0), c(1.0, 0.0)));
        let stack = Stack::builder(clad(), clad())
            .layer(
                uniaxial(
                    constant(c(ordinary, 0.0), c(1.0, 0.0)),
                    constant(c(extraordinary, 0.0), c(1.0, 0.0)),
                    axis,
                ),
                Length::centimetres(2.0 * half),
            )
            .finalise();
        for polarisation in POLARISATIONS {
            let (p, w, r, p_c) = match polarisation {
                Polarisation::TransverseElectric => (1.0, epsilon[1], 1.0, 1.0),
                Polarisation::TransverseMagnetic => (epsilon[0], 1.0, epsilon[2], cladding),
            };
            let q = |beta: f64| (p * (w * k0 * k0 - beta * beta / r)).sqrt();
            let alpha = |beta: f64| (beta * beta - cladding * k0 * k0).sqrt();
            let conditions: [&dyn Fn(f64) -> f64; 2] = [
                &|beta| {
                    q(beta) * (q(beta) * half).sin() / p
                        - alpha(beta) / p_c * (q(beta) * half).cos()
                },
                &|beta| {
                    q(beta) * (q(beta) * half).cos() / p
                        + alpha(beta) / p_c * (q(beta) * half).sin()
                },
            ];
            let (lo, hi) = (k0 * cladding.sqrt(), k0 * (w * r).sqrt());
            let grid: Vec<f64> = (1..4000)
                .map(|i| lo + (hi - lo) * i as f64 / 4000.0)
                .collect();
            for condition in conditions {
                for pair in grid.windows(2) {
                    let (mut a, mut b) = (pair[0], pair[1]);
                    if condition(a).signum() == condition(b).signum() {
                        continue;
                    }
                    for _ in 0..200 {
                        let m = 0.5 * (a + b);
                        if condition(m).signum() == condition(a).signum() {
                            a = m;
                        } else {
                            b = m;
                        }
                    }
                    let beta = 0.5 * (a + b);
                    modes += 1;

                    let coordinates = CanonicalCoordinates::new(
                        HoloJ1::constant(arr0(c(k0, 0.0))),
                        HoloJ1::variable(arr0(c(beta, 0.0))),
                    );
                    macro_rules! newton_step {
                        ($backend:expr) => {{
                            let evaluator =
                                ComplexPlaneEvaluator::<HoloJ1, _, _>::compile(&stack, $backend)
                                    .unwrap();
                            let exterior =
                                principal_exterior_wavevectors(evaluator.stack(), &coordinates);
                            let determinant = evaluator
                                .checked_determinant(&coordinates, &exterior, polarisation)
                                .unwrap()
                                .into_inner();
                            (determinant.value()[()] / determinant.first()[()]).norm() / beta
                        }};
                    }
                    for error in [
                        newton_step!(Transfer2::new()),
                        newton_step!(Scatter2::new()),
                    ] {
                        assert!(
                            error < 1e-10,
                            "axis {axis} {polarisation:?} β = {beta}: {error:e}"
                        );
                        worst = worst.max(error);
                    }
                }
            }
        }
    }
    assert!(modes >= 12, "{modes} modes, worst {worst:e}");
    println!("uniaxial slab: {modes} modes, worst relative β error {worst:e}");
}

#[test]
fn driven_fields_satisfy_maxwell_and_interface_conditions() {
    // ∇×E = i k0 B and ∇×H = −i k0 D, with ∂x = i k∥, D = εE and B = μH
    // entry by entry (H scaled by Z0).
    let (k0, beta, h) = (2.4, 1.3, 1e-5);
    let stack = biaxial_stack();
    let thickness: Vec<f64> = stack
        .layers_left_to_right()
        .iter()
        .map(|layer| layer.thickness().as_centimetres())
        .collect();
    let triple = |centre: f64| [centre - h, centre, centre + h].map(Length::centimetres);
    let mut interior = FieldSampling::new()
        .left_exterior(ExteriorSampling::distances(triple(0.1)))
        .right_exterior(ExteriorSampling::distances(triple(0.1)));
    for (index, d) in thickness.iter().enumerate() {
        interior = interior.layer(index, LayerSampling::offsets(triple(0.4 * d)));
    }
    let mut boundaries =
        FieldSampling::new().left_exterior(ExteriorSampling::point(Length::zero()));
    for index in 0..thickness.len() {
        boundaries = boundaries.layer(index, LayerSampling::boundaries());
    }
    boundaries = boundaries.right_exterior(ExteriorSampling::point(Length::zero()));
    let i = c(0.0, 1.0);

    for polarisation in POLARISATIONS {
        for_each_backend!(evaluator, {
            let state = evaluator
                .retain(scalar_real_input(k0, beta), &stack, polarisation)
                .unwrap();
            let point = state.project_point(&()).unwrap();
            for side in [IncidentSide::Left, IncidentSide::Right] {
                let excitation = point.excitation(side).unwrap();
                let fields = excitation.evaluate_fields(&interior).unwrap();
                let constitutive = excitation.evaluate_constitutive_fields(&interior).unwrap();
                let (e, m) = (fields.quantity().electric(), fields.quantity().magnetic());
                let (d, b) = (
                    constitutive.quantity().electric_displacement(),
                    constitutive.quantity().magnetic_induction(),
                );
                // Triples: left exterior (outward, so z decreases), layers, right exterior.
                for triple in 0..thickness.len() + 2 {
                    let sign = if triple == 0 { -1.0 } else { 1.0 };
                    let n = 3 * triple + 1;
                    let dz = |f: &Array1<C>| (f[n + 1] - f[n - 1]) * sign / (2.0 * h);
                    let scale = k0
                        * (e.x()[n].norm() + e.y()[n].norm() + m.x()[n].norm() + m.y()[n].norm());
                    let residuals = [
                        -dz(e.y()) - i * k0 * b.x()[n],
                        dz(e.x()) - i * beta * e.z()[n] - i * k0 * b.y()[n],
                        i * beta * e.y()[n] - i * k0 * b.z()[n],
                        -dz(m.y()) + i * k0 * d.x()[n],
                        dz(m.x()) - i * beta * m.z()[n] + i * k0 * d.y()[n],
                        i * beta * m.y()[n] + i * k0 * d.z()[n],
                    ];
                    for (component, residual) in residuals.iter().enumerate() {
                        assert!(
                            residual.norm() < 1e-7 * scale,
                            "{polarisation:?} region {triple} component {component}: {residual:e} vs {scale:e}"
                        );
                    }
                }

                // Tangential E, H and normal D, B are continuous across every interface.
                let fields = excitation.evaluate_fields(&boundaries).unwrap();
                let constitutive = excitation
                    .evaluate_constitutive_fields(&boundaries)
                    .unwrap();
                let (e, m) = (fields.quantity().electric(), fields.quantity().magnetic());
                let (d, b) = (
                    constitutive.quantity().electric_displacement(),
                    constitutive.quantity().magnetic_induction(),
                );
                for interface in 0..=thickness.len() {
                    let (a, z) = (2 * interface, 2 * interface + 1);
                    for (label, f) in [
                        ("Ex", e.x()),
                        ("Ey", e.y()),
                        ("Hx", m.x()),
                        ("Hy", m.y()),
                        ("Dz", d.z()),
                        ("Bz", b.z()),
                    ] {
                        assert_close(
                            f[z],
                            f[a],
                            1e-12,
                            &format!("{polarisation:?} interface {interface} {label}"),
                        );
                    }
                }
            }
        });
    }
}

#[test]
fn lossy_diagonal_stack_energy_and_flux_identities() {
    const POINTS: usize = 2001;
    let integrate = |density: &[f64], thickness: f64| {
        let interior: f64 = density[1..density.len() - 1].iter().sum();
        thickness / (density.len() - 1) as f64
            * (0.5 * (density[0] + density[density.len() - 1]) + interior)
    };
    let stack = biaxial_stack();
    let thickness: Vec<f64> = stack
        .layers_left_to_right()
        .iter()
        .map(|layer| layer.thickness().as_centimetres())
        .collect();
    let mut sampling = FieldSampling::new();
    for index in 0..thickness.len() {
        sampling = sampling.layer(
            index,
            LayerSampling::uniform_with_boundaries(POINTS, true, true),
        );
    }

    for polarisation in POLARISATIONS {
        for_each_backend!(evaluator, {
            let state = evaluator
                .retain(scalar_real_input(2.4, 1.3), &stack, polarisation)
                .unwrap();
            let point = state.project_point(&()).unwrap();
            for side in [IncidentSide::Left, IncidentSide::Right] {
                let excitation = point.excitation(side).unwrap();
                let flux_loss = excitation.layer_power().unwrap();
                let dissipation = excitation.layer_dissipation().unwrap();
                let energy = excitation.layer_energy_dispersive().unwrap();
                let dissipation_density =
                    excitation.evaluate_dissipation_density(&sampling).unwrap();
                let energy_density = excitation.evaluate_energy_density(&sampling).unwrap();
                let dissipation_density =
                    dissipation_density.quantity().total().as_slice().unwrap();
                let energy_density = energy_density.quantity();

                let mut absorbed = 0.0;
                for (index, &d) in thickness.iter().enumerate() {
                    let layer = FiniteLayerIndex::new(index);
                    let range = index * POINTS..(index + 1) * POINTS;
                    let loss = flux_loss.get(layer).unwrap().absorbed()[()];
                    absorbed += loss;
                    let label = format!("{polarisation:?} {side:?} layer {index}");

                    let analytic = dissipation.get(layer).unwrap().total()[()];
                    assert!(
                        (analytic - loss).abs() < 1e-12,
                        "{label}: dissipation {analytic} vs flux loss {loss}"
                    );
                    let integrated = integrate(&dissipation_density[range.clone()], d);
                    assert!(
                        (integrated - loss).abs() < 1e-6 * loss.abs().max(1e-3),
                        "{label}: {integrated} vs {loss}"
                    );

                    let layer_energy = energy.get(layer).unwrap();
                    for (density, analytic) in [
                        (energy_density.electric(), layer_energy.electric()[()]),
                        (energy_density.magnetic(), layer_energy.magnetic()[()]),
                    ] {
                        let integrated = integrate(&density.as_slice().unwrap()[range.clone()], d);
                        assert!(
                            (integrated - analytic).abs() < 1e-6 * analytic.abs(),
                            "{label}: energy {integrated} vs {analytic}"
                        );
                    }
                }
                let power = excitation.power();
                assert!(absorbed > 1e-3);
                assert!(
                    (absorbed - power.absorptance()[()]).abs() < 1e-12,
                    "{absorbed} vs {}",
                    power.absorptance()[()]
                );
            }
        });
    }
}

#[test]
fn full_regions_are_refused_by_name() {
    let full = || {
        principal_region(
            [0, 1, 2].map(|i| constant(c(2.0 + i as f64, 0.0), c(1.0, 0.0))),
            Rotation::from_euler(0.3, 0.6, 0.0),
        )
    };
    let vacuum = || isotropic(constant(c(1.0, 0.0), c(1.0, 0.0)));
    let in_layer = Stack::builder(vacuum(), vacuum())
        .layer(full(), Length::centimetres(0.2))
        .finalise();
    let in_exterior = Stack::builder(vacuum(), full()).finalise();
    let unsupported = |region| ConstitutiveError::Unsupported {
        region,
        class: ConstitutiveClass::Full,
    };
    let input = || scalar_real_input(2.0, 0.5);

    for polarisation in POLARISATIONS {
        match RealAxisEvaluator::new(Scatter2::new()).evaluate(input(), &in_layer, polarisation) {
            Err(RealAxisEvaluationError::Backend { source }) => assert_eq!(
                source,
                crate::backend::Scatter2Error::Constitutive(unsupported(StackRegion::Layer(
                    FiniteLayerIndex::new(0)
                )))
            ),
            other => panic!("unexpected {:?}", other.map(|_| ())),
        }
        match RealAxisEvaluator::new(Transfer2::new()).evaluate(input(), &in_layer, polarisation) {
            Err(RealAxisEvaluationError::Backend { source }) => assert_eq!(
                source,
                crate::backend::Transfer2Error::Constitutive(unsupported(StackRegion::Layer(
                    FiniteLayerIndex::new(0)
                )))
            ),
            other => panic!("unexpected {:?}", other.map(|_| ())),
        }
        match RealAxisEvaluator::new(Scatter2::new()).evaluate(input(), &in_exterior, polarisation)
        {
            Err(RealAxisEvaluationError::Constitutive(error)) => {
                assert_eq!(error, unsupported(StackRegion::RightExterior))
            }
            other => panic!("unexpected {:?}", other.map(|_| ())),
        }
    }
}

#[test]
fn walls_and_pmls_in_diagonal_exteriors() {
    let (k0, beta) = (2.1, 0.8);
    let layer = || isotropic(dispersive(c(2.4, 0.05), c(1.0, 0.0)));
    let exterior = constant(c(1.7, 0.02), c(1.2, 0.01));
    let terminations: [crate::Termination<f64>; 3] = [
        Wall::new(WallKind::PerfectElectric, Length::centimetres(0.45)).into(),
        Wall::new(WallKind::PerfectMagnetic, Length::centimetres(0.45)).into(),
        Pml::new(
            Length::centimetres(0.4),
            Length::centimetres(0.3),
            StretchProfile::polynomial(2.0, c(1.5, 2.0)),
            WallKind::PerfectElectric,
        )
        .into(),
    ];
    for termination in terminations {
        // An isotropic tensor exterior against the scalar one.
        let build = |right: Region| {
            Stack::builder(isotropic(constant(c(1.0, 0.0), c(1.0, 0.0))), right)
                .layer(layer(), Length::centimetres(0.25))
                .finalise()
                .with_right_termination(termination)
        };
        let (scalar, tensor) = (build(isotropic(exterior)), build(diagonal([exterior; 3])));
        for polarisation in POLARISATIONS {
            for_each_backend!(evaluator, {
                let reflection = |stack: &AnalyticalTensorStack<C>| {
                    evaluator
                        .evaluate(scalar_real_input(k0, beta), stack, polarisation)
                        .unwrap()
                        .amplitudes(IncidentSide::Left)
                        .unwrap()
                        .reflection()[()]
                };
                assert_close(
                    reflection(&tensor),
                    reflection(&scalar),
                    1e-12,
                    "terminated reflection",
                );
            });
            let coordinates = CanonicalCoordinates::new(
                HoloJ0::constant(arr0(c(k0, -0.05))),
                HoloJ0::constant(arr0(c(beta, 0.1))),
            );
            let determinant = |stack: &AnalyticalTensorStack<C>| {
                let evaluator =
                    ComplexPlaneEvaluator::<HoloJ0, _, _>::compile(stack, Transfer2::new())
                        .unwrap();
                let exterior = ExteriorWavevectors::closed_right(
                    principal_exterior_wavevectors(evaluator.stack(), &coordinates)
                        .left()
                        .unwrap()
                        .scalar()
                        .unwrap()
                        .clone(),
                );
                evaluator
                    .checked_determinant(&coordinates, &exterior, polarisation)
                    .unwrap()
                    .value()[()]
            };
            assert_close(
                determinant(&tensor),
                determinant(&scalar),
                1e-12,
                "terminated characteristic",
            );
        }
    }

    // A wall behind a gap of uniaxial exterior is the wall behind an explicit
    // layer of that material: the gap takes the exterior's per-polarisation
    // κ and admittance.
    for axis in 0..3 {
        let uniaxial = || {
            uniaxial(
                constant(c(2.25, 0.01), c(1.1, 0.0)),
                constant(c(4.0, 0.02), c(0.9, 0.0)),
                axis,
            )
        };
        let wall = Wall::new(WallKind::PerfectElectric, Length::centimetres(0.4));
        let gap = Stack::builder(isotropic(constant(c(1.0, 0.0), c(1.0, 0.0))), uniaxial())
            .layer(layer(), Length::centimetres(0.25))
            .finalise()
            .with_right_termination(wall);
        let explicit = Stack::builder(
            isotropic(constant(c(1.0, 0.0), c(1.0, 0.0))),
            isotropic(constant(c(1.0, 0.0), c(1.0, 0.0))),
        )
        .layer(layer(), Length::centimetres(0.25))
        .layer(uniaxial(), Length::centimetres(0.15))
        .finalise()
        .with_right_termination(wall);
        for polarisation in POLARISATIONS {
            for_each_backend!(evaluator, {
                let reflection = |stack: &AnalyticalTensorStack<C>| {
                    evaluator
                        .evaluate(scalar_real_input(k0, beta), stack, polarisation)
                        .unwrap()
                        .amplitudes(IncidentSide::Left)
                        .unwrap()
                        .reflection()[()]
                };
                assert_close(
                    reflection(&gap),
                    reflection(&explicit),
                    1e-12,
                    "gap against explicit layer",
                );
            });
        }
    }
}

#[test]
fn transfer_and_scatter_agree_on_diagonal_stacks() {
    let uniaxial_stack = Stack::builder(
        isotropic(constant(c(1.0, 0.0), c(1.0, 0.0))),
        uniaxial(
            constant(c(1.6, 0.0), c(1.0, 0.0)),
            constant(c(2.4, 0.0), c(1.0, 0.0)),
            2,
        ),
    )
    .layer(
        uniaxial(
            dispersive(c(2.25, 0.0), c(1.0, 0.0)),
            dispersive(c(4.0, 0.0), c(1.0, 0.0)),
            1,
        ),
        Length::centimetres(0.3),
    )
    .layer(
        isotropic(dispersive(c(1.9, 0.0), c(1.0, 0.0))),
        Length::centimetres(0.2),
    )
    .finalise();
    let (value, first, second) = (1e-11, 1e-10, 1e-9);

    for stack in [uniaxial_stack, biaxial_stack()] {
        for polarisation in POLARISATIONS {
            for axes in [
                (
                    Parameter::Spectral,
                    Parameter::LayerThickness(FiniteLayerIndex::new(0)),
                ),
                (Parameter::InPlane, Parameter::Spectral),
            ] {
                macro_rules! evaluate {
                    ($evaluator:expr) => {
                        $evaluator
                            .evaluate_bivariate_second(
                                scalar_real_input(2.2, 1.1),
                                &stack,
                                polarisation,
                                axes.0,
                                axes.1,
                            )
                            .unwrap()
                    };
                }
                let (t, s) = (
                    evaluate!(RealAxisEvaluator::new(Transfer2::new())),
                    evaluate!(RealAxisEvaluator::new(Scatter2::new())),
                );
                for side in [IncidentSide::Left, IncidentSide::Right] {
                    let (ta, sa) = (t.amplitudes(side).unwrap(), s.amplitudes(side).unwrap());
                    let (tp, sp) = (t.power(side).unwrap(), s.power(side).unwrap());
                    for (tolerance, a, b, p, q) in [
                        (value, ta.value(), sa.value(), tp.value(), sp.value()),
                        (
                            first,
                            ta.gradient().axis0(),
                            sa.gradient().axis0(),
                            tp.gradient().axis0(),
                            sp.gradient().axis0(),
                        ),
                        (
                            first,
                            ta.gradient().axis1(),
                            sa.gradient().axis1(),
                            tp.gradient().axis1(),
                            sp.gradient().axis1(),
                        ),
                        (
                            second,
                            ta.hessian().axis0_axis0(),
                            sa.hessian().axis0_axis0(),
                            tp.hessian().axis0_axis0(),
                            sp.hessian().axis0_axis0(),
                        ),
                        (
                            second,
                            ta.hessian().axis0_axis1(),
                            sa.hessian().axis0_axis1(),
                            tp.hessian().axis0_axis1(),
                            sp.hessian().axis0_axis1(),
                        ),
                        (
                            second,
                            ta.hessian().axis1_axis1(),
                            sa.hessian().axis1_axis1(),
                            tp.hessian().axis1_axis1(),
                            sp.hessian().axis1_axis1(),
                        ),
                    ] {
                        assert_close(
                            a.reflection()[()],
                            b.reflection()[()],
                            tolerance,
                            "reflection",
                        );
                        assert_close(
                            a.transmission()[()],
                            b.transmission()[()],
                            tolerance,
                            "transmission",
                        );
                        for (x, y) in [
                            (p.reflectance()[()], q.reflectance()[()]),
                            (p.transmittance()[()], q.transmittance()[()]),
                            (p.absorptance()[()], q.absorptance()[()]),
                        ] {
                            assert_close(c(x, 0.0), c(y, 0.0), tolerance, "power");
                        }
                    }
                }
            }

            // Retained boundary waves.
            macro_rules! waves {
                ($evaluator:expr, $side:expr) => {
                    $evaluator
                        .retain(scalar_real_input(2.2, 1.1), &stack, polarisation)
                        .unwrap()
                        .project_point(&())
                        .unwrap()
                        .excitation($side)
                        .unwrap()
                        .boundary_waves()
                        .unwrap()
                };
            }
            for side in [IncidentSide::Left, IncidentSide::Right] {
                let (t, s) = (
                    waves!(RealAxisEvaluator::new(Transfer2::new()), side),
                    waves!(RealAxisEvaluator::new(Scatter2::new()), side),
                );
                for (t, s) in t.iter().zip(s.iter()) {
                    for (a, b) in [(t.left(), s.left()), (t.right(), s.right())] {
                        assert_close(a.forward()[()], b.forward()[()], value, "forward wave");
                        assert_close(a.backward()[()], b.backward()[()], value, "backward wave");
                    }
                }
            }

            // The characteristic and its bivariate derivatives in (k0, k∥).
            let coordinates = || {
                CanonicalCoordinates::new(
                    <HoloJB2 as SeedJet>::variable(arr0(c(2.2, -0.1)), 0).unwrap(),
                    <HoloJB2 as SeedJet>::variable(arr0(c(1.4, 0.15)), 1).unwrap(),
                )
            };
            macro_rules! determinant {
                ($backend:expr) => {{
                    let evaluator =
                        ComplexPlaneEvaluator::<HoloJB2, _, _>::compile(&stack, $backend).unwrap();
                    let coordinates = coordinates();
                    let exterior = evaluate_stack_exterior_wavevectors::<ComplexPlane, _, _>(
                        &coordinates,
                        evaluator.stack(),
                        polarisation.into(),
                    )
                    .unwrap();
                    evaluator
                        .checked_determinant(&coordinates, &exterior, polarisation)
                        .unwrap()
                        .into_inner()
                }};
            }
            let (t, s) = (
                determinant!(Transfer2::new()),
                determinant!(Scatter2::new()),
            );
            for (tolerance, a, b) in [
                (value, t.value(), s.value()),
                (first, t.axis0(), s.axis0()),
                (first, t.axis1(), s.axis1()),
                (second, t.axis0_axis0(), s.axis0_axis0()),
                (second, t.axis0_axis1(), s.axis0_axis1()),
                (second, t.axis1_axis1(), s.axis1_axis1()),
            ] {
                assert_close(a[()], b[()], tolerance, "characteristic");
            }
        }
    }
}

#[test]
fn scattering_backend_resolves_thick_evanescent_diagonal_layers() {
    // Total internal reflection at a uniaxial layer far thicker than its
    // decay length: the reflection is that of a uniaxial half-space.
    let (k0, beta) = (2.0, 2.0 * 5.0_f64.sqrt());
    let dense = || isotropic(constant(c(6.0, 0.0), c(1.0, 0.0)));
    for axis in 0..3 {
        let layer = || {
            uniaxial(
                constant(c(2.25, 0.0), c(1.0, 0.0)),
                constant(c(4.0, 0.0), c(1.0, 0.0)),
                axis,
            )
        };
        let thick = Stack::builder(dense(), dense())
            .layer(layer(), Length::centimetres(400.0))
            .finalise();
        let half_space = Stack::builder(dense(), layer()).finalise();
        for polarisation in POLARISATIONS {
            let evaluator = RealAxisEvaluator::new(Scatter2::new());
            let thick = evaluator
                .evaluate(scalar_real_input(k0, beta), &thick, polarisation)
                .unwrap();
            let half_space = evaluator
                .evaluate(scalar_real_input(k0, beta), &half_space, polarisation)
                .unwrap();
            let reflection = half_space
                .amplitudes(IncidentSide::Left)
                .unwrap()
                .reflection()[()];
            assert!((reflection.norm() - 1.0).abs() < 1e-14);
            let amplitudes = thick.amplitudes(IncidentSide::Left).unwrap();
            assert_close(
                amplitudes.reflection()[()],
                reflection,
                1e-13,
                "thick reflection",
            );
            assert!(amplitudes.transmission()[()].norm() < 1e-100);
        }
    }
}

/// diag(o, o, e) with the extraordinary ε and μ linearly interpolated between
/// tabulated samples, and not finite outside them: real-axis values only.
struct TabulatedExtraordinary {
    ordinary: ComplexQuadratic,
    wavenumbers: Vec<f64>,
    extraordinary: Vec<[C; 2]>,
}

impl TabulatedExtraordinary {
    /// Tabulate `extraordinary` at `wavenumbers`.
    fn new(
        ordinary: ComplexQuadratic,
        extraordinary: ComplexQuadratic,
        wavenumbers: &[f64],
    ) -> Self {
        use crate::Material;
        let extraordinary = wavenumbers
            .iter()
            .map(|&k| {
                [
                    extraordinary.relative_permittivity(crate::Scalar(k)),
                    extraordinary.relative_permeability(crate::Scalar(k)),
                ]
            })
            .collect();
        Self {
            ordinary,
            wavenumbers: wavenumbers.to_vec(),
            extraordinary,
        }
    }

    fn tensor<X>(&self, k: f64, quantity: usize) -> crate::Tensor3<X>
    where
        X: crate::ComplexScalar<RealField = f64>,
    {
        use crate::Material;
        let ordinary: X = [
            self.ordinary.relative_permittivity(crate::Scalar(k)),
            self.ordinary.relative_permeability(crate::Scalar(k)),
        ][quantity];
        let extraordinary = self
            .wavenumbers
            .windows(2)
            .zip(self.extraordinary.windows(2))
            .find(|(nodes, _)| nodes[0] <= k && k <= nodes[1])
            .map_or(C::new(f64::NAN, f64::NAN), |(nodes, values)| {
                let t = (k - nodes[0]) / (nodes[1] - nodes[0]);
                values[0][quantity] * (1.0 - t) + values[1][quantity] * t
            });
        let zero = X::zero();
        [
            [ordinary, zero, zero],
            [zero, ordinary, zero],
            [
                zero,
                zero,
                X::from_parts(extraordinary.re, extraordinary.im),
            ],
        ]
    }
}

impl crate::TensorMaterial for TabulatedExtraordinary {
    type Real = f64;
    const PRINCIPAL_AXES: bool = true;

    fn relative_permittivity_tensor<I, X>(&self, k: I) -> I::Mapped<crate::Tensor3<X>>
    where
        X: crate::ComplexScalar<RealField = f64>,
        I: crate::Sampled<Elem = f64>,
    {
        k.map(|k| self.tensor(k, 0))
    }

    fn relative_permeability_tensor<I, X>(&self, k: I) -> I::Mapped<crate::Tensor3<X>>
    where
        X: crate::ComplexScalar<RealField = f64>,
        I: crate::Sampled<Elem = f64>,
    {
        k.map(|k| self.tensor(k, 1))
    }
}

#[test]
fn tabulated_tensor_regions_give_the_analytical_real_axis_response_at_the_samples() {
    use crate::{MaterialRegion, TensorMaterialHandle, TensorStack, material::MaterialHandle};

    let wavenumbers = [1.5, 1.8, 2.1, 2.4, 2.7];
    let (ordinary, extraordinary) = (
        dispersive(c(2.25, 0.05), c(1.05, 0.01)),
        dispersive(c(4.0, 0.2), c(0.95, 0.02)),
    );
    let (vacuum, cladding) = (
        constant(c(1.0, 0.0), c(1.0, 0.0)),
        dispersive(c(1.9, 0.08), c(1.05, 0.02)),
    );
    // Optic axis along stack x and along stack z: both diagonal classes.
    let rotations = [
        Rotation::from_optic_axis(FRAC_PI_2, 0.0),
        Rotation::identity(),
    ];

    let tabulated = |rotation| {
        MaterialRegion::oriented(
            TensorMaterialHandle::new(TabulatedExtraordinary::new(
                ordinary,
                extraordinary,
                &wavenumbers,
            )),
            rotation,
        )
    };
    let scalar = |material| MaterialRegion::from(MaterialHandle::new(material));
    let real_axis: TensorStack<C> = Stack::builder(scalar(vacuum), tabulated(rotations[1]))
        .layer(tabulated(rotations[0]), Length::centimetres(0.21))
        .layer(scalar(cladding), Length::centimetres(0.17))
        .finalise();
    let analytical_region =
        |rotation| principal_region([ordinary, ordinary, extraordinary], rotation);
    let analytical: AnalyticalTensorStack<C> =
        Stack::builder(isotropic(vacuum), analytical_region(rotations[1]))
            .layer(analytical_region(rotations[0]), Length::centimetres(0.21))
            .layer(isotropic(cladding), Length::centimetres(0.17))
            .finalise();
    assert_eq!(
        real_axis.layers_left_to_right()[0].material().class(),
        ConstitutiveClass::Diagonal
    );
    assert_eq!(
        real_axis.right_exterior().class(),
        ConstitutiveClass::Diagonal
    );

    let input = || sampled_real_input(&wavenumbers, &[0.9; 5]);
    for polarisation in POLARISATIONS {
        for_each_backend!(evaluator, {
            let tabulated = evaluator
                .evaluate(input(), &real_axis, polarisation)
                .unwrap();
            let expected = evaluator
                .evaluate(input(), &analytical, polarisation)
                .unwrap();
            for side in [IncidentSide::Left, IncidentSide::Right] {
                let (a, b) = (
                    tabulated.amplitudes(side).unwrap(),
                    expected.amplitudes(side).unwrap(),
                );
                assert_arrays_close(a.reflection(), b.reflection(), 1e-14, "reflection");
                assert_arrays_close(a.transmission(), b.transmission(), 1e-14, "transmission");
                let (p, q) = (
                    tabulated.power(side).unwrap(),
                    expected.power(side).unwrap(),
                );
                for (x, y) in [
                    (p.reflectance(), q.reflectance()),
                    (p.transmittance(), q.transmittance()),
                    (p.absorptance(), q.absorptance()),
                ] {
                    assert_arrays_close(&x.mapv(C::from), &y.mapv(C::from), 1e-14, "power");
                }
            }
        });
    }

    // Outside the table the tabulated entries are not finite: a typed refusal.
    let non_finite = ConstitutiveError::NonFinite {
        region: StackRegion::RightExterior,
    };
    for polarisation in POLARISATIONS {
        for_each_backend!(evaluator, {
            match evaluator.evaluate(scalar_real_input(3.0, 0.9), &real_axis, polarisation) {
                Err(RealAxisEvaluationError::Constitutive(error)) => assert_eq!(error, non_finite),
                other => panic!("unexpected {:?}", other.map(|_| ())),
            }
        });
    }
}
