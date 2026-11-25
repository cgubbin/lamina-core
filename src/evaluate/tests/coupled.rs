//! The 4×4 backends against each other, against the 2×2 backends, and
//! against Maxwell's equations.

use lamina_units::Length;
use ndarray::arr0;

use crate::{
    AnalyticalMaterialHandle, AnalyticalMaterialRegion, AnalyticalTensorStack, Backend,
    CanonicalCoordinates, Channel, ChannelError, ComplexPlane, ComplexPlaneEvaluator,
    ExteriorWavevectors, FiniteLayerIndex, IncidentSide, JonesAmplitudes, PlaneWaveAmplitudes, Pml,
    Polarisation, ProjectJonesAmplitudes, RealAxisEvaluationError, RealAxisEvaluator, Rotation,
    Scatter2, Scatter4, Scatter4Error, SeedJet, Stack, StackRegion, StretchProfile, Transfer2,
    Transfer4, Transfer4Error, Wall, WallKind,
    algebra::ScalarAlgebra,
    backend::{
        CoupledEntries, CoupledWorkspace, PlaneWaveSolution, evaluate_stack_exterior_wavevectors,
    },
    material::TensorConstitutiveLift,
    test_support::{
        jet::{HoloJ0, RealJ0, RealJ2},
        materials::{ComplexQuadratic, principal_region},
        planar::scalar_real_input,
    },
};

pub(super) type C = num_complex::Complex64;
pub(super) type Region = AnalyticalMaterialRegion<C>;

const POLARISATIONS: [Polarisation; 2] = [
    Polarisation::TransverseElectric,
    Polarisation::TransverseMagnetic,
];
pub(super) const SIDES: [IncidentSide; 2] = [IncidentSide::Left, IncidentSide::Right];

pub(super) fn c(re: f64, im: f64) -> C {
    C::new(re, im)
}

pub(super) fn assert_close(actual: C, expected: C, tolerance: f64, label: &str) {
    let error = (actual - expected).norm();
    assert!(
        error <= tolerance * expected.norm().max(1.0),
        "{label}: expected {expected}, got {actual}, error {error:e}"
    );
}

pub(super) fn material(epsilon: C, mu: C) -> ComplexQuadratic {
    ComplexQuadratic {
        epsilon: [epsilon, c(0.03, 0.01), c(-0.002, 0.001)],
        mu: [mu, c(-0.01, 0.004), C::default()],
    }
}

pub(super) fn lossless(epsilon: f64) -> ComplexQuadratic {
    ComplexQuadratic {
        epsilon: [c(epsilon, 0.0), C::default(), C::default()],
        mu: [c(1.0, 0.0), C::default(), C::default()],
    }
}

pub(super) fn isotropic(material: ComplexQuadratic) -> Region {
    AnalyticalMaterialHandle::new(material).into()
}

pub(super) fn oriented(axes: [ComplexQuadratic; 3], rotation: Rotation<f64>) -> Region {
    principal_region(axes, rotation)
}

/// A lossy, dispersive, magnetic stack of a biaxial, an isotropic and a
/// uniaxial layer (optic axis along z) between isotropic exteriors.
pub(super) fn diagonal_stack() -> AnalyticalTensorStack<C> {
    Stack::builder(
        isotropic(material(c(1.2, 0.0), c(1.0, 0.0))),
        isotropic(material(c(2.1, 0.02), c(1.1, 0.0))),
    )
    .layer(
        oriented(
            [
                material(c(2.2, 0.15), c(1.1, 0.03)),
                material(c(3.1, 0.05), c(0.9, 0.02)),
                material(c(4.0, 0.3), c(1.3, 0.01)),
            ],
            Rotation::identity(),
        ),
        Length::centimetres(0.21),
    )
    .layer(
        isotropic(material(c(1.7, 0.1), c(1.0, 0.0))),
        Length::centimetres(0.1),
    )
    .layer(
        oriented(
            [
                material(c(2.5, 0.1), c(1.0, 0.0)),
                material(c(2.5, 0.1), c(1.0, 0.0)),
                material(c(5.0, 0.2), c(1.0, 0.0)),
            ],
            Rotation::identity(),
        ),
        Length::centimetres(0.13),
    )
    .finalise()
}

/// A lossy stack of general orientations: no channel is invariant.
fn full_stack() -> AnalyticalTensorStack<C> {
    full_stack_with(0.13)
}

/// [`full_stack`] with a second layer of thickness `thickness` cm.
pub(super) fn full_stack_with(thickness: f64) -> AnalyticalTensorStack<C> {
    Stack::builder(
        isotropic(material(c(1.0, 0.0), c(1.0, 0.0))),
        oriented(
            [lossless(2.1), lossless(2.1), lossless(2.9)],
            Rotation::from_euler(0.4, 0.9, 0.2),
        ),
    )
    .layer(
        oriented(
            [
                material(c(2.2, 0.15), c(1.1, 0.03)),
                material(c(3.1, 0.05), c(0.9, 0.02)),
                material(c(4.0, 0.3), c(1.3, 0.01)),
            ],
            Rotation::from_euler(0.3, 0.7, -0.5),
        ),
        Length::centimetres(0.21),
    )
    .layer(
        oriented(
            [
                material(c(2.5, 0.1), c(1.0, 0.0)),
                material(c(2.5, 0.1), c(1.0, 0.0)),
                material(c(5.0, 0.2), c(1.0, 0.0)),
            ],
            Rotation::from_optic_axis(0.8, 0.6),
        ),
        Length::centimetres(thickness),
    )
    .finalise()
}

/// The real-axis solution of `backend` at one point.
type Solved<B> = PlaneWaveSolution<<B as Backend<RealJ0, crate::RealAxis>>::Entries>;

fn solve<B>(
    backend: B,
    stack: &AnalyticalTensorStack<C>,
    k0: f64,
    beta: f64,
    channel: Channel,
) -> Solved<B>
where
    B: Backend<RealJ0, crate::RealAxis>,
{
    RealAxisEvaluator::new(backend)
        .evaluate(scalar_real_input(k0, beta), stack, channel)
        .unwrap()
        .into_parts()
        .0
}

fn jones(
    solution: &PlaneWaveSolution<CoupledEntries<RealJ0>>,
    side: IncidentSide,
) -> JonesAmplitudes<RealJ0> {
    solution
        .entries()
        .project_jones_amplitudes(solution.context(), side)
}

fn entry(value: Option<&RealJ0>) -> C {
    value.map_or(c(f64::NAN, f64::NAN), |value| value.value()[()])
}

fn amplitudes(
    solution: &PlaneWaveSolution<crate::backend::scatter2::Scatter2ProjectiveEntries<RealJ0>>,
    side: IncidentSide,
) -> (C, C) {
    let amplitudes: PlaneWaveAmplitudes<RealJ0> = solution.amplitudes(side);
    let (r, t) = amplitudes.into_parts();
    (r.value()[()], t.value()[()])
}

type Coupled<J> = PlaneWaveSolution<CoupledEntries<J>>;

/// The Transfer4 and Scatter4 solutions at one point.
fn coupled(
    stack: &AnalyticalTensorStack<C>,
    k0: f64,
    beta: f64,
    channel: Channel,
) -> [Coupled<RealJ0>; 2] {
    [
        solve(Transfer4::new(), stack, k0, beta, channel),
        solve(Scatter4::new(), stack, k0, beta, channel),
    ]
}

/// The Transfer4 and Scatter4 retained workspaces of the coupled channel.
fn retained(stack: &AnalyticalTensorStack<C>, k0: f64, beta: f64) -> [CoupledWorkspace<RealJ0>; 2] {
    let input = || scalar_real_input(k0, beta);
    [
        RealAxisEvaluator::new(Transfer4::new())
            .retain(input(), stack, Channel::Coupled)
            .unwrap()
            .workspace()
            .clone(),
        RealAxisEvaluator::new(Scatter4::new())
            .retain(input(), stack, Channel::Coupled)
            .unwrap()
            .workspace()
            .clone(),
    ]
}

/// Every entry of two Jones matrices agrees.
fn assert_jones_close(
    actual: &JonesAmplitudes<RealJ0>,
    expected: &JonesAmplitudes<RealJ0>,
    tolerance: f64,
    label: &str,
) {
    for row in 0..expected.size() {
        for column in 0..expected.size() {
            let label = format!("{label} ({row}, {column})");
            assert_close(
                entry(actual.reflection(row, column)),
                entry(expected.reflection(row, column)),
                tolerance,
                &label,
            );
            assert_close(
                entry(actual.transmission(row, column)),
                entry(expected.transmission(row, column)),
                tolerance,
                &label,
            );
        }
    }
}

const POINTS: [(f64, f64); 3] = [(2.0, 0.3), (2.6, 1.9), (3.1, 4.2)];

/// Jones matrices that are diagonal with entries `(r, t)` of `expected`
/// within `tolerance`, off-diagonal ones below `1e-12`.
fn assert_diagonal(
    jones: &JonesAmplitudes<RealJ0>,
    expected: &[(C, C)],
    tolerance: f64,
    label: &str,
) {
    assert_eq!(jones.size(), expected.len(), "{label}");
    for row in 0..expected.len() {
        for column in 0..expected.len() {
            let (r, t, tolerance) = if row == column {
                (expected[row].0, expected[row].1, tolerance)
            } else {
                (C::default(), C::default(), 1e-12)
            };
            let label = format!("{label} ({row}, {column})");
            assert_close(entry(jones.reflection(row, column)), r, tolerance, &label);
            assert_close(entry(jones.transmission(row, column)), t, tolerance, &label);
        }
    }
}

#[test]
fn scalar_and_coupled_channels_match_the_2x2_backends_on_separable_stacks() {
    // A uniaxial crystal with its optic axis on y, turned about y: a full
    // (separable) region whose tensor is the diagonal one to round-off.
    let axes = [
        material(c(2.2, 0.1), c(1.1, 0.0)),
        material(c(3.4, 0.2), c(0.9, 0.01)),
        material(c(2.2, 0.1), c(1.1, 0.0)),
    ];
    let uniaxial = |rotation| {
        Stack::builder(isotropic(lossless(1.0)), isotropic(lossless(1.8)))
            .layer(oriented(axes, rotation), Length::centimetres(0.25))
            .finalise()
    };
    let full = uniaxial(Rotation::from_euler(0.0, 0.7, 0.0));
    assert_eq!(
        full.layers_left_to_right()[0].material().class(),
        crate::ConstitutiveClass::Full
    );

    // (4×4 stack, 2×2 stack, tolerance)
    for (four, two, tolerance) in [
        (diagonal_stack(), diagonal_stack(), 1e-10),
        (full, uniaxial(Rotation::identity()), 1e-12),
    ] {
        for ((k0, beta), side) in POINTS
            .into_iter()
            .flat_map(|point| SIDES.map(|side| (point, side)))
        {
            let label = format!("{side:?} k0 {k0} β {beta}");
            let expected = POLARISATIONS.map(|polarisation| {
                amplitudes(
                    &solve(Scatter2::new(), &two, k0, beta, polarisation.into()),
                    side,
                )
            });
            for (polarisation, expected) in POLARISATIONS.iter().zip(expected) {
                for solution in coupled(&four, k0, beta, (*polarisation).into()) {
                    assert_diagonal(&jones(&solution, side), &[expected], tolerance, &label);
                }
            }
            // The coupled channel: TE then TM on the diagonal.
            for solution in coupled(&four, k0, beta, Channel::Coupled) {
                assert_diagonal(&jones(&solution, side), &expected, tolerance, &label);
            }
        }
    }
}

#[test]
fn power_matches_the_2x2_backends_and_conserves_flux_without_loss() {
    let stack = diagonal_stack();
    for (k0, beta) in POINTS {
        for polarisation in POLARISATIONS {
            let reference = solve(Transfer2::new(), &stack, k0, beta, polarisation.into());
            for (side, solution) in SIDES.into_iter().flat_map(|side| {
                coupled(&stack, k0, beta, polarisation.into()).map(|solution| (side, solution))
            }) {
                let (expected, actual) = (reference.power(side), solution.power(side));
                for (actual, expected) in [
                    (&actual.reflectance()[0], expected.reflectance()),
                    (&actual.transmittance()[0], expected.transmittance()),
                ] {
                    assert!(
                        (actual.value()[()] - expected.value()[()]).abs() < 1e-10,
                        "{side:?}"
                    );
                }
            }
        }
    }

    // A lossless stack of Hermitian (real symmetric) tensors: R + T = 1 for
    // every incident eigen-wave, propagating on both sides.
    let stack = Stack::builder(
        isotropic(lossless(1.0)),
        oriented(
            [lossless(2.1), lossless(2.1), lossless(2.9)],
            Rotation::from_euler(0.4, 0.9, 0.2),
        ),
    )
    .layer(
        oriented(
            [lossless(2.2), lossless(3.1), lossless(4.0)],
            Rotation::from_euler(0.3, 0.7, -0.5),
        ),
        Length::centimetres(0.21),
    )
    .finalise();
    for solution in coupled(&stack, 3.0, 1.1, Channel::Coupled) {
        for side in SIDES {
            let power = solution.power(side);
            for wave in 0..2 {
                assert!(
                    power.absorptance()[wave].value()[()].abs() < 1e-12,
                    "{side:?} {wave}"
                );
                assert!(power.reflectance()[wave].value()[()] > 0.0);
            }
        }
    }
}

#[test]
fn transfer4_and_scatter4_agree_on_full_tensor_stacks() {
    let stack = full_stack();
    for (k0, beta) in POINTS {
        let [transfer, scatter] = coupled(&stack, k0, beta, Channel::Coupled);
        for side in [crate::StackSide::Left, crate::StackSide::Right] {
            let conditioning = transfer.context().basis_conditioning(side).unwrap();
            assert!(conditioning > 0.1 && conditioning <= 1.0, "{conditioning}");
        }
        for side in SIDES {
            let label = format!("{side:?} k0 {k0} β {beta}");
            assert_jones_close(
                &jones(&scatter, side),
                &jones(&transfer, side),
                1e-10,
                &label,
            );
            // Off-axis tensors couple TE and TM.
            assert!(entry(jones(&transfer, side).reflection(0, 1)).norm() > 1e-4);
            let (tp, sp) = (transfer.power(side), scatter.power(side));
            for wave in 0..2 {
                let difference =
                    tp.absorptance()[wave].value()[()] - sp.absorptance()[wave].value()[()];
                assert!(difference.abs() < 1e-10, "{label}");
            }
        }
    }

    // Retained fields in every region.
    let [transfer, scatter] = retained(&stack, 2.6, 1.9);
    let depth = RealJ0::lift_constant(arr0(c(0.05, 0.0)));
    let regions = [
        StackRegion::LeftExterior,
        StackRegion::Layer(FiniteLayerIndex::new(0)),
        StackRegion::Layer(FiniteLayerIndex::new(1)),
        StackRegion::RightExterior,
    ];
    for (side, wave, region) in SIDES
        .into_iter()
        .flat_map(|side| (0..2).flat_map(move |wave| regions.map(|region| (side, wave, region))))
    {
        let t = transfer.fields(side, wave, region, &depth).unwrap();
        let s = scatter.fields(side, wave, region, &depth).unwrap();
        for (a, b) in t.iter().flatten().zip(s.iter().flatten()) {
            assert_close(
                b.value()[()],
                a.value()[()],
                1e-10,
                &format!("{side:?} {wave} {region:?}"),
            );
        }
    }
}

/// The winding number of `f` around the rectangle `[lo, hi]`.
fn winding(f: impl Fn(C) -> C, lo: C, hi: C) -> f64 {
    let corners = [lo, c(hi.re, lo.im), hi, c(lo.re, hi.im), lo];
    let steps = 4000;
    let mut phase = 0.0;
    let mut previous = f(lo);
    for pair in corners.windows(2) {
        for step in 1..=steps {
            let z = pair[0] + (pair[1] - pair[0]) * (step as f64 / steps as f64);
            let value = f(z);
            phase += (value / previous).arg();
            previous = value;
        }
    }
    phase / std::f64::consts::TAU
}

/// The decaying normal wavenumber `i√(β² − εμ k₀²)` of an isotropic
/// exterior, continuous above its light line.
fn decaying(epsilon_mu: f64, k0: f64, beta: C) -> HoloJ0 {
    <HoloJ0 as SeedJet>::constant(arr0(
        c(0.0, 1.0) * (beta * beta - epsilon_mu * k0 * k0).sqrt(),
    ))
}

/// The exterior branches at complex `β` of decaying isotropic exteriors
/// `(ε_L, ε_R)` on the open sides.
fn exterior(
    k0: f64,
    beta: C,
    exteriors: [Option<f64>; 2],
    channel: Channel,
) -> ExteriorWavevectors<HoloJ0> {
    let branch = |epsilon: Option<f64>| {
        epsilon.map(|epsilon| {
            let kappa = decaying(epsilon, k0, beta);
            match channel {
                Channel::Coupled => crate::ExteriorBranch::Pair([kappa.clone(), kappa]),
                _ => crate::ExteriorBranch::Scalar(kappa),
            }
        })
    };
    ExteriorWavevectors::from_branches(branch(exteriors[0]), branch(exteriors[1]))
}

/// The characteristic of `backend` at complex `β`.
macro_rules! characteristic {
    ($backend:expr, $stack:expr, $k0:expr, $exteriors:expr, $channel:expr) => {{
        let evaluator = ComplexPlaneEvaluator::<HoloJ0, _, _>::compile($stack, $backend).unwrap();
        let (k0, exteriors, channel) = ($k0, $exteriors, $channel);
        move |beta: C| {
            let coordinates = CanonicalCoordinates::new(
                <HoloJ0 as SeedJet>::constant(arr0(c(k0, 0.0))),
                <HoloJ0 as SeedJet>::constant(arr0(beta)),
            );
            evaluator
                .checked_determinant(
                    &coordinates,
                    &exterior(k0, beta, exteriors, channel),
                    channel,
                )
                .unwrap()
                .value()[()]
        }
    }};
}

/// The decaying exteriors of [`waveguide`] on its open sides.
fn waveguide_exteriors(closed: bool) -> [Option<f64>; 2] {
    if closed {
        [None, None]
    } else {
        [Some(1.0), Some(2.0736)]
    }
}

/// A lossless slab waveguide with a biaxial core between isotropic
/// claddings, optionally closed by a PML on the left and a wall on the right.
fn waveguide(core: Region, closed: bool) -> AnalyticalTensorStack<C> {
    let stack = Stack::builder(isotropic(lossless(1.0)), isotropic(lossless(2.0736)))
        .layer(core, Length::centimetres(1.0))
        .finalise();
    if !closed {
        return stack;
    }
    stack
        .with_left_termination(Pml::new(
            Length::centimetres(-0.9),
            Length::centimetres(0.6),
            StretchProfile::polynomial(2.0, c(2.0, 3.0)),
            WallKind::PerfectElectric,
        ))
        .with_right_termination(Wall::new(
            WallKind::PerfectMagnetic,
            Length::centimetres(1.7),
        ))
}

#[test]
fn coupled_characteristic_factorises_and_counts_te_plus_tm_modes() {
    let k0 = 4.0;
    let (lo, hi) = (c(5.85, -0.05), c(8.6, 0.05));
    let core = || {
        oriented(
            [lossless(4.0), lossless(4.4), lossless(3.6)],
            Rotation::identity(),
        )
    };

    for closed in [false, true] {
        let stack = waveguide(core(), closed);
        let exteriors = waveguide_exteriors(closed);
        let te = characteristic!(
            Transfer2::new(),
            &stack,
            k0,
            exteriors,
            Channel::TransverseElectric
        );
        let tm = characteristic!(
            Transfer2::new(),
            &stack,
            k0,
            exteriors,
            Channel::TransverseMagnetic
        );
        let coupled = characteristic!(Transfer4::new(), &stack, k0, exteriors, Channel::Coupled);
        let scattered = characteristic!(Scatter4::new(), &stack, k0, exteriors, Channel::Coupled);

        let (n_te, n_tm) = (winding(&te, lo, hi), winding(&tm, lo, hi));
        let n_coupled = winding(&coupled, lo, hi);
        let label = format!("closed {closed}: TE {n_te}, TM {n_tm}, coupled {n_coupled}");
        assert!(n_te.round() >= 2.0 && n_tm.round() >= 2.0, "{label}");
        assert!((n_coupled - n_te - n_tm).abs() < 1e-6, "{label}");
        assert!(
            (winding(&scattered, lo, hi) - n_coupled).abs() < 1e-6,
            "{label}"
        );
        // The gauge g = F / (F_TE F_TM) and the Scatter4 / Transfer4 ratio are
        // zero-free on the region.
        assert!(
            winding(|beta| coupled(beta) / (te(beta) * tm(beta)), lo, hi).abs() < 1e-6,
            "{label}"
        );
        assert!(
            winding(|beta| scattered(beta) / coupled(beta), lo, hi).abs() < 1e-6,
            "{label}"
        );

        // The 4×4 scalar channels have the 2×2 zeros.
        let scalar: [(Channel, &dyn Fn(C) -> C); 2] = [
            (Channel::TransverseElectric, &te),
            (Channel::TransverseMagnetic, &tm),
        ];
        for (channel, two) in scalar {
            let four = characteristic!(Transfer4::new(), &stack, k0, exteriors, channel);
            assert!(
                winding(|beta| four(beta) / two(beta), lo, hi).abs() < 1e-6,
                "{label} {channel:?}"
            );
        }
    }

    // Rotating the core's optic axis off the principal axes couples TE and
    // TM; Transfer4 and Scatter4 count the same hybrid modes.
    for closed in [false, true] {
        let stack = waveguide(
            oriented(
                [lossless(4.0), lossless(4.4), lossless(3.6)],
                Rotation::from_euler(0.5, 0.6, 0.3),
            ),
            closed,
        );
        let exteriors = waveguide_exteriors(closed);
        let coupled = characteristic!(Transfer4::new(), &stack, k0, exteriors, Channel::Coupled);
        let scattered = characteristic!(Scatter4::new(), &stack, k0, exteriors, Channel::Coupled);
        let n = winding(&coupled, lo, hi);
        assert!(
            n.round() >= 4.0 && (n - n.round()).abs() < 1e-6,
            "closed {closed}: {n}"
        );
        assert!(winding(|beta| scattered(beta) / coupled(beta), lo, hi).abs() < 1e-6);
    }
}

#[test]
fn fields_satisfy_maxwell_equations_for_full_complex_tensors() {
    let tensor = |values: [(f64, f64); 3], rotation| {
        oriented(
            values.map(|(epsilon, mu)| material(c(epsilon, 0.2 * epsilon), c(mu, 0.05))),
            rotation,
        )
    };
    let thicknesses = [0.2, 0.15, 0.1];
    let layers = [
        (
            [(2.2, 1.0), (3.1, 1.3), (4.0, 0.8)],
            Rotation::from_euler(1.1, 0.4, -0.7),
        ),
        (
            [(1.5, 1.2), (2.8, 1.0), (2.0, 1.1)],
            Rotation::from_euler(-0.6, 1.2, 0.3),
        ),
        (
            [(3.5, 0.9), (2.1, 1.0), (2.6, 1.3)],
            Rotation::from_optic_axis(0.7, 2.1),
        ),
    ];
    let left = isotropic(material(c(1.3, 0.0), c(1.0, 0.0)));
    let right = tensor(
        [(2.0, 1.1), (2.4, 0.9), (3.0, 1.2)],
        Rotation::from_euler(0.2, 0.5, 0.9),
    );
    let stack = layers
        .into_iter()
        .zip(thicknesses)
        .fold(
            Stack::builder(left, right),
            |builder, ((values, rotation), thickness)| {
                builder.layer(tensor(values, rotation), Length::centimetres(thickness))
            },
        )
        .finalise();
    let (k0, beta) = (2.7, 1.3);
    let tensors = |index: usize| {
        let layer = &stack.layers_left_to_right()[index];
        let [epsilon, mu] =
            <RealJ0 as TensorConstitutiveLift<crate::RealAxis, _>>::constitutive_tensors(
                StackRegion::Layer(FiniteLayerIndex::new(index)),
                layer.material(),
                &RealJ0::lift_constant(arr0(c(k0, 0.0))),
            )
            .unwrap();
        [epsilon.full(), mu.full()]
            .map(|tensor| tensor.map(|row| row.map(|entry| entry.value()[()])))
    };

    for (backend, workspace) in retained(&stack, k0, beta).iter().enumerate() {
        for side in SIDES {
            for wave in 0..2 {
                let fields = |layer: usize, depth: f64| {
                    let [e, eta] = workspace
                        .fields(
                            side,
                            wave,
                            StackRegion::Layer(FiniteLayerIndex::new(layer)),
                            &RealJ0::lift_constant(arr0(c(depth, 0.0))),
                        )
                        .unwrap();
                    [e.map(|v| v.value()[()]), eta.map(|v| v.value()[()])]
                };
                for (layer, thickness) in thicknesses.iter().enumerate() {
                    let [epsilon, mu] = tensors(layer);
                    let apply = |tensor: [[C; 3]; 3], field: [C; 3]| {
                        std::array::from_fn::<C, 3, _>(|row| {
                            (0..3)
                                .map(|column| tensor[row][column] * field[column])
                                .sum()
                        })
                    };
                    let z = 0.4 * thickness;
                    let h = 1e-4;
                    let [e, eta] = fields(layer, z);
                    let [e_plus, eta_plus] = fields(layer, z + h);
                    let [e_minus, eta_minus] = fields(layer, z - h);
                    let dz = |plus: C, minus: C| (plus - minus) / (2.0 * h);
                    let ikx = c(0.0, beta);
                    let ik0 = c(0.0, k0);
                    let (mu_eta, epsilon_e) = (apply(mu, eta), apply(epsilon, e));
                    // ∇×E = i k₀ μ η, ∇×η = −i k₀ ε E with ∂_x = iβ, ∂_y = 0.
                    let residuals = [
                        -dz(e_plus[1], e_minus[1]) - ik0 * mu_eta[0],
                        dz(e_plus[0], e_minus[0]) - ikx * e[2] - ik0 * mu_eta[1],
                        ikx * e[1] - ik0 * mu_eta[2],
                        -dz(eta_plus[1], eta_minus[1]) + ik0 * epsilon_e[0],
                        dz(eta_plus[0], eta_minus[0]) - ikx * eta[2] + ik0 * epsilon_e[1],
                        ikx * eta[1] + ik0 * epsilon_e[2],
                    ];
                    let scale = e.iter().chain(&eta).map(|v| v.norm()).fold(0.0, f64::max) * k0;
                    for (index, residual) in residuals.iter().enumerate() {
                        assert!(
                            residual.norm() < 1e-6 * scale,
                            "backend {backend} {side:?} wave {wave} layer {layer} equation {index}: {residual}"
                        );
                    }
                }
                // Tangential E, H and normal D, B are continuous.
                for (layer, thickness) in thicknesses.iter().take(2).enumerate() {
                    let [e_left, eta_left] = fields(layer, *thickness);
                    let [e_right, eta_right] = fields(layer + 1, 0.0);
                    let ([epsilon_l, mu_l], [epsilon_r, mu_r]) =
                        (tensors(layer), tensors(layer + 1));
                    let normal = |tensor: [[C; 3]; 3], field: [C; 3]| {
                        (0..3)
                            .map(|column| tensor[2][column] * field[column])
                            .sum::<C>()
                    };
                    for (left, right) in [
                        (e_left[0], e_right[0]),
                        (e_left[1], e_right[1]),
                        (eta_left[0], eta_right[0]),
                        (eta_left[1], eta_right[1]),
                        (normal(epsilon_l, e_left), normal(epsilon_r, e_right)),
                        (normal(mu_l, eta_left), normal(mu_r, eta_right)),
                    ] {
                        assert_close(left, right, 1e-10, "continuity");
                    }
                }
            }
        }
    }
}

#[test]
fn layer_cutoff_refuses_scatter4_and_keeps_transfer4_smooth() {
    // At β = 1.5 k₀ the layer's forward and backward eigenvalues meet at 0.
    let stack = Stack::builder(isotropic(lossless(1.0)), isotropic(lossless(3.0)))
        .layer(
            oriented(
                [lossless(2.25), lossless(2.25), lossless(2.25)],
                Rotation::from_euler(0.3, 0.8, 0.1),
            ),
            Length::centimetres(0.3),
        )
        .finalise();
    let (k0, cutoff) = (2.0, 3.0);
    match RealAxisEvaluator::new(Scatter4::new()).evaluate(
        scalar_real_input(k0, cutoff),
        &stack,
        Channel::Coupled,
    ) {
        Err(RealAxisEvaluationError::Backend { source }) => assert_eq!(
            source,
            Scatter4Error::DefectiveLayer {
                layer: FiniteLayerIndex::new(0)
            }
        ),
        other => panic!("unexpected {:?}", other.map(|_| ())),
    }

    // Transfer4's characteristic and its first and second k∥ derivatives at
    // the cutoff agree with central differences.
    let evaluator =
        ComplexPlaneEvaluator::<RealJ2, _, _>::compile(&stack, Transfer4::new()).unwrap();
    let at = |beta: f64| {
        let coordinates = CanonicalCoordinates::new(
            <RealJ2 as SeedJet>::constant(arr0(c(k0, 0.0))),
            <RealJ2 as SeedJet>::variable(arr0(c(beta, 0.0)), 0).unwrap(),
        );
        let exterior = evaluate_stack_exterior_wavevectors::<ComplexPlane, _, _>(
            &coordinates,
            evaluator.stack(),
            Channel::Coupled,
        )
        .unwrap();
        evaluator
            .checked_determinant(&coordinates, &exterior, Channel::Coupled)
            .unwrap()
            .into_inner()
    };
    let h = 1e-4;
    let (centre, plus, minus) = (at(cutoff), at(cutoff + h), at(cutoff - h));
    let value = |jet: &RealJ2| jet.value()[()];
    assert_close(
        centre.first()[()],
        (value(&plus) - value(&minus)) / (2.0 * h),
        1e-6,
        "first",
    );
    assert_close(
        centre.second()[()],
        (value(&plus) - 2.0 * value(&centre) + value(&minus)) / (h * h),
        1e-4,
        "second",
    );
    assert_close(plus.second()[()], minus.second()[()], 1e-2, "continuity");
}

#[test]
fn thick_evanescent_layer_stays_finite_in_scatter4() {
    let stack = |thickness: f64| {
        Stack::builder(isotropic(lossless(2.5)), isotropic(lossless(2.5)))
            .layer(
                oriented(
                    [lossless(1.2), lossless(1.4), lossless(1.1)],
                    Rotation::from_euler(0.3, 0.8, 0.1),
                ),
                Length::centimetres(thickness),
            )
            .finalise()
    };
    let (k0, beta) = (10.0, 15.0);
    // Moderate: both agree.
    let [transfer, scatter] = coupled(&stack(1.0), k0, beta, Channel::Coupled);
    for side in SIDES {
        for (row, column) in [(0, 0), (0, 1), (1, 1)] {
            assert_close(
                entry(jones(&scatter, side).transmission(row, column)),
                entry(jones(&transfer, side).transmission(row, column)),
                1e-8,
                "moderate",
            );
        }
    }
    // Optically very thick: Transfer4 overflows, Scatter4 stays finite with
    // vanishing transmission.
    let thick = stack(200.0);
    match RealAxisEvaluator::new(Transfer4::new()).evaluate(
        scalar_real_input(k0, beta),
        &thick,
        Channel::Coupled,
    ) {
        Err(RealAxisEvaluationError::Backend { source }) => {
            assert_eq!(source, Transfer4Error::NonFinite)
        }
        other => panic!("unexpected {:?}", other.map(|_| ())),
    }
    let scatter = solve(Scatter4::new(), &thick, k0, beta, Channel::Coupled);
    assert!(scatter.determinant().value().all_finite());
    for side in SIDES {
        let jones = jones(&scatter, side);
        for row in 0..2 {
            for column in 0..2 {
                assert!(entry(jones.reflection(row, column)).is_finite());
                assert!(entry(jones.transmission(row, column)).norm() < 1e-100);
            }
        }
    }
}

#[test]
fn unsupported_and_non_invariant_channels_are_refused() {
    let stack = full_stack();
    let input = || scalar_real_input(2.0, 0.5);
    match RealAxisEvaluator::new(Scatter2::new()).evaluate(
        input(),
        &diagonal_stack(),
        Channel::Coupled,
    ) {
        Err(RealAxisEvaluationError::Backend { source }) => assert_eq!(
            source,
            crate::backend::Scatter2Error::Channel(ChannelError::Unsupported {
                channel: Channel::Coupled
            })
        ),
        other => panic!("unexpected {:?}", other.map(|_| ())),
    }
    let thin = Stack::builder(isotropic(lossless(1.0)), isotropic(lossless(1.0)))
        .layer(
            stack.layers_left_to_right()[0].material().clone(),
            Length::centimetres(0.1),
        )
        .finalise();
    match RealAxisEvaluator::new(Transfer4::new()).evaluate(
        input(),
        &thin,
        Polarisation::TransverseElectric,
    ) {
        Err(RealAxisEvaluationError::Backend { source }) => assert_eq!(
            source,
            Transfer4Error::Channel(ChannelError::NotInvariant {
                region: StackRegion::Layer(FiniteLayerIndex::new(0)),
                channel: Channel::TransverseElectric,
            })
        ),
        other => panic!("unexpected {:?}", other.map(|_| ())),
    }
}

#[test]
fn derivatives_agree_between_backends_and_with_central_differences() {
    use crate::Parameter;
    let stack = full_stack();
    let (k0, beta) = (2.6, 1.9);
    for parameter in [
        Parameter::Spectral,
        Parameter::InPlane,
        Parameter::LayerThickness(FiniteLayerIndex::new(1)),
    ] {
        let input = || scalar_real_input(k0, beta);
        let reflection = |solution: Coupled<RealJ2>| {
            let jones = solution
                .entries()
                .project_jones_amplitudes(solution.context(), IncidentSide::Left);
            [
                jones.reflection(0, 1).unwrap().clone(),
                jones.transmission(1, 0).unwrap().clone(),
            ]
        };
        let transfer = reflection(
            RealAxisEvaluator::new(Transfer4::new())
                .evaluate_second(input(), &stack, Channel::Coupled, parameter)
                .unwrap()
                .into_parts()
                .0,
        );
        let scatter = reflection(
            RealAxisEvaluator::new(Scatter4::new())
                .evaluate_second(input(), &stack, Channel::Coupled, parameter)
                .unwrap()
                .into_parts()
                .0,
        );
        for (t, s) in transfer.iter().zip(&scatter) {
            assert_close(
                s.first()[()],
                t.first()[()],
                1e-8,
                &format!("{parameter:?} first"),
            );
            assert_close(
                s.second()[()],
                t.second()[()],
                1e-7,
                &format!("{parameter:?} second"),
            );
        }

        // Central differences of the values.
        let h = 1e-4;
        let shifted = |offset: f64| {
            let (k, b, thickness) = match parameter {
                Parameter::Spectral => (k0 + offset, beta, 0.13),
                Parameter::InPlane => (k0, beta + offset, 0.13),
                Parameter::LayerThickness(_) => (k0, beta, 0.13 + offset),
            };
            let solution = solve(
                Transfer4::new(),
                &full_stack_with(thickness),
                k,
                b,
                Channel::Coupled,
            );
            entry(jones(&solution, IncidentSide::Left).reflection(0, 1))
        };
        let derivative = (shifted(h) - shifted(-h)) / (2.0 * h);
        assert_close(
            transfer[0].first()[()],
            derivative,
            1e-6,
            &format!("{parameter:?} difference"),
        );
    }
}

#[test]
fn continued_scatter4_split_counts_like_transfer4_away_from_the_real_axis() {
    // A lossy rotated-uniaxial core: its hybrid modes lie well off the real
    // axis. The rectangle stops short of the core's own branch points near
    // β = k₀ √ε (about 7.65 + 0.84i and beyond).
    let k0 = 4.0;
    let lossy = |epsilon: f64| ComplexQuadratic {
        epsilon: [c(epsilon, 0.8), C::default(), C::default()],
        mu: [c(1.0, 0.0), C::default(), C::default()],
    };
    let stack = waveguide(
        oriented(
            [lossy(4.0), lossy(4.4), lossy(3.6)],
            Rotation::from_euler(0.5, 0.6, 0.3),
        ),
        false,
    );
    let exteriors = waveguide_exteriors(false);
    let (lo, hi) = (c(5.85, 0.1), c(7.3, 3.0));
    let centre = (lo + hi) * 0.5;
    let transfer = characteristic!(Transfer4::new(), &stack, k0, exteriors, Channel::Coupled);
    let continued = characteristic!(
        Scatter4::continued_from(c(k0, 0.0), centre),
        &stack,
        k0,
        exteriors,
        Channel::Coupled
    );
    let n = winding(&transfer, lo, hi);
    assert!(
        n.round() >= 2.0 && (n - n.round()).abs() < 1e-6,
        "Transfer4 {n}"
    );
    assert!((winding(&continued, lo, hi) - n).abs() < 1e-6);
    assert!(winding(|beta| continued(beta) / transfer(beta), lo, hi).abs() < 1e-6);
}

#[test]
fn continuation_through_a_layer_cutoff_is_refused() {
    // The layer's forward and backward eigenvalues meet at β = 1.5 k₀.
    let stack = Stack::builder(isotropic(lossless(1.0)), isotropic(lossless(3.0)))
        .layer(
            oriented(
                [lossless(2.25), lossless(2.25), lossless(2.25)],
                Rotation::from_euler(0.3, 0.8, 0.1),
            ),
            Length::centimetres(0.3),
        )
        .finalise();
    let k0 = 2.0;
    let evaluator = ComplexPlaneEvaluator::<HoloJ0, _, _>::compile(
        &stack,
        Scatter4::continued_from(c(k0, 0.0), c(2.4, 0.0)),
    )
    .unwrap();
    let at = |beta: C| {
        let coordinates = CanonicalCoordinates::new(
            <HoloJ0 as SeedJet>::constant(arr0(c(k0, 0.0))),
            <HoloJ0 as SeedJet>::constant(arr0(beta)),
        );
        evaluator.checked_determinant(
            &coordinates,
            &evaluate_stack_exterior_wavevectors::<ComplexPlane, _, _>(
                &coordinates,
                evaluator.stack(),
                Channel::Coupled,
            )
            .unwrap(),
            Channel::Coupled,
        )
    };
    assert!(at(c(2.6, 0.0)).is_ok());
    match at(c(3.3, 0.0)) {
        Err(crate::DeterminantEvaluationError::Backend { source }) => assert_eq!(
            source,
            Scatter4Error::Continuation {
                layer: FiniteLayerIndex::new(0)
            }
        ),
        other => panic!("unexpected {:?}", other.map(|_| ())),
    }
}
