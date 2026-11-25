//! Stacks closed by walls and PMLs: characteristic, zeros, fields and
//! derivatives.
use ndarray::{Ix0, arr0};
use num_complex::Complex64 as C;

use crate::{
    CanonicalCoordinates, ComplexPlane, ComplexPlaneEvaluator, Constant, ExteriorSampling,
    ExteriorWavevectors, FieldSampling, FieldSamplingError, FiniteLayerIndex, LayerSampling,
    ModalFieldBackend, Polarisation, ScalarAlgebra, Scatter2, Stack, StackCompileError, StackSide,
    Transfer2, ValidationConfig,
    backend::{TerminationError, evaluate_stack_exterior_wavevectors},
    compile_canonical_constant_stack, compile_canonical_thickness_stack,
    input::canonical::{CanonicalLayer, CanonicalTermination},
    stack::{Pml, StretchProfile, Termination, Wall, WallKind},
    test_support::jet::HoloJ2,
    units::Length,
};

type J = HoloJ2;
type Material = Constant<f64>;

const TE: Polarisation = Polarisation::TransverseElectric;
const TM: Polarisation = Polarisation::TransverseMagnetic;
const POLARISATIONS: [Polarisation; 2] = [TE, TM];
const KINDS: [WallKind; 2] = [WallKind::PerfectElectric, WallKind::PerfectMagnetic];

#[derive(Copy, Clone, Debug)]
enum Backend {
    Transfer,
    Scatter,
}

const BACKENDS: [Backend; 2] = [Backend::Transfer, Backend::Scatter];

fn c(re: f64) -> C {
    C::new(re, 0.0)
}

fn relative(a: C, b: C) -> f64 {
    (a - b).norm() / a.norm().max(b.norm()).max(1e-300)
}

fn dirichlet(kind: WallKind, polarisation: Polarisation) -> bool {
    matches!(
        (kind, polarisation),
        (WallKind::PerfectElectric, Polarisation::TransverseElectric)
            | (WallKind::PerfectMagnetic, Polarisation::TransverseMagnetic)
    )
}

fn wall(kind: WallKind, z_cm: f64) -> Termination<f64> {
    Wall::new(kind, Length::centimetres(z_cm)).into()
}

fn pml(start: f64, thickness: f64, order: f64, strength: C, backing: WallKind) -> Termination<f64> {
    Pml::new(
        Length::centimetres(start),
        Length::centimetres(thickness),
        StretchProfile::polynomial(order, strength),
        backing,
    )
    .into()
}

/// The characteristic on a canonical stack, with the open sides' principal
/// exterior wavevectors.
fn determinant(
    stack: crate::input::CanonicalStack<Material, J>,
    backend: Backend,
    k0: J,
    beta: J,
    polarisation: Polarisation,
) -> Result<J, String> {
    let coordinates = CanonicalCoordinates::new(k0, beta);
    let exterior = evaluate_stack_exterior_wavevectors::<ComplexPlane, _, J>(
        &coordinates,
        &stack,
        polarisation,
    )
    .unwrap();
    match backend {
        Backend::Transfer => ComplexPlaneEvaluator::from_canonical_stack(stack, Transfer2::new())
            .checked_determinant(&coordinates, &exterior, polarisation)
            .map(|d| d.into_inner())
            .map_err(|e| format!("{e:?}")),
        Backend::Scatter => ComplexPlaneEvaluator::from_canonical_stack(stack, Scatter2::new())
            .checked_determinant(&coordinates, &exterior, polarisation)
            .map(|d| d.into_inner())
            .map_err(|e| format!("{e:?}")),
    }
}

/// The same stack with materials behind analytical handles.
fn analytical(stack: &Stack<Material, f64>) -> crate::AnalyticalMaterialStack<C> {
    let handle = |material: &Material| crate::AnalyticalMaterialHandle::new(*material);
    let mut builder = Stack::builder(
        handle(stack.left_exterior()),
        handle(stack.right_exterior()),
    );
    for layer in stack.layers_left_to_right() {
        builder = builder.layer(handle(layer.material()), layer.thickness());
    }
    builder
        .left_termination(*stack.left_termination())
        .right_termination(*stack.right_termination())
        .finalise()
}

fn compile(stack: &Stack<Material, f64>) -> crate::input::CanonicalStack<Material, J> {
    compile_canonical_constant_stack::<_, J>(stack, Ix0(), &ValidationConfig::permissive()).unwrap()
}

/// The characteristic as a jet in β.
fn characteristic(
    stack: &Stack<Material, f64>,
    backend: Backend,
    k0: C,
    beta: C,
    polarisation: Polarisation,
) -> J {
    let (k0, beta) = (J::constant(arr0(k0)), J::variable(arr0(beta)));
    determinant(compile(stack), backend, k0, beta, polarisation).unwrap()
}

/// Newton refinement in β from `start`.
fn refine(f: impl Fn(C) -> J, start: C) -> C {
    let mut beta = start;
    for _ in 0..60 {
        let d = f(beta);
        let step = d.value()[()] / d.first()[()];
        beta -= step;
        if step.norm() <= 1e-15 * beta.norm().max(1.0) {
            break;
        }
    }
    beta
}

/// A homogeneous medium of relative permittivity `EPS` between walls at
/// −0.3 cm and 0.9 cm, with one finite layer of the same medium on [0, 0.5].
const EPS: f64 = 2.25;
const K0: f64 = 3.0;
const WIDTH: f64 = 1.2;

/// Two layers, one lossy, between different exteriors.
fn lossy_stack(left: Termination<f64>, right: Termination<f64>) -> Stack<Material, f64> {
    Stack::builder(Constant::dielectric(1.5), Constant::dielectric(2.0))
        .layer(Constant::dielectric(4.0), Length::centimetres(0.4))
        .layer(Constant::new(2.0, 1.3), Length::centimetres(0.3))
        .left_termination(left)
        .right_termination(right)
        .finalise()
}

fn box_stack(left: Termination<f64>, right: Termination<f64>) -> Stack<Material, f64> {
    Stack::builder(Constant::dielectric(EPS), Constant::dielectric(EPS))
        .layer(Constant::dielectric(EPS), Length::centimetres(0.5))
        .left_termination(left)
        .right_termination(right)
        .finalise()
}

/// The parallel-plate ladder of normal wavenumbers κ_m for two wall
/// conditions on a box of (possibly complex) width `width`.
fn ladder(left_dirichlet: bool, right_dirichlet: bool, width: C, count: usize) -> Vec<C> {
    let pi = std::f64::consts::PI;
    (0..count)
        .filter_map(|m| {
            let m = m as f64;
            match (left_dirichlet, right_dirichlet) {
                (true, true) => (m >= 1.0).then(|| c(m * pi) / width),
                (false, false) => Some(c(m * pi) / width),
                _ => Some(c((m + 0.5) * pi) / width),
            }
        })
        .collect()
}

fn beta_of(kappa: C) -> C {
    (c(EPS * K0 * K0) - kappa * kappa).sqrt()
}

#[test]
fn parallel_plate_zeros_follow_the_ladder() {
    for polarisation in POLARISATIONS {
        for left in KINDS {
            for right in KINDS {
                let stack = box_stack(wall(left, -0.3), wall(right, 0.9));
                let kappas = ladder(
                    dirichlet(left, polarisation),
                    dirichlet(right, polarisation),
                    c(WIDTH),
                    5,
                );
                for kappa in kappas {
                    let expected = beta_of(kappa);
                    for backend in BACKENDS {
                        // Scatter2's directional basis degenerates at κ = 0.
                        if kappa.norm() == 0.0 && matches!(backend, Backend::Scatter) {
                            continue;
                        }
                        let start = expected * C::new(1.0 + 2e-3, 1e-3) + C::new(0.0, 1e-3);
                        let found = refine(
                            |beta| characteristic(&stack, backend, c(K0), beta, polarisation),
                            start,
                        );
                        assert!(
                            relative(found, expected) < 1e-10,
                            "{polarisation:?} {left:?}/{right:?} {backend:?} κ={kappa}: \
                             found {found}, expected {expected}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn closed_characteristic_is_even_and_continuous_across_cuts() {
    let stack = lossy_stack(
        wall(WallKind::PerfectMagnetic, -0.2),
        pml(0.8, 0.6, 2.0, C::new(0.5, 3.0), WallKind::PerfectElectric),
    );
    let k0 = c(K0);
    for polarisation in POLARISATIONS {
        for backend in BACKENDS {
            for beta in [C::new(2.0, 0.3), C::new(5.0, -0.7), C::new(0.4, 3.0)] {
                let plus = characteristic(&stack, backend, k0, beta, polarisation);
                let minus = characteristic(&stack, backend, k0, -beta, polarisation);
                assert!(
                    relative(plus.value()[()], minus.value()[()]) < 1e-10,
                    "{polarisation:?} {backend:?} {beta}"
                );
            }
            // Across the principal cut of each material's κ²: εμk0² − β² on the
            // negative real axis, approached from both sides.
            for eps_mu in [1.5, 2.0, 4.0, 2.6] {
                let beta_squared = eps_mu * K0 * K0 + 3.0;
                for tiny in [1e-9, -1e-9] {
                    let beta = C::new(beta_squared, tiny).sqrt();
                    let above = characteristic(&stack, backend, k0, beta, polarisation);
                    let below = characteristic(
                        &stack,
                        backend,
                        k0,
                        C::new(beta_squared, -tiny).sqrt(),
                        polarisation,
                    );
                    // A step across a branch cut would be of the order of the
                    // value itself; compare on the local scale |D| + |β D'|.
                    let scale = above.value()[()].norm() + (beta * above.first()[()]).norm();
                    assert!(
                        (above.value()[()] - below.value()[()]).norm() < 1e-7 * scale,
                        "{polarisation:?} {backend:?} εμ={eps_mu}"
                    );
                }
            }
        }
    }
}

#[test]
fn backends_agree_on_closed_stacks() {
    let terminations = [
        (wall(WallKind::PerfectElectric, -0.1), Termination::Open),
        (Termination::Open, wall(WallKind::PerfectMagnetic, 0.95)),
        (
            pml(-0.05, 0.7, 1.0, C::new(0.0, 2.0), WallKind::PerfectMagnetic),
            wall(WallKind::PerfectElectric, 0.7),
        ),
        (
            pml(-0.2, 0.5, 2.0, C::new(1.0, 3.0), WallKind::PerfectElectric),
            pml(0.9, 0.4, 3.0, C::new(0.5, 2.5), WallKind::PerfectMagnetic),
        ),
    ];
    for (left, right) in terminations {
        let stack = Stack::builder(Constant::dielectric(1.5), Constant::new(2.0, 1.1))
            .layer(Constant::dielectric(4.0), Length::centimetres(0.4))
            .layer(Constant::new(2.0, 1.3), Length::centimetres(0.3))
            .left_termination(left)
            .right_termination(right)
            .finalise();
        for polarisation in POLARISATIONS {
            for beta in [C::new(2.0, 0.3), C::new(5.5, -0.2), C::new(1.0, 2.0)] {
                let transfer = characteristic(&stack, Backend::Transfer, c(K0), beta, polarisation);
                let scatter = characteristic(&stack, Backend::Scatter, c(K0), beta, polarisation);
                assert!(relative(transfer.value()[()], scatter.value()[()]) < 1e-10);
                assert!(relative(transfer.first()[()], scatter.first()[()]) < 1e-9);
                assert!(relative(transfer.second()[()], scatter.second()[()]) < 1e-8);
            }
        }
    }
}

/// The canonical stack in which the PML of `stack` is replaced by a layer of
/// the exterior material with complex thickness τ, both sides closed by walls
/// on the outermost interfaces.
fn complex_thickness_box(
    left_gap: f64,
    layers: &[(Material, f64)],
    right_gap: f64,
    tau: C,
    left_kind: WallKind,
    right_kind: WallKind,
    exterior: Material,
) -> crate::input::CanonicalStack<Material, J> {
    let constant = |value: C| J::constant(arr0(value));
    let mut canonical = vec![CanonicalLayer::new(exterior, constant(c(left_gap)))];
    canonical.extend(
        layers
            .iter()
            .map(|(material, thickness)| CanonicalLayer::new(*material, constant(c(*thickness)))),
    );
    canonical.push(CanonicalLayer::new(exterior, constant(c(right_gap))));
    canonical.push(CanonicalLayer::new(exterior, constant(tau)));
    let wall = |kind| CanonicalTermination::Wall {
        kind,
        position_cm: constant(c(0.0)),
        gap_cm: constant(c(0.0)),
    };
    crate::input::CanonicalStack::new(exterior, exterior, canonical)
        .with_terminations(wall(left_kind), wall(right_kind))
}

#[test]
fn pml_equals_complex_thickness_layer_up_to_constant_normalisation() {
    let exterior = Constant::dielectric(1.5);
    let layers = [
        (Constant::dielectric(4.0), 0.4),
        (Constant::new(2.0, 1.3), 0.3),
    ];
    let (thickness, order, strength) = (0.5, 2.0, C::new(1.0, 3.0));
    let tau = c(thickness) * (c(1.0) + strength / (order + 1.0));
    let k0 = c(K0);
    let reference = (k0 * k0 * 1.5).sqrt();
    let normalisation = (C::i() * reference * tau).exp();

    for left_kind in KINDS {
        for backing in KINDS {
            let stack = Stack::builder(exterior, exterior)
                .layer(layers[0].0, Length::centimetres(0.4))
                .layer(layers[1].0, Length::centimetres(0.3))
                .left_termination(wall(left_kind, -0.2))
                .right_termination(pml(0.8, thickness, order, strength, backing))
                .finalise();
            for polarisation in POLARISATIONS {
                for backend in BACKENDS {
                    for beta in [C::new(2.0, 0.3), C::new(5.5, -0.2), C::new(0.3, 2.0)] {
                        let closed = characteristic(&stack, backend, k0, beta, polarisation);
                        let equivalent = determinant(
                            complex_thickness_box(
                                0.2, &layers, 0.1, tau, left_kind, backing, exterior,
                            ),
                            backend,
                            J::constant(arr0(k0)),
                            J::variable(arr0(beta)),
                            polarisation,
                        )
                        .unwrap();
                        let ratio = closed.value()[()] / equivalent.value()[()];
                        assert!(
                            relative(ratio, normalisation) < 1e-10,
                            "{polarisation:?} {backend:?} {beta}: ratio {ratio}"
                        );
                        assert!(
                            relative(closed.first()[()], normalisation * equivalent.first()[()])
                                < 1e-9
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn equal_optical_thickness_profiles_give_equal_characteristics() {
    // τ = d (1 + c/(m+1)) = d (1 + i) for both profiles.
    let profiles = [(2.0, C::new(0.0, 3.0)), (0.0, C::new(0.0, 1.0))];
    let stacks = profiles.map(|(order, strength)| {
        Stack::builder(Constant::dielectric(1.0), Constant::dielectric(2.0))
            .layer(Constant::dielectric(4.0), Length::centimetres(0.4))
            .right_termination(pml(0.6, 0.5, order, strength, WallKind::PerfectElectric))
            .finalise()
    });
    for polarisation in POLARISATIONS {
        for backend in BACKENDS {
            for beta in [C::new(2.0, 0.3), C::new(5.5, -0.2)] {
                let a = characteristic(&stacks[0], backend, c(K0), beta, polarisation);
                let b = characteristic(&stacks[1], backend, c(K0), beta, polarisation);
                assert!(relative(a.value()[()], b.value()[()]) < 1e-12);
            }
        }
    }
}

#[test]
fn thick_pmls_stay_finite_agree_and_keep_the_complex_box_ladder() {
    let strength = C::new(0.5, 4.0);
    let order = 2.0;
    for thickness in [2.0, 10.0, 25.0, 50.0] {
        let tau = c(thickness) * (c(1.0) + strength / (order + 1.0));
        for polarisation in POLARISATIONS {
            for left in KINDS {
                for backing in KINDS {
                    let stack = box_stack(
                        wall(left, -0.3),
                        pml(0.9, thickness, order, strength, backing),
                    );
                    // |Im κ_P τ| spans tens to hundreds over these points.
                    for beta in [C::new(0.5, 0.1), C::new(4.0, -0.3), C::new(1.0, 2.5)] {
                        let transfer =
                            characteristic(&stack, Backend::Transfer, c(K0), beta, polarisation);
                        let scatter =
                            characteristic(&stack, Backend::Scatter, c(K0), beta, polarisation);
                        assert!(transfer.all_finite() && scatter.all_finite());
                        assert!(
                            relative(transfer.value()[()], scatter.value()[()]) < 1e-8,
                            "d={thickness} {polarisation:?} {beta}: {} vs {}",
                            transfer.value()[()],
                            scatter.value()[()]
                        );
                    }
                    let width = c(WIDTH) + tau;
                    let spacing = c(std::f64::consts::PI) / width;
                    for kappa in ladder(
                        dirichlet(left, polarisation),
                        dirichlet(backing, polarisation),
                        width,
                        4,
                    ) {
                        let expected = beta_of(kappa);
                        for backend in BACKENDS {
                            if kappa.norm() == 0.0 && matches!(backend, Backend::Scatter) {
                                continue;
                            }
                            // Start a small fraction of the ladder spacing away.
                            let found = refine(
                                |beta| characteristic(&stack, backend, c(K0), beta, polarisation),
                                beta_of(kappa + spacing * C::new(1e-2, 1e-2)),
                            );
                            assert!(
                                relative(found, expected) < 1e-9,
                                "d={thickness} {polarisation:?} {backend:?}: {found} vs {expected}"
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn out_of_range_pml_scale_is_a_typed_failure() {
    let thickness = 200.0;
    let strength = C::new(0.0, 4.0);
    let stack = box_stack(
        wall(WallKind::PerfectElectric, -0.3),
        pml(0.9, thickness, 2.0, strength, WallKind::PerfectElectric),
    );
    let tau = c(thickness) * (c(1.0) + strength / 3.0);
    let at = |backend, beta: C| {
        let (k0, beta) = (J::constant(arr0(c(K0))), J::variable(arr0(beta)));
        determinant(compile(&stack), backend, k0, beta, TE)
    };
    // Overflow far from β = 0; underflow where κ_P τ is real.
    let kappa_real_product = (-C::i() * tau.arg()).exp() * 2.0;
    for beta in [C::new(0.0, 20.0), beta_of(kappa_real_product)] {
        for backend in BACKENDS {
            let error = at(backend, beta).unwrap_err();
            assert!(
                error.contains("PmlScaleOutOfRange") && error.contains("Right"),
                "{backend:?} {beta}: {error}"
            );
        }
    }
    // Near β = 0 the normalised scale is of order one.
    for backend in BACKENDS {
        at(backend, c(0.1)).unwrap();
    }
}

#[test]
fn exterior_wavevector_pattern_must_match_the_terminations() {
    let walled = box_stack(wall(WallKind::PerfectElectric, -0.3), Termination::Open);
    let pml_closed = box_stack(
        Termination::Open,
        pml(0.9, 0.5, 2.0, C::new(0.0, 2.0), WallKind::PerfectElectric),
    );
    let coordinates =
        CanonicalCoordinates::new(J::constant(arr0(c(K0))), J::constant(arr0(c(1.0))));
    let kappa = || J::constant(arr0(c(4.0)));
    // The correct patterns evaluate.
    for (stack, exterior) in [
        (&walled, ExteriorWavevectors::closed_left(kappa())),
        (&pml_closed, ExteriorWavevectors::closed_right(kappa())),
    ] {
        ComplexPlaneEvaluator::from_canonical_stack(compile(stack), Transfer2::new())
            .checked_determinant(&coordinates, &exterior, TE)
            .unwrap();
        ComplexPlaneEvaluator::from_canonical_stack(compile(stack), Scatter2::new())
            .checked_determinant(&coordinates, &exterior, TE)
            .unwrap();
    }
    let cases = [
        // A wavevector for a wall side.
        (
            &walled,
            ExteriorWavevectors::new(kappa(), kappa()),
            StackSide::Left,
            false,
        ),
        // No wavevector for an open side.
        (
            &walled,
            ExteriorWavevectors::closed(),
            StackSide::Right,
            true,
        ),
        // A wavevector for a PML side.
        (
            &pml_closed,
            ExteriorWavevectors::new(kappa(), kappa()),
            StackSide::Right,
            false,
        ),
    ];
    for (stack, exterior, side, open) in cases {
        let expected = TerminationError::ExteriorBranchMismatch { side, open };
        let transfer =
            ComplexPlaneEvaluator::from_canonical_stack(compile(stack), Transfer2::new())
                .determinant(&coordinates, &exterior, TE)
                .map(|_| ())
                .unwrap_err();
        assert_eq!(
            transfer,
            crate::backend::Transfer2Error::Termination(expected.clone())
        );
        let scatter = ComplexPlaneEvaluator::from_canonical_stack(compile(stack), Scatter2::new())
            .determinant(&coordinates, &exterior, TE)
            .map(|_| ())
            .unwrap_err();
        assert_eq!(
            scatter,
            crate::backend::Scatter2Error::Termination(expected)
        );
    }
}

#[test]
fn lossless_walled_stack_has_a_real_spectrum() {
    let stack = Stack::builder(Constant::dielectric(1.0), Constant::dielectric(1.5))
        .layer(Constant::dielectric(4.0), Length::centimetres(0.4))
        .layer(Constant::dielectric(2.0), Length::centimetres(0.3))
        .left_termination(wall(WallKind::PerfectMagnetic, -0.5))
        .right_termination(wall(WallKind::PerfectElectric, 1.1))
        .finalise();
    for polarisation in POLARISATIONS {
        let mut found = Vec::new();
        for re in [0.5, 2.0, 3.5, 5.0, 6.5] {
            for im in [-1.5, 0.3, 1.5, 4.0] {
                let beta = refine(
                    |beta| characteristic(&stack, Backend::Transfer, c(K0), beta, polarisation),
                    C::new(re, im),
                );
                let residual = characteristic(&stack, Backend::Transfer, c(K0), beta, polarisation);
                if residual.value()[()].norm()
                    < 1e-8 * residual.first()[()].norm() * beta.norm().max(1.0)
                {
                    found.push(beta);
                }
            }
        }
        assert!(found.len() > 5, "{polarisation:?}: {found:?}");
        for beta in found {
            let squared = beta * beta;
            assert!(
                squared.im.abs() < 1e-8 * squared.norm().max(1.0),
                "{polarisation:?}: β² = {squared}"
            );
        }
    }
}

/// Real roots of `f` on (lo, hi), by scanning for sign changes and bisecting.
fn real_roots(f: impl Fn(f64) -> f64, lo: f64, hi: f64) -> Vec<f64> {
    let steps = 4000;
    let mut roots = Vec::new();
    let at = |i: usize| lo + (hi - lo) * i as f64 / steps as f64;
    for i in 1..steps - 1 {
        let (mut a, mut b) = (at(i), at(i + 1));
        if f(a).signum() == f(b).signum() {
            continue;
        }
        for _ in 0..200 {
            let mid = 0.5 * (a + b);
            if mid == a || mid == b {
                break;
            }
            if f(mid).signum() == f(a).signum() {
                a = mid;
            } else {
                b = mid;
            }
        }
        roots.push(0.5 * (a + b));
    }
    roots
}

#[test]
fn wall_reproduces_odd_or_even_modes_of_the_mirrored_stack() {
    // A slab of thickness a on a wall at z = 0, open on the right, against the
    // symmetric slab of thickness 2a in the same cladding. The guided modes of
    // the symmetric slab solve, with q² = 4k0² − β², p² = β² − k0² and weights
    // w = 1/μ (TE) or 1/ε (TM):
    //   even: w1 q sin(qa) − w2 p cos(qa) = 0,   odd: w1 q cos(qa) + w2 p sin(qa) = 0.
    let a = 0.35;
    for polarisation in POLARISATIONS {
        let w1 = match polarisation {
            TE => 1.0,
            TM => 0.25,
        };
        let q = |beta: f64| (4.0 * K0 * K0 - beta * beta).sqrt();
        let p = |beta: f64| (beta * beta - K0 * K0).sqrt();
        let even = real_roots(
            |b| w1 * q(b) * (q(b) * a).sin() - p(b) * (q(b) * a).cos(),
            K0,
            2.0 * K0,
        );
        let odd = real_roots(
            |b| w1 * q(b) * (q(b) * a).cos() + p(b) * (q(b) * a).sin(),
            K0,
            2.0 * K0,
        );
        assert!(!even.is_empty() && !odd.is_empty(), "{polarisation:?}");
        for kind in KINDS {
            let walled = Stack::builder(Constant::dielectric(1.0), Constant::dielectric(1.0))
                .layer(Constant::dielectric(4.0), Length::centimetres(a))
                .left_termination(wall(kind, 0.0))
                .finalise();
            let expected_odd = dirichlet(kind, polarisation);
            for (modes, odd_modes) in [(&even, false), (&odd, true)] {
                for &beta in modes {
                    for backend in BACKENDS {
                        let value = characteristic(&walled, backend, c(K0), c(beta), polarisation);
                        let scale = value.first()[()].norm() * beta;
                        let is_mode = value.value()[()].norm() < 1e-9 * scale;
                        assert_eq!(
                            is_mode,
                            odd_modes == expected_odd,
                            "{polarisation:?} {kind:?} {backend:?} β={beta}: D={}",
                            value.value()[()]
                        );
                    }
                }
            }
        }
    }
}

fn component(
    fields: &crate::OwnedModeFields,
    polarisation: Polarisation,
    dirichlet: bool,
    index: usize,
) -> C {
    let e = fields.quantity().electric();
    let h = fields.quantity().magnetic();
    match (polarisation, dirichlet) {
        (TE, true) => e.y()[index],
        (TE, false) => h.x()[index],
        (TM, true) => h.y()[index],
        (TM, false) => e.x()[index],
    }
}

fn normal_flux(fields: &crate::OwnedModeFields, index: usize) -> f64 {
    let e = fields.quantity().electric();
    let h = fields.quantity().magnetic();
    0.5 * (e.x()[index] * h.y()[index].conj() - e.y()[index] * h.x()[index].conj()).re
}

fn largest(fields: &crate::OwnedModeFields) -> f64 {
    let e = fields.quantity().electric();
    let h = fields.quantity().magnetic();
    [e.x(), e.y(), e.z(), h.x(), h.y(), h.z()]
        .into_iter()
        .flatten()
        .map(|v| v.norm())
        .fold(0.0, f64::max)
}

#[test]
fn reconstructed_modes_satisfy_the_wall_conditions() {
    for polarisation in POLARISATIONS {
        for left in KINDS {
            for right in KINDS {
                let stack = Stack::builder(Constant::dielectric(1.5), Constant::dielectric(2.0))
                    .layer(Constant::dielectric(4.0), Length::centimetres(0.4))
                    .layer(Constant::dielectric(2.2), Length::centimetres(0.3))
                    .left_termination(wall(left, -0.25))
                    .right_termination(wall(right, 0.95))
                    .finalise();
                // On the real axis the characteristic of a lossless walled
                // stack has a constant phase; take the highest real zero.
                let real = |beta: f64| {
                    let d = characteristic(&stack, Backend::Transfer, c(K0), c(beta), polarisation)
                        .value()[()];
                    d.re + d.im
                };
                let beta = c(*real_roots(real, 0.1, 2.0 * K0).last().unwrap());
                let d = characteristic(&stack, Backend::Transfer, c(K0), beta, polarisation);
                assert!(d.value()[()].norm() < 1e-9 * d.first()[()].norm() * beta.norm());
                let sampling = FieldSampling::new()
                    .left_exterior(ExteriorSampling::distances([
                        Length::centimetres(0.25),
                        Length::centimetres(0.1),
                    ]))
                    .layer(0, LayerSampling::uniform_with_boundaries(5, true, true))
                    .layer(1, LayerSampling::uniform_with_boundaries(5, true, true))
                    .right_exterior(ExteriorSampling::distances([
                        Length::centimetres(0.1),
                        Length::centimetres(0.25),
                    ]));
                for backend in [
                    &Transfer2::new() as &dyn ModalFieldBackend,
                    &Scatter2::new(),
                ] {
                    let fields = backend
                        .unscaled_fields(
                            &analytical(&stack),
                            c(K0),
                            beta,
                            [None, None],
                            polarisation,
                            &sampling,
                            c(1.0),
                        )
                        .unwrap();
                    let fields = &fields.fields;
                    let scale = largest(fields);
                    let last = fields.sampling().len() - 1;
                    let at_left = component(fields, polarisation, dirichlet(left, polarisation), 0);
                    let at_right =
                        component(fields, polarisation, dirichlet(right, polarisation), last);
                    assert!(
                        at_left.norm() < 1e-9 * scale && at_right.norm() < 1e-9 * scale,
                        "{polarisation:?} {left:?}/{right:?}: {at_left} {at_right} scale {scale}"
                    );
                    // A real-β mode of a lossless closed stack carries no
                    // normal flux anywhere, in particular not into the walls.
                    for index in 0..=last {
                        assert!(normal_flux(fields, index).abs() < 1e-9 * scale * scale);
                    }
                }
            }
        }
    }
}

#[test]
fn pml_fields_follow_the_stretched_coordinate() {
    // A PEC-backed PML closes the right side of a homogeneous box; for TE the
    // standing wave vanishes at the backing wall, so beyond the PML start
    // E_y(ζ) / E_y(0) = sin(κ(τ − ζ̃(ζ))) / sin(κτ).
    let thickness = 0.5;
    let profiles = [(2.0, C::new(0.0, 3.0)), (0.0, C::new(0.0, 1.0))];
    let zetas = [0.0, 0.1, 0.25, 0.4, 0.5];
    let mut in_pml = Vec::new();
    let mut outside = Vec::new();
    for (order, strength) in profiles {
        let stack = box_stack(
            wall(WallKind::PerfectElectric, -0.3),
            pml(0.9, thickness, order, strength, WallKind::PerfectElectric),
        );
        let tau = c(thickness) * (c(1.0) + strength / (order + 1.0));
        let beta = C::new(3.7, 0.2);
        let kappa = (c(EPS * K0 * K0) - beta * beta).sqrt();
        let sampling = FieldSampling::new()
            .layer(0, LayerSampling::uniform_with_boundaries(3, true, true))
            .right_exterior(ExteriorSampling::distances(
                zetas.map(|zeta| Length::centimetres(0.4 + zeta)),
            ));
        for backend in [
            &Transfer2::new() as &dyn ModalFieldBackend,
            &Scatter2::new(),
        ] {
            let fields = backend
                .unscaled_fields(
                    &analytical(&stack),
                    c(K0),
                    beta,
                    [None, None],
                    TE,
                    &sampling,
                    c(1.0),
                )
                .unwrap()
                .fields;
            let ey = fields.quantity().electric().y();
            let start = ey[3];
            let reference = ey[0];
            for (index, zeta) in zetas.iter().enumerate() {
                let fraction = zeta / thickness;
                let stretched =
                    c(*zeta) + strength * thickness / (order + 1.0) * fraction.powf(order + 1.0);
                let expected = (kappa * (tau - stretched)).sin() / (kappa * tau).sin();
                assert!(
                    (ey[3 + index] / start - expected).norm() < 1e-9,
                    "m={order} ζ={zeta}: {} vs {expected}",
                    ey[3 + index] / start
                );
            }
            outside.push(ey.iter().take(4).map(|v| v / reference).collect::<Vec<_>>());
            in_pml.push(ey[5] / reference);
        }
    }
    // Equal τ: the same fields outside the PML, different ones inside.
    for (a, b) in outside[0].iter().zip(&outside[2]) {
        assert!((a - b).norm() < 1e-10);
    }
    assert!((in_pml[0] - in_pml[2]).norm() > 1e-3);
}

#[test]
fn sampling_beyond_the_wall_is_refused() {
    let stack = box_stack(
        wall(WallKind::PerfectElectric, -0.3),
        pml(0.9, 0.5, 2.0, C::new(0.0, 2.0), WallKind::PerfectMagnetic),
    );
    let at_walls = FieldSampling::<f64>::new()
        .left_exterior(ExteriorSampling::point(Length::centimetres(0.3)))
        .right_exterior(ExteriorSampling::point(Length::centimetres(0.9)));
    at_walls.resolve(&stack).unwrap();
    for (sampling, side) in [
        (
            FieldSampling::new().left_exterior(ExteriorSampling::point(Length::centimetres(0.31))),
            StackSide::Left,
        ),
        (
            FieldSampling::new().right_exterior(ExteriorSampling::point(Length::millimetres(9.1))),
            StackSide::Right,
        ),
    ] {
        assert!(matches!(
            sampling.resolve(&stack),
            Err(FieldSamplingError::BeyondTermination { side: found, .. }) if found == side
        ));
        assert!(
            Transfer2::new()
                .unscaled_fields(
                    &analytical(&stack),
                    c(K0),
                    c(1.0),
                    [None, None],
                    TE,
                    &sampling,
                    c(1.0)
                )
                .is_err()
        );
    }
}

/// The characteristic of `stack` with walls fixed in z, as a function of the
/// thickness of layer 0.
fn with_thickness(stack: &Stack<Material, f64>, thickness: f64) -> Stack<Material, f64> {
    let mut builder = Stack::builder(*stack.left_exterior(), *stack.right_exterior());
    for (index, layer) in stack.layers_left_to_right().iter().enumerate() {
        let d = if index == 0 {
            Length::centimetres(thickness)
        } else {
            layer.thickness()
        };
        builder = builder.layer(*layer.material(), d);
    }
    builder
        .left_termination(*stack.left_termination())
        .right_termination(*stack.right_termination())
        .finalise()
}

#[test]
fn derivatives_through_closed_sides_match_central_differences() {
    let stacks = [
        lossy_stack(
            wall(WallKind::PerfectMagnetic, -0.2),
            wall(WallKind::PerfectElectric, 0.9),
        ),
        lossy_stack(
            Termination::Open,
            pml(0.8, 0.6, 2.0, C::new(0.5, 3.0), WallKind::PerfectMagnetic),
        ),
    ];
    let h = 1e-4;
    for stack in &stacks {
        for polarisation in POLARISATIONS {
            for backend in BACKENDS {
                let (k0, beta) = (c(K0), C::new(3.3, 0.2));
                let value = |stack: &Stack<Material, f64>, k0: C, beta: C| {
                    characteristic(stack, backend, k0, beta, polarisation).value()[()]
                };
                let base = value(stack, k0, beta);
                let check = |jet: J, plus: C, minus: C| {
                    let first = (plus - minus) / (2.0 * h);
                    let second = (plus - 2.0 * base + minus) / (h * h);
                    assert!(relative(jet.value()[()], base) < 1e-12);
                    assert!(
                        relative(jet.first()[()], first) < 1e-6,
                        "{} vs {first}",
                        jet.first()[()]
                    );
                    assert!(
                        relative(jet.second()[()], second) < 1e-4,
                        "{} vs {second}",
                        jet.second()[()]
                    );
                };
                let determinant =
                    |stack, k0, beta| determinant(stack, backend, k0, beta, polarisation).unwrap();
                let constant = |value: C| J::constant(arr0(value));
                // β
                check(
                    characteristic(stack, backend, k0, beta, polarisation),
                    value(stack, k0, beta + h),
                    value(stack, k0, beta - h),
                );
                // k0
                check(
                    determinant(compile(stack), J::variable(arr0(k0)), constant(beta)),
                    value(stack, k0 + h, beta),
                    value(stack, k0 - h, beta),
                );
                // Thickness of layer 0 with the walls fixed in z; the thickness
                // is the only seeded variable when k0 and β are constant.
                let canonical = compile_canonical_thickness_stack::<_, J>(
                    stack,
                    Ix0(),
                    &ValidationConfig::permissive(),
                    FiniteLayerIndex::new(0),
                )
                .unwrap();
                check(
                    determinant(canonical, constant(k0), constant(beta)),
                    value(&with_thickness(stack, 0.4 + h), k0, beta),
                    value(&with_thickness(stack, 0.4 - h), k0, beta),
                );
            }
        }
    }
}

#[test]
fn outward_thickness_derivative_at_a_wall_on_the_last_interface_is_unavailable() {
    let stack = Stack::builder(Constant::dielectric(1.5), Constant::dielectric(2.0))
        .layer(Constant::dielectric(4.0), Length::centimetres(0.4))
        .layer(Constant::dielectric(2.0), Length::centimetres(0.3))
        .right_termination(wall(WallKind::PerfectElectric, 0.4 + 0.3))
        .finalise();
    let error = compile_canonical_thickness_stack::<_, J>(
        &stack,
        Ix0(),
        &ValidationConfig::permissive(),
        FiniteLayerIndex::new(1),
    )
    .unwrap_err();
    assert_eq!(
        error,
        StackCompileError::OutwardThicknessDerivative {
            layer: FiniteLayerIndex::new(1),
            side: StackSide::Right,
        }
    );
    // The value itself compiles, with a zero gap.
    let canonical = compile(&stack);
    assert_eq!(
        canonical.right_termination().gap_cm().unwrap().value()[()],
        c(0.0)
    );
}

#[test]
fn canonical_terminations_carry_positions_gaps_and_optical_thickness() {
    let stack = Stack::builder(Constant::dielectric(1.5), Constant::dielectric(2.0))
        .layer(Constant::dielectric(4.0), Length::millimetres(4.0))
        .layer(Constant::dielectric(2.0), Length::millimetres(3.0))
        .left_termination(wall(WallKind::PerfectElectric, -0.2))
        .right_termination(pml(
            0.9,
            0.5,
            2.0,
            C::new(1.0, 3.0),
            WallKind::PerfectMagnetic,
        ))
        .finalise();
    let canonical = compile_canonical_thickness_stack::<_, J>(
        &stack,
        Ix0(),
        &ValidationConfig::permissive(),
        FiniteLayerIndex::new(0),
    )
    .unwrap();
    let CanonicalTermination::Wall {
        kind: WallKind::PerfectElectric,
        position_cm,
        gap_cm,
    } = canonical.left_termination()
    else {
        panic!("left PEC wall expected");
    };
    assert!(relative(position_cm.value()[()], c(-0.2)) < 1e-15);
    assert!(relative(gap_cm.value()[()], c(0.2)) < 1e-15);
    assert_eq!(gap_cm.first()[()], c(0.0));
    let CanonicalTermination::Pml(right) = canonical.right_termination() else {
        panic!("right PML expected");
    };
    assert_eq!(right.backing, WallKind::PerfectMagnetic);
    assert!(relative(right.start_cm.value()[()], c(0.9)) < 1e-15);
    assert!(relative(right.gap_cm.value()[()], c(0.2)) < 1e-14);
    // The layer is given in millimetres: d gap / d(thickness in mm) = −0.1 cm/mm.
    assert!(relative(right.gap_cm.first()[()], c(-0.1)) < 1e-15);
    let tau = c(0.5) * (c(1.0) + C::new(1.0, 3.0) / 3.0);
    assert!(relative(right.optical_thickness_cm.value()[()], tau) < 1e-15);
    assert_eq!(right.optical_thickness_cm.first()[()], c(0.0));
    // Finite layers are unchanged.
    assert_eq!(canonical.layer_count(), 2);
}

#[test]
fn closed_modes_refuse_quasinormal_normalisation() {
    let stack = Stack::from_analytical_materials(Constant::vacuum(), Constant::vacuum())
        .analytical_layer(Constant::dielectric(4.0), Length::centimetres(0.3))
        .finalise()
        .with_right_termination(wall(WallKind::PerfectElectric, 0.5));
    let evaluator = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Transfer2::new()).unwrap();
    let state = evaluator
        .retain(
            CanonicalCoordinates::new(J::constant(arr0(c(K0))), J::constant(arr0(c(4.0)))),
            ExteriorWavevectors::closed_right(J::constant(arr0((c(K0 * K0) - 16.0).sqrt()))),
            TE,
        )
        .unwrap();
    let error = state.mode().map(|_| ()).unwrap_err();
    assert_eq!(format!("{error:?}"), "ClosedStack");
}

#[test]
fn plane_waves_come_only_from_an_open_side() {
    use crate::{
        CoordinateInput, Coordinates, InPlaneCoordinate, IncidentSide, RealAxisEvaluator,
        SolveRequestError, SpectralCoordinate,
        units::{InverseLengthUnit, LengthUnit},
    };
    use ndarray::array;

    let coordinates = Coordinates::new(
        SpectralCoordinate::VacuumWavelength(LengthUnit::Nanometre),
        InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre),
    );
    let input = || {
        CoordinateInput::intrinsic(coordinates, array![600.0, 900.0], array![0.0, 4.0e4]).unwrap()
    };
    let sampling = FieldSampling::new().layer(0, LayerSampling::uniform(3));
    let pml = |start| {
        pml(
            start,
            1.0e-4,
            2.0,
            C::new(0.0, 3.0),
            WallKind::PerfectElectric,
        )
    };
    // Each closure on the right and on the left of the layers, and whether a
    // lossless stack closed by it reflects everything.
    for (right, left, total_reflection) in [
        (
            wall(WallKind::PerfectElectric, 1.0e-4),
            wall(WallKind::PerfectElectric, -0.1e-4),
            true,
        ),
        (
            wall(WallKind::PerfectMagnetic, 0.9e-4),
            wall(WallKind::PerfectMagnetic, -0.1e-4),
            true,
        ),
        (pml(0.9e-4), pml(-0.1e-4), false),
    ] {
        for (open, closed) in [
            (IncidentSide::Left, IncidentSide::Right),
            (IncidentSide::Right, IncidentSide::Left),
        ] {
            let stack = Stack::builder(Constant::dielectric(1.0), Constant::dielectric(1.0))
                .layer(Constant::dielectric(2.25), Length::nanometres(500.0))
                .layer(Constant::dielectric(1.7), Length::nanometres(250.0))
                .finalise();
            let stack = match closed {
                IncidentSide::Right => stack.with_right_termination(right),
                IncidentSide::Left => stack.with_left_termination(left),
            };
            for polarisation in POLARISATIONS {
                macro_rules! check {
                    ($backend:expr) => {
                        let evaluator = RealAxisEvaluator::new($backend);
                        let result = evaluator.evaluate(input(), &stack, polarisation).unwrap();
                        assert_eq!(
                            result.amplitudes(closed).map(|_| ()).unwrap_err(),
                            SolveRequestError::ClosedIncidence { side: closed }
                        );
                        assert!(result.power(closed).is_err());
                        let amplitudes = result.amplitudes(open).unwrap();
                        let power = result.power(open).unwrap();
                        let (reflection, transmission) = amplitudes.into_parts();
                        for (r, t) in reflection.iter().zip(transmission.iter()) {
                            assert_eq!(*t, c(0.0));
                            if total_reflection {
                                assert!((r.norm() - 1.0).abs() < 1e-12, "|r| = {}", r.norm());
                            } else {
                                assert!(r.norm() < 1.0);
                            }
                        }
                        let (_, transmittance, _) = power.into_parts();
                        assert!(transmittance.iter().all(|t| *t == 0.0));

                        let retained = evaluator.retain(input(), &stack, polarisation).unwrap();
                        let point = retained.project_point(&0).unwrap();
                        assert_eq!(
                            point.excitation(closed).map(|_| ()).unwrap_err(),
                            SolveRequestError::ClosedIncidence { side: closed }
                        );
                        let error = point
                            .excitation(open)
                            .unwrap()
                            .evaluate_fields(&sampling)
                            .map(|_| ())
                            .unwrap_err();
                        assert!(format!("{error:?}").contains("DrivenFieldsOfClosedStack"));
                    };
                }
                check!(Transfer2::new());
                check!(Scatter2::new());
            }
        }
    }
}
