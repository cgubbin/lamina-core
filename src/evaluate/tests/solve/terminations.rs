//! Stacks closed by walls and PMLs: characteristic, zeros, fields and
//! derivatives.
use ndarray::{Ix0, arr0};
use num_complex::Complex64 as C;

use crate::{
    CanonicalCoordinates, ComplexPlane, ComplexPlaneEvaluator, Constant, ExteriorSampling,
    ExteriorWavevectors, FieldSampling, FieldSamplingError, FiniteLayerIndex, LayerSampling,
    ModalFieldBackend, Polarisation, Scatter2, Stack, StackCompileError, StackSide, Transfer2,
    ValidationConfig,
    backend::{IsotropicMediumQuantities, TerminationError},
    compile_canonical_constant_stack, compile_canonical_thickness_stack,
    input::canonical::CanonicalTermination,
    stack::{Termination, Wall, WallKind},
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

/// Principal exterior wavevectors for the open sides of `stack`.
fn exterior(
    stack: &crate::input::CanonicalStack<Material, J>,
    coordinates: &CanonicalCoordinates<J>,
) -> ExteriorWavevectors<J> {
    let kappa = |material| {
        IsotropicMediumQuantities::evaluate::<ComplexPlane, _>(material, coordinates).into_kappa()
    };
    ExteriorWavevectors::from_options(
        stack
            .left_termination()
            .is_open()
            .then(|| kappa(stack.left_exterior())),
        stack
            .right_termination()
            .is_open()
            .then(|| kappa(stack.right_exterior())),
    )
}

/// The characteristic on a canonical stack, as a jet in β (`beta_variable`)
/// or in k0.
fn characteristic_canonical(
    stack: crate::input::CanonicalStack<Material, J>,
    backend: Backend,
    k0: C,
    beta: C,
    polarisation: Polarisation,
    beta_variable: bool,
) -> Result<J, String> {
    let (k0, beta) = if beta_variable {
        (J::constant(arr0(k0)), J::variable(arr0(beta)))
    } else {
        (J::variable(arr0(k0)), J::constant(arr0(beta)))
    };
    let coordinates = CanonicalCoordinates::new(k0, beta);
    let exterior = exterior(&stack, &coordinates);
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

fn characteristic(
    stack: &Stack<Material, f64>,
    backend: Backend,
    k0: C,
    beta: C,
    polarisation: Polarisation,
) -> J {
    characteristic_canonical(compile(stack), backend, k0, beta, polarisation, true).unwrap()
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
    let stack = Stack::builder(Constant::dielectric(1.5), Constant::dielectric(2.0))
        .layer(Constant::dielectric(4.0), Length::centimetres(0.4))
        .layer(Constant::new(2.0, 1.3), Length::centimetres(0.3))
        .left_termination(wall(WallKind::PerfectMagnetic, -0.2))
        .right_termination(wall(WallKind::PerfectElectric, 1.4))
        .finalise();
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
            wall(WallKind::PerfectMagnetic, -0.75),
            wall(WallKind::PerfectElectric, 0.7),
        ),
        (
            wall(WallKind::PerfectElectric, -0.7),
            wall(WallKind::PerfectMagnetic, 1.3),
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

#[test]
fn exterior_wavevector_pattern_must_match_the_terminations() {
    let walled = box_stack(wall(WallKind::PerfectElectric, -0.3), Termination::Open);
    let right_walled = box_stack(Termination::Open, wall(WallKind::PerfectMagnetic, 0.9));
    let coordinates =
        CanonicalCoordinates::new(J::constant(arr0(c(K0))), J::constant(arr0(c(1.0))));
    let kappa = || J::constant(arr0(c(4.0)));
    // The correct patterns evaluate.
    for (stack, exterior) in [
        (&walled, ExteriorWavevectors::closed_left(kappa())),
        (&right_walled, ExteriorWavevectors::closed_right(kappa())),
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
        // A wavevector for the other wall side.
        (
            &right_walled,
            ExteriorWavevectors::new(kappa(), kappa()),
            StackSide::Right,
            false,
        ),
    ];
    for (stack, exterior, side, open) in cases {
        let expected = TerminationError::ExteriorBranchMismatch {
            side,
            open,
            supplied: !open,
        };
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
fn sampling_beyond_the_wall_is_refused() {
    let stack = box_stack(
        wall(WallKind::PerfectElectric, -0.3),
        wall(WallKind::PerfectMagnetic, 1.4),
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
        Stack::builder(Constant::dielectric(1.5), Constant::dielectric(2.0))
            .layer(Constant::dielectric(4.0), Length::centimetres(0.4))
            .layer(Constant::new(2.0, 1.3), Length::centimetres(0.3))
            .left_termination(wall(WallKind::PerfectMagnetic, -0.2))
            .right_termination(wall(WallKind::PerfectElectric, 0.9))
            .finalise(),
        Stack::builder(Constant::dielectric(1.5), Constant::dielectric(2.0))
            .layer(Constant::dielectric(4.0), Length::centimetres(0.4))
            .layer(Constant::new(2.0, 1.3), Length::centimetres(0.3))
            .left_termination(Termination::Open)
            .right_termination(wall(WallKind::PerfectMagnetic, 1.4))
            .finalise(),
    ];
    let value = |stack: &Stack<Material, f64>, backend, k0: C, beta: C, polarisation| {
        characteristic_canonical(compile(stack), backend, k0, beta, polarisation, true)
            .unwrap()
            .value()[()]
    };
    let h = 1e-4;
    for stack in &stacks {
        for polarisation in POLARISATIONS {
            for backend in BACKENDS {
                let (k0, beta) = (c(K0), C::new(3.3, 0.2));
                let base = value(stack, backend, k0, beta, polarisation);
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
                // β
                check(
                    characteristic_canonical(compile(stack), backend, k0, beta, polarisation, true)
                        .unwrap(),
                    value(stack, backend, k0, beta + h, polarisation),
                    value(stack, backend, k0, beta - h, polarisation),
                );
                // k0
                check(
                    characteristic_canonical(
                        compile(stack),
                        backend,
                        k0,
                        beta,
                        polarisation,
                        false,
                    )
                    .unwrap(),
                    value(stack, backend, k0 + h, beta, polarisation),
                    value(stack, backend, k0 - h, beta, polarisation),
                );
                // Thickness of layer 0 with the walls fixed in z; the thickness
                // is the only seeded variable when k0 and β are constant.
                let coordinates =
                    CanonicalCoordinates::new(J::constant(arr0(k0)), J::constant(arr0(beta)));
                let canonical = compile_canonical_thickness_stack::<_, J>(
                    stack,
                    Ix0(),
                    &ValidationConfig::permissive(),
                    FiniteLayerIndex::new(0),
                )
                .unwrap();
                let exterior = exterior(&canonical, &coordinates);
                let thickness_jet = match backend {
                    Backend::Transfer => {
                        ComplexPlaneEvaluator::from_canonical_stack(canonical, Transfer2::new())
                            .checked_determinant(&coordinates, &exterior, polarisation)
                            .unwrap()
                            .into_inner()
                    }
                    Backend::Scatter => {
                        ComplexPlaneEvaluator::from_canonical_stack(canonical, Scatter2::new())
                            .checked_determinant(&coordinates, &exterior, polarisation)
                            .unwrap()
                            .into_inner()
                    }
                };
                check(
                    thickness_jet,
                    value(
                        &with_thickness(stack, 0.4 + h),
                        backend,
                        k0,
                        beta,
                        polarisation,
                    ),
                    value(
                        &with_thickness(stack, 0.4 - h),
                        backend,
                        k0,
                        beta,
                        polarisation,
                    ),
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
fn canonical_terminations_carry_positions_and_gaps() {
    let stack = Stack::builder(Constant::dielectric(1.5), Constant::dielectric(2.0))
        .layer(Constant::dielectric(4.0), Length::millimetres(4.0))
        .layer(Constant::dielectric(2.0), Length::millimetres(3.0))
        .left_termination(wall(WallKind::PerfectElectric, -0.2))
        .right_termination(wall(WallKind::PerfectMagnetic, 0.9))
        .finalise();
    let canonical = compile_canonical_thickness_stack::<_, J>(
        &stack,
        Ix0(),
        &ValidationConfig::permissive(),
        FiniteLayerIndex::new(0),
    )
    .unwrap();
    let CanonicalTermination::Wall(left) = canonical.left_termination() else {
        panic!("left wall expected");
    };
    assert_eq!(left.kind(), WallKind::PerfectElectric);
    assert!(relative(left.position_cm().value()[()], c(-0.2)) < 1e-15);
    assert!(relative(left.gap_cm().value()[()], c(0.2)) < 1e-15);
    assert_eq!(left.gap_cm().first()[()], c(0.0));
    let CanonicalTermination::Wall(right) = canonical.right_termination() else {
        panic!("right wall expected");
    };
    assert_eq!(right.kind(), WallKind::PerfectMagnetic);
    assert!(relative(right.position_cm().value()[()], c(0.9)) < 1e-15);
    assert!(relative(right.gap_cm().value()[()], c(0.2)) < 1e-14);
    // The layer is given in millimetres: d gap / d(thickness in mm) = −0.1 cm/mm.
    assert!(relative(right.gap_cm().first()[()], c(-0.1)) < 1e-15);
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
    let closures = [
        wall(WallKind::PerfectElectric, 1.0e-4),
        wall(WallKind::PerfectMagnetic, 0.9e-4),
    ];
    for closure in closures {
        for (open, closed) in [
            (IncidentSide::Left, IncidentSide::Right),
            (IncidentSide::Right, IncidentSide::Left),
        ] {
            let mut stack = Stack::builder(Constant::dielectric(1.0), Constant::dielectric(1.0))
                .layer(Constant::dielectric(2.25), Length::nanometres(500.0))
                .layer(Constant::dielectric(1.7), Length::nanometres(250.0))
                .finalise();
            stack = match closed {
                IncidentSide::Right => stack.with_right_termination(closure),
                IncidentSide::Left => stack.with_left_termination(match closure {
                    Termination::Wall(wall) => {
                        Wall::new(wall.kind(), Length::centimetres(-0.1e-4)).into()
                    }
                    other => other,
                }),
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
                            // A lossless stack closed by a wall reflects everything.
                            assert!((r.norm() - 1.0).abs() < 1e-12, "|r| = {}", r.norm());
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
