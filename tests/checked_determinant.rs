use lamina_core::{
    AnalyticalMaterialStack, CanonicalCoordinates, ComplexPlaneEvaluator, Constant,
    DeterminantEvaluationError, ExteriorWavevectors, ModeJet1, Polarisation, ScalarAlgebra,
    Scatter2, SeedJet, Stack, Transfer2, units::Length,
};
use ndarray::{Ix0, arr0};
use num_complex::Complex64 as C;

type J = ModeJet1<C, Ix0>;

fn constant(value: C) -> J {
    J::lift_constant(arr0(value))
}

fn stack(thickness: Option<f64>) -> AnalyticalMaterialStack<C> {
    let mut builder = Stack::from_analytical_materials(Constant::vacuum(), Constant::vacuum());
    if let Some(d) = thickness {
        builder = builder.analytical_layer(Constant::dielectric(4.0), Length::centimetres(d));
    }
    builder.finalise()
}

fn coordinates() -> CanonicalCoordinates<J> {
    CanonicalCoordinates::new(
        J::variable(arr0(C::new(2.0, 0.1)), 0).unwrap(),
        constant(C::new(0.3, 0.0)),
    )
}

fn exterior() -> ExteriorWavevectors<J> {
    let coordinates = coordinates();
    let k = coordinates
        .vacuum_angular_wavenumber()
        .square()
        .subtract(&coordinates.parallel_angular_wavenumber().square())
        .sqrt();
    ExteriorWavevectors::new(k.clone(), k)
}

#[test]
fn checked_public_api_preserves_regular_values_and_derivatives() {
    let stack = stack(Some(0.3));
    let evaluator = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Scatter2::new()).unwrap();
    for pol in [
        Polarisation::TransverseElectric,
        Polarisation::TransverseMagnetic,
    ] {
        let raw = evaluator
            .determinant(&coordinates(), &exterior(), pol)
            .unwrap()
            .into_inner();
        let checked = evaluator
            .checked_determinant(&coordinates(), &exterior(), pol)
            .unwrap()
            .into_inner();
        assert_eq!(raw.value(), checked.value());
        assert_eq!(raw.first(), checked.first());
    }
}

#[test]
fn checked_rejects_each_nonfinite_input_including_derivatives() {
    let stack = stack(None);
    let evaluator = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Scatter2::new()).unwrap();
    for component in 0..4 {
        for derivative_only in [false, true] {
            let bad = if derivative_only {
                // The primal is finite; only the carried derivative fails.
                J::variable(arr0(C::new(0.0, 0.0)), 0)
                    .unwrap()
                    .sqrt()
                    .add(&constant(C::new(2.0, 0.0)))
            } else {
                constant(C::new(f64::NAN, 0.0))
            };
            let mut inputs = [
                constant(C::new(2.0, 0.0)),
                constant(C::new(0.3, 0.0)),
                constant(C::new(1.0, 0.0)),
                constant(C::new(1.0, 0.0)),
            ];
            inputs[component] = bad;
            let [k0, kp, left, right] = inputs;
            assert!(matches!(
                evaluator.checked_determinant(
                    &CanonicalCoordinates::new(k0, kp),
                    &ExteriorWavevectors::new(left, right),
                    Polarisation::TransverseElectric
                ),
                Err(DeterminantEvaluationError::NonFiniteInput)
            ));
        }
    }
}

#[test]
fn checked_accepts_a_finite_zero_characteristic() {
    let stack = stack(None);
    let evaluator = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Transfer2::new()).unwrap();
    // A homogeneous interface with opposite selected sheets has D=-i(YL+YR)=0.
    let k = constant(C::new(1.0, 0.0));
    let result = evaluator
        .checked_determinant(
            &coordinates(),
            &ExteriorWavevectors::new(k.clone(), k.negate()),
            Polarisation::TransverseElectric,
        )
        .unwrap();
    assert_eq!(result.value().value()[()], C::new(0.0, 0.0));
    assert!(result.value().all_finite());
}

#[test]
fn checked_detects_output_derivative_overflow_with_finite_primal() {
    let stack = stack(Some(0.3));
    let evaluator = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Transfer2::new()).unwrap();
    let coordinates = CanonicalCoordinates::new(
        J::variable(arr0(C::new(0.0, 0.0)), 0)
            .unwrap()
            .scale(C::new(1e308, 0.0))
            .add(&constant(C::new(2.0, 0.0))),
        constant(C::new(0.3, 0.0)),
    );
    let exterior = ExteriorWavevectors::new(constant(C::new(1.0, 0.0)), constant(C::new(1.0, 0.0)));
    let raw = evaluator
        .determinant(&coordinates, &exterior, Polarisation::TransverseElectric)
        .unwrap();
    assert!(raw.value().value()[()].is_finite());
    assert!(!raw.value().all_finite());
    assert!(matches!(
        evaluator.checked_determinant(&coordinates, &exterior, Polarisation::TransverseElectric),
        Err(DeterminantEvaluationError::NonFiniteOutput)
    ));
}

#[test]
fn checked_preserves_transfer_backend_error_source() {
    use std::error::Error;
    let stack = stack(Some(1000.0));
    let evaluator = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Transfer2::new()).unwrap();
    let coordinates =
        CanonicalCoordinates::new(constant(C::new(1.0, 0.0)), constant(C::new(3.0, 0.0)));
    let exterior = ExteriorWavevectors::new(constant(C::new(0.0, 2.0)), constant(C::new(0.0, 2.0)));
    let raw = evaluator
        .determinant(&coordinates, &exterior, Polarisation::TransverseElectric)
        .err()
        .unwrap();
    let error = evaluator
        .checked_determinant(&coordinates, &exterior, Polarisation::TransverseElectric)
        .err()
        .unwrap();
    assert!(error.source().is_some());
    match error {
        DeterminantEvaluationError::Backend { source } => assert_eq!(source, raw),
        other => panic!("expected preserved backend failure, got {other:?}"),
    }
}

#[test]
fn opposite_admittance_interface_has_finite_scattering_continuation() {
    // Both constitutive factors are finite and nonzero. Principal interior
    // kappa=+1 gives Y=-1, opposite the exterior Y=+1 for either TE or TM.
    let d = 0.3_f64;
    let stack: AnalyticalMaterialStack<C> =
        Stack::from_analytical_materials(Constant::vacuum(), Constant::vacuum())
            .analytical_layer(Constant::new(-1.0, -1.0), Length::centimetres(d))
            .finalise();
    let coordinates =
        CanonicalCoordinates::new(constant(C::new(1.0, 0.0)), constant(C::new(0.0, 0.0)));
    let exterior = ExteriorWavevectors::new(constant(C::new(1.0, 0.0)), constant(C::new(1.0, 0.0)));
    let transfer = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Transfer2::new()).unwrap();
    let scatter = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Scatter2::new()).unwrap();
    let expected = C::new(2.0 * d.sin(), -2.0 * d.cos());
    for pol in [
        Polarisation::TransverseElectric,
        Polarisation::TransverseMagnetic,
    ] {
        let regular = transfer
            .checked_determinant(&coordinates, &exterior, pol)
            .unwrap();
        assert!((regular.value().value()[()] - expected).norm() < 1e-13);
        let raw = scatter.determinant(&coordinates, &exterior, pol).unwrap();
        assert!(raw.value().all_finite());
        let checked = scatter
            .checked_determinant(&coordinates, &exterior, pol)
            .unwrap();
        assert!((checked.value().value()[()] - expected).norm() < 1e-13);
        assert_eq!(checked.value().value(), raw.value().value());
        let retained = scatter
            .retain(coordinates.clone(), exterior.clone(), pol)
            .unwrap();
        assert!((retained.determinant().value().value()[()] - expected).norm() < 1e-13);
    }
}

#[test]
fn thick_evanescent_characteristic_is_correct_or_explicitly_fails() {
    let coordinates =
        CanonicalCoordinates::new(constant(C::new(1.0, 0.0)), constant(C::new(3.0, 0.0)));
    let exterior = ExteriorWavevectors::new(
        constant(C::new(0.0, 8.0_f64.sqrt())),
        constant(C::new(0.0, 8.0_f64.sqrt())),
    );
    for d in [0.3_f64, 5.0, 1000.0] {
        let stack = stack(Some(d));
        let transfer = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Transfer2::new()).unwrap();
        let scatter = ComplexPlaneEvaluator::<J, _, _>::compile(&stack, Scatter2::new()).unwrap();
        for (pol, factor) in [
            (Polarisation::TransverseElectric, 1.0),
            (Polarisation::TransverseMagnetic, 4.0),
        ] {
            let t = transfer.checked_determinant(&coordinates, &exterior, pol);
            let s = scatter.checked_determinant(&coordinates, &exterior, pol);
            if d < 10.0 {
                // Independent hyperbolic slab equation with q=-5, xi=sqrt(8).
                let phase = 5.0_f64.sqrt() * d;
                let expected = 2.0 * 8.0_f64.sqrt() * phase.cosh()
                    + (8.0 * factor + 5.0 / factor) * phase.sinh() / 5.0_f64.sqrt();
                for value in [t.unwrap(), s.unwrap()] {
                    let relative =
                        (value.value().value()[()] - C::new(expected, 0.0)).norm() / expected.abs();
                    assert!(relative < 1e-12, "relative error {relative}");
                }
            } else {
                assert!(matches!(t, Err(DeterminantEvaluationError::Backend { .. })));
                assert!(matches!(
                    s,
                    Err(DeterminantEvaluationError::NonFiniteOutput)
                ));
            }
        }
    }
}
