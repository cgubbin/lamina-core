use crate::{
    CanonicalCoordinates, ComplexPlaneEvaluator, Constant, ExteriorWavevectors, Polarisation,
    ScalarAlgebra, Scatter2, Stack, Transfer2,
    test_support::{C, assertions::assert_complex_close, c, jet::HoloJ2},
    units::Length,
};
use ndarray::arr0;

fn evaluate(
    k0: C,
    kp: C,
    spectral: bool,
    pol: Polarisation,
    multilayer: bool,
    delta: f64,
    scatter: bool,
) -> HoloJ2 {
    let mut builder = Stack::builder(Constant::vacuum(), Constant::vacuum())
        .layer(Constant::new(-1.0 + delta, -1.0), Length::centimetres(0.3));
    if multilayer {
        builder = builder.layer(Constant::vacuum(), Length::centimetres(0.2));
    }
    let stack = builder.finalise();
    let (k0, kp) = if spectral {
        (HoloJ2::variable(arr0(k0)), HoloJ2::constant(arr0(kp)))
    } else {
        (HoloJ2::constant(arr0(k0)), HoloJ2::variable(arr0(kp)))
    };
    let k = k0.square().subtract(&kp.square()).sqrt();
    let exterior = ExteriorWavevectors::new(k.clone(), k);
    let coordinates = CanonicalCoordinates::new(k0, kp);
    if scatter {
        let evaluator =
            ComplexPlaneEvaluator::<HoloJ2, _, _>::compile(&stack, Scatter2::new()).unwrap();
        let result = evaluator
            .checked_determinant(&coordinates, &exterior, pol)
            .unwrap()
            .into_inner();
        let retained = evaluator.retain(coordinates, exterior, pol).unwrap();
        let rd = retained.determinant().into_inner();
        assert_complex_close(rd.value()[()], result.value()[()], 1e-12);
        assert_complex_close(rd.first()[()], result.first()[()], 1e-12);
        assert_complex_close(rd.second()[()], result.second()[()], 1e-12);
        result
    } else {
        ComplexPlaneEvaluator::<HoloJ2, _, _>::compile(&stack, Transfer2::new())
            .unwrap()
            .checked_determinant(&coordinates, &exterior, pol)
            .unwrap()
            .into_inner()
    }
}

#[test]
fn negative_layer_and_multilayer_continue_with_analytic_derivatives() {
    for multilayer in [false, true] {
        let d = if multilayer { 0.1 } else { 0.3 };
        for spectral in [false, true] {
            for pol in [
                Polarisation::TransverseElectric,
                Polarisation::TransverseMagnetic,
            ] {
                let k0 = c(1.0);
                let kp = c(0.3);
                let k = (k0 * k0 - kp * kp).sqrt();
                let dk = if spectral { k0 / k } else { -kp / k };
                let ddk = if spectral {
                    -kp * kp / (k * k * k)
                } else {
                    -k0 * k0 / (k * k * k)
                };
                let phase = (C::i() * k * d).exp();
                let value = -2.0 * C::i() * k * phase;
                let first_k = -2.0 * C::i() * phase * (c(1.0) + C::i() * k * d);
                let second_k = phase * (c(4.0 * d) + 2.0 * C::i() * k * d * d);
                let expected = [value, first_k * dk, second_k * dk * dk + first_k * ddk];
                for scatter in [false, true] {
                    let actual = evaluate(k0, kp, spectral, pol, multilayer, 0.0, scatter);
                    assert!(actual.all_finite());
                    for (actual, expected) in
                        [actual.value()[()], actual.first()[()], actual.second()[()]]
                            .into_iter()
                            .zip(expected)
                    {
                        assert_complex_close(actual, expected, 2e-12);
                    }
                }
            }
        }
    }
}

fn independent_multilayer(k0: C, kp: C, pol: Polarisation, delta: f64) -> C {
    let k = (k0 * k0 - kp * kp).sqrt();
    let q = (1.0 - delta) * k0 * k0 - kp * kp;
    let kn = q.sqrt();
    let f = if pol == Polarisation::TransverseElectric {
        -1.0
    } else {
        -1.0 + delta
    };
    let a = (kn * 0.3).cos();
    let b = -(kn * 0.3).sin() * f / kn;
    let c0 = (kn * 0.3).sin() * kn / f;
    let e = (k * 0.2).cos();
    let g = -(k * 0.2).sin() / k;
    let h = (k * 0.2).sin() * k;
    // Independent scalar state matching across both finite layers.
    let field = (a * e + b * h) + (a * g + b * e) * C::i() * k;
    let slope = (c0 * e + a * h) + (c0 * g + a * e) * C::i() * k;
    -C::i() * k * field - slope
}

#[test]
fn physical_intermediate_interface_pole_is_regular_in_complete_multilayer() {
    use crate::backend::scatter2::{Scatter2ProjectiveEntries, cascade_projection};
    let k = (1.0_f64 - 0.3 * 0.3).sqrt();
    let constant = |v: C| HoloJ2::constant(arr0(v));
    // The physical vacuum -> negative-index interface, followed by propagation,
    // has d=0 and n21!=0. Later interfaces remove the subnetwork pole.
    let interface = Scatter2ProjectiveEntries::from_parts(
        constant(c(0.0)),
        constant(c(2.0 * k)),
        constant(c(-2.0 * k)),
        constant(c(2.0 * k)),
        constant(c(-2.0 * k)),
        constant(c(0.0)),
    );
    let p = constant((C::i() * k * 0.3).exp());
    let propagation = Scatter2ProjectiveEntries::from_parts(
        constant(c(1.0)),
        constant(c(0.0)),
        p.clone(),
        p.clone(),
        constant(c(0.0)),
        p.square(),
    );
    let prefix = cascade_projection(&interface, &propagation);
    assert_eq!(prefix.denominator().value()[()], c(0.0));
    assert!(prefix.n21().value()[()].norm() > 1.0);
    for delta in [-1e-6, 0.0, 1e-6] {
        for pol in [
            Polarisation::TransverseElectric,
            Polarisation::TransverseMagnetic,
        ] {
            for spectral in [false, true] {
                let k0 = c(1.0);
                let kp = c(0.3);
                let step = 1e-4;
                let actual = evaluate(k0, kp, spectral, pol, true, delta, true);
                let reference = independent_multilayer(k0, kp, pol, delta);
                let (plus, minus) = if spectral {
                    (
                        independent_multilayer(k0 + step, kp, pol, delta),
                        independent_multilayer(k0 - step, kp, pol, delta),
                    )
                } else {
                    (
                        independent_multilayer(k0, kp + step, pol, delta),
                        independent_multilayer(k0, kp - step, pol, delta),
                    )
                };
                assert_complex_close(actual.value()[()], reference, 1e-12);
                assert_complex_close(actual.first()[()], (plus - minus) / (2.0 * step), 1e-7);
                assert_complex_close(
                    actual.second()[()],
                    (plus - 2.0 * reference + minus) / (step * step),
                    2e-7,
                );
            }
        }
    }
}
