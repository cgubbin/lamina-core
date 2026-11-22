//! Independent slab-characteristic checks at regular interior cutoffs.
use crate::{
    CanonicalCoordinates, ComplexPlaneEvaluator, Constant, ExteriorWavevectors, Polarisation,
    ScalarAlgebra, Scatter2, Stack, Transfer2,
    test_support::{C, assertions::assert_complex_close, c, jet::HoloJ2},
    units::Length,
};
use ndarray::arr0;

const D: f64 = 0.3;

// Independent power series, used without a square root on either side of q=0.
fn reference(k0: C, kp: C, pol: Polarisation) -> C {
    let q = 4.0 * k0 * k0 - kp * kp;
    let mut cosine = c(1.0);
    let mut sinc = c(1.0);
    let mut ct = c(1.0);
    let mut st = c(1.0);
    for n in 1..30 {
        ct *= -q * D * D / ((2 * n - 1) * (2 * n)) as f64;
        st *= -q * D * D / ((2 * n) * (2 * n + 1)) as f64;
        cosine += ct;
        sinc += st;
    }
    let (f, exterior_factor) = match pol {
        Polarisation::TransverseElectric => (1.0, 1.0),
        Polarisation::TransverseMagnetic => (4.0, 9.0),
    };
    let xi = -C::i() * (9.0 * k0 * k0 - kp * kp).sqrt() / exterior_factor;
    2.0 * xi * cosine + (f * xi * xi - q / f) * D * sinc
}

fn solve(k0: C, kp: C, spectral: bool, pol: Polarisation, scatter: bool) -> HoloJ2 {
    let stack = Stack::builder(Constant::dielectric(9.0), Constant::dielectric(9.0))
        .layer(Constant::dielectric(4.0), Length::centimetres(D))
        .finalise();
    let (k0, kp) = if spectral {
        (HoloJ2::variable(arr0(k0)), HoloJ2::constant(arr0(kp)))
    } else {
        (HoloJ2::constant(arr0(k0)), HoloJ2::variable(arr0(kp)))
    };
    let k = k0.square().scale(c(9.0)).subtract(&kp.square()).sqrt();
    let exterior = ExteriorWavevectors::new(k.clone(), k);
    let coordinates = CanonicalCoordinates::new(k0, kp);
    if scatter {
        ComplexPlaneEvaluator::<HoloJ2, _, _>::compile(&stack, Scatter2::new())
            .unwrap()
            .checked_determinant(&coordinates, &exterior, pol)
            .unwrap()
            .into_inner()
    } else {
        ComplexPlaneEvaluator::<HoloJ2, _, _>::compile(&stack, Transfer2::new())
            .unwrap()
            .checked_determinant(&coordinates, &exterior, pol)
            .unwrap()
            .into_inner()
    }
}

#[test]
fn transfer_slab_cutoff_and_nearby_derivatives_match_independent_characteristic() {
    for pol in [
        Polarisation::TransverseElectric,
        Polarisation::TransverseMagnetic,
    ] {
        for spectral in [false, true] {
            for delta in [0.0, -1e-7, 1e-7, 0.03] {
                let k0 = c(1.0 + delta);
                let kp = c(2.0);
                let actual = solve(k0, kp, spectral, pol, false);
                assert!(actual.all_finite());
                assert_complex_close(actual.value()[()], reference(k0, kp, pol), 2e-12);
                let h = 2e-4;
                let (plus, minus) = if spectral {
                    (reference(k0 + h, kp, pol), reference(k0 - h, kp, pol))
                } else {
                    (reference(k0, kp + h, pol), reference(k0, kp - h, pol))
                };
                assert_complex_close(actual.first()[()], (plus - minus) / (2.0 * h), 2e-6);
                assert_complex_close(
                    actual.second()[()],
                    (plus - 2.0 * reference(k0, kp, pol) + minus) / (h * h),
                    3e-6,
                );
            }
        }
    }
}

#[test]
fn transfer_continues_across_interior_square_root_cut() {
    for pol in [
        Polarisation::TransverseElectric,
        Polarisation::TransverseMagnetic,
    ] {
        for spectral in [false, true] {
            // Interior q lies just above/below the negative real axis; exterior
            // q stays near the positive real axis and is not the cut being tested.
            for im in [-1e-8, 0.0, 1e-8] {
                let k0 = C::new(1.0, im);
                let kp = c(2.1);
                let transfer = solve(k0, kp, spectral, pol, false);
                let scatter = solve(k0, kp, spectral, pol, true);
                assert_complex_close(transfer.value()[()], reference(k0, kp, pol), 1e-12);
                assert_complex_close(transfer.value()[()], scatter.value()[()], 1e-12);
                assert_complex_close(transfer.first()[()], scatter.first()[()], 2e-11);
                assert_complex_close(transfer.second()[()], scatter.second()[()], 2e-10);
                let h = 1e-5;
                let (plus, minus) = if spectral {
                    (reference(k0 + h, kp, pol), reference(k0 - h, kp, pol))
                } else {
                    (reference(k0, kp + h, pol), reference(k0, kp - h, pol))
                };
                assert_complex_close(transfer.first()[()], (plus - minus) / (2.0 * h), 1e-7);
            }
        }
    }
}
