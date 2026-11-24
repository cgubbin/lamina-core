//! Guided-mode fields decay monotonically deep into both exteriors.
//!
//! An accepted root carries a round-off residual, so a reconstruction that
//! keeps the incoming exterior wave grows as exp(+kappa d) and overtakes the
//! decaying outgoing wave some tens of decay lengths out. Mode fields keep only
//! the outgoing wave in each exterior.
use lamina_core::{
    Constant, ExteriorSampling, FieldSampling, LayerSampling, ModalFieldBackend, Polarisation,
    Scatter2, Stack, Transfer2,
};
use lamina_units::{Length, LengthUnit};
use num_complex::Complex64 as C;

const K0: f64 = std::f64::consts::TAU; // cm^-1, a 1 cm wavelength
const CORE: f64 = 4.0;
const D: f64 = 1.0; // cm

/// Fundamental TE mode of a symmetric vacuum / eps 4 / vacuum slab, from
/// tan(q d / 2) = kappa / q, bisected to machine precision.
fn fundamental_beta() -> f64 {
    let f = |b: f64| {
        let q = (CORE * K0 * K0 - b * b).sqrt();
        let kappa = (b * b - K0 * K0).sqrt();
        (q * D / 2.).tan() - kappa / q
    };
    // The fundamental lies where q d / 2 is below pi / 2.
    let q_max = std::f64::consts::PI / D;
    let (mut lo, mut hi) = (
        (CORE * K0 * K0 - q_max * q_max).sqrt() + 1e-12,
        CORE.sqrt() * K0,
    );
    loop {
        let mid = 0.5 * (lo + hi);
        if mid == lo || mid == hi {
            return mid;
        }
        if f(mid) > 0. { lo = mid } else { hi = mid }
    }
}

fn check(backend: &impl ModalFieldBackend) {
    let beta = fundamental_beta();
    let kappa = (beta * beta - K0 * K0).sqrt();
    let stack = Stack::from_analytical_materials(Constant::vacuum(), Constant::vacuum())
        .analytical_layer(Constant::dielectric(CORE), Length::centimetres(D))
        .finalise();
    // Out to 60 decay lengths: the field falls by e^-60.
    let distances: Vec<f64> = (0..=600).map(|i| i as f64 * 0.1 / kappa).collect();
    let sampling = FieldSampling::new()
        .left_exterior(ExteriorSampling::distances_in(
            distances.clone(),
            LengthUnit::Centimetre,
        ))
        .layer(0, LayerSampling::boundaries())
        .right_exterior(ExteriorSampling::distances_in(
            distances.clone(),
            LengthUnit::Centimetre,
        ));
    let result = backend
        .unscaled_fields(
            &stack,
            C::new(K0, 0.),
            C::new(beta, 0.),
            [Some(C::new(0., kappa)); 2],
            Polarisation::TransverseElectric,
            &sampling,
            C::new(1., 0.),
        )
        .unwrap();
    let ey = result.fields.quantity().electric().y();
    let n = distances.len();
    let (left, rest) = ey.as_slice().unwrap().split_at(n);
    let (core, right) = rest.split_at(2);
    for side in [left, right] {
        for (i, pair) in side.windows(2).enumerate() {
            assert!(
                pair[1].norm() < pair[0].norm(),
                "|Ey| rises at {} decay lengths",
                (i + 1) as f64 * 0.1
            );
        }
        // And at the analytic rate.
        let ratio = side[n - 1] / side[0];
        assert!((ratio - C::new((-60.0_f64).exp(), 0.)).norm() < 1e-10 * (-60.0_f64).exp());
    }
    // Tangential Ey is continuous at both interfaces, so the exterior field at
    // the interface is unchanged by dropping the round-off incoming wave.
    assert!((left[0] - core[0]).norm() < 1e-12 * core[0].norm());
    assert!((right[0] - core[1]).norm() < 1e-12 * core[1].norm());
}

#[test]
fn scatter2_guided_mode_decays_into_both_exteriors() {
    check(&Scatter2::new());
}

#[test]
fn transfer2_guided_mode_decays_into_both_exteriors() {
    check(&Transfer2::new());
}
