//! Entire functions of a squared phase, including derivatives at cutoff.
//!
//! C(x)=cos(sqrt(x)), S(x)=sin(sqrt(x))/sqrt(x). Both are single-valued
//! entire functions, even where a pointwise square-root jet is singular.

use nalgebra::ComplexField;
use num_traits::{One, Zero};

fn integer<C: ComplexField + Copy>(mut n: usize) -> C {
    let mut result = C::zero();
    let mut bit = C::one();
    while n != 0 {
        if n & 1 != 0 {
            result += bit;
        }
        bit += bit;
        n >>= 1;
    }
    result
}

// Twenty terms suffice for f32/f64 at |x| <= 1, including two derivatives.
// Sum each derivative's own series to avoid cancellation and division by x.
fn series<C: ComplexField + Copy>(x: C, sinc: bool, derivative: usize) -> C {
    let offset = usize::from(sinc);
    let mut term = match (sinc, derivative) {
        (_, 0) => C::one(),
        (false, 1) => -C::one() / integer::<C>(2),
        (true, 1) => -C::one() / integer::<C>(6),
        (false, 2) => C::one() / integer::<C>(12),
        (true, 2) => C::one() / integer::<C>(60),
        _ => unreachable!(),
    };
    let mut sum = term;
    for n in 1..20 {
        let k = n + derivative;
        term = -term * x * integer::<C>(k)
            / integer::<C>(n)
            / integer::<C>(2 * k - 1 + offset)
            / integer::<C>(2 * k + offset);
        sum += term;
    }
    sum
}

/// Value and first two scalar derivatives of C (sinc=false) or S (sinc=true).
pub(crate) fn expansion<C: ComplexField + Copy>(x: C, sinc: bool) -> [C; 3] {
    if x.modulus() <= C::RealField::one() {
        return [series(x, sinc, 0), series(x, sinc, 1), series(x, sinc, 2)];
    }
    // The hyperbolic form also supports negative real inputs for real jets.
    let (c, s) = if x.real() < C::RealField::zero() {
        let root = (-x).sqrt();
        (root.cosh(), root.sinh() / root)
    } else {
        let root = x.sqrt();
        (root.cos(), root.sin() / root)
    };
    let two = integer::<C>(2);
    let ds = (c - s) / (two * x);
    if sinc {
        [s, ds, (-s / two - integer::<C>(3) * ds) / (two * x)]
    } else {
        [c, -s / two, -ds / two]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64 as C;

    #[test]
    fn exact_cutoff_coefficients() {
        assert_eq!(expansion(0.0_f64, false), [1.0, -0.5, 1.0 / 12.0]);
        assert_eq!(expansion(0.0_f64, true), [1.0, -1.0 / 6.0, 1.0 / 60.0]);
    }

    #[test]
    fn regular_forms_and_series_agree_across_switch() {
        for angle in [0.0_f64, 0.7, 2.5, std::f64::consts::PI] {
            for radius in [0.999999, 1.0, 1.000001, 2.0] {
                let x = C::from_polar(radius, angle);
                for sinc in [false, true] {
                    let actual = expansion(x, sinc);
                    for (order, value) in actual.iter().enumerate() {
                        assert!((*value - series(x, sinc, order)).norm() < 2e-14);
                    }
                }
            }
        }
        assert!((expansion(-4.0_f64, false)[0] - 2.0_f64.cosh()).abs() < 1e-14);
    }

    #[test]
    fn derivatives_agree_with_independent_finite_differences() {
        for x in [C::new(0.3, 0.7), C::new(-3.0, 0.1), C::new(4.0, -2.0)] {
            for sinc in [false, true] {
                let [v, d, dd] = expansion(x, sinc);
                let h = 1e-4;
                let plus = expansion(x + h, sinc)[0];
                let minus = expansion(x - h, sinc)[0];
                assert!(((plus - minus) / (2.0 * h) - d).norm() < 1e-8);
                assert!(((plus - 2.0 * v + minus) / (h * h) - dd).norm() < 1e-7);
            }
        }
    }

    #[test]
    fn every_jet_family_continues_mixed_arrays_through_zero() {
        use crate::algebra::{
            ArrayJet0, ArrayJet1, ArrayJet2, ArrayJetBivariate1, ArrayJetBivariate2,
            HolomorphicParameter as P, ScalarAlgebra,
        };
        use crate::differential::{BivariateGradient, BivariateHessian};
        use ndarray::{Array1, Ix1, array};
        let values = array![
            C::new(0.0, 0.0),
            C::new(1e-12, -1e-12),
            C::new(-2.0, 0.3),
            C::new(0.4, -0.2)
        ];
        let fill = |v: f64| Array1::from_elem(values.len(), C::new(v, 0.0));
        for sinc in [false, true] {
            let j0 = ArrayJet0::<C, Ix1, P>::new(values.clone()).squared_phase_function(sinc);
            let j1 = ArrayJet1::<C, Ix1, P>::from_parts(values.clone(), fill(2.0))
                .squared_phase_function(sinc);
            let j2 = ArrayJet2::<C, Ix1, P>::from_parts(values.clone(), fill(2.0), fill(3.0))
                .squared_phase_function(sinc);
            let b1 = ArrayJetBivariate1::<C, Ix1, P>::from_parts(
                values.clone(),
                BivariateGradient::new(fill(2.0), fill(4.0)),
            )
            .squared_phase_function(sinc);
            let b2 = ArrayJetBivariate2::<C, Ix1, P>::from_parts(
                values.clone(),
                BivariateGradient::new(fill(2.0), fill(4.0)),
                BivariateHessian::new(fill(3.0), fill(5.0), fill(7.0)),
            )
            .squared_phase_function(sinc);
            assert!(
                j0.all_finite()
                    && j1.all_finite()
                    && j2.all_finite()
                    && b1.all_finite()
                    && b2.all_finite()
            );
            for (i, x) in values.iter().enumerate() {
                let [v, d, dd] = expansion(*x, sinc);
                for actual in [
                    j0.value()[i],
                    j1.value()[i],
                    j2.value()[i],
                    b1.value()[i],
                    b2.value()[i],
                ] {
                    assert!((actual - v).norm() < 1e-14);
                }
                for actual in [j1.first()[i], j2.first()[i], b1.axis0()[i], b2.axis0()[i]] {
                    assert!((actual - 2.0 * d).norm() < 1e-14);
                }
                assert!((b1.axis1()[i] - 4.0 * d).norm() < 1e-14);
                assert!((b2.axis1()[i] - 4.0 * d).norm() < 1e-14);
                assert!((j2.second()[i] - (4.0 * dd + 3.0 * d)).norm() < 1e-14);
                assert!((b2.axis0_axis0()[i] - (4.0 * dd + 3.0 * d)).norm() < 1e-14);
                assert!((b2.axis0_axis1()[i] - (8.0 * dd + 5.0 * d)).norm() < 1e-14);
                assert!((b2.axis1_axis1()[i] - (16.0 * dd + 7.0 * d)).norm() < 1e-14);
            }
        }
    }
}
