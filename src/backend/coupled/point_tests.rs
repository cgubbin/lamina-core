//! The coupled algebra on point jets against the array jets it replaces for
//! one-sample evaluations: every value and derivative agrees bit for bit.

use ndarray::{Ix0, arr0};
use num_complex::Complex64 as C;

use super::{
    berreman::{
        Region,
        tests::{random, random_tensor},
    },
    exponential::exponential,
    matrix::Matrix,
    split::{Direction, Split},
};
use crate::{
    algebra::{
        ArrayJet0, ArrayJet1, ArrayJet2, HolomorphicParameter, InlineJet, RealParameter,
        ScalarAlgebra, SeedJet,
    },
    backend::Channel,
    material::ConstitutiveTensor,
};

/// Every component of a one-sample jet, as bit patterns.
trait Bits {
    fn bits(&self) -> Vec<u64>;
}

fn complex_bits(value: C) -> [u64; 2] {
    [value.re.to_bits(), value.im.to_bits()]
}

impl<P> Bits for ArrayJet0<C, Ix0, P> {
    fn bits(&self) -> Vec<u64> {
        complex_bits(self.value()[()]).to_vec()
    }
}

impl<P> Bits for ArrayJet1<C, Ix0, P> {
    fn bits(&self) -> Vec<u64> {
        [self.value()[()], self.first()[()]]
            .into_iter()
            .flat_map(complex_bits)
            .collect()
    }
}

impl<P> Bits for ArrayJet2<C, Ix0, P> {
    fn bits(&self) -> Vec<u64> {
        [self.value()[()], self.first()[()], self.second()[()]]
            .into_iter()
            .flat_map(complex_bits)
            .collect()
    }
}

fn assert_same<J>(array: &Matrix<J>, point: &Matrix<J::Inline>, label: &str)
where
    J: InlineJet<Scalar = C> + Bits,
{
    let point = point.clone().into_map(J::from_inline);
    assert_eq!(
        (array.rows(), array.columns()),
        (point.rows(), point.columns())
    );
    for (index, (a, p)) in array.entries().iter().zip(point.entries()).enumerate() {
        assert_eq!(a.bits(), p.bits(), "{label}: entry {index}");
    }
}

/// A jet whose value is `value` and whose derivatives are those of `value`
/// times a seeded variable at `1 + 0.3i` (a constant `value` for jets
/// without derivatives).
fn jet<J: SeedJet<Scalar = C, Dimension = Ix0> + ScalarAlgebra>(value: C) -> J {
    match J::variable(arr0(C::new(1.0, 0.3)), 0) {
        Ok(variable) => ScalarAlgebra::multiply(&variable, &J::lift_constant(arr0(value))),
        Err(_) => J::lift_constant(arr0(value)),
    }
}

/// Generators, splits, exponentials, solves and determinants of random
/// tensors, on the array jets `J` and on their point form.
fn compare<J>()
where
    J: InlineJet<Scalar = C, Dimension = Ix0> + SeedJet + Bits,
{
    let mut seed = 17;
    let mut splits = 0;
    for _ in 0..4 {
        let mut tensor = |conjugate| {
            ConstitutiveTensor::Full(
                random_tensor(&mut seed, conjugate).map(|row| row.map(jet::<J>)),
            )
        };
        let tensors = [tensor(false), tensor(true)];
        let (k0, xi) = (jet::<J>(C::new(1.3, 0.0)), jet::<J>(C::new(0.4, 0.05)));
        let array = Region::new(&tensors, &k0, &xi, Channel::Coupled);
        let point = Region::new(
            &tensors.each_ref().map(|tensor| tensor.map(J::to_inline)),
            &k0.to_inline(),
            &xi.to_inline(),
            Channel::Coupled,
        );
        assert_same(&array.generator, &point.generator, "generator");

        let thickness = C::new(0.0, random(&mut seed).re + 1.0);
        assert_same(
            &exponential(&array.generator.scale(thickness)),
            &exponential(&point.generator.scale(thickness)),
            "exponential",
        );

        let rhs = Matrix::unit_columns(4, &[2, 1], array.generator.get(0, 0));
        let rhs_point = Matrix::unit_columns(4, &[2, 1], point.generator.get(0, 0));
        assert_same(
            &array.generator.solve(&rhs),
            &point.generator.solve(&rhs_point),
            "solve",
        );
        assert_same(
            &Matrix::from_fn(1, 1, |_, _| array.generator.determinant()),
            &Matrix::from_fn(1, 1, |_, _| point.generator.determinant()),
            "determinant",
        );

        let split = Split::select(&array.generator, &array.generator, Direction::Forward);
        let split_point = Split::select(&point.generator, &point.generator, Direction::Forward);
        match (split, split_point) {
            (Some(array), Some(point)) => {
                assert_same(&array.projector, &point.projector, "projector");
                splits += 1;
            }
            (None, None) => {}
            _ => panic!("the split exists on one form only"),
        }
    }
    assert!(splits > 0, "no generator split");
}

#[test]
fn first_order_point_jets_match_the_array_jets() {
    compare::<ArrayJet1<C, Ix0, HolomorphicParameter>>();
}

#[test]
fn second_order_point_jets_match_the_array_jets() {
    compare::<ArrayJet2<C, Ix0, RealParameter>>();
}

#[test]
fn value_point_jets_match_the_array_jets() {
    compare::<ArrayJet0<C, Ix0, RealParameter>>();
}
