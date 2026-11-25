//! A scalar material, or an isotropic region, whose ε or μ is not finite at the
//! requested wavenumber is refused by name through both 2×2 backends, on the
//! real axis and in the complex plane.

use lamina_units::Length;
use ndarray::arr0;

use crate::{
    AnalyticalMaterialHandle, AnalyticalMaterialRegion, CanonicalCoordinates,
    ComplexPlaneEvaluator, ConstitutiveError, DeterminantEvaluationError, ExteriorWavevectors,
    FiniteLayerIndex, MaterialRegion, Parameter, Polarisation, RealAxisEvaluationError,
    RealAxisEvaluator, Scatter2, Stack, StackRegion, Transfer2,
    backend::{Scatter2Error, Transfer2Error},
    material::MaterialHandle,
    test_support::{
        jet::{HoloJ0, HoloJ1, HoloJ2},
        materials::ComplexQuadratic,
        planar::scalar_real_input,
    },
};

type C = num_complex::Complex64;

const POLARISATIONS: [Polarisation; 2] = [
    Polarisation::TransverseElectric,
    Polarisation::TransverseMagnetic,
];

const LAYER: StackRegion = StackRegion::Layer(FiniteLayerIndex::new(1));

fn constant(epsilon: C, mu: C) -> ComplexQuadratic {
    ComplexQuadratic {
        epsilon: [epsilon, C::default(), C::default()],
        mu: [mu, C::default(), C::default()],
    }
}

fn vacuum() -> ComplexQuadratic {
    constant(C::from(1.0), C::from(1.0))
}

/// Non-finite ε, then non-finite μ: the value of a tabulated material
/// evaluated off its table.
fn non_finite() -> [ComplexQuadratic; 2] {
    [
        constant(C::new(f64::NAN, 0.0), C::from(1.0)),
        constant(C::from(2.25), C::new(0.0, f64::INFINITY)),
    ]
}

/// vacuum | 2.25 | `bad` | vacuum: the non-finite material is layer 1.
fn stack<M>(region: impl Fn(ComplexQuadratic) -> M, bad: ComplexQuadratic) -> Stack<M, f64> {
    Stack::builder(region(vacuum()), region(vacuum()))
        .layer(
            region(constant(C::from(2.25), C::from(1.0))),
            Length::centimetres(0.2),
        )
        .layer(region(bad), Length::centimetres(0.3))
        .finalise()
}

trait BackendConstitutive {
    fn constitutive(&self) -> Option<ConstitutiveError>;
}

impl BackendConstitutive for Scatter2Error {
    fn constitutive(&self) -> Option<ConstitutiveError> {
        match self {
            Self::Constitutive(error) => Some(*error),
            _ => None,
        }
    }
}

impl BackendConstitutive for Transfer2Error {
    fn constitutive(&self) -> Option<ConstitutiveError> {
        match self {
            Self::Constitutive(error) => Some(*error),
            _ => None,
        }
    }
}

fn real_axis<T, X, B: BackendConstitutive>(
    result: Result<T, RealAxisEvaluationError<X, B>>,
    label: &str,
) -> ConstitutiveError {
    match result {
        Err(RealAxisEvaluationError::Backend { source }) => source
            .constitutive()
            .unwrap_or_else(|| panic!("{label}: backend error is not constitutive")),
        Err(RealAxisEvaluationError::Constitutive(error)) => error,
        Err(_) => panic!("{label}: unexpected error"),
        Ok(_) => panic!("{label}: a non-finite material was accepted"),
    }
}

fn complex_plane<T, B: BackendConstitutive>(
    result: Result<T, DeterminantEvaluationError<B>>,
    label: &str,
) -> ConstitutiveError {
    match result {
        Err(DeterminantEvaluationError::Backend { source }) => source
            .constitutive()
            .unwrap_or_else(|| panic!("{label}: backend error is not constitutive")),
        Err(_) => panic!("{label}: unexpected error"),
        Ok(_) => panic!("{label}: a non-finite material was accepted"),
    }
}

/// `$body` once with `$evaluator` a real-axis evaluator on each backend.
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

/// Each listed real-axis evaluation, `evaluate` for values and
/// `evaluate_first` and `evaluate_second` for spectral jets.
macro_rules! check_real_axis {
    ($stack:expr) => {
        check_real_axis!($stack; evaluate_first, evaluate_second)
    };
    ($stack:expr; $($derivative:ident),*) => {{
        let stack = $stack;
        let expected = ConstitutiveError::NonFinite { region: LAYER };
        let input = || scalar_real_input(2.0, 0.5);
        for polarisation in POLARISATIONS {
            for_each_backend!(evaluator, {
                let label = format!("{polarisation:?}");
                assert_eq!(
                    real_axis(evaluator.evaluate(input(), &stack, polarisation), &label),
                    expected,
                    "{label} evaluate"
                );
                $(
                    assert_eq!(
                        real_axis(
                            evaluator.$derivative(input(), &stack, polarisation, Parameter::Spectral),
                            &label
                        ),
                        expected,
                        "{label} {}",
                        stringify!($derivative)
                    );
                )*
            });
        }
    }};
}

#[test]
fn non_finite_scalar_and_isotropic_layers_are_refused_on_the_real_axis() {
    for bad in non_finite() {
        check_real_axis!(stack(AnalyticalMaterialHandle::<C>::new, bad));
        check_real_axis!(stack(MaterialHandle::<C>::new, bad););
        check_real_axis!(stack(
            |m| AnalyticalMaterialRegion::<C>::from(AnalyticalMaterialHandle::new(m)),
            bad
        ));
        check_real_axis!(stack(|m| MaterialRegion::<C>::from(MaterialHandle::new(m)), bad););
    }
}

#[test]
fn non_finite_scalar_and_isotropic_exteriors_are_refused_on_the_real_axis() {
    for bad in non_finite() {
        let scalar = Stack::builder(
            AnalyticalMaterialHandle::<C>::new(vacuum()),
            AnalyticalMaterialHandle::new(bad),
        )
        .finalise();
        let region = |m| AnalyticalMaterialRegion::<C>::from(AnalyticalMaterialHandle::new(m));
        let tensor = Stack::builder(region(vacuum()), region(bad)).finalise();
        let expected = ConstitutiveError::NonFinite {
            region: StackRegion::RightExterior,
        };
        for polarisation in POLARISATIONS {
            for_each_backend!(evaluator, {
                let input = || scalar_real_input(2.0, 0.5);
                let label = format!("{polarisation:?}");
                assert_eq!(
                    real_axis(evaluator.evaluate(input(), &scalar, polarisation), &label),
                    expected
                );
                assert_eq!(
                    real_axis(evaluator.evaluate(input(), &tensor, polarisation), &label),
                    expected
                );
            });
        }
    }
}

/// The checked characteristic at jet orders 0, 1 and 2 on both backends, with
/// vacuum exteriors' outgoing wavevectors.
macro_rules! check_complex_plane {
    ($stack:expr) => {{
        let stack = $stack;
        let expected = ConstitutiveError::NonFinite { region: LAYER };
        let (k0, kp) = (C::new(2.0, -0.1), C::new(0.5, 0.05));
        macro_rules! at_order {
            ($jet:ty, $seed:ident) => {{
                let coordinates =
                    CanonicalCoordinates::new(<$jet>::$seed(arr0(k0)), <$jet>::constant(arr0(kp)));
                let k0 = coordinates.vacuum_angular_wavenumber();
                let kp = coordinates.parallel_angular_wavenumber();
                let kappa = k0.multiply(k0).subtract(&kp.multiply(kp)).sqrt();
                let exterior = ExteriorWavevectors::new(kappa.clone(), kappa);
                for polarisation in POLARISATIONS {
                    let label = format!("{polarisation:?} {}", stringify!($jet));
                    let scatter =
                        ComplexPlaneEvaluator::<$jet, _, _>::compile(&stack, Scatter2::new())
                            .unwrap();
                    assert_eq!(
                        complex_plane(
                            scatter.checked_determinant(&coordinates, &exterior, polarisation),
                            &label
                        ),
                        expected,
                        "{label} Scatter2"
                    );
                    let transfer =
                        ComplexPlaneEvaluator::<$jet, _, _>::compile(&stack, Transfer2::new())
                            .unwrap();
                    assert_eq!(
                        complex_plane(
                            transfer.checked_determinant(&coordinates, &exterior, polarisation),
                            &label
                        ),
                        expected,
                        "{label} Transfer2"
                    );
                }
            }};
        }
        at_order!(HoloJ0, constant);
        at_order!(HoloJ1, variable);
        at_order!(HoloJ2, variable);
    }};
}

#[test]
fn non_finite_scalar_and_isotropic_layers_are_refused_in_the_complex_plane() {
    for bad in non_finite() {
        check_complex_plane!(stack(AnalyticalMaterialHandle::<C>::new, bad));
        check_complex_plane!(stack(
            |m| AnalyticalMaterialRegion::<C>::from(AnalyticalMaterialHandle::new(m)),
            bad
        ));
    }
}
