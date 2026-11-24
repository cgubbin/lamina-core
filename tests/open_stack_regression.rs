//! Open stacks keep their numerical results bit for bit.
//!
//! Each line of the golden file pairs a label with the 64-bit FNV-1a hash of
//! the `Debug` rendering of a plane-wave response, derivative, characteristic
//! or set of mode fields of a stack open on both sides. `Debug` prints every
//! `f64` with its shortest round-trip representation, so an equal hash means
//! equal bits. Set `LAMINA_BLESS=1` to rewrite the file.
use lamina_core::{
    CanonicalCoordinates, Constant, CoordinateInput, Coordinates, DifferentiableMaterialStack,
    ExteriorSampling, ExteriorWavevectors, FieldSampling, FiniteLayerIndex, InPlaneCoordinate,
    IncidentSide, LayerSampling, ModalCharacteristicDirection, ModalDerivativeBackend,
    ModalFieldBackend, ModeJet1, Parameter, Polarisation, RealAxisEvaluator, ScalarAlgebra,
    Scatter2, SeedJet, SpectralCoordinate, Stack, Transfer2,
    units::{InverseLengthUnit, Length, LengthUnit},
};
use ndarray::{Ix0, arr0, array};
use num_complex::Complex64 as C;

const GOLDEN: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/open_stack_regression.txt"
);

const POLARISATIONS: [Polarisation; 2] = [
    Polarisation::TransverseElectric,
    Polarisation::TransverseMagnetic,
];

fn kappas(left: C, right: C) -> [Option<C>; 2] {
    [Some(left), Some(right)]
}

fn real_axis_lines(lines: &mut Vec<String>) {
    let stack: DifferentiableMaterialStack<C> =
        Stack::from_differentiable_materials(Constant::vacuum(), Constant::new(2.1, 1.0))
            .differentiable_layer(Constant::dielectric(2.25), Length::nanometres(500.0))
            .differentiable_layer(Constant::new(3.0, 1.2), Length::nanometres(120.0))
            .differentiable_layer(Constant::dielectric(1.7), Length::micrometres(0.3))
            .finalise();
    let coordinates = Coordinates::new(
        SpectralCoordinate::VacuumWavelength(LengthUnit::Nanometre),
        InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre),
    );
    let input = || {
        CoordinateInput::intrinsic(
            coordinates,
            array![600.0, 800.0, 1000.0],
            array![0.0, 5.0e3, 1.2e5],
        )
        .unwrap()
    };
    let sampling = FieldSampling::new()
        .left_exterior(ExteriorSampling::uniform(3, Length::nanometres(200.0)))
        .layer(0, LayerSampling::uniform_with_boundaries(4, true, true))
        .layer(1, LayerSampling::uniform(3))
        .layer(2, LayerSampling::boundaries())
        .right_exterior(ExteriorSampling::uniform(3, Length::nanometres(300.0)));

    macro_rules! with_backend {
        ($backend:expr, $name:literal) => {
            let evaluator = RealAxisEvaluator::new($backend);
            for polarisation in POLARISATIONS {
                let result = evaluator.evaluate(input(), &stack, polarisation).unwrap();
                let first = evaluator
                    .evaluate_first(
                        input(),
                        &stack,
                        polarisation,
                        Parameter::LayerThickness(FiniteLayerIndex::new(1)),
                    )
                    .unwrap();
                let retained = evaluator.retain(input(), &stack, polarisation).unwrap();
                for side in [IncidentSide::Left, IncidentSide::Right] {
                    lines.push(format!(
                        "{} {:?} {:?} amplitudes {:?}",
                        $name,
                        polarisation,
                        side,
                        result.amplitudes(side).unwrap()
                    ));
                    lines.push(format!(
                        "{} {:?} {:?} power {:?}",
                        $name,
                        polarisation,
                        side,
                        result.power(side).unwrap()
                    ));
                    lines.push(format!(
                        "{} {:?} {:?} thickness derivative {:?}",
                        $name,
                        polarisation,
                        side,
                        first.amplitudes(side).unwrap()
                    ));
                    for index in 0..3 {
                        let point = retained.project_point(&index).unwrap();
                        let fields = point
                            .excitation(side)
                            .unwrap()
                            .evaluate_fields(&sampling)
                            .unwrap();
                        lines.push(format!(
                            "{} {:?} {:?} {} fields {:?}",
                            $name,
                            polarisation,
                            side,
                            index,
                            fields.quantity()
                        ));
                    }
                }
            }
        };
    }

    with_backend!(Scatter2::new(), "scatter");
    with_backend!(Transfer2::new(), "transfer");
}

fn modal_lines(lines: &mut Vec<String>) {
    let stack = Stack::from_analytical_materials(Constant::vacuum(), Constant::new(2.0, 1.0))
        .analytical_layer(Constant::dielectric(4.0), Length::centimetres(0.3))
        .analytical_layer(Constant::new(2.5, 1.1), Length::centimetres(0.2))
        .finalise();
    let k0 = C::new(7.0, -0.02);
    let points = [
        C::new(9.0, 0.01),
        C::new(11.5, -0.2),
        C::new(3.0, 0.5),
        C::new(0.0, 4.0),
    ];
    let sampling = FieldSampling::new()
        .left_exterior(ExteriorSampling::uniform(3, Length::centimetres(0.4)))
        .layer(0, LayerSampling::uniform_with_boundaries(4, true, true))
        .layer(1, LayerSampling::uniform(3))
        .right_exterior(ExteriorSampling::uniform(3, Length::centimetres(0.5)));

    macro_rules! with_backend {
        ($backend:expr, $name:literal) => {
            let backend = $backend;
            for polarisation in POLARISATIONS {
                for beta in points {
                    let left = (k0 * k0 - beta * beta).sqrt();
                    let right = (2.0 * k0 * k0 - beta * beta).sqrt();
                    for (label, branches) in
                        [("principal", (left, right)), ("flipped", (-left, right))]
                    {
                        let fields = backend.unscaled_fields(
                            &stack,
                            k0,
                            beta,
                            kappas(branches.0, branches.1),
                            polarisation,
                            &sampling,
                            C::new(0.7, 0.2),
                        );
                        match fields {
                            Ok(fields) => lines.push(format!(
                                "{} {:?} {} {} fields {:?} {:?} {:?} {:?}",
                                $name,
                                polarisation,
                                beta,
                                label,
                                fields.fields.quantity(),
                                fields.seed_residual,
                                fields.characteristic,
                                fields.seed_scale
                            )),
                            Err(error) => lines.push(format!(
                                "{} {:?} {} {} fields error {}",
                                $name, polarisation, beta, label, error
                            )),
                        }

                        type J = ModeJet1<C, Ix0>;
                        let jet = |value: C, rate: C| {
                            <J as SeedJet>::constant(arr0(value)).add(
                                &<J as SeedJet>::variable(arr0(C::new(0.0, 0.0)), 0)
                                    .unwrap()
                                    .scale(rate),
                            )
                        };
                        let coordinates = CanonicalCoordinates::new(
                            jet(k0, C::new(0.0, 0.0)),
                            jet(beta, C::new(1.0, 0.0)),
                        );
                        let exterior = ExteriorWavevectors::new(
                            jet(branches.0, -beta / branches.0),
                            jet(branches.1, -beta / branches.1),
                        );
                        let direction = ModalCharacteristicDirection {
                            stack: &stack,
                            coordinates: &coordinates,
                            exterior: &exterior,
                            polarisation,
                            thickness_rates: &[0.3, -0.1],
                            constitutive_rates: &[(C::new(0.0, 0.0), C::new(0.0, 0.0)); 4],
                        };
                        lines.push(format!(
                            "{} {:?} {} {} characteristic {:?}",
                            $name,
                            polarisation,
                            beta,
                            label,
                            backend
                                .directional_characteristic(direction)
                                .map_err(|e| e.to_string())
                        ));
                        lines.push(format!(
                            "{} {:?} {} {} directional fields {:?}",
                            $name,
                            polarisation,
                            beta,
                            label,
                            backend
                                .directional_fields(direction, &sampling, C::new(1.0, 0.0))
                                .map(|fields| (
                                    fields.fields.quantity().clone(),
                                    fields.tangent.quantity().clone()
                                ))
                                .map_err(|e| e.to_string())
                        ));
                    }
                }
            }
        };
    }

    with_backend!(Scatter2::new(), "scatter");
    with_backend!(Transfer2::new(), "transfer");
}

fn fnv1a(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

#[test]
fn open_stack_results_are_unchanged_bit_for_bit() {
    let mut lines = Vec::new();
    real_axis_lines(&mut lines);
    modal_lines(&mut lines);
    let actual = lines
        .iter()
        .map(|line| {
            let label = line.split(" fields ").next().unwrap_or(line);
            let label = label.split(" {").next().unwrap_or(label);
            let label = label.split(" Ok(").next().unwrap_or(label);
            let label = label.split(" Err(").next().unwrap_or(label);
            format!("{:016x} {}", fnv1a(line), label)
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";

    if std::env::var_os("LAMINA_BLESS").is_some() {
        std::fs::create_dir_all(std::path::Path::new(GOLDEN).parent().unwrap()).unwrap();
        std::fs::write(GOLDEN, &actual).unwrap();
        return;
    }

    let expected = std::fs::read_to_string(GOLDEN).unwrap();
    for (index, (actual, expected)) in actual.lines().zip(expected.lines()).enumerate() {
        assert_eq!(actual, expected, "line {index} differs");
    }
    assert_eq!(actual.lines().count(), expected.lines().count());
}
