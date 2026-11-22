use lamina_core::{
    CoordinateInput, Coordinates, ElectromagneticFields, ExteriorSampling, FieldPosition,
    FieldReconstructionError, FieldSampling, FieldSamplingError, InPlaneCoordinate, IncidentSide,
    LayerSampling, MaterialStack, PointProjectionError, Polarisation, RealAxisEvaluator,
    ResolvedFieldSampling, Scatter2, SpatialResponse, SpectralCoordinate, Stack, VectorField,
    material::Constant,
    units::{InverseLengthUnit, Length, LengthUnit},
};
use ndarray::{Ix1, arr0, array};
use num_complex::Complex64;

#[test]
fn projected_fields_match_scalar_solve_and_obey_interface_continuity() {
    let stack: MaterialStack<Complex64> =
        Stack::from_materials(Constant::vacuum(), Constant::vacuum())
            .material_layer(Constant::dielectric(2.25), Length::nanometres(500.0))
            .finalise();
    let coordinates = Coordinates::new(
        SpectralCoordinate::VacuumWavelength(LengthUnit::Nanometre),
        InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre),
    );
    let evaluator = RealAxisEvaluator::new(Scatter2::new());
    let sampling = FieldSampling::new()
        .left_exterior(ExteriorSampling::point(Length::zero()))
        .layer(0, LayerSampling::uniform_with_boundaries(5, true, true))
        .right_exterior(ExteriorSampling::point(Length::zero()));

    for polarisation in [
        Polarisation::TransverseElectric,
        Polarisation::TransverseMagnetic,
    ] {
        let retained = evaluator
            .retain(
                CoordinateInput::intrinsic(
                    coordinates,
                    array![600.0, 800.0, 1000.0],
                    array![0.0, 0.0, 0.0],
                )
                .unwrap(),
                &stack,
                polarisation,
            )
            .unwrap();
        assert!(matches!(
            retained.project_point(&3),
            Err(PointProjectionError::OutOfBounds { .. })
        ));
        let point = retained.project_point(&1).unwrap();
        let direct = evaluator
            .retain(
                CoordinateInput::intrinsic(coordinates, arr0(800.0), arr0(0.0)).unwrap(),
                &stack,
                polarisation,
            )
            .unwrap();

        for side in [IncidentSide::Left, IncidentSide::Right] {
            let excitation = point.excitation(side).unwrap();
            let response: SpatialResponse<ElectromagneticFields<VectorField<Complex64, Ix1>>, f64> =
                excitation.evaluate_fields(&sampling).unwrap();
            let expected = direct
                .excitation(side)
                .unwrap()
                .evaluate_fields(&sampling)
                .unwrap();
            let resolved: &ResolvedFieldSampling<f64> = response.sampling();
            assert_eq!(resolved, expected.sampling());
            assert_eq!(resolved.len(), 7);
            assert!(matches!(
                resolved.positions()[0],
                FieldPosition::LeftExterior { .. }
            ));
            let fields = response.quantity();
            let reference = expected.quantity();
            for (actual, expected) in [
                (fields.electric().x(), reference.electric().x()),
                (fields.electric().y(), reference.electric().y()),
                (fields.electric().z(), reference.electric().z()),
                (fields.magnetic().x(), reference.magnetic().x()),
                (fields.magnetic().y(), reference.magnetic().y()),
                (fields.magnetic().z(), reference.magnetic().z()),
            ] {
                for (a, b) in actual.iter().zip(expected) {
                    assert!((*a - *b).norm() < 1e-11);
                }
            }
            // Independently check tangential continuity on opposite sides of
            // both interfaces: exterior/film and film/exterior.
            for component in [
                fields.electric().x(),
                fields.electric().y(),
                fields.magnetic().x(),
                fields.magnetic().y(),
            ] {
                for (left, right) in [(0, 1), (5, 6)] {
                    assert!((component[left] - component[right]).norm() < 1e-11);
                }
            }
            let invalid = FieldSampling::new().layer(1, LayerSampling::uniform(3));
            assert!(matches!(
                excitation.evaluate_fields(&invalid),
                Err(FieldReconstructionError::FieldSampling(
                    FieldSamplingError::LayerOutOfBounds { .. }
                ))
            ));
        }
    }
}
