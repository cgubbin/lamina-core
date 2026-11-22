use lamina_core::{
    CoordinateInput, Coordinates, ExteriorSampling, FieldSampling, InPlaneCoordinate, IncidentSide,
    LayerSampling, MaterialStack, Polarisation, RealAxisEvaluator, Scatter2, SpectralCoordinate,
    Stack,
    material::Constant,
    units::{InverseLengthUnit, Length, LengthUnit},
};
use ndarray::array;
use num_complex::Complex64;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stack: MaterialStack<Complex64> =
        Stack::from_materials(Constant::vacuum(), Constant::vacuum())
            .material_layer(Constant::dielectric(2.25), Length::nanometres(500.0))
            .finalise();
    let coordinates = Coordinates::new(
        SpectralCoordinate::VacuumWavelength(LengthUnit::Nanometre),
        InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre),
    );
    let input = CoordinateInput::intrinsic(
        coordinates,
        array![600.0, 800.0, 1000.0],
        array![0.0, 0.0, 0.0],
    )?;
    let evaluator = RealAxisEvaluator::new(Scatter2::new());
    let retained = evaluator.retain(input, &stack, Polarisation::TransverseElectric)?;

    // Project the 800 nm excitation before reconstructing spatial fields.
    let point = retained.project_point(&1)?;
    let excitation = point.excitation(IncidentSide::Left)?;
    let sampling = FieldSampling::new()
        .left_exterior(ExteriorSampling::point(Length::zero()))
        .layer(0, LayerSampling::uniform_with_boundaries(5, true, true))
        .right_exterior(ExteriorSampling::point(Length::zero()));
    let response = excitation.evaluate_fields(&sampling)?;

    // Fields use the excitation's unit incident-flux normalization. Samples
    // retain region identity, including opposite sides of each interface.
    println!("Fields at 800 nm (unit incident flux)");
    for (position, ey) in response
        .sampling()
        .positions()
        .iter()
        .zip(response.quantity().electric().y())
    {
        println!("{position:?}: Ey = {ey}");
    }
    Ok(())
}
