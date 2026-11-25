use crate::{
    Parameter, Polarisation, RealAxisEvaluator,
    backend::{ExteriorContextProvider, Scatter2, Transfer2},
    test_support::{
        assertions::assert_first_jet_close,
        planar::{scalar_real_input, two_layer_stack},
    },
};

#[test]
fn retained_backends_agree_on_exterior_constitutive_quantities() {
    let stack = two_layer_stack();

    let scatter = RealAxisEvaluator::new(Scatter2::new())
        .retain_first(
            scalar_real_input(2.5, 0.31),
            &stack,
            Polarisation::TransverseElectric,
            Parameter::Spectral,
        )
        .unwrap();

    let transfer = RealAxisEvaluator::new(Transfer2::new())
        .retain_first(
            scalar_real_input(2.5, 0.31),
            &stack,
            Polarisation::TransverseElectric,
            Parameter::Spectral,
        )
        .unwrap();

    let scatter = scatter.project_point(&()).unwrap();
    let transfer = transfer.project_point(&()).unwrap();

    let scatter = scatter.workspace().solution().context();
    let transfer = transfer.workspace().solution().context();

    assert_first_jet_close(scatter.left().epsilon(), transfer.left().epsilon());
    assert_first_jet_close(scatter.left().mu(), transfer.left().mu());
    assert_first_jet_close(scatter.right().epsilon(), transfer.right().epsilon());
    assert_first_jet_close(scatter.right().mu(), transfer.right().mu());
}
