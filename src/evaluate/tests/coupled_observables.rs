//! Layer observables and closed-form integrals of the 4×4 backends: against
//! the 2×2 backends on separable stacks, against quadrature of the sampled
//! fields, and against the flux balance.

use lamina_units::Length;
use ndarray::arr0;

use crate::{
    AnalyticalTensorStack, Backend, Channel, FiniteLayerIndex, IncidentSide, Polarisation,
    RealAxisEvaluator, Rotation, Scatter2, Scatter4, Stack, StackRegion, Transfer4,
    algebra::ScalarAlgebra,
    backend::{CoupledSolution, CoupledWorkspace, PlaneWaveSolutionSource, RetainedLayers},
    material::ConstitutiveTensor,
    observable::ProjectPower,
    test_support::{jet::RealJ0, planar::scalar_real_input},
};

use super::coupled::{
    C, SIDES, assert_close, c, diagonal_stack, isotropic, lossless, material, oriented,
};

type Real = <RealJ0 as crate::algebra::ComplexJet>::RealJet;

fn real(value: &Real) -> f64 {
    value.value()[()]
}

fn retain<B>(
    backend: B,
    stack: &AnalyticalTensorStack<C>,
    k0: f64,
    beta: f64,
    channel: Channel,
) -> CoupledWorkspace<RealJ0>
where
    B: Backend<RealJ0, crate::RealAxis, Workspace = CoupledWorkspace<RealJ0>>,
{
    RealAxisEvaluator::new(backend)
        .retain(scalar_real_input(k0, beta), stack, channel)
        .unwrap()
        .workspace()
        .clone()
}

fn tensors(
    _: StackRegion,
    tensors: &[ConstitutiveTensor<RealJ0>; 2],
) -> [ConstitutiveTensor<RealJ0>; 2] {
    tensors.clone()
}

fn unit(
    _: StackRegion,
    tensors: &[ConstitutiveTensor<RealJ0>; 2],
) -> [ConstitutiveTensor<RealJ0>; 2] {
    let one = |tensor: &ConstitutiveTensor<RealJ0>| {
        let like = tensor.entries().next().unwrap();
        ConstitutiveTensor::Isotropic(like.one_like())
    };
    [one(&tensors[0]), one(&tensors[1])]
}

/// A lossy stack of general orientations between lossless exteriors.
fn lossy_full_stack() -> AnalyticalTensorStack<C> {
    let tensor = |values: [(f64, f64); 3], rotation| {
        oriented(
            values.map(|(epsilon, mu)| material(c(epsilon, 0.2 * epsilon), c(mu, 0.05))),
            rotation,
        )
    };
    Stack::builder(
        isotropic(lossless(1.3)),
        oriented(
            [lossless(2.0), lossless(2.4), lossless(3.0)],
            Rotation::from_euler(0.2, 0.5, 0.9),
        ),
    )
    .layer(
        tensor(
            [(2.2, 1.0), (3.1, 1.3), (4.0, 0.8)],
            Rotation::from_euler(1.1, 0.4, -0.7),
        ),
        Length::centimetres(0.2),
    )
    .layer(
        tensor(
            [(1.5, 1.2), (2.8, 1.0), (2.0, 1.1)],
            Rotation::from_euler(-0.6, 1.2, 0.3),
        ),
        Length::centimetres(0.15),
    )
    .finalise()
}

/// The lossless counterpart of [`lossy_full_stack`].
fn lossless_full_stack() -> AnalyticalTensorStack<C> {
    Stack::builder(
        isotropic(lossless(1.3)),
        oriented(
            [lossless(2.0), lossless(2.4), lossless(3.0)],
            Rotation::from_euler(0.2, 0.5, 0.9),
        ),
    )
    .layer(
        oriented(
            [lossless(2.2), lossless(3.1), lossless(4.0)],
            Rotation::from_euler(1.1, 0.4, -0.7),
        ),
        Length::centimetres(0.2),
    )
    .finalise()
}

/// The 2×2 layer observables of one excitation, values only.
struct Scalar {
    power: Vec<[f64; 3]>,
    energy: Vec<[f64; 3]>,
    dissipation: Vec<[f64; 3]>,
    participation: Vec<[f64; 3]>,
    confinement: [f64; 3],
}

fn three(parts: (impl AsRef<[f64]>, impl AsRef<[f64]>, impl AsRef<[f64]>)) -> [f64; 3] {
    [
        parts.0.as_ref()[0],
        parts.1.as_ref()[0],
        parts.2.as_ref()[0],
    ]
}

fn first_layer(index: FiniteLayerIndex) -> bool {
    index.get() == 0
}

fn scalar(
    stack: &AnalyticalTensorStack<C>,
    k0: f64,
    beta: f64,
    polarisation: Polarisation,
    side: IncidentSide,
) -> Scalar {
    let state = RealAxisEvaluator::new(Scatter2::new())
        .retain(scalar_real_input(k0, beta), stack, polarisation)
        .unwrap();
    let excitation = state.excitation(side).unwrap();
    let array = |value: &ndarray::Array0<f64>| [value[()]];
    let power = excitation.layer_power().unwrap();
    let energy = excitation.layer_energy_nondispersive().unwrap();
    let dissipation = excitation.layer_dissipation().unwrap();
    let participation = excitation.layer_participation_nondispersive().unwrap();
    let confinement = excitation
        .layer_confinement_by_nondispersive(first_layer)
        .unwrap();
    Scalar {
        power: power
            .iter()
            .map(|p| {
                three((
                    array(p.left_flux()),
                    array(p.right_flux()),
                    array(p.absorbed()),
                ))
            })
            .collect(),
        energy: energy
            .iter()
            .map(|e| three((array(e.electric()), array(e.magnetic()), array(e.total()))))
            .collect(),
        dissipation: dissipation
            .iter()
            .map(|e| three((array(e.electric()), array(e.magnetic()), array(e.total()))))
            .collect(),
        participation: participation
            .iter()
            .map(|e| three((array(e.electric()), array(e.magnetic()), array(e.total()))))
            .collect(),
        confinement: three((
            array(confinement.electric()),
            array(confinement.magnetic()),
            array(confinement.total()),
        )),
    }
}

/// The 4×4 layer observables of one excitation, values only.
fn coupled(workspace: &CoupledWorkspace<RealJ0>, side: IncidentSide, wave: usize) -> Scalar {
    let values = |parts: (Real, Real, Real)| [real(&parts.0), real(&parts.1), real(&parts.2)];
    let energy = workspace.layer_energy(side, wave, tensors).unwrap();
    let participation = energy.participation().unwrap();
    let confinement = energy
        .confinement_by(|index, _| first_layer(index))
        .unwrap();
    let collect = |layers: Vec<[f64; 3]>| layers;
    Scalar {
        power: collect(
            workspace
                .layer_power(side, wave)
                .unwrap()
                .into_inner()
                .into_iter()
                .map(|p| values(p.into_parts()))
                .collect(),
        ),
        energy: energy
            .iter()
            .map(|e| values(e.clone().into_parts()))
            .collect(),
        dissipation: workspace
            .layer_dissipation(side, wave)
            .unwrap()
            .into_inner()
            .into_iter()
            .map(|e| values(e.into_parts()))
            .collect(),
        participation: participation
            .into_inner()
            .into_iter()
            .map(|e| values(e.into_parts()))
            .collect(),
        confinement: values(confinement.into_parts()),
    }
}

fn assert_observables_close(actual: &Scalar, expected: &Scalar, tolerance: f64, label: &str) {
    let pairs = [
        ("power", &actual.power, &expected.power),
        ("energy", &actual.energy, &expected.energy),
        ("dissipation", &actual.dissipation, &expected.dissipation),
        (
            "participation",
            &actual.participation,
            &expected.participation,
        ),
    ];
    for (name, actual, expected) in pairs {
        assert_eq!(actual.len(), expected.len(), "{label} {name}");
        for (layer, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
            for part in 0..3 {
                assert_close(
                    c(a[part], 0.0),
                    c(e[part], 0.0),
                    tolerance,
                    &format!("{label} {name} layer {layer} part {part}"),
                );
            }
        }
    }
    for part in 0..3 {
        assert_close(
            c(actual.confinement[part], 0.0),
            c(expected.confinement[part], 0.0),
            tolerance,
            &format!("{label} confinement part {part}"),
        );
    }
}

#[test]
fn separable_observables_match_the_2x2_backends() {
    let stack = diagonal_stack();
    let polarisations = [
        Polarisation::TransverseElectric,
        Polarisation::TransverseMagnetic,
    ];
    for (k0, beta) in [(2.0, 0.3), (2.6, 1.9)] {
        let coupled_workspaces = [
            retain(Scatter4::new(), &stack, k0, beta, Channel::Coupled),
            retain(Transfer4::new(), &stack, k0, beta, Channel::Coupled),
        ];
        for (wave, polarisation) in polarisations.into_iter().enumerate() {
            for side in SIDES {
                let label = format!("{k0} {beta} {polarisation:?} {side:?}");
                let expected = scalar(&stack, k0, beta, polarisation, side);
                // Scatter4: the scalar channel and the coupled channel's
                // eigen-wave of the same polarisation (TE first).
                let channel = retain(Scatter4::new(), &stack, k0, beta, polarisation.into());
                assert_observables_close(
                    &coupled(&channel, side, 0),
                    &expected,
                    1e-9,
                    &format!("Scatter4 scalar {label}"),
                );
                assert_observables_close(
                    &coupled(&coupled_workspaces[0], side, wave),
                    &expected,
                    1e-9,
                    &format!("Scatter4 coupled {label}"),
                );
                // Transfer4 smoke test on the coupled channel.
                assert_observables_close(
                    &coupled(&coupled_workspaces[1], side, wave),
                    &expected,
                    1e-8,
                    &format!("Transfer4 coupled {label}"),
                );
            }
        }
    }
}

#[test]
fn dissipation_equals_flux_loss_and_layers_telescope_to_the_absorptance() {
    let stack = lossy_full_stack();
    let (k0, beta) = (2.7, 1.3);
    for (name, workspace) in [
        (
            "Scatter4",
            retain(Scatter4::new(), &stack, k0, beta, Channel::Coupled),
        ),
        (
            "Transfer4",
            retain(Transfer4::new(), &stack, k0, beta, Channel::Coupled),
        ),
    ] {
        for side in SIDES {
            let solution = workspace.solution();
            let power = solution.entries().project_power(solution.context(), side);
            for wave in 0..2 {
                let label = format!("{name} {side:?} {wave}");
                let layers = workspace.layer_power(side, wave).unwrap();
                let dissipation = workspace.layer_dissipation(side, wave).unwrap();
                let mut absorbed = 0.0;
                for (layer, loss) in layers.iter().zip(dissipation.iter()) {
                    assert!(real(loss.total()) > 1e-3, "{label}: lossy layers dissipate");
                    assert_close(
                        c(real(loss.total()), 0.0),
                        c(real(layer.absorbed()), 0.0),
                        1e-9,
                        &format!("{label} dissipation"),
                    );
                    absorbed += real(layer.absorbed());
                }
                let absorptance = real(&power.absorptance()[wave]);
                assert_close(c(absorbed, 0.0), c(absorptance, 0.0), 1e-9, &label);
            }
        }
    }
}

#[test]
fn lossless_stacks_balance_energy() {
    let stack = lossless_full_stack();
    let (k0, beta) = (2.7, 1.3);
    let workspace = retain(Scatter4::new(), &stack, k0, beta, Channel::Coupled);
    let solution = workspace.solution();
    for side in SIDES {
        let power = solution.entries().project_power(solution.context(), side);
        for wave in 0..2 {
            let total = real(&power.reflectance()[wave]) + real(&power.transmittance()[wave]);
            assert_close(c(total, 0.0), c(1.0, 0.0), 1e-10, "R + T");
            for layer in workspace.layer_power(side, wave).unwrap().iter() {
                assert!(real(layer.absorbed()).abs() < 1e-10);
            }
            for layer in workspace.layer_dissipation(side, wave).unwrap().iter() {
                assert!(real(layer.total()).abs() < 1e-10);
            }
            for layer in workspace.layer_energy(side, wave, tensors).unwrap().iter() {
                assert!(real(layer.electric()) > 0.0 && real(layer.magnetic()) > 0.0);
            }
        }
    }
}

/// Composite five-point Gauss–Legendre quadrature of `f` on `[0, length]`.
fn quadrature(length: f64, panels: usize, f: impl Fn(f64) -> C) -> C {
    let nodes = [
        (-0.906_179_845_938_664, 0.236_926_885_056_189),
        (-0.538_469_310_105_683, 0.478_628_670_499_366),
        (0.0, 0.568_888_888_888_889),
        (0.538_469_310_105_683, 0.478_628_670_499_366),
        (0.906_179_845_938_664, 0.236_926_885_056_189),
    ];
    let h = length / panels as f64;
    (0..panels)
        .flat_map(|panel| {
            nodes.map(|(node, weight)| (h * (panel as f64 + 0.5 + 0.5 * node), 0.5 * h * weight))
        })
        .map(|(z, weight)| f(z) * weight)
        .sum()
}

type Fields = [[C; 3]; 2];

fn sampled(solution: &CoupledSolution<'_, RealJ0>, region: StackRegion, depth: f64) -> Fields {
    solution
        .fields(region, &RealJ0::lift_constant(arr0(c(depth, 0.0))))
        .unwrap()
        .map(|field| field.map(|entry| entry.value()[()]))
}

fn apply(tensor: &ConstitutiveTensor<RealJ0>, field: &[C; 3]) -> [C; 3] {
    let full = tensor.full();
    std::array::from_fn(|row| {
        (0..3)
            .map(|column| full[row][column].value()[()] * field[column])
            .sum()
    })
}

/// Every closed-form integral of `region` against quadrature on `[0, length]`.
fn assert_integrals_match(
    a: &CoupledSolution<'_, RealJ0>,
    b: &CoupledSolution<'_, RealJ0>,
    region: StackRegion,
    length: f64,
    weights: Option<[ConstitutiveTensor<RealJ0>; 2]>,
    tolerance: f64,
    label: &str,
) {
    let weight = |part: usize, field: &[C; 3]| match &weights {
        Some(tensors) => apply(&tensors[part], field),
        None => *field,
    };
    let given = |_: StackRegion, tensors: &[ConstitutiveTensor<RealJ0>; 2]| {
        weights
            .clone()
            .unwrap_or_else(|| unit(StackRegion::LeftExterior, tensors))
    };
    let panels = 200;
    let pair = |hermitian: bool, part: usize| {
        quadrature(length, panels, |z| {
            let (fa, fb) = (sampled(a, region, z), sampled(b, region, z));
            let wb = weight(part, &fb[part]);
            (0..3)
                .map(|k| if hermitian { fa[part][k].conj() } else { fa[part][k] } * wb[k])
                .sum()
        })
    };
    let hermitian = a.hermitian_overlap(b, region, given).unwrap();
    let bilinear = a.bilinear_overlap(b, region, given).unwrap();
    for (part, (h, l)) in [
        (hermitian.electric(), bilinear.electric()),
        (hermitian.magnetic(), bilinear.magnetic()),
    ]
    .into_iter()
    .enumerate()
    {
        assert_close(
            h.value()[()],
            pair(true, part),
            tolerance,
            &format!("{label} Hermitian {part}"),
        );
        assert_close(
            l.value()[()],
            pair(false, part),
            tolerance,
            &format!("{label} bilinear {part}"),
        );
    }
    let power = quadrature(length, panels, |z| {
        let [e, eta] = sampled(a, region, z);
        (e[1] * eta[2].conj() - e[2] * eta[1].conj()) * 0.5
    });
    assert_close(
        c(real(&a.longitudinal_power(region).unwrap()), 0.0),
        c(power.re, 0.0),
        tolerance,
        &format!("{label} longitudinal power"),
    );
}

#[test]
fn layer_integrals_match_quadrature_of_sampled_fields() {
    let stack = lossy_full_stack();
    let (k0, beta) = (2.7, 1.3);
    for (name, workspace) in [
        (
            "Scatter4",
            retain(Scatter4::new(), &stack, k0, beta, Channel::Coupled),
        ),
        (
            "Transfer4",
            retain(Transfer4::new(), &stack, k0, beta, Channel::Coupled),
        ),
    ] {
        let a = workspace.excitation(IncidentSide::Left, 0).unwrap();
        let b = workspace.excitation(IncidentSide::Right, 1).unwrap();
        for index in 0..2 {
            let region = StackRegion::Layer(FiniteLayerIndex::new(index));
            let thickness = workspace.layer_thickness(index).unwrap().value()[()].re;
            let tensors = workspace.layer_tensors(index);
            for (pair, (x, y)) in [(&a, &a), (&a, &b), (&b, &a)].into_iter().enumerate() {
                for weights in [None, tensors.clone()] {
                    assert_integrals_match(
                        x,
                        y,
                        region,
                        thickness,
                        weights,
                        1e-10,
                        &format!("{name} layer {index} pair {pair}"),
                    );
                }
            }
        }
    }
}

#[test]
fn integrals_hold_at_a_layer_cutoff_and_over_decaying_exteriors() {
    // At β = 1.5 k₀ the layer's forward and backward eigenvalues meet at 0
    // (a defective generator); the left exterior is evanescent, the right
    // one propagating.
    let stack = Stack::builder(isotropic(lossless(1.0)), isotropic(lossless(3.0)))
        .layer(
            oriented(
                [lossless(2.25), lossless(2.25), lossless(2.25)],
                Rotation::from_euler(0.3, 0.8, 0.1),
            ),
            Length::centimetres(0.3),
        )
        .layer(
            oriented(
                [lossless(2.0), lossless(2.6), lossless(3.3)],
                Rotation::from_euler(0.7, 0.2, 1.1),
            ),
            Length::centimetres(0.25),
        )
        .finalise();
    let (k0, beta) = (2.0, 3.0);
    let workspace = retain(Transfer4::new(), &stack, k0, beta, Channel::Coupled);
    let a = workspace.excitation(IncidentSide::Right, 0).unwrap();
    let b = workspace.excitation(IncidentSide::Right, 1).unwrap();
    for (index, thickness) in [0.3, 0.25].into_iter().enumerate() {
        let region = StackRegion::Layer(FiniteLayerIndex::new(index));
        assert_integrals_match(
            &a,
            &b,
            region,
            thickness,
            None,
            1e-10,
            &format!("layer {index}"),
        );
        assert_integrals_match(
            &a,
            &b,
            region,
            thickness,
            workspace.layer_tensors(index),
            1e-10,
            &format!("weighted layer {index}"),
        );
    }
    // Decay rate k₀ √(ξ² − 1) ≈ 2.2: the tail beyond 12 cm is below 1e-20.
    for (x, y) in [(&a, &a), (&a, &b)] {
        assert_integrals_match(
            x,
            y,
            StackRegion::LeftExterior,
            12.0,
            None,
            1e-10,
            "left tail",
        );
    }
    // The incident side holds incoming waves: no tail.
    assert!(a.longitudinal_power(StackRegion::RightExterior).is_none());
    // A left excitation's own side is incident.
    let left = workspace.excitation(IncidentSide::Left, 0).unwrap();
    assert!(left.longitudinal_power(StackRegion::LeftExterior).is_none());
}
