# Proposal

## Why

The solver now returns accepted modal points and validated local implicit slopes, but users still have to rebuild searches to follow a mode as frequency or in-plane wavenumber changes. A seed-based continuation API will turn isolated roots into useful dispersion branches without repeating a global search at every sample.

## What Changes

- Add an owned branch trace starting from one accepted `ModePoint` and advancing over an ordered real control-parameter schedule.
- Use physical local derivatives to predict the dependent complex coordinate, then correct from that seed in a bounded local neighborhood while preserving the source strategy and exterior branch labels.
- Support propagation branches (`k0` varied, complex `beta` solved) and spectral branches (`beta` varied, complex `k0` solved) for direct-frequency and supported conformal-spectral sources.
- Retain every accepted intermediate point, correction evidence, backend and source provenance; return typed partial progress if a step cannot be continued.
- Validate schedules and settings, and stop explicitly at branch cutoffs, degenerate slopes, correction failure or exhausted step budgets. Do not silently launch a global search or claim a global branch identity.
- Extend the physical benchmark suite with branch reports for plasmon modes, an avoided crossing, and the stopped-light planar waveguide studied in Adam Freddie Page's Imperial thesis. Reports SHALL include reproducible dispersion plots and transverse field profiles at points around the avoided crossing and stopped-light points.

## Capabilities

### New Capabilities

- `modal-branch-continuation`: seed-based, evidence-preserving continuation of accepted modal points over real physical control parameters.

### Modified Capabilities

- `modal-physical-benchmarks`: add physical branch-tracking fixtures and reproducible dispersion/field reports for the avoided-crossing and stopped-light cases.

## Impact

Adds a public continuation API, owned branch and step outcomes, validated predictor/corrector settings, tests, a runnable branch example and reproducible report artifacts in lamina-mode. It builds on the physical derivatives, field recovery and QNM normalization already merged. Seeded correction will reuse lamina-mode's checked evaluation and refinement machinery with existing core stack and chart primitives; no parallel device representation is introduced. The first version follows real control-parameter paths; complex-valued parameter paths, automatic ZGV optimization and design optimization remain future work. Benchmark plots report measured local slopes and fields without treating the solver as a ZGV certificate.
