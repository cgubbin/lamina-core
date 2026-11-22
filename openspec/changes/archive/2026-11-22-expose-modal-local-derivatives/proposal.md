# Proposal

## Why

The public evaluator returns a derivative in its search chart coordinate, which is insufficient for physical dispersion slopes and device sensitivities. Branch continuation needs separately identified physical characteristic derivatives, with material dispersion and exterior-sheet dependence retained.

## What Changes

- Expose a separately fallible local evaluation of the checked characteristic and first derivatives with respect to physical k0, physical beta and selected finite-layer thicknesses at an accepted point.
- Include actual permittivity/permeability dispersion and branch-preserving exterior derivatives, with a required local analyticity declaration covering the requested physical variations.
- Reuse core's first-order jets, parameter selectors and thickness-unit conventions. Add a narrowly scoped public directional stack-compilation entry point because the existing public compiler seeds all thicknesses as constants.
- Return owned values, parameter/unit labels, source and backend provenance, checked derivative diagnostics and guarded simple-root implicit slopes.
- Verify physical gradients against native derivatives, independent equations and finite-difference step sweeps on slabs, interface plasmons and the finite-film field fixtures.

This first increment covers core's existing spectral, in-plane and layer-thickness parameters. Material spectral derivatives are included; arbitrary material-coefficient sensitivities, field derivatives, Hessians, ZGV/exceptional-point solvers, branch tracking and optimization remain outside scope.

## Capabilities

### New Capabilities

- `modal-local-derivatives`: Physical first-order characteristic gradients and conditional simple-root slopes on accepted exterior branches, with explicit parameter units and failure evidence.

### Modified Capabilities

None. Existing native-coordinate derivatives, accepted points and search reports retain their meaning.

## Impact

Primarily sibling lamina-mode's accepted-point evaluation and public diagnostics, with a small lamina-core extension around the existing directional stack compiler. No new material/parameter framework, numerical differentiation production path or runtime dependencies are planned. Reuse the physical benchmarks and field-stage fixtures.

## Sequence

Third planned change, after `preserve-plasmon-benchmarks` and in the working sequence after `recover-modal-fields`. This completes the fields/local-derivatives milestone and hands off to branch continuation; QNM normalization, features and lamina-design optimization follow separately.
