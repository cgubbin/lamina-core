# Proposal

## Why

Complex-frequency solving is currently restricted by the constant-exterior conformal transform, although lamina-core already supports dispersive analytical stacks and material derivatives. Add a local frequency-chart path before freezing the public solver API, so that the API accepts core stacks without a parallel stack model or silently frozen exterior dispersion.

## What Changes

- Accept an existing lamina-core analytical stack for direct complex-k0 searches at fixed real in-plane wavenumber. Retain the current constant-exterior conformal path as a separate specialised strategy.
- Evaluate actual exterior constitutive values and derivatives throughout each frequency chart. Continue each exterior wavevector analytically from an explicit reference rather than selecting principal square roots independently at evaluation points.
- Require explicit material/characteristic and exterior-branch domain assumptions covering the requested rectangle and its numerical evaluation halo. Record caller-declared pole/threshold exclusions and reject intersecting charts without claiming automatic singularity discovery or certification from samples.
- Generalise the internal fallible chart boundary so discovery, refinement, physical descriptions and finite chart orchestration can share the existing evidence-preserving pipeline. Record native chart coordinates and units explicitly.
- Support explicit complex exterior branch references and real outgoing/decaying reference verification under its existing positive-dielectric assumptions, now applying to material values at the anchor rather than exterior constancy.
- Add constant-limit comparisons, dispersive TE/TM device regressions, material-derivative checks, branch-cut continuation and failure-propagation tests.
- Update the public guide and roadmap: local dispersive charts precede the public point-solver API; fields, branches, QNM normalisation, features and design optimisation retain their order.

## Capabilities

### New Capabilities

- `dispersive-spectral-search`: finite complex-frequency search charts with actual dispersive exterior materials, explicit analytic branches, domain obligations and traceable search evidence.

### Modified Capabilities

- `modal-search-problem`: scope the constant-exterior requirement to its conformal strategy, add a coherent stack-based direct-frequency path, and distinguish native chart coordinates in provenance.

## Impact

Primary target is sibling lamina-mode: chart interfaces, dispersion/search adapters, problem/domain construction, radiation provenance, refinement/orchestration/physical-description integration, tests and documentation. Reuse lamina-core stack, layer, material, polarization, unit and constitutive-lifting primitives. No new stack/material framework, dependency or upstream implementation change is expected. Interfaces remain private. Automatic full-surface uniformisation, pole/threshold discovery, contour-with-holes integration, homotopy, disjoint-chart radiation transport, public exports, broader physical labels, fields and QNM normalisation are outside this change.
