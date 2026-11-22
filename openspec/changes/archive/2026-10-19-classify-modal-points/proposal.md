# Proposal

## Why

The solver now finds and numerically validates roots while preserving their exterior sheets, but it cannot yet say what a root represents physically. We need explicit radiation references and traceable physical classifications before exposing a public point-solver interface.

## What Changes

- Add an opt-in, problem-owned radiation-reference construction path, distinct from the signed-coordinate sheet. At a supplied real physical reference within the requested domain, verify requested outgoing/decaying or incoming/growing exterior choices for supported positive-index dielectric exteriors; continue the existing analytic map without reselecting roots pointwise.
- Introduce physical mode-point descriptions for accepted roots from either single-chart or batch refinement. Retain physical coordinates, polarization, ordered exterior wavevectors, numerical evidence, domain identity and radiation-reference provenance.
- Report exterior decay/growth and spatial or temporal attenuation independently, with configured numerical boundary bands. Report exterior flux direction only where the formulation supports it.
- Implement the paper's guided/leaky/source/absorptive/gain interpretation for supported fixed-real-frequency dielectric propagation modes. Report boundary, unsupported and unavailable outcomes explicitly instead of guessing labels.
- Keep classification separate from numerical acceptance and duplicate grouping. Unsupported interpretation or failed physical-data evaluation must leave accepted roots, failures and coverage evidence intact.
- Extend the public mathematical guide with classification assumptions, radiation-reference examples and the handoff to a public solver API.

## Capabilities

### New Capabilities

- `modal-point-classification`: traceable physical descriptions and assumption-aware classification of numerically accepted mode points.

### Modified Capabilities

- `modal-search-problem`: add explicit validated radiation-reference selection and retain its provenance alongside the existing coordinate branch and analytic-domain assumption.

## Impact

Implementation belongs in sibling lamina-mode: a physical interpretation module, problem/domain construction and provenance, adapters over existing single-chart and batch results, regression tests, README and docs/modal-transform.md. Core primitives and dependencies are expected to remain unchanged. Existing arbitrary-sheet explicit domains remain supported with unspecified radiation provenance. Modules remain private: public exports and user-facing solver examples form the next change. This change does not add automatic chart coverage, a causal radiation-sheet proof for arbitrary complex exteriors, field reconstruction, QNM normalisation, continuation or optimisation.
