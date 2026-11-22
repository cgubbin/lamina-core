# Proposal

## Why

The conformal map formulas retain the exterior square-root sheets, but their documentation incorrectly identifies a width-pi strip as the full sheet-resolved surface. Before exposing a public solver or interpreting roots physically, we must reconcile our conventions with Bakhtazad et al. (1997), correct misleading coverage claims, and preserve a roadmap toward fields, branches and optimisation.

## What Changes

- Establish a documented correspondence between the paper's fixed-frequency transform and Lamina's coordinates, ordered exterior wavevectors, units and time/space conventions.
- Correct the distinction between the squared-coordinate period i*pi and the exterior-sheet period 2*i*pi; audit both in-plane and spectral geometry and correct any implementation that treats opposite exterior sheets as identical.
- Retain current canonical charts as explicitly limited finite selections, reporting omitted sheet sectors as well as margins and truncations. Explicit domains remain the route to other sectors; automatic exhaustive coverage is deferred.
- Publish an accessible mathematical guide linked from the lamina-mode README, with corrected rustdoc, source attribution, classification assumptions and differences from the paper.
- Record the longer-term roadmap in that guide and this design, including point solving, classification, fields, continuation, QNM normalisation, feature tracking and lamina-design integration.
- Add independent mathematical and physical regressions that distinguish opposite exterior sheets. Do not infer classification or all-outgoing QNM status from numerical root acceptance.

## Capabilities

### New Capabilities

None. This change clarifies the existing chart contract rather than introducing a classification API.

### Modified Capabilities

- `modal-chart-orchestration`: distinguish full exterior-sheet coverage from squared-coordinate repetition and expose canonical sheet-sector omissions in the coverage contract.

## Impact

Primary target: sibling lamina-mode, especially src/conformal/in_plane, src/conformal/wavenumber, src/problem/orchestration, their tests, README and a linked public mathematical guide. OpenSpec planning remains in lamina-core. Existing explicit-domain searches and conservative duplicate grouping remain supported. The modules are private; no public solver API or new dependency is introduced. Core changes are expected only if a demonstrated convention defect requires them. Runtime classification, automatic full-surface chart generation, fields, continuation, normalisation and optimisation are subsequent changes.
