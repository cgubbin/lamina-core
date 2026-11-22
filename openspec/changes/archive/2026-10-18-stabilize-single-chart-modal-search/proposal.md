# Proposal

## Why

The next solver milestone needs a dependable single-chart search foundation. `lamina-mode` compiles against the current dependencies, but its test baseline has two failures, some derivative tests discard evaluation errors, and its determinant adapter bypasses the checked core API; core scattering composition also still fails at removable interface and intermediate-pole charts.

## What Changes

- Establish a passing lamina-mode baseline: fix the residual comparison, repair the adaptive contour/chart boundary mismatch, assert derivative-test Results, and retire unused prototype code and debug output.
- Carry a complementary scattering determinant in core and construct interface projective entries directly, avoiding divisions by opposite-admittance sums and common intermediate projective denominators.
- Preserve the characteristic gauge and jet derivatives; validate repaired cases against independent slab equations and the transfer reference where well-conditioned.
- Route mode dispersion evaluation through `checked_determinant`, preserving failures through argloc and modal-search errors. Preserve existing successful-result termination and unresolved-region diagnostics without claiming root acceptance or complete coverage.
- Keep the requested search rectangle separate from a validated evaluation halo; report an invalid domain rather than expanding a chart across a branch point.
- Preserve the longer roadmap and document the remaining cutoff, material-domain, and representability restrictions.

## Capabilities

### New Capabilities

- `single-chart-modal-search`: branch-preserving checked evaluation, contour-domain handling, error propagation, and diagnostic preservation for one analytic chart.

### Modified Capabilities

- `modal-search-domain`: add requirements for removable scattering-chart continuation and compatibility with the existing characteristic contract; retain the six established requirements.

## Impact

The coordinated implementation targets core scattering backend/projection/workspace/tests and the sibling `/home/kit/Repos/lamina-mode` dispersion, conformal chart, search, tests, and documentation. No changes to argloc, lamina-materials, or lamina-units are planned, and no new numerical dependency is required. Planning artifacts live in core's existing OpenSpec root. Its current apply scope permits core edits only: before applying sibling tasks, the execution scope must explicitly include lamina-mode; do not initialize or alter its workflow configuration as a side effect.

The plan preserves existing raw core APIs and valid ordinary scattering behavior. It does not deliver a public solver API, candidate refinement or multiplicity resolution, automatic material-domain certification, multi-chart orchestration, exact-cutoff scattering fields, continuation, or a universal overflow cure. Those remain later roadmap stages. This turn writes planning artifacts only.
