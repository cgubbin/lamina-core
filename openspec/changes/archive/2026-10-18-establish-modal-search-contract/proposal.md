# Proposal

## Why

The rough `lamina-mode` search combines validated exterior conformal charts with a determinant whose full analytic and numerical domain is not yet established. Reliable root counts require distinguishing physical mode zeros, material singularities, and removable representation failures before expanding the solver API.

## What Changes

- Add an opt-in checked complex-plane determinant entry point that rejects non-finite inputs and any non-finite output jet component while preserving backend errors.
- Evaluate the transfer characteristic analytically through regular interior cutoffs, including its derivatives, using squared-wavevector entire functions.
- Document the characteristic gauge, its relation to inverse transmission, caller responsibility for analytic pole-free domains, and the difference between a numerical check and domain certification.
- Add physical and algebraic regressions that expose scattering chart failures; report remaining failures through the checked API rather than claiming complete coverage.
- Preserve the longer roadmap: dependency baseline; determinant contract; bound search problem; validated single-chart results; chart orchestration and common geometries; public API; fields; continuation and performance.
- Neither fixed-frequency nor fixed-in-plane-wavenumber searches have priority; both should use the shared pipeline.

## Capabilities

### New Capabilities

- `modal-search-domain`: checked point evaluation, characteristic semantics, removable transfer cutoff continuation, and explicit domain responsibilities for downstream mode searches.

### Modified Capabilities

None. The current OpenSpec capability inventory is empty.

## Impact

Implementation is limited to `lamina-core`: complex-plane evaluation, scalar jet algebra, transfer layer propagation, tests, exports, and documentation. Preserve the existing raw determinant API and backend error types. No new production dependency is required. The sibling `lamina-mode` and `argloc` were inspected as consumers; migrating them, automatic domain certification, a full scattering representation redesign, candidate acceptance, and root searches are follow-up changes. This foundation explicitly does not make the existing scattering determinant globally suitable for zero counting.
