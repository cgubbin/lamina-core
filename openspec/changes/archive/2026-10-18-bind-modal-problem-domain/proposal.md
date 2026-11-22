# Proposal

## Why

The checked single-chart kernel now passes its numerical regressions, but callers can still combine a solver with an unrelated chart or incompatible exterior material assumptions. Bind the physical problem and its declared domain before adding candidate acceptance, so subsequent results have a coherent physical interpretation.

## What Changes

- Add an internal problem construction path binding one stack, polarization, fixed physical coordinate, and search direction to stack-consistent conformal charts and checked evaluation.
- Derive complex-in-plane map data from the stack at its fixed spectral coordinate. Construct complex-spectral problems from explicit constant exterior materials shared by the stack and map, allowing dispersive finite layers without freezing dispersive exteriors.
- Bind one requested conformal rectangle, explicit reference/sheet selection, and validated evaluation halo to a problem-owned domain.
- Require an explicit caller declaration of characteristic analyticity and absence of poles over the closed domain and a surrounding neighborhood, including the evaluation halo. Preserve this declaration as an assumption, never a numerical certificate.
- Carry problem/domain provenance alongside existing estimates and diagnostics, preventing domains from unrelated problems being submitted through this path.
- Keep both search directions equally supported within current map restrictions; retain raw internal test seams and checked evaluator semantics.

Candidate refinement, mode acceptance, automatic material-domain certification, dispersive spectral exteriors, automatic chart orchestration, fields, and a public solver facade are excluded.

## Capabilities

### New Capabilities

- `modal-search-problem`: Coherent physical problem construction, compatible exterior-map data, explicit domain assumptions, and bound single-chart execution provenance.

### Modified Capabilities

None. Existing `modal-search-domain` and `single-chart-modal-search` requirements remain authoritative and unchanged; the new capability constrains their composition.

## Impact

Implementation primarily affects sibling `lamina-mode`: a new internal problem/domain layer, integration with dispersion and search, tests, and developer documentation. Planning stays in lamina-core's OpenSpec root. Existing stack accessors and material evaluation interfaces appear sufficient; no upstream dependency or lockfile migration is planned. Apply must reconcile the core-only CLI edit scope with the already authorized sibling executor before writing mode source. No implementation is included in this proposal.
