# Proposal

## Why

The solver can now discover and numerically accept roots within one owned chart, but a physical search can require several charts with overlaps and excluded regions. Coordinate those searches without hiding failed coverage, counting overlapping observations twice, or merging roots from different radiation sheets.

## What Changes

- Add an internal search plan containing an ordered set of domains belonging to one physical problem, with explicit branch selections and analytic-domain assumptions.
- Connect the existing canonical geometry builders: one in-plane chart or the lower/upper spectral charts where supported. Record finite truncations, excluded strips, and unsupported geometry explicitly.
- Support explicit finite chart plans for symmetric exteriors, normal incidence, and other cases outside canonical geometry support; require caller-selected bounds instead of inventing an automatic truncation.
- Execute and refine each chart independently, retaining every chart's discovery result, accepted roots, unresolved candidates, budgets, and typed failures.
- Produce a coverage ledger distinguishing executed charts, failed charts, unresolved search/refinement evidence, and declared exclusions. Coverage remains relative to the finite plan and caller assumptions.
- Group duplicate accepted observations only where overlapping charts supply compatible enclosure and physical branch evidence. Preserve all contributing observations; leave ambiguous matches separate and visible.

The first version is sequential and internal. It does not add automatic cluster splitting, seam recovery, global completeness certification, arbitrary physical-window mapping, dispersive spectral exteriors, fields, continuation, or the public facade.

## Capabilities

### New Capabilities

- `modal-chart-orchestration`: Problem-owned finite chart plans, independent execution, coverage accounting, and conservative reconciliation of overlapping root observations.

### Modified Capabilities

None. Existing single-chart search, problem/domain, and numerical acceptance requirements remain unchanged; this capability composes them.

## Impact

Implementation primarily affects lamina-mode's private problem, conformal geometry, and refinement integration, plus a new orchestration layer, tests, and developer documentation. Core's characteristic evaluation and argloc's search contracts remain unchanged. Planning stays in core; apply must reconcile sibling execution scope using the authorized executor. No upstream edits, dependency migrations, or lockfile changes are planned.
