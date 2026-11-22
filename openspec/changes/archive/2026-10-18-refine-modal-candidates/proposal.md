# Proposal

## Why

Bound searches now preserve physical and domain provenance, but their returned coordinates remain contour-moment estimates. Add an explicit numerical refinement and acceptance stage so downstream work can distinguish converged isolated roots from unresolved clusters, stalled candidates, and evaluation failures.

## What Changes

- Add an optional post-search refinement operation for results belonging to the same problem and exact domain.
- Refine only zero estimates representing one root, using checked conformal value/derivative evaluation on the existing branch and within the candidate enclosure and requested domain.
- Require documented residual and undamped correction criteria before producing a numerically accepted root; record configuration and convergence evidence.
- Return explicit outcomes for clusters, invalid candidates, boundary limits, derivative failures, exhausted budgets, and checked evaluation errors, retaining their causes.
- Preserve the complete original bound search result and one outcome per estimate, keeping local root acceptance separate from global coverage and material-domain assumptions.

Cluster splitting, repeated-root refinement, completeness certification, cross-chart deduplication, fields, continuation, and the public facade remain later work. Both search directions receive equal physical regression coverage.

## Capabilities

### New Capabilities

- `modal-candidate-validation`: Branch-preserving bounded refinement, numerical acceptance criteria, candidate outcomes, and preserved discovery evidence.

### Modified Capabilities

None. Existing search contracts continue to return estimates; refinement is an additive stage that consumes them without relabeling the original result.

## Impact

Implementation primarily affects lamina-mode's private problem/domain layer and a new refinement module, plus tests and developer documentation. The core checked determinant and argloc search contracts remain unchanged. Planning lives in core; sibling execution scope must use the already authorized executor or an explicitly reconciled equivalent. No upstream edits, dependency migrations, or lockfile changes are expected. This proposal contains planning artifacts only.
