# Tasks

## 1. Execution scope and baseline

- [x] 1.1 Reconcile core-local CLI scope with an authorized sibling mode executor, preserve starting diffs, and run locked/offline baselines; verify core's 1481 and mode's 216 passing tests or record concrete differences, with no workflow configuration or lockfile changes.

## 2. Refinement contract and provenance

- [x] 2.1 Add validated refinement configuration, numerical evidence, explicit candidate outcomes, and a privately constructed accepted-root type; test invalid/non-finite tolerances and scales, overflowing thresholds, zero budgets, and inaccessible acceptance construction, and document gauge/units and proposed default criteria.
- [x] 2.2 Add immutable exact-domain identity shared by domain clones and originating search provenance; verify clones can refine their result while foreign problems and separately constructed same-problem domains fail before evaluation, preserving existing problem/domain tests.
- [x] 2.3 Implement candidate eligibility and enclosure validation; verify clusters and repeated multiplicities remain unresolved, unsupported targets and contradictory metadata are non-accepted, invalid/outside enclosures are rejected, and no ineligible candidate invokes optical evaluation; document numerical isolation evidence versus certification.

## 3. Bounded local refinement

- [x] 3.1 Implement checked first-derivative Newton updates with bounded residual-decreasing backtracking; verify simple complex polynomial roots converge, undamped outward steps are reduced before evaluation, and every evaluated point remains strictly inside its enclosure and requested domain without clamping or halo enlargement.
- [x] 3.2 Enforce both residual and undamped correction criteria with a fresh final checked bundle and branch-consistent coordinate projection; test exact simple zeros, tiny scaled functions away from roots, misleading small damped steps, and zero-derivative zeros, ensuring only valid simple-root evidence constructs accepted outcomes.
- [x] 3.3 Add deterministic iteration, trial, and total-evaluation budgets and failure diagnostics; test zero/exhausted budgets, boundary limitations, lack of residual progress, unchanged points, and non-finite correction arithmetic, verifying evaluation counters include initial/trial/final checks and cannot exceed configured limits.
- [x] 3.4 Retain checked input/backend/output failures and projection errors as typed candidate outcomes; inject failures at initial, trial, and final checks, verify source traversal and last-valid evidence, and document that optical failures are not silently hidden by backtracking.

## 4. Bound reports and preserved discovery

- [x] 4.1 Expose optional bound refinement over a borrowed successful search result and return one outcome per original estimate; verify original estimate indices, metadata, branch provenance, complete run summary, leaves, termination, and caller assumption remain unchanged, including empty discovery.
- [x] 4.2 Test mixed accepted, failed, and unresolved candidates plus successful budget-limited discovery; verify candidate failures do not omit estimates or suppress independent outcomes and local acceptance never upgrades global coverage or invents a partial result from failed discovery, with matching module documentation.

## 5. Physical integration and handoff

- [x] 5.1 Extend the existing TE/TM physical guided-slab and spectral outgoing-cavity fixtures through bound refinement; verify independent root agreement, checked residual/correction acceptance, branch-consistent coordinates, recorded enclosure containment, and both reference-selection forms while preserving discovery-only tests.
- [x] 5.2 Update the mode README with runnable refinement regression commands, configuration/acceptance semantics, and a construction-to-search-to-refinement flow; verify the commands execute and explain unresolved clusters, gauge-sensitive residuals, numerical acceptance, and remaining coverage limitations.
- [x] 5.3 Run all-target tests, doctests, Clippy, documentation, changed-file formatting, and whitespace checks in both affected repositories; verify no new failures/warnings, record established unrelated baseline exceptions, and confirm no unrequested upstream source or lockfile changes.
- [x] 5.4 Run `openspec validate refine-modal-candidates --strict`, review the implementation against new and existing requirements, and record final counts, exclusions, and chart orchestration as the next roadmap phase in design.md before completing the change.

## Workflow follow-up

- Review the complete result, sync its delta spec, and archive when requested.
- Continue with chart orchestration and overlap/coverage handling, then the public physical-input API, fields, continuation, and profiling-led performance work. Preserve cluster/repeated-root work and numerical exclusions in the roadmap.
