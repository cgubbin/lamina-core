# Tasks

## 1. Execution scope and preserved baseline

- [x] 1.1 Reconcile the core-local CLI scope with an authorized sibling lamina-mode executor before writing source; verify actual writable/executor roots, record both starting diffs, and preserve existing upstream work without creating workflow configuration.
- [x] 1.2 Run locked/offline all-target baselines in core and mode and compare with the archived 1481/206 passing tests; record current results and known warnings in design.md and verify no unrequested lockfile changes.

## 2. Coherent physical problem construction

- [x] 2.1 Add the internal immutable problem context and in-plane constructor using one analytical stack, polarization, and fixed real spectral coordinate; derive both exterior constitutive pairs from that stack and verify TE/TM coordinate/exterior jets and checked values against direct evaluation, including dispersive exterior evaluation at fixed k0 and invalid-input rejection.
- [x] 2.2 Add explicit constant exterior descriptors and a mode-local analytical constant adapter for spectral stack assembly; verify finite complex epsilon/mu, real/complex evaluation agreement, zero derivatives, lossy/magnetic constants, and rejection of non-finite/zero map factors, and document why arbitrary erased handles do not establish constancy.
- [x] 2.3 Add the spectral constructor assembling its optical stack and map from the same descriptors plus ordered analytical finite layers; verify dispersive finite-layer values and conformal derivatives against direct checked evaluation, reject or make unavailable arbitrary dispersive-exterior construction, and test both fixed-coordinate validation and epsilon/mu consistency individually.
- [x] 2.4 Document the two supported construction forms and immutability obligations alongside the module; verify the bound interface has no independent map/exterior replacement path and no mutable context replacement, retaining low-level mathematical test seams without exposing a public facade.

## 3. Problem-owned domains and explicit assumptions

- [x] 3.1 Implement problem-owned domain creation for explicit rectangles with existing sheet and explicit reference-coordinate selection; verify both map families preserve reference coordinates, exterior jets, original contour bounds, and a single validated halo, with unsafe chart/halo failures occurring before search.
- [x] 3.2 Bind domains to immutable problem ownership; verify a foreign problem's domain is rejected before any optical evaluation, clones of the same problem can reuse it, and repeated execution does not expand its halo or mutate geometry.
- [x] 3.3 Require and retain an explicit caller-assumed holomorphic/pole-free declaration with no default; verify missing declarations cannot create executable domains and document fixed-spectral versus varying-spectral material obligations over a neighborhood containing the entire evaluation region, without claiming certification from finite samples.

## 4. Bound execution and provenance

- [x] 4.1 Delegate bound execution to the existing checked single-chart kernel and wrap its complete result with read-only problem/domain provenance; verify polarization, direction, fixed coordinate, owner context, requested/evaluation bounds, reference, and assumption survive while estimates, multiplicities, leaves, run summary, and termination remain intact.
- [x] 4.2 Preserve construction and optical error sources through the bound path; verify invalid stack/chart construction, foreign-domain attempts, and injected checked input/backend/output failures return identifiable errors without panics or successful empty results, and document construction errors separately from root acceptance.
- [x] 4.3 Exercise budget-limited and cluster results through the bound path; verify unresolved leaves and multiplicity metadata remain available and document that provenance and successful execution do not certify accepted roots or complete coverage.

## 5. Physical integration and roadmap handoff

- [x] 5.1 Migrate or extend the existing physical TE/TM slab regressions in both map families to use bound problem/domain construction; verify independent guided and outgoing-cavity roots, conformal derivatives, and branch consistency while preserving focused low-level tests and current explicit symmetric/normal-incidence chart support.
- [x] 5.2 Update the mode README with runnable focused commands and both construction flows; verify commands execute and clearly describe constant spectral exteriors, analytic caller assumptions, internal-only API status, and the continuing refinement/multi-chart roadmap.
- [x] 5.3 Run all-target tests, doctests, Clippy, documentation builds, changed-file formatting, and whitespace checks in the affected repositories; verify no new failures or warnings, record unrelated baseline exceptions, and confirm dependency lockfiles and upstream source remain unchanged unless a separately justified need was approved.
- [x] 5.4 Run `openspec validate bind-modal-problem-domain --strict`, review delivered behavior against the new and existing main requirements, and update design.md with final counts, limitations, and candidate refinement/validation as the next bounded phase before marking implementation complete.

## Workflow follow-up

- Review the completed implementation, sync its delta spec, and archive when requested.
- Continue with candidate refinement/validation, then chart orchestration, public physical-input API, fields, continuation, and profiling-led performance work; keep the numerical exclusions in design.md visible.
