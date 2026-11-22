# Design

## Context

See proposal.md for motivation and specs/modal-candidate-validation/spec.md for behavior. The archived `2026-10-18-bind-modal-problem-domain` handoff identifies this phase as next. Current verified baseline is core 1481 and mode 216 passing all-target tests, with established formatting/rustdoc/Clippy exceptions documented in that archive.

Inspection of mode and the path-dependent argloc source establishes:

- `SingularPointEstimate` includes location, multiplicity, enclosure, kind (`SinglePoint` or `ClusterCentroid`), and target. ArgLoc's localisation computes a contour-moment centroid and sets kind from enclosed count; no local root polishing is implemented there.
- `ModalEstimate` pairs that metadata with a `ConformalMode`, which currently records a projected coordinate and residual without acceptance evidence. Existing search APIs explicitly document these as estimates.
- `ModalProblem::evaluate` already returns checked F and dF/dz on an owned domain. Its chart retains a validated evaluation halo, reference coordinate, and problem ownership. Neither a second derivative nor a new core evaluation API is required.
- `ProblemSearchResult` stores the complete engine-derived result and immutable domain provenance. Provenance currently checks problem identity but has no identity token distinguishing domains of the same problem.
- Existing tests include injected checked failures, a two-root cluster with budget termination, and independently engineered TE/TM guided and spectral outgoing slab roots. They provide integration fixtures for this phase.

## Goals / Non-Goals

**Goals:** numerical acceptance for isolated simple zero estimates; bounded work and branch-preserving evaluation; explicit outcomes preserving original evidence; equal treatment of both map families.

**Non-Goals:** rigorous root/isolation certificates, cluster splitting, multiple-root Newton methods, contour recounting, duplicate merging, automatic retries of discovery, scale-management changes to the characteristic, multi-chart coverage, new public exports, fields, or continuation. Existing raw search semantics remain untouched.

## Decisions

### 1. Make refinement an additive bound operation

Add a private refinement module and a method on `ModalProblem` taking the originating domain, a borrowed `ProblemSearchResult`, and validated refinement configuration. Return a report borrowing the original result and owning one candidate outcome per original index. This preserves every original field, including elapsed time, without reconstructing or cloning engine results. Existing `search` remains discovery-only; no accepted-root claims are added to `ConformalMode` or existing iterators.

Add a private shared identity token per constructed domain, cloned into its provenance and search result. Check both problem ownership and exact domain identity before refinement. Cloned domains share identity; independently constructed equal-valued domains do not. Rectangle equality alone is insufficient because reference/configuration provenance can differ. The token requires no global counter or persistence ID.

Use a separate accepted-root type with a private constructor, containing final z, physical coordinate, residual, derivative, correction magnitude, numerical thresholds, source estimate index, and domain provenance. An estimate cannot become accepted simply by being renamed. The report retains its configuration and iteration/evaluation evidence. No results are deduplicated in this phase.

### 2. Treat engine isolation metadata conservatively

Eligibility requires zero target, kind SinglePoint, and multiplicity exactly one. Clusters remain `UnresolvedCluster`, including a single coincident repeated root represented by multiplicity greater than one. Contradictory/unsupported metadata is `IneligibleEstimate`. Require finite, ordered enclosure bounds contained in the requested rectangle, and a finite starting z strictly inside their intersection. Reject invalid metadata before optical evaluation; starting points on the boundary receive a boundary outcome.

The enclosure and count are discovery evidence under the existing caller assumption, not a new mathematical isolation certificate. Do not resolve a cluster with ordinary Newton and assign its full multiplicity to whichever root attracts it.

### 3. Use bounded damped Newton in conformal coordinates

At each current point, request one checked bundle through the existing domain. Compute the undamped correction `delta=F/F_prime`. Require a finite nonzero derivative and finite correction before testing acceptance. For a non-accepted point, try `z_trial=z-lambda*delta`, starting with lambda=1 and halving on rejection. Require strict interior membership in both the original enclosure and requested rectangle before evaluating; do not use the halo as accepted geometry, clamp, move contours, or build another chart.

For evaluated trials require strict residual reduction; a finite but unproductive trial is reduced again within the configured budget. Reuse its checked bundle when it becomes the next current point. If no admissible trial exists, return BoundaryLimited; if admissible trials cannot reduce the residual, return NoProgress. Classify non-finite step arithmetic and representational stagnation explicitly. A checked evaluation failure terminates that candidate immediately with its typed source; do not silently hide it by backtracking.

Limit accepted Newton updates, backtracking trials per update, and total bundle evaluations per candidate, including the initial and final acceptance checks. Suggested initial defaults are 32 updates, 12 reductions per update, and 256 total evaluations. Zero work budgets are valid and yield BudgetExceeded unless no evaluation is needed for an ineligible/cluster outcome; there is no acceptance without a checked point. Counters and deterministic reason precedence must be documented and tested. Configuration rejects non-finite/negative tolerances, zero/non-finite residual scales, and tolerance combinations with no positive bound before any evaluation.

A finite nonzero derivative avoids a gauge-dependent hard derivative floor. Unrepresentable division still yields an explicit derivative/correction failure. Multiple-root conditioning and extra derivative orders remain separate work.

### 4. Define acceptance using both residual and undamped correction

At a checked point strictly inside both bounds require

```
abs(F) <= residual_atol + residual_rtol * residual_scale
abs(F/F_prime) <= z_atol + z_rtol * max(1, abs(z))
```

The configured residual scale is fixed for the run, positive, in the same units/gauge as the core characteristic, and never estimated from a near-zero initial residual. Preserve raw F and its gauge. Document defaults (proposed residual_atol=1e-10, residual_rtol=0, residual_scale=1, z_atol=1e-10, z_rtol=1e-10) as numerical conventions, not universal physical tolerances; expose configuration for callers using different characteristic scales. Reject non-finite computed thresholds too.

Evaluate the proposed final point afresh before acceptance, charging that evaluation to the same budget. Recheck F, derivative, correction, both bounds, and physical coordinate projection on the original chart; record the fresh bundle as evidence. A zero F with zero derivative is not accepted as an isolated simple root. A small damped move does not satisfy the correction test: always use the undamped F/F_prime. This prevents tiny damping or a globally tiny characteristic from masquerading as convergence. Acceptance is numerical conditional on caller assumptions, with no root uncertainty or rigorous error-bound claim.

### 5. Separate per-candidate failure from discovery failure

Global input errors are invalid configuration or foreign problem/domain provenance and occur before candidate evaluation. Every eligible candidate independently yields Accepted or a diagnostic non-accepted status; evaluation failures retain the existing typed error source plus last valid point/evidence where available. Continue to other candidates after a candidate failure. This new report is an explicit collection of post-search outcomes, not a fabricated successful partial result from a failed argloc run: refinement only consumes an already successful bound discovery result.

Keep original multiplicities, leaves, termination, summary, and analytic assumption accessible unchanged through the report. A root may be accepted from budget-limited discovery, but coverage remains incomplete. Record candidate failures rather than filtering them out; empty discovery yields an empty outcome list with its original diagnostics retained.

## Risks / Trade-offs

- Small raw residuals depend on characteristic scaling → require the correction criterion and explicit gauge-aware residual configuration; test tiny scaled functions away from roots.
- Enclosure bounds can stop otherwise convergent Newton steps → return BoundaryLimited, preserving the enclosure; splitting/recounting belongs to later work.
- Finite arithmetic can stagnate or overflow → detect non-finite steps, unchanged trial points, and exhausted budgets before acceptance.
- A numerical root does not prove material holomorphy or uniqueness → retain the caller assumption and original engine isolation evidence without certification wording.
- Borrowed reports constrain result lifetimes → acceptable for the internal API; avoid premature cloning/serialization architecture.
- Root duplicates may remain → preserve one outcome per estimate and defer overlap/deduplication to chart orchestration.

## Migration Plan

Planning remains in core; source changes primarily target sibling mode. Before apply, reconcile CLI scope with the authorized sibling executor and preserve working-tree diffs. No upstream or dependency changes are planned. Implement configuration/outcomes and domain identity, then the local refinement routine with synthetic tests, then bound report integration and physical regressions. Keep low-level discovery and core APIs intact; rollback removes the additive operation without reverting existing problem/search contracts.

Run locked/offline all-target tests, doctests, Clippy, docs, changed-file formatting, and strict OpenSpec validation. Record existing baseline exceptions separately. Update the handoff with actual counts and delivered limitations. Review, sync, and archive only after implementation is verified.

## Roadmap and Resume Point

After this phase, address chart orchestration, coverage, overlaps, and deliberate symmetric-exterior/normal-incidence support, including duplicate reconciliation by branch identity. Then expose a physical-input public API and executable examples, connect accepted roots to fields/normalization, add continuation and mode matching, and profile repeated solves before performance changes. Cluster splitting and repeated-root treatment remain explicit follow-up numerical work; exact scattering cutoffs, exceptional-chart fields, finite homogeneous dynamic range, and dispersive spectral exterior maps remain exclusions. Neither search direction takes priority.

## Execution record

Implementation uses the previously approved sibling `git apply` executor; the
CLI planning root remains core and workflow configuration is unchanged. Starting
diffs are preserved in /tmp/refine-start-core.diff and /tmp/refine-start-mode.diff.
Locked/offline baselines reproduce core's 1481 and mode's 216 passing tests.
No dependency or lockfile changes were needed.

## Delivered verification and handoff

The private mode refinement layer is complete. Problem-owned execution validates
configuration, problem ownership, and an exact shared domain token before any
candidate evaluation. Reports borrow the complete original discovery and retain
one indexed outcome per estimate. Accepted roots have private construction and
record the fresh checked value/derivative, undamped correction, thresholds,
physical coordinate, configuration, counters, and domain provenance.

Eligibility treats canonical multiplicity-greater-than-one cluster estimates as
unresolved and contradictory single-point metadata as ineligible. Bounds remain
strictly inside the original enclosure and requested rectangle; the existing
validated halo is unchanged. Full steps and at most twelve halvings are attempted
by default. Zero work limits disable candidate evaluations. Update/evaluation
budget exhaustion is distinct from bounded backtracking failure: the latter
reports BoundaryLimited when every proposal is exterior, NoProgress for finite
unproductive evaluated trials, or Stagnated for unchanged floating-point points.
All trial proposals count, including unevaluated exterior proposals; the trial
counter saturates rather than overflowing on extreme configured limits.

New regression coverage includes invalid configuration and dynamic threshold
overflow; cluster/inconsistent metadata without evaluations; polynomial
convergence and outward-step damping; tiny residuals with large undamped
corrections; exact simple and derivative-degenerate zeros; work limits and
representational failure; errors at initial, trial, fresh final, and projection
checks; fresh finite checks vetoing acceptance; exact-domain rejection and clone
reuse; empty and mixed reports; successful incomplete discovery retaining its
leaves and accepting an isolated root; and opposite-coordinate-sheet sign
preservation. Existing TE/TM guided and outgoing spectral cavity fixtures now
also exercise numerical acceptance with both sheet and explicit-reference
construction, while retaining their discovery-only assertions.

Separate temporary compile probes confirmed that a caller cannot invoke an
accepted-root constructor (E0599) or populate its private fields (E0451). Probe
sources and manifest path adjustments exist only in /tmp; neither repository
contains diagnostic patches or dependency overrides.

Final locked/offline verification on 2026-10-18:

- Core all-target tests: 1473 unit tests plus eight integration tests pass (1481
  total). Doctests: one passes, with the same two existing ignored examples.
- Mode all-target tests: 228 pass, zero failures or ignored tests, up from the
  216-test baseline. Its doctest command succeeds with zero runnable tests.
- Clippy and no-dependency documentation builds succeed in both repositories.
  Core Clippy is clean. Mode retains exactly the eight existing warnings about
  type complexity, partial-order comparisons, and unit-struct default calls;
  none originate in the new refinement code. Core rustdoc retains its two known
  broken links (TransferStabilityCheck in src/lib.rs and Parameter::InPlane in
  src/input/coordinate.rs); mode docs remain warning-free.
- Every changed Rust file passes rustfmt, and both diffs pass whitespace checks.
  Repository-wide formatting checks still report only the established core
  exceptions: src/evaluate/real_axis/excitation.rs,
  src/evaluate/tests/retain/fields/electromagnetic_modal.rs, and
  src/observable/mod.rs. Mode's workspace-wide check traverses that dependency.
- Core's tracked working-tree diff matches its pre-apply snapshot exactly;
  only this change's planning artifacts were authored there. Earlier mode work
  was preserved. Neither lockfile changed, and no upstream source was edited.

Acceptance remains numerical and conditional on the caller's analytic-domain
assumption and discovery metadata. It certifies neither root uniqueness nor
complete coverage. Clusters, repeated-root treatment, exceptional scattering
charts, finite homogeneous dynamic range, and dispersive spectral exteriors
remain excluded as described above.

The next bounded phase is chart orchestration with explicit coverage and overlap
handling, preserving branch identity and reconciling duplicate accepted roots
without dropping unresolved regions. Public API/examples, fields, continuation,
and profiling-led performance follow in that order. Review the implementation,
then sync this capability and archive when requested.
