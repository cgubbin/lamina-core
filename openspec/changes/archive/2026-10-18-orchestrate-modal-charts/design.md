# Design

## Context

See proposal.md for motivation and specs/modal-chart-orchestration/spec.md for behavior. The archived `2026-10-18-refine-modal-candidates` handoff identifies this as the next phase and records core 1481 and mode 228 passing all-target tests.

Read-only inspection finds existing `ComplexInPlaneMap::search_region` and `ComplexWavenumberMap::search_regions` builders. The spectral builder returns lower/upper rectangles separated from the real-z branch-point line and periodic boundaries by an explicit margin. It rejects nonpositive/complex exterior products for canonical geometry and reports degeneracy when r=0, including equal exterior products and fixed k_parallel=0. Explicit charts already support useful cases outside those canonical restrictions; the spectral normal-incidence cavity regression uses one.

`ModalProblem` creates immutable owned domains. `DomainProvenance` now carries exact domain identity; `ProblemSearchResult` owns complete discovery evidence; `RefinementReport` borrows its discovery and owns indexed candidate outcomes. Accepted roots retain fresh checked numerical evidence and branch provenance. There is no chart-plan orchestration or duplicate reconciliation today. Existing main requirements remain authoritative and unchanged.

## Goals / Non-Goals

**Goals:** finite coordinated plans in both search directions, explicit coverage limitations, independent typed chart outcomes, and conservative duplicate grouping for supported overlaps.

**Non-Goals:** a public facade, arbitrary physical-window preimage tessellation, automatic contour deformation/seam recovery, global completeness or material-domain certificates, cluster splitting, periodic-preimage equivalence outside overlapping conformal domains, parallel execution, or changes to the underlying maps and characteristic gauge.

## Decisions

### 1. Own ordered domains in an immutable plan

Add a private orchestration module with a plan holding a nonempty vector of domains for one cloned problem context. Validate problem ownership before execution; reject repeated exact domain identities as accidental duplicate jobs, while allowing intentionally distinct overlapping domains. Assign stable plan-local chart indices and a private shared plan identity so discovery batches cannot be submitted against another plan. Clones share identity. Keep domain identity intact for refinement; do not rebuild domains or their halos.

Expose an explicit plan constructor from already validated domains and optional declared exclusion records. Its coverage scope is exactly the listed finite domains on their branches; arbitrary omitted physical space is outside the declared scope, not automatically inferred or certified as a geometric complement.

### 2. Wrap existing canonical builders without changing their geometry

Add direction-specific problem methods that derive geometry from the problem-owned map, then create domains using explicit per-chart reference selections, configuration, and caller assumptions. The in-plane method wraps one existing canonical region; the spectral method wraps the existing lower/upper regions in their stable order. Do not silently choose radiation-sheet semantics from a physical maximum: require explicit selection for each returned chart. Bounds, extent inputs, and imaginary margins remain inspectable.

Canonical coverage metadata records the finite conformal real bounds, the excluded strips near the real-z branch-point line and periodically identified boundaries, and that a physical maximum sets truncation rather than defining a rectangular complex physical window. In-plane metadata preserves its own builder's margins and fundamental-strip conventions. Do not invent a global physical-region certificate from these convenience builders.

For degenerate/unsupported canonical inputs retain the original typed reason and provide a clear explicit-domain route. Symmetric exteriors and normal incidence are supported through explicit finite plans, tested physically. The existing single preimage escaping to Re(z)=-infinity means a maximum alone cannot define all finite bounds; silently padding that infinity would change the requested physics. New automatic degenerate geometry is not part of this phase.

### 3. Stage discovery and refinement to avoid self-referential reports

First execute the plan sequentially into an owned discovery batch: one indexed entry per chart, storing either `ProblemSearchResult` or the original typed discovery error. Preserve failures and continue independent chart jobs. There are no cancellations, automatic retries, or fabricated empty results. Validate all global configuration before any chart executes.

Then refine successful entries against the exact plan-owned domains into a batch report borrowing the discovery batch. Each entry retains its local `RefinementReport`; failed discoveries are represented by references to their existing error entries and do not receive manufactured refinement reports. The aggregate can borrow plan/batch data as needed. This staged API avoids a struct owning discoveries while containing references into itself, and preserves elapsed time and all metadata without reconstruction. The report is an explicit set of independent chart outcomes, not a successful partial return from a failed single-chart argloc call.

### 4. Reconcile only supported overlap evidence

Collect accepted observations indexed by chart and candidate. Never merge unresolved clusters or failure outcomes. Validate finite nonnegative matching tolerances before execution/reconciliation. Use separate absolute/relative tolerances for conformal positions, physical varying coordinates, and exterior wavevectors, with documented units and inspectable values; defaults should be conservative numerical conventions tested against the fixtures, not a uniqueness guarantee.

For a merge require different charts with positive-area conformal overlap, accepted z values close under conformal tolerance, each accepted point strictly inside both source enclosures, and agreement of full canonical coordinates and ordered (left,right) exterior kappa values evaluated on their original charts. These checks make source isolation evidence relevant and prevent close but disjointly enclosed distinct modes being collapsed. Reference point equality is unnecessary: two references can continue the same branch. Coordinate proximity alone is insufficient. Never swap left/right exteriors.

When both same-sign and sign-reversed exterior comparisons fall within tolerance, including near a vanishing exterior channel, classify the pair as ambiguous and keep it separate. Preserve signature-evaluation failures as reconciliation diagnostics, not discovery failures or lost roots. Require all members of a duplicate group to pass pairwise compatibility; use stable chart/candidate ordering for deterministic greedy grouping and prevent transitive distance chains. Keep a representative observation, chosen deterministically, and all contributing accepted evidence. Do not average roots, sum multiplicities, change acceptance evidence, or merge two candidates from one chart.

Disjoint periodic representations are deliberately not merged in this phase. Possible coordinate-close observations lacking the supported overlap/enclosure evidence can remain separate with an ambiguity diagnostic. Call the output duplicate groups/observations, not a certified globally unique mode count.

### 5. Make coverage evidence additive and inspectable

Provide a ledger per chart containing requested/evaluation bounds, reference/assumption provenance, discovery success/failure, original termination and leaves, and every refinement status. Include canonical exclusion/truncation metadata or explicit caller-declared exclusions. Preserve exact source references rather than reducing all leaf data to a guessed boolean resolved flag; existing result leaves are not uniformly equivalent to unresolved cells.

Summaries can count executed/failed/budget-terminated charts and accepted/non-accepted observations, but should not expose an unconditional `complete` flag. Duplicate grouping changes presentation only. A successful overlapping chart does not erase another chart's failure, budget termination, clusters, or excluded regions; proving substitute coverage would require additional geometric/counting work. Empty accepted groups with failed charts is distinguishable from a successful zero-root search.

## Risks / Trade-offs

- Independent execution can retain useful roots after another chart fails → require a visible ledger and original typed error for every failed job.
- Near-root tolerances can merge different modes → require mutually compatible source enclosures, ordered branch signatures, complete pairwise compatibility, and preserve ambiguity rather than forcing merges.
- Canonical margins omit physical roots near seams → report exclusions and leave seam recovery explicit; no completeness claim.
- Multiple returned charts may cover different physical branches → require explicit reference selections and preserve all identities.
- Borrowed reports impose lifetimes → keep staged discovery/refinement and defer public API ergonomics to the next phase.
- Distinct periodic preimages may remain duplicate observations → intentional conservative limit; do not claim global uniqueness.

## Migration Plan

Planning stays in core; implementation primarily targets sibling mode using the authorized executor after reconciling CLI scope and preserving existing diffs. No upstream, dependency, or lockfile changes are planned. Add plans and geometry tests first, staged execution and coverage tests next, then branch-aware grouping with synthetic and physical fixtures. Keep existing single-chart APIs and reports intact; rollback removes the additive orchestration layer.

Verification includes both canonical geometry families; explicit symmetric and normal-incidence finite plans; whole-chart failures and budget outcomes; unresolved clusters; overlapping observations on the same and different sheets; near-cutoff ambiguity; disjoint close roots; nontransitive matching; and physical TE/TM roots. Run locked tests, docs, Clippy, changed-file formatting, whitespace checks, and strict OpenSpec validation, recording inherited exceptions separately.

## Roadmap and Resume Point

Next expose the physical-input public solver API with sensible configuration, result ergonomics, and executable examples based on this staged internal pipeline. Then connect accepted roots to retained fields/normalization, continuation and mode matching, and profiling-led performance work. Automatic seam recovery, cluster splitting/repeated-root handling, disjoint periodic duplicate equivalence, dispersive spectral exterior maps, exact scattering cutoffs, exceptional-point reconstruction, and arbitrary homogeneous dynamic range remain visible numerical work. Neither search direction takes priority.

## Execution record

The previously approved sibling git-apply executor authorizes mode source edits;
core remains the CLI planning root. Starting diffs are preserved in
/tmp/chart-start-core.diff and /tmp/chart-start-mode.diff. Locked/offline baselines
reproduce core 1481 and mode 228 passing tests, with inherited warnings unchanged.
No workflow configuration, dependency versions, or lockfiles changed.

The orchestration discovery wrapper exposes validated numeric settings and builds
an upstream ArgumentConfig from them. The integrator's configuration fields are
private and lack validation getters; accepting arbitrary opaque configuration
would defeat pre-execution validation. This wrapper retains default quadrature
and undeformed contours while exposing validated tolerances/evaluation limits.
Each stage validates its own global settings before executing that stage's jobs.

## Delivered verification and handoff

The private `problem::orchestration` layer now owns nonempty ordered plans,
canonical and explicit geometry metadata, validated staged discovery settings,
independent indexed discovery outcomes, exact-domain batch refinement, a coverage
ledger, physical branch signatures, and conservative duplicate observation groups.
Plan and domain clones preserve identity; foreign batches/domains and duplicate
exact jobs are rejected. Successful and failed chart entries retain original
results or typed errors, including source traversal. Borrowed refinement reports
preserve summaries, elapsed time, leaves, termination, and caller assumptions.

Ten new regressions cover plan ownership, invalid global settings before work,
canonical bounds/exclusions, symmetric/normal-incidence explicit domains, canonical
unsupported lossy exteriors with valid explicit plans, independent success/failure/
success execution, overlapping failure evidence, budget-limited clusters and
accepted roots, candidate failures, signature errors, ordered radiation channels,
near-zero sign ambiguity, distinct enclosures, same-chart exclusions, deterministic
representatives, and nontransitive matching. Physical TE/TM tests check repeated
observations on overlapping in-plane guided-slab and explicit normal-incidence
spectral charts. Both canonical spectral lower/upper charts participate in an
asymmetric physical slab test with an independent scalar round-trip reference.

That canonical regression explicitly pins its exterior signs: for its selected
references the lower-chart reference uses negative left and positive right
longitudinal roots; the conjugate upper-chart solution reverses those signs.
The independently computed TE reference is approximately
`k0 = 5.29160883243551 - 0.71000693473403 i cm^-1` on the lower chart.
An all-outgoing scalar reference belongs to a different exterior sheet and is not
interchangeable with those selected canonical charts. The all-outgoing cavity
reference remains separately enforced by explicit normal-incidence plans. No
radiation-sheet coverage guarantee was inferred from canonical geometry.

Matching uses fixed configurable absolute/relative tolerances (defaults 1e-8),
mutual enclosure containment, positive-area chart overlap, full canonical
coordinates, and ordered exterior values. All group members must be pairwise
compatible. Every original observation remains available, including signature
failures and ambiguous matches; representatives are not averaged and multiplicity
is not summed. Disjoint periodic representations and same-chart candidates are
not merged. Coverage is relative to the finite plan: no completeness flag,
automatic recovery of failed overlap, or seam coverage claim is introduced.

Final locked/offline verification on 2026-10-18:

- Core all-target tests: 1473 unit tests plus eight integration tests pass
  (1481 total). Doctests: one passes; the same two examples remain ignored.
- Mode all-target tests: 238 pass, zero failures or ignored tests, up from 228.
  The ten orchestration tests and README-focused physical/error commands run
  successfully. Mode's doctest command succeeds with zero runnable tests.
- Clippy and no-dependency documentation builds succeed in both repositories.
  Core Clippy is clean; mode retains the same eight existing warnings about
  type complexity, partial-order comparisons, and unit-struct default calls.
  No warning originates in orchestration. Core rustdoc retains its two known
  broken links, TransferStabilityCheck in src/lib.rs and Parameter::InPlane in
  src/input/coordinate.rs. Mode docs remain warning-free.
- Changed Rust files pass rustfmt; both diffs pass whitespace checks.
  Workspace-wide formatting still reports only the established core exceptions:
  src/evaluate/real_axis/excitation.rs,
  src/evaluate/tests/retain/fields/electromagnetic_modal.rs, and
  src/observable/mod.rs. Mode's all-workspace check traverses that dependency.
- Core's tracked working-tree diff matches its pre-apply snapshot exactly;
  only this change's planning artifacts were authored there. Previous mode work
  remains intact. No upstream source, workflow configuration, dependencies, or
  lockfiles were changed, and no temporary override or diagnostic patch remains
  in production source.

The next bounded phase is the public physical-input solver API with documented
branch choices, finite search scope, staged result ergonomics, sensible settings,
and executable examples. It must not describe canonical regions as universally
all-outgoing or exhaustive. Fields/normalization, continuation, and profiling-led
performance follow. Remaining seam recovery, cluster/repeated-root work, periodic
reconciliation, dispersive spectral exteriors, scattering cutoffs, exceptional
reconstruction, and dynamic-range limitations stay visible in the roadmap above.
Review the implementation, then sync this capability and archive when requested.
