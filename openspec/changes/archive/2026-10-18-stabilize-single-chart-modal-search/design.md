# Design

## Context

See proposal.md and the two delta specs for scope and behavior. The previous core foundation is archived at `../archive/2026-10-18-establish-modal-search-contract/`; it delivered checked point evaluation and analytic transfer cutoff continuation. The six main `modal-search-domain` requirements remain authoritative.

Planning checks on 2026-10-18 used the existing Nix environment, locked/offline dependencies, and a temporary target directory. `lamina-mode` compiles; tests report 196 passed and two failed. The residual/coordinate equality failure differs by approximately 4e-15; the single-chart slab search fails during refinement with `Chart(CoordinateOutsideRegion)`. Three derivative-helper calls in dispersion.rs discard a Result. The obsolete function.rs prototype is not declared in lib.rs and is not a current compiler failure. The sibling working tree was clean after the checks. No source or lockfile was changed during planning.

Core currently builds ordinary interface matrices before projecting them, so division by YL+YR can fail before projective accumulation. The current cascade inserts d_left*d_right into its output homogeneous coordinates and can collapse to all zeros at an intermediate pole. The final characteristic is `2*xi_left*d/n21`, not simply d. Retained reconstruction still uses ordinary component matrices and must remain compatible at regular points.

## Goals / Non-Goals

**Goals:** remove the identified interface and intermediate-pole representation failures; preserve analytic derivatives and the characteristic gauge; establish trustworthy existing single-chart regressions; carry checked failures into the search; preserve contour coverage and diagnostic metadata.

**Non-Goals:** a new public facade, candidate refinement/acceptance, cluster resolution, automatic material singularity certification, multi-chart orchestration, continuation, general scattering cutoff regularization, exceptional-point field reconstruction, or guaranteed representability at arbitrary optical thickness. Neither complex-frequency nor complex-in-plane search takes priority. No production changes to argloc, materials, or units are planned.

## Decisions

### 1. Carry the complementary determinant instead of dividing out common factors

Extend internal homogeneous entries `(d,a,b,c,e)` with h satisfying

```
S = [a b; c e] / d
b*c - a*e = d*h
```

Carry h as independent algebraic data through composition; do not reconstruct it as `(b*c-a*e)/d` at poles. For a left tuple `(dl,a,b,c,e,hl)` and right tuple `(dr,u,v,w,x,hr)`, use

```
d_new = dl*dr - e*u
a_new = a*dr + hl*u
b_new = b*v
c_new = c*w
e_new = x*dl + hr*e
h_new = hl*hr - a*x
```

Substituting the invariant into the old rational cascade cancels its redundant dl/dr factors before numerical evaluation. All formulas use additions and products, so they carry full jets through an intermediate pole. A common rescaling must scale h along with all other entries to preserve the invariant. Do not introduce modulus-based pointwise rescaling or drop derivative factors.

Construct interface tuples directly:

```
d = YL+YR; a = YL-YR; b = 2*YR; c = 2*YL; e = YR-YL; h = YL+YR
```

Propagation has `(1,0,p,p,0,p*p)`, p=exp(i*kappa*thickness); identity has `(1,0,1,1,0,1)`. Arbitrary regular ordinary entries use d=1 and h=b*c-a*e. Pole fixtures must supply h from their analytic parameterization; a five-component pole tuple does not determine h uniquely.

For the archived counterexample, left `(0,1,2,3,6,0)` and right `(1,.2,.7,.8,.1,.54)` now produce `(-1.2,1,1.4,2.4,3.24,-.1)`, a usable finite response. A planning-only algebraic check matched the reduced recurrence against ordinary composition on 100 regular random examples to approximately 1.2e-16 absolute error and checked the invariant. The completed implementation adds independent Rust value/derivative and physical-stack regressions (results below).

This representation is preferred over evaluating Transfer2 as a fallback, because fallback does not cure exponential transfer instability. It is preferred over searching raw d, because raw homogeneous zeros do not by themselves establish characteristic zeros. Preserve `D=2*xi_left*d/n21` and the existing boundary gauge wherever defined; d=0 with n21 nonzero is an intended modal zero.

### 2. Separate total accumulation from optional ordinary reconstruction data

The total response and determinant accumulate direct projective components, never passing through ordinary interface division. Extend the workspace append boundary to accept those components. Convert to ordinary components only for existing retained reconstruction where requested. Maintain ordinary reconstruction semantics on regular stacks and both incident sides; extend its exceptional-chart guarantees only in a separate change.

The response-only and retained determinant must agree, including at the repaired physical opposite-admittance fixture. Non-finite optional directional reconstruction data must not contaminate the total determinant. This is an internal representation refactor: do not widen public backend types or replace the raw evaluator API. Update the archived characterization regression from expected failure to expected continuation, preserving independent reference evidence.

Exact interior kappa=0 still degenerates the scattering directional basis, and zero polarization factors or material poles are not covered by this repair. Do not claim otherwise; keep the transfer cutoff reference and checked numerical failures. Large finite homogeneous coordinates can still overflow, and transmission can underflow. No silently altered gauge or universal scale-management policy is part of this change.

### 3. Preserve contour geometry with a validated numerical evaluation halo

Keep the caller's search rectangle R as independent search state; do not infer it later from an enlarged chart. Build an evaluation chart on R_plus containing R and a small roundoff margin per axis. Default margin is `64*epsilon*max(1,abs(min),abs(max),width)` with outward expansion verified to be finite and strict. Reuse the original map, reference point, and reference coordinate via the existing with_reference_coordinate constructor, which validates the enlarged region and preserves the sheet. Pass the original R to argloc.

The margin is a numerical guard, not a search tolerance or a branch-crossing permission. If bounds cannot expand finitely or map validation rejects R_plus, return a domain-construction error before search. Points outside R_plus still fail the existing chart checks; no clamping, broad contains tolerance, contour shrinking, or principal-sheet reselection is allowed. Store/report R as the searched region; document the evaluation-only R_plus.

Reproduce the existing failure with an instrumented test recorder to establish actual evaluation coordinates. The proposed guard handles floating-point overshoot; it must not conceal an incorrect coordinate transform or a point materially outside R. If that trace establishes a different defect requiring argloc changes, stop and revise scope before editing a third repository. Add endpoint/corner, affine roundtrip, branchpoint-rejection, and genuine outside-point tests in lamina-mode.

### 4. Use the checked core adapter and preserve existing error ownership

Change CompiledPlaneWaveModeSolver's associated error to DeterminantEvaluationError<B::Error> and forward to checked_determinant. Keep the chart-provided coordinate/exterior jets intact, and keep one seeded solve for F and F_prime. DispersionError::Solver wraps that error. The locked argloc/quadtree stack retains it structurally, but its oracle variant omits Error::source. A small tracked ComplexFunction adapter therefore preserves the most recent evaluation failure as the typed source of ModalSearchError::Search alongside the complete argloc error. The observer shares state across engine clones and does not evaluate twice. Successful evaluation clears a stale observed failure. This compatibility adapter can be removed when upstream source traversal is complete. A post-estimate projection failure remains ModalSearchError::Mode. No panic, fabricated value, empty success, or automatic backend substitution replaces an error.

Successful results retain the complete SingularPointEstimate, leaves, summary, and termination. An error is not promised to carry a successful partial result that argloc did not return; preserve any engine state already in its error rather than inventing a recovery contract. Use a tight-budget test for successful incomplete termination and an injected evaluation failure for the separate error path. Existing ConformalMode naming is not a new acceptance guarantee; document that these remain estimates and multiplicities can describe clusters.

### 5. Make the baseline tests enforce their evidence

Replace the exact residual/coordinate comparison with the fixture's explicit numerical tolerance, retaining a sign-flip negative check. Assert or propagate every Result from derivative comparisons; never use `let _` to silence them. Remove dbg output and unused imports, and remove the uncompiled obsolete function.rs and unused mapping prototype after confirming no consumers. Replace the template README with focused internal-workflow and limitations documentation; the library remains private until the later public API change. Do not invent dependency migrations when the current compiled path already builds.

## Implementation order and verification

1. Establish the two-repository execution scope and reproduce current baselines.
2. Restore lamina-mode test validity and the contour/evaluation-region distinction with colocated tests and docs.
3. Add the invariant-carrying reduced composition, direct interface components, and workspace plumbing with algebraic, derivative, and regular retained compatibility tests.
4. Turn the physical opposite-admittance and intermediate-pole cases into continuation successes, with TE/TM and both varied coordinates; retain thick-layer failure honesty.
5. Integrate the checked adapter and source-preserving search tests. Run known slab searches in both conformal map families and preserve successful incomplete results.
6. Run all-target tests, doctests, Clippy, documentation, touched-file formatting, and strict OpenSpec validation. Keep independent scalar slab/ordinary Redheffer references in addition to cross-backend agreement. Treat known unrelated core formatting and rustdoc warnings separately; do not silently weaken acceptance to make new tests pass.

## Migration Plan

Artifacts remain in core's OpenSpec root. The CLI currently restricts implementation edits to core, and the session filesystem does not grant sibling writes. Before apply begins sibling tasks, explicitly reconcile the execution scope for `/home/kit/Repos/lamina-mode`; this proposal does not alter workflow configuration or grant filesystem access. No source changes occurred in planning. Do not mark the coordinated change complete after only the core half lands.

Land internal core fixes before switching the consumer to checked evaluation, then verify both dependency-linked crates together. Keep current path dependency versions unless a real locked resolution change requires updating mode's lockfile. Public raw determinant/backend error types remain intact. Roll back the consumer adapter separately from core internals if required; retain the new regression evidence and roadmap. Deployment, publishing, and third-party dependency releases are outside this change.

## Risks / Trade-offs

- Complement data can be initialized or scaled incorrectly: verify `b*c-a*e=d*h`, regular composition, associativity, derivatives, and explicit pole limits.
- Homogeneous dynamic range remains finite: retain typed checked failures and avoid claims of global numerical stability.
- Ordinary retained reconstruction remains singular at exceptional charts: promise determinant continuation there and regular reconstruction compatibility only.
- Enlarging the chart can encounter a coordinate branch point: validate the halo and return a domain error rather than changing the contour or sheet.
- Existing discarded derivative errors may reveal more failures: fix in-scope mathematics or report the needed scope revision; do not ignore them.
- Generic materials have no domain certificate: continue requiring explicit caller analytic-domain responsibility.

## Roadmap and resume point

This change delivers the next part of the trustworthy search kernel: a passing baseline, repaired removable scattering charts, and checked single-chart integration. It does not complete the usable solver milestone.

Continue afterward with:

1. Bind stack-derived exterior data, charts, polarization, and fixed coordinate into a coherent problem with explicit domain assumptions.
2. Refine/validate isolated candidates, distinguish multiplicities and unresolved clusters, and preserve coverage evidence.
3. Orchestrate charts for both search directions; reconcile overlaps by sheet and support symmetric exteriors and normal incidence deliberately.
4. Expose the physical-input public API, sensible defaults, executable examples, and complete documentation.
5. Connect accepted modes to retained fields, separating root acceptance, reconstruction, and normalization failures.
6. Add continuation, sheet-preserving chart deformation, parameter-sweep mode matching, and light-line transitions; treat dispersive complex-frequency exteriors as a distinct map problem.
7. Profile repeated solves before caching, scale management, and parallel performance work.

Neither search direction takes priority. Keep exact-cutoff scattering regularization and robust exceptional-chart field reconstruction visible as remaining numerical work before any broader coverage claim. Upon completion, record final test counts, remaining exclusions, and the next bounded change here before archiving.


## Delivered verification and handoff (2026-10-18)

All 19 implementation tasks are complete across core and mode. Sibling edits were
applied through the explicitly approved `git -C /home/kit/Repos/lamina-mode apply`
executor; workflow configuration was not changed. The later user authorization
also permitted necessary upstream edits, but none were required: no source in
argloc, quadtree, materials, or units was edited by this change. Their pre-existing
working tree changes were preserved. Neither lockfile changed, and no temporary
dependency override or instrumentation patch remains in production paths.

The reproduced core baseline was 1468 unit tests plus eight integration tests;
mode was 196 passing and two failing tests. The recorded adaptive point had real
coordinate `0.3299293863210311` against upper bound `0.329929386321031`, confirming
a roughly 1.1e-16 roundoff overshoot rather than a materially exterior contour.
The retained test-only recorder supports boundary rejection assertions.

Delivered guarantees are the invariant-carrying reduced cascade, polynomial
interfaces, total/retained determinant isolation, checked one-seed optical
integration, unchanged contours with validated halos, preserved physical branch
coordinates, identifiable optical failures, and unmodified successful search
diagnostics. Independent references cover TE/TM opposite-admittance continuation
and its first/second derivatives, a physical intermediate-pole multilayer,
regular Redheffer composition, guided asymmetric slab roots, and complex spectral
outgoing cavity poles. The latter uses an explicit normal-incidence chart and
makes no automatic symmetric-exterior orchestration claim. Both search directions
are tested without prioritizing one over the other.

Final locked/offline checks:

- Core all-target tests: 1473 unit tests and eight integration tests pass (1481
  total). Doctests: one passes, two existing ignored examples remain ignored.
- Mode all-target tests: 206 pass, zero failures or ignored tests. Doctests:
  zero runnable tests, successful command.
- Clippy and no-dependency documentation builds succeed in both crates. Core
  Clippy is clean. Mode retains eight existing warnings (type complexity,
  partial-order comparisons, and unit-struct default construction); new warnings
  were corrected. Core rustdoc retains two existing broken links:
  `TransferStabilityCheck` in src/lib.rs and `Parameter::InPlane` in
  src/input/coordinate.rs. Mode documentation emits no warnings.
- Changed Rust files pass rustfmt. Repository-wide formatting checks still
  report only the pre-existing core files src/evaluate/real_axis/excitation.rs,
  src/evaluate/tests/retain/fields/electromagnetic_modal.rs, and
  src/observable/mod.rs; mode's `--all` check also traverses this path dependency.
  Both diffs pass whitespace checks.

Remaining exclusions are exact scattering directional cutoffs, arbitrary
homogeneous dynamic range, exceptional-point ordinary field reconstruction,
material singularity certification, and candidate acceptance/completeness beyond
the engine's retained diagnostics. No fallback backend or characteristic gauge
change was introduced.

The next bounded change should bind stack-derived exterior data, polarization,
fixed physical coordinate, and validated chart/domain assumptions into one
coherent problem contract. Its acceptance evidence should reject mismatched
stack/map inputs, preserve branch identity in both map families, and carry the
explicit material-domain obligation. Candidate refinement/validation follows
that contract; the seven-phase roadmap above remains the continuation plan.
Review the coordinated diff, sync the two delta specs, then archive this change
when requested.
