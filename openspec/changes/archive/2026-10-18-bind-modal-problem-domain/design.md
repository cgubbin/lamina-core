# Design

## Context

See proposal.md for motivation and specs/modal-search-problem/spec.md for the construction contract. The archived `2026-10-18-stabilize-single-chart-modal-search` change supplies the checked single-chart kernel and continuing roadmap. Existing main specifications remain unchanged.

Read-only inspection finds:

- `lamina-mode/src/dispersion.rs` accepts an arbitrary solver/chart pair; `CompiledPlaneWaveModeSolver` separately stores an evaluator and polarization. Its checked adapter preserves exterior jets.
- `search.rs` owns the requested rectangle separately from the chart's validated halo and retains engine diagnostics. Its error observer preserves optical causes omitted by the locked upstream oracle source chain.
- Both map constructors receive exterior constitutive values explicitly. Their fixed coordinates are currently real. The spectral map documents constant exterior constitutive products and is not a dispersive-exterior map.
- Core `Stack` exposes exterior materials and ordered layers. `AnalyticalMaterialHandle` supports complex material evaluation and derivatives but erases model capabilities; handle equality is identity, and a type name or a zero sampled derivative proves no constancy.
- Core's built-in `Constant<R>` holds real epsilon and mu; raw maps also support complex constant exteriors. A mode-local constant material adapter can preserve that support without changing upstream material interfaces.
- Mode modules are private. This change can establish an internal boundary before designing a stable public API. Last verified baseline: core 1481 all-target tests and mode 206; baseline formatting/rustdoc/Clippy exceptions are recorded in the archived design.

## Goals / Non-Goals

**Goals:** one construction authority for physical inputs and maps; explicit material compatibility for both directions; immutable domain and result provenance; composition with the existing checked kernel.

**Non-Goals:** general material certificates, acceptance/refinement, searching multiple charts, physical-range-to-chart orchestration, a public facade, automatic symmetric-exterior handling, new radiation-sheet choices, or changes to the core characteristic gauge. Preserve current real fixed-coordinate restrictions. Explicit charts already supported for symmetric/normal-incidence cases remain usable; no new automatic construction policy is implied.

## Decisions

### 1. Own a problem definition and construct its solver and map together

Add an internal `problem` module with direction-specific constructors and a shared problem context. Names are provisional; the ownership contract is not. Retain an owned stack snapshot, polarization, fixed-coordinate value with its direction and current units convention, compiled checked solver, and map definition. Do not expose mutable replacement accessors. Clones share the immutable context, and changing physical inputs constructs a new context.

For in-plane searches, accept an analytical stack and fixed real k0; evaluate each exterior's epsilon and mu at that exact point through existing core interfaces and construct `ComplexInPlaneMap` from those values. Evaluation still uses the original materials. Validate the same finite/nonzero and fixed-coordinate restrictions as the raw map; do not invent a broader physical-input API or new automatic chart selection.

For spectral searches, use a constructor accepting explicit constant exterior descriptors and ordered analytical finite layers. Assemble the actual optical stack and `ComplexWavenumberMap` from the same descriptors. A small private constant-material adapter implementing existing real/complex material and derivative interfaces supports finite complex epsilon/mu with zero spectral derivatives. Erase it into the same analytical handles used for the finite layers. Do not accept an unrelated complete stack plus replacement constants, or offer an arbitrary-handle route that implies verified constancy. Existing constant real stacks can be expressed through this path; dispersive finite layers remain fully evaluated, with their derivatives carried by the core jet.

This constructive approach avoids adding global constancy metadata to erased handles or incorrectly inferring it from samples. Raw solver/chart seams remain available for colocated mathematical tests but are not the bound production construction path. Numerical identity checks on epsilon*mu alone are insufficient because polarization admittance depends on the individual factors.

### 2. Domain construction belongs to the problem

Expose a problem method creating one domain from an explicit conformal rectangle, chart configuration, and reference selection. Support the existing coordinate-sheet and explicit reference-coordinate forms, preserving their existing validation. Build the raw map internally; callers do not provide another map. Reuse `ConformalChart` and `SingleChartModalSearch` to construct the evaluation margin once, preserving the original rectangle and analytic reference.

Use an immutable owner token shared with the problem (for example an Arc identity) and check identity before execution, or an equivalent ownership interface preventing foreign-domain submission. Independently constructed equal-valued problems do not share domain ownership; clones do. Avoid process-global counters or approximate physical equality as identity. Keep requested rectangle, evaluation bounds, reference point/coordinate, and direction-specific context available through read-only provenance accessors. Do not expose operations that accidentally add a second halo on reuse.

### 3. Require an assumption, without implying certification

Domain creation requires an explicit named value such as `CallerAssumedHolomorphicPoleFree`; it has no Default and is not created implicitly by finite evaluations. Store that value in domain/result provenance. Its documentation defines the obligation over a neighborhood of the entire closed searched domain containing the evaluation halo, including material validity, exterior branch consistency, and representability of the chosen characteristic. The in-plane obligation applies at fixed k0; the spectral obligation applies across mapped k0 values and their neighborhood.

This is a deliberate caller assertion, not an unsafe Rust precondition, mathematical proof, or root acceptance certificate. A boolean silently defaulted to true or a `Certified` name would misrepresent what the software knows. Coordinate validation and checked evaluation continue to establish only their existing guarantees.

### 4. Wrap the existing execution result with provenance

A bound search delegates to the existing checked single-chart machinery and returns a wrapper containing its unmodified `ModalSearchResult` plus immutable problem/domain context. Preserve elapsed time and other summary fields by moving the result, not reconstructing it from selected fields. Domain ownership and construction errors occur before adaptive evaluation. Optical failures retain their existing source chain; construction errors retain relevant map, chart, and stack causes. No unwrap is used to convert fallible construction into a panic.

Provenance records polarization, direction and fixed coordinate, stack context identity, requested/evaluation bounds, reference selection, and the explicit domain assumption. It should support later refinement without prematurely designing serialization, persistence IDs, or accepted-mode objects.

## Risks / Trade-offs

- Erased models cannot prove spectral exterior constancy → construct spectral stacks from known constant descriptors; test that no arbitrary-handle shortcut freezes an exterior.
- Complex constant adapter could diverge between real/complex evaluation paths → test epsilon, mu, and derivatives in both paths, plus lossy/magnetic examples supported by raw maps.
- Domain identity may reject separately built equivalent problems → intentional conservative behavior; shared clones reuse the same context.
- Custom materials can use interior mutability → document that material behavior must remain unchanged over a problem lifetime; owning/cloning a handle is not a deep immutable snapshot of arbitrary external state.
- Domain declaration could be mistaken for validation → use assumption language throughout and keep finite-sample checks explicitly separate.
- Halo reuse could change geometry → store the already validated search/domain once and assert unchanged bounds across reuse.
- Existing low-level constructors remain permissive → migrate the real physical-search fixtures to the bound path while retaining focused low-level tests.

## Migration Plan

Planning lives in core; implementation primarily lives in sibling mode. Before applying, reconcile the CLI's core-only edit scope using the previously approved sibling patch executor or an explicitly authorized equivalent, without changing workflow configuration. Preserve current uncommitted work in every crate. No dependency/lockfile changes or upstream edits are expected; report a concrete need before broadening implementation.

Add the problem and domain layer, then connect result provenance and migrate physical fixtures/developer examples. Keep low-level tests proving the kernel in isolation. Verify TE/TM in both directions, dispersive-layer derivatives, mismatch prevention, foreign-domain rejection, explicit assumptions, and checked errors. Run locked tests, doctests, Clippy, docs, touched-file formatting, and strict spec validation. Treat the documented baseline warnings separately. Rollback removes the new internal layer and restores fixture entry points without reverting prior kernel repairs.

## Roadmap and Resume Point

This change ends with a coherent internal physical problem and explicitly assumed domain; results remain estimates. Next refine and validate isolated candidates with branch-preserving steps, residual/convergence criteria, multiplicity handling, and unresolved-region evidence. Then orchestrate charts and overlaps, including deliberate symmetric-exterior/normal-incidence coverage; expose the public physical-input API; connect accepted modes to fields and normalization; add continuation and mode matching; profile before performance work. Dispersive spectral exteriors, exact scattering cutoffs, exceptional-chart reconstruction, and arbitrary homogeneous dynamic range remain visible numerical limitations. Neither search direction takes priority.

## Execution record

The selected planning root is core. Implementation uses the previously authorized
`git -C /home/kit/Repos/lamina-mode apply` executor for sibling writes; no workflow
configuration is changed. Starting diffs are preserved in /tmp/bind-start-core.diff
and /tmp/bind-start-mode.diff. The locked/offline baselines reproduce core's 1481
and mode's 206 passing tests. Existing formatting, Clippy, and rustdoc exceptions
remain as recorded in the archived predecessor. No lockfiles changed.

## Delivered verification and next phase

The implementation adds mode's private `problem` module, an explicit complex
constant-exterior material adapter, immutable context/domain provenance, and a
wrapper around the complete existing single-chart result. It uses complex f64
scalar jets and real fixed coordinates in cm^-1, matching the current fixtures;
backend and conformal map types remain generic internally. No public facade or
new upstream material metadata was introduced. Chart invalid-region errors now
expose their original region cause through Error::source; bound optical errors
are boxed to preserve their complete source information without large error
return values.

Ten new regressions cover stack-derived fixed-spectral exterior values, TE/TM
checked values and exterior jets, real/complex constant values and zero
derivatives, dispersive finite-layer derivatives, both map families and reference
forms, early invalid construction, foreign-domain rejection before evaluation,
checked input/backend/output source chains, budget-limited clusters, stable halo
reuse, result provenance, and physical guided/outgoing slab roots in both search
directions. Existing low-level kernel regressions remain intact. The spectral
normal-incidence cavity continues to use an explicit chart; no automatic coverage
or candidate acceptance claim is added.

Temporary compile probes in /tmp, using a separate copied crate and the same
locked dependencies, confirmed the construction restrictions: omitting the
analytic assumption fails with E0061, supplying arbitrary analytical exterior
handles fails with E0308, and replacing context/map fields fails with E0616. No
probe source or temporary dependency override was added to either working tree.
The explicit enum has no Default; its caller-assumption wording and immutable
provenance distinguish it from coordinate validation and numerical sampling.

Final locked/offline verification:

- Core: 1473 unit tests plus eight integration tests pass (1481 total).
  Doctests: one passes and the same two existing examples remain ignored.
- Mode: 216 all-target tests pass, zero failures or ignored tests. All ten bound
  problem tests pass, including the physical searches and README's focused
  regression commands. Its doctest command succeeds with zero runnable tests.
- Clippy and no-dependency documentation builds succeed in both repositories.
  Core Clippy is clean; mode retains the same eight existing warnings concerning
  type complexity, partial-order comparisons, and unit-struct default calls.
  Core documentation retains its two existing broken links, TransferStabilityCheck
  in src/lib.rs and Parameter::InPlane in src/input/coordinate.rs. Mode docs are
  warning-free. No new warnings remain.
- All changed Rust files pass rustfmt and both diffs pass whitespace checks.
  Repository-wide formatting still reports only the three established core
  exceptions: src/evaluate/real_axis/excitation.rs,
  src/evaluate/tests/retain/fields/electromagnetic_modal.rs, and
  src/observable/mod.rs. Mode's all-workspace formatting traverses that dependency.
- Core's tracked working-tree diff matches the saved pre-apply diff exactly;
  only this change's planning artifacts were authored there. Previous mode work
  remains intact. Both lockfiles are unchanged, and no upstream source was edited.

The next bounded proposal is candidate refinement and validation on these owned
domains: preserve the selected branch during refinement, define residual and
convergence criteria, separate individual roots from unresolved multiplicities,
and retain coverage evidence when refinement fails. The remaining roadmap and
numerical exclusions above remain authoritative. Review the coordinated result,
then sync this new capability and archive when requested.
