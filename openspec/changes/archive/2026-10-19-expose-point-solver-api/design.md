# Design

## Context

See proposal.md for motivation. Planning lives in lamina-core; implementation is in sibling lamina-mode. The current baseline is 260 tests, six existing clippy warnings and no public solver exports. `src/lib.rs` includes the mathematical guide but keeps `problem`, `refinement`, `search`, `dispersion` and chart implementations private.

Observed implementation constraints:

- `ModalProblem<M,S>` has working fixed-frequency in-plane, constant-exterior conformal spectral and supplied-stack direct-frequency constructors. The first and third accept core analytical stacks; the constant spectral constructor builds its coherent stack from constant exterior evidence and core analytical layers.
- `ProblemDomain<M,S,H>` and `ChartPlan<M,S,H>` carry raw map/chart/compiled-solver parameters that callers should not have to name. Domains already retain immutable context, opaque identity, native-coordinate kind, requested/evaluation bounds, assumptions and references.
- Discovery owns its results. `RefinementReport` borrows discovery, and `BatchReport` borrows the discovery batch and plan. `PointDescription` borrows an accepted root. These types work internally but cannot simply be returned together inside a self-contained public convenience report.
- `AcceptedRoot` already owns cloneable numerical evidence and domain provenance. Reports already preserve failures, unresolved outcomes and overlap groups. The public layer must retain those contracts rather than implement a second solver.
- Core supplies stacks, layers, analytical material handles, polarization, backends and units. Mode's `ConstantExterior` is existing constancy construction evidence, not a new stack representation. Erased handles cannot prove exterior constancy.

The existing modal-search-problem, modal-chart-orchestration, modal-candidate-validation, modal-point-classification and dispersive-spectral-search specifications remain unchanged and authoritative.

## Goals / Non-Goals

**Goals:** A downstream caller can construct and run each supported strategy, inspect staged discovery/refinement or use a convenience solve, return an owned report from a helper, retain a point value for later work, and diagnose failures using public types. Public examples prove this boundary from outside the library.

**Non-Goals:** No new numerical method, generic user-defined chart plugin API, automatic material classification/strategy switching, mixed-strategy plan, serialization format, asynchronous execution, CLI, publication, fields, continuation or normalization. This interface does not promise arbitrary physical-window inversion into conformal charts. Scalar support stays complex f64 with real fixed k0 or beta, matching the implemented solvers.

## Decisions

### 1. Add a small facade over the three existing problem strategies

Add a public API module, with selected crate-root re-exports, rather than making every current module public. Use distinct public problem/domain/plan wrappers for propagation, constant-exterior conformal spectral and direct-frequency spectral paths. Their private fields contain the current implementation types; raw map and jet parameters do not appear in required caller annotations. Selectively re-export existing value/configuration types when their complete signature graph is public and their terminology is appropriate. Wrapping is for ownership and API clarity, not duplicate numerics.

Use explicit strategy names such as `PropagationProblem`, `ConformalSpectralProblem` and `FrequencyProblem` (final spelling can be refined consistently during implementation). There is no single spectral constructor that silently freezes a supplied stack or switches strategy based on sampled derivatives. Fixed-frequency and direct-frequency constructors accept existing core analytical stacks. The constant conformal constructor accepts the existing validated constant exterior evidence plus core analytical layers; it must never accept an independent full stack together with replacement constants.

Use core backends directly, with Scatter2 as the ergonomic default and explicit backend construction available. Preserve supported checked backend selection, including Transfer2 where its existing bounds permit it. Backend type parameters are acceptable; requiring callers to spell internal scalar jets, raw maps or compiled evaluator types is not. Compile public-only clients using both built-in backends to verify the signature boundary. Do not introduce a backend framework or silent fallback. Expose a public complex-f64 alias so clients need not add a separate algebra crate merely to construct a reference.

Alternative rejected: blanket `pub mod` or re-exporting raw generic aliases. That would accidentally commit the public API to internal chart algebra and lifetimes before it provides a usable device workflow.

### 2. Make coordinates and declarations visible at the boundary

Keep separate conformal and frequency region inputs backed by existing rectangles. Present a neutral `native_coordinate` name in public evidence and corrections rather than exporting `z`/`conformal` names as if frequency values were dimensionless. Provenance identifies native kind and units; physical k0/beta and exterior kappa use cm^-1, frequency-chart derivatives are dF/dk0, and conformal derivatives are dF/dz. No new units system is introduced. Public settings use documented native absolute and dimensionless relative tolerances, mapping exactly to existing settings.

Expose coordinate-sheet/reference controls and per-exterior real radiation choices. Explicit complex references retain unspecified physical radiation; verified real anchors retain the existing conditional qualification. Characteristic/material analyticity, additional direct-frequency nonzero-q assumptions and known exclusions must be supplied explicitly. Convenience execution cannot insert a declaration. Expose continuation cutoff, integral/consistency settings and finite budgets without a new numerical convention.

Existing canonical conformal geometry and explicit finite plans are available. Direct frequency uses explicit rectangles only. Retain invalid/foreign/repeated-domain and foreign-discovery checks, stable domain identities through clones, original halo semantics and exact plan ordering. No inferred global coverage, periodic frequency surface or automatic threshold avoidance.

Alternative rejected: a single unlabelled complex rectangle and generic tolerance field. It makes otherwise valid-looking calls mix dimensionless conformal bounds with physical frequency values.

### 3. Offer staged execution and an owned convenience report

Public plans support the conceptual flow:

```text
physical problem -> validated domains -> finite plan
plan.discover(discovery_settings) -> owned discovery batch
plan.refine(&discovery, refinement_settings, matching_settings) -> inspectable refinement
plan.solve(discovery_settings, refinement_settings, matching_settings) -> owned solve report
```

A single-domain solve is the same finite-plan pipeline with one domain, avoiding a second acceptance algorithm. Optional descriptions are requested separately on the accepted observations/report; they are not a prerequisite for numerical success. The public staged result can borrow the supplied discovery where documented, but the convenience `solve` report must not borrow a local plan, discovery variable or problem.

Introduce an owned representation of report data using owned discovery, candidate outcomes and provenance, plus indices for chart/candidate/group relationships. Move data after existing computation, or share immutable discovery through Arc where appropriate; avoid cloning whole discovery trees per point. Factor an owned internal report payload out of the borrowed presentation if that is the simplest way to avoid double execution. Preserve settings, per-chart failures, optional refinement for successful charts, all candidate statuses and work/evidence, original estimates/leaves/summaries/terminations, enclosures, accepted branch signatures, groups and reconciliation diagnostics. No self-references, leaked allocations or unsafe lifetime extension. Run discovery/refinement once; do not rerun numerics to manufacture an owned report.

A public accepted point owns its accepted root/evidence/provenance, so it can be retained independently of an iterator/report view. It remains a point on an exact source domain, not a globally unique branch identity. Keep acceptance construction controlled: callers cannot manufacture an accepted root by filling public fields. An immutable source context can remain shared as it already is. Dropping execution objects does not imply fields or normalization were computed. Returning the report from a helper is a required external compile/runtime test.

Alternative rejected: return only a `Vec<ModePoint>` or accepted representatives. That discards coverage and would make empty, failed and incomplete solves look alike. Another rejected alternative is owning a discovery batch next to a refinement report borrowing it.

### 4. Keep failures and optional interpretation structured

Expose public construction, domain/configuration, execution, refinement-input and interpretation diagnostics. Adapt private wrappers at the boundary while preserving public failure category and source chains. Backend-specific original errors may be retained behind a public source object when necessary, provided their causes remain inspectable and known checked input/backend/output and exterior side/kind/work diagnostics have structured public access. No string-only normalization, swallowed source or fabricated success. Do not require arbitrary backend error types to implement new traits solely to format a public report.

An invalid configuration or foreign association fails the stage before optical work. A valid multi-chart execution retains per-chart discovery errors, and successful charts still run. Candidate errors remain non-accepted outcomes. Optional physical interpretation returns owned details or a public diagnostic attached to the same accepted point; neither changes acceptance or groups. Batch interpretation uses retained signature evidence under the existing semantics rather than silently recomputing roots.

Expose the current dielectric support matrix and component data. For spectral points, show temporal/exterior behaviour and radiation provenance; fixed-frequency paper labels and instantaneous complex-frequency flux remain explicitly unsupported. Keep arbitrary-sheet, verified-other and outgoing-continuation cases distinguishable. Do not name every temporally decaying root a QNM or supply normalization/mode volume fields with placeholder values.

Alternative rejected: one boxed error message for every stage, or interpreting a failure to classify as failure to find a mode. Both lose distinctions needed by downstream design/continuation code.

### 5. Prove usability from outside the crate

Add integration tests under lamina-mode/tests using only public library imports and existing core/material/unit APIs. Cover all three strategies, real and complex reference paths, native units, configurable controls, foreign ownership, returned owned reports/points, failed charts, continuation causes, unresolved clusters, overlap contributors/opposite branches and unsupported/unavailable interpretation. Compare representative public roots and physical derivatives with the existing independent fixtures, retaining the conformal chain rule. Internal tests still exercise deeper numerical fault injection; do not publicly export synthetic test solvers for that purpose.

Add runnable examples, with expected commands documented:

- `cargo run --example guided_mode --locked --offline`: analytical dielectric slab at fixed real k0, accepted propagation point and supported guided classification.
- `cargo run --example conformal_cavity --locked --offline`: constant-exterior normal-incidence cavity pole on explicit finite conformal bounds and selected radiation sheet.
- `cargo run --example dispersive_substrate --locked --offline`: actual analytical Lorentz exterior from lamina-materials, direct-frequency chart, real outgoing anchor, required assumptions and known guard regions.

Each executable checks an expected root/result and prints concise physical coordinates, numerical evidence and relevant coverage/radiation limitations. Reuse fixture mathematics, not private test modules or include-path tricks. Run the examples; all-target compilation alone is insufficient. Use a small valid search region to keep runtime practical without altering tolerances just to force success.

Replace the README's private-only opening with a public quickstart that is compiled by rustdoc or directly links to the tested example source. At least one complete external workflow must be a non-ignored doctest. Keep deeper private mathematical snippets labelled as illustrative when they remain ignored. Rustdoc should expose a concise public type graph and link units, assumptions, reference modes, errors and coverage. Correct the template package description; no package publication or unrelated dependency cleanup.

### 6. Preserve the roadmap and ownership boundaries

After this change, the next milestone is mode fields and local derivatives. Later work remains branch continuation, QNM normalization and derived observables, feature sensitivities, and objectives/optimization in lamina-design. Public reports retain evidence suitable for those later stages without implementing their operations now. Core continues to own reusable device/material/field primitives; mode owns modal search and interpretation; design consumes points/branches/features. Mode must not depend on design.

## Risks / Trade-offs

- [Facade grows into a second framework] -> Keep wrappers thin, reuse the one numerical pipeline and core primitives, and validate against existing independent regressions.
- [Ownership conversion loses evidence or duplicates expensive work] -> Use owned payloads and source indices; compare discovery/candidate/group counts and errors, and test returning reports after locals drop.
- [Private aliases leak through associated errors or signatures] -> Compile external clients and inspect generated rustdoc with private-interface lint checks; expose only the necessary public value/error types.
- [Convenient defaults are mistaken for physical assumptions] -> Numerical defaults are allowed, analytic declarations and radiation choices stay explicit; document units and unsupported interpretations.
- [Public API appears to promise global coverage or normalized QNMs] -> Keep coverage ledgers and conditional radiation evidence mandatory report data; examples display the limits.
- [Public report snapshots consume memory] -> Share discovery once and retain point-level evidence/provenance without copying entire trees. Serialization and large-sweep storage policy are later work.

## Migration Plan

Add the facade while keeping internal constructors and all 260 existing tests working. Adapt report ownership only where needed, then add external integration tests, examples and compiled quickstart. Run locked offline all-target tests, doctests, clippy, no-dependency rustdoc, changed-file formatting and all example executions. Review public signatures and every new spec scenario; run strict OpenSpec validation. No core changes are expected; any demonstrated missing primitive must be justified and kept scoped. There is no existing public solver API to migrate and no deployment or publication step. If the facade needs revision before adoption, internal numerical behavior and archived specifications remain the reference.

## Concrete facade representation

The implementation uses three named public problem wrappers with shared opaque
`Domain<B>` and `Plan<B>` handles. Their private variants retain the three distinct
internal strategies; cross-strategy associations are rejected before execution,
and provenance exposes `Strategy`. This avoids duplicating the public report and
execution surface while keeping raw chart and compiled evaluator types private.
The owned public report is `SolveReport`, with cloneable `ModePoint` observations,
shared `Discovery`, chart/candidate indices, native-coordinate evidence and
structured `Diagnostic` sources. Reports and points also retain the selected
backend's diagnostic type name, not a serialization identifier or retained field
workspace. Both built-in backends are exercised by public-only tests. Implementing
new external core backends is not added here: core's current backend trait uses
private response types, and expanding that contract is outside this facade change.

See [verification.md](verification.md) for final scenario coverage, exported API,
checks and the fields/local-derivatives handoff.
