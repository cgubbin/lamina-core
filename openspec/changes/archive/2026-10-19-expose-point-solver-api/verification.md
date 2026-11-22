# Public point API review — 2026-10-19

Planning remains in lamina-core; implementation is in sibling lamina-mode.
Pre-existing uncommitted work in both repositories is preserved.

## Exported surface

- Problems: `PropagationProblem`, `ConformalSpectralProblem`, `FrequencyProblem`;
  Scatter2 defaults and explicit supported core backend construction.
- Execution: opaque strategy-checked `Domain<B>` and `Plan<B>`, explicit finite
  plans, existing canonical conformal builders, checked `evaluate`, staged
  `discover`/`refine` and owned convenience `solve`.
- Ownership/evidence: `Discovery`, `ChartDiscovery`, `Estimate`, `SolveReport`,
  `ModePoint`, `Provenance`, `Strategy`, `NumericalEvidence`, `ChartOutcomes`,
  `CandidateOutcome` and `CandidateStatus`. Discovery is shared once via Arc;
  accepted points retain exact root evidence and source context independently.
- Controls: distinct conformal/frequency bounds, required analytic declarations,
  constant-exterior evidence, explicit and verified radiation references,
  exclusions, discovery/continuation settings and public native-unit refinement
  and matching settings. Core device/material/polarization/backend/unit primitives
  remain the inputs. `Complex` is a convenience alias for complex f64.
- Diagnostics/interpretation: `Diagnostic`, `Stage`, `OpticalFailure`, selected
  public cause types, source traversal and foreign-association identification;
  separately fallible point/report descriptions with the existing physical
  support matrix and applied tolerances.

The raw analytic chart trait, raw maps, compiled solver and seeded jet parameters
are not required public inputs. There is no second device model or root algorithm.
Owned reports use original chart/candidate/group indices, retain each accepted
observation and all outcomes, and do not rerun numerics to obtain ownership.

## Scenario coverage

| Requirement | Evidence |
| --- | --- |
| Physical construction | `external_construction_all_strategies_and_backends` constructs all three strategies with core primitives, checks built-in backend agreement and actual reference provenance. The conformal constructor requires constant descriptors rather than arbitrary exterior handles. |
| Domain/reference controls | `external_domain_units_derivatives_and_ownership` compares native units and dF/dz / k0 against dF/dk0; frequency error tests cover exclusions, invalid continuation settings and side-specific work. A compile-fail doctest demonstrates the mandatory declaration parameter. |
| Finite execution | Public tests construct explicit and both canonical conformal plans; empty/repeated/foreign/mixed strategy associations and foreign discovery are rejected, and unsupported normal-incidence canonical geometry requires explicit bounds. Domain clone identity and requested scope are retained. |
| Owned reports/points | `owned_reports_points_and_overlap_survive_execution_locals` returns a report from a helper after local execution objects drop, checks discovery/candidate/group evidence and retains a point after dropping the report. Backend/strategy provenance is inspectable. |
| Evidence-preserving inspection | `frequency_owned_coverage_keeps_overlapping_failure_and_all_outcomes` keeps overlapping success and continuation failure with original causes and coverage. `canonical_plans_and_budget_clusters_keep_public_diagnostics` retains cluster multiplicity, leaves, budget termination and non-accepted outcomes. Opposite-branch observations remain separate in a public frequency plan. Existing internal same-coordinate sheet-reconciliation regressions remain in the baseline. |
| Structured diagnostics | Public tests retain actual Transfer2 backend errors and checked output causes; continuation diagnostics preserve side/kind/work in discovery. The adapter unit regression checks all three checked input/backend/output categories and original sources, without exporting a synthetic solver. |
| Optional interpretation | `accepted_guided_point_survives_unavailable_interpretation_and_checked_failure` checks supported guided classification, injects a model fault after acceptance and retains the accepted point, groups and discovery alongside unavailable interpretation and later refinement failure. This deliberate mutable model is a diagnostic fixture; stable behavior remains a real caller obligation. Spectral tests retain unspecified/outgoing provenance and unsupported fixed-frequency labels. |
| Executable usage | `guided_mode`, `conformal_cavity` and `dispersive_substrate` use only public imports and real core/material primitives. Each checks an independent expected result; the dispersive example also rejects the frozen-exterior comparison. A complete public quickstart compiles and executes as a non-ignored doctest. |

## Verification commands

Run from lamina-mode with its Rust/Nix environment and locked offline dependencies:

- `cargo test --all-targets --locked --offline`
- `cargo test --doc --locked --offline`
- `cargo clippy --all-targets --locked --offline`
- `RUSTDOCFLAGS='-D rustdoc::broken_intra_doc_links' cargo doc --no-deps --locked --offline`
- Changed-file rustfmt, local README/guide/source links and heading-target checks.
- `cargo run --example guided_mode --locked --offline`
- `cargo run --example conformal_cavity --locked --offline`
- `cargo run --example dispersive_substrate --locked --offline`
- `openspec validate expose-point-solver-api --strict --json`

Final results: **270 tests pass** (261 library tests and nine external integration
tests), preserving all 260 baseline tests and adding one diagnostic-adapter test
plus nine external API regressions. All three example targets compile and their
executables run successfully. The guided mode is beta=3.5 cm^-1; the constant
cavity pole is k0≈5.235987755983−2.537102031436i cm^-1; the Lorentz substrate pole
is k0≈5.222211243166−2.711901290352i cm^-1, agreeing with independent checks.

Doctests pass: one complete public workflow executes and one missing-declaration
example fails compilation as intended; one pre-existing private illustration
remains intentionally ignored. Clippy succeeds with the same six existing warnings
(two partial-order comparisons, two complex types and two unit-struct defaults),
with none added. Strict-link rustdoc succeeds and all public facade types render.
Formatting passes for 14 changed Rust files; local links and heading targets
resolve. Strict OpenSpec validation passes. No core implementation changes were
needed, so core tests were not rerun. Cargo.toml changed only its template package
description; no dependency or lockfile change was introduced.

## Limits and next milestone

The public interface exposes the existing local, conditional numerical contracts.
It does not certify unknown poles/thresholds, extend acceptance into halos, recover
omitted sheets or excluded scope, split clusters, reconcile mixed strategies or
establish a globally unique mode count. Explicit complex references do not acquire
physical outgoing qualification. Spectral named fixed-frequency labels and
instantaneous complex-frequency flux remain unsupported. No fields, QNM
normalization, mode volumes, sweeps, feature sensitivities or optimization are
computed by this change.

Reports and points own/share evidence and source context, not retained optical
field workspaces. Backend type names are diagnostics rather than a serialization
contract. Built-in core backends are supported; creating new external core backend
implementations would require a separate review of core's currently private
response-type boundary. No core implementation, new dependency, unrelated upstream
repair, CLI or publication step was required here.

The next roadmap milestone is **mode fields and local derivatives**, with fallible
field reconstruction separate from root acceptance. Then come adaptive branch
continuation, QNM normalization and derived observables, features/sensitivities and
optimization in lamina-design. Core owns reusable device/material/field primitives;
mode owns modal search and interpretation; design owns objectives and algorithms.
The public and mathematical guides retain that order.
