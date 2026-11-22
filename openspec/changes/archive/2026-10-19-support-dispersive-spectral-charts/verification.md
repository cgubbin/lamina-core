# Implementation review — 2026-10-19

Implementation is in sibling lamina-mode. Planning and this review remain in
lamina-core. Existing uncommitted work in both repositories was preserved.

## Scenario review

| Contract | Implementation and runnable evidence |
| --- | --- |
| Actual dispersive stack | `FrequencyProblem::frequency` compiles the supplied core analytical stack and retains its actual exterior handles; `dispersive_te_tm_slab_roots_and_derivatives_are_not_frozen` independently checks TE/TM roots, characteristic values and derivatives at/away from a root. The frozen-reference root differs by more than 1e-3 cm^-1. |
| Analytical material capabilities | The constructor requires core `AnalyticalMaterialStack<C>`; its handles require analytical values/derivatives. There is no real-sample conversion, fitted model, replacement exterior or constancy inference. |
| Fixed-reference branches | `dispersive_branches_cross_principal_cut_and_keep_both_material_derivatives` checks q=exp(k), shuffled/repeated samples, independent side signs, both epsilon and mu derivatives, and finite differences across a principal cut. |
| Domain declarations and exclusions | `frequency_domain_exclusions_cover_halo_and_never_claim_certification` checks native units, immutable reference identity, explicit exterior assumptions, requested/halo boundary exclusions, invalid bounds/settings/references and genuine outside-region failure. |
| Controlled fallibility | `exterior_budgets_depth_and_nonfinite_derivatives_remain_typed` and `sampled_cutoff_and_right_side_budget_preserve_identity_and_work` check per-side budget/work evidence, recursion exhaustion, derivative failure and sampled cutoff without a fallback strategy. |
| Reference versus radiation | `actual_anchor_radiation_is_separate_from_explicit_complex_branches` checks actual lossless dispersive anchor values, outgoing and incoming choices, unsupported lossy/magnetic anchors, cutoff and usable complex explicit branches with unspecified radiation. Zero frequency introduces no artificial coordinate branch point. |
| Local search and ownership | The physical TE/TM slab regression checks discovery/refinement, accepted native coordinates, exact source domains, foreign-problem rejection and supported spectral descriptions without a fixed-frequency named label. Shared pre-existing ownership/refinement regressions remain intact. |
| Multiple frequency charts | `frequency_overlap_keeps_native_evidence_and_finite_coverage` retains both overlap observations and their enclosures, reconciles them, keeps finite/excluded scope and compares single/batch descriptions. `frequency_plans_keep_opposite_sheets_failures_and_clusters` separates opposite exterior signatures and keeps successful/failed charts, original causes, terminations, leaves and multiplicities. |
| Later-stage failures | `frequency_refinement_and_description_retain_original_material_failures` injects model faults after discovery: candidate failures retain typed causes and discovery, and unavailable physical data retain the accepted root. Fault injection deliberately violates stable-model lifetime assumptions solely to test diagnostics. |
| Constant-exterior compatibility | `constant_te_tm_limit_matches_conformal_chain_rule_and_roots` checks both polarizations on matching sheets, values, dF/dk0=(dF/dz)/(dk0/dz), direct/conformal accepted physical roots and their distinct native coordinate units. |

## Checks

Commands run from lamina-mode with its Rust/Nix environment and locked offline
dependencies:

- `cargo test --all-targets --locked --offline`: 260 passing tests versus the
  250-test baseline (10 new frequency regressions).
- `cargo test --doc --locked --offline`: succeeds; one intentionally ignored
  private-API illustration, with no compiled public examples promised here.
- `cargo clippy --all-targets --locked --offline`: succeeds with six remaining
  pre-existing warnings, versus eight at baseline; no new warnings. Remaining
  categories are two partial-order comparisons, two complex types and two
  unit-struct defaults. The shared search signature warning persists with its
  new generic spelling; this is the same existing warning site.
- `cargo doc --no-deps --locked --offline`: succeeds; rendered guide content and
  anchors were inspected. Local README/guide file links and heading targets resolve.
- Changed-file rustfmt checks pass for 13 Rust files.
- `openspec validate support-dispersive-spectral-charts --strict --json`: succeeds.

No core implementation or dependency change was needed for this change, so core
checks were not rerun. No stack, layer, material, polarization or unit framework
was added. Frequency bounds wrap the existing rectangle machinery; all device
primitives come from core. No public exports, automatic chart strategy selection,
backend substitution, unrelated upstream edits or dependency changes were made.

## Limits and handoff

This implements local rectangular frequency charts, conditional on caller
knowledge of characteristic/material analyticity and holomorphic nonzero exterior
q on a neighbourhood of the closed evaluation region. Known exclusions must miss
that whole region and its halo. Empty exclusions and successful samples do not
certify it; even removable q zeros are conservatively excluded. There is no
contour-with-holes support, automatic threshold/pole discovery, global causal-sheet
coverage or cross-strategy reconciliation. Continuation is tolerance limited and
adds material-jet/quadrature cost; explicit finite work controls are retained.

Explicit complex references establish branch consistency only. Supported real
positive dielectric anchors additionally verify radiation at that anchor.
Spectral temporal/exterior behaviour and radiation provenance are available;
complex-frequency instantaneous flux and the paper's fixed-frequency named types
remain unsupported. Root acceptance, coverage, physical descriptions, later field
recovery and QNM normalization remain separate evidence/failure stages.

The public handoff is **point-solver API and executable examples for both spectral
strategies**, exposing native units, exact domain provenance, continuation controls,
assumptions, failures and coverage. Later milestones are fields/local derivatives,
branch continuation, QNM normalization/derived observables, features/sensitivities,
then optimization in lamina-design. Core owns reusable stack/material/field
primitives; mode owns search and modal interpretation; design owns objectives and
optimization. The public guide's roadmap records that order.
