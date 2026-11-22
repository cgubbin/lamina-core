# Design

## Context

See `proposal.md` for motivation and `specs/` for the observable contract. The existing `lamina-mode` API returns owned accepted `ModePoint`s with per-point problem provenance and local implicit derivatives. Its refinement path performs checked local correction, while physical problem construction binds a fixed real control coordinate and source strategy. Current plasmon and field examples generate reproducible point data, but there is no reusable branch trace or maintained physical report for avoided crossings and stopped-light points.

## Goals / Non-Goals

**Goals:**

- Add a mode-level continuation abstraction that can reconstruct each next physical problem from the original source description and requested real control coordinate.
- Keep continuation state independent of discovery state and preserve valid partial traces when a step stops.
- Deliver branch samples and field observations as machine-readable report inputs, with plots generated from those inputs.

**Non-Goals:**

- Modify `lamina-core` stack primitives or add a competing device/stack model.
- Certify completeness, global mode identity, exceptional points, or ZGVs.
- Normalize every field profile automatically when the selected mode or boundary conditions do not support the existing QNM normalization contract.
- Optimize device parameters in this change.

## Decisions

1. **Keep continuation in `lamina-mode`, using a stored source descriptor.** A source descriptor records the physical stack, polarization, source strategy, fixed/varied coordinate convention, and explicit backend. For each station, it constructs the corresponding problem using existing lamina-core stack and chart types. This avoids retaining a `Plan` whose domain identity is intentionally fixed, and avoids adding stack primitives upstream. Rebuilding the physical problem is distinct from repeating global root discovery.

2. **Predict with the physical implicit slope and correct locally.** The predictor advances the solved complex coordinate from the previous accepted point's local derivative. The corrector evaluates/refines in a bounded neighborhood of that prediction and produces a fresh accepted `ModePoint` with the new problem's provenance. Reusing the ordinary global `Plan::refine` is not sufficient because it binds a discovery to one plan/domain identity; instead factor or expose a narrow internal checked-seed correction path that accepts the newly constructed problem and predicted root. Do not silently fall back to unrestricted global discovery.

3. **Represent adaptivity in the trace.** The caller supplies ordered target stations and bounded retry/minimum-step settings. On a failed interval, the controller bisects toward the target and records any accepted inserted stations in order. A typed stop outcome records the requested target, last accepted point, attempted step, and cause. The result owns the seed, accepted points, and outcomes even when incomplete.

4. **Preserve exterior branch choices at each physical station.** The source descriptor must retain explicit anchors or the branch-selection policy used by the seed. At each new physical frequency, reevaluate material response and carry the corresponding branch choices through correction; do not replace them with principal square roots after a failure. If the source cannot provide branch-consistent exteriors, return a typed unsupported/branch failure.

5. **Separate numerical traces from the physical report pipeline.** Keep existing point-search plasmon validation intact. Add fixtures and a report runner that produce versioned machine-readable branch samples and field samples; the plotter consumes only the freshly generated report directory and rejects incomplete data. Include the maintained interface plasmon in both complex cases, a two-root planar avoided crossing, and the planar stopped-light structure documented by Adam Freddie Page's Imperial thesis. The thesis configuration and reference values must be extracted and checked against the source before they become fixture constants; the report records citation, units, parameter provenance, and any differences between thesis and implementation.

6. **Make field gauge and normalization explicit per sample.** Recover fields from the exact accepted point and record physical coordinate, exterior branches, gauge, and normalization status. Use the existing QNM normalization only where its assumptions are met; otherwise report the declared amplitude gauge and mark QNM normalization unavailable. Select samples before, near, and after the avoided crossing and around each reported stopped-light point.

7. **Call ZGV evidence a measured local slope.** Report the finite/local derivative and its convergence or comparison to a source value where available. A small derivative is a diagnostic observation, not a solver certificate or a globally optimized ZGV.

## Risks / Trade-offs

- [The local corrector may converge to a nearby root near an avoided crossing] → retain both independently seeded traces, record correction evidence, and avoid assigning global identity from coordinate ordering alone.
- [Near a ZGV the predictor slope can become singular or poorly conditioned] → use typed slope outcomes, bounded step reduction, and explicit stop diagnostics; do not divide through a degenerate slope.
- [The thesis fixture may not map exactly to available material models or cited parameter conventions] → record the source and conversion trail, compare against published reference points, and expose unsupported material assumptions instead of silently fitting them.
- [Field amplitudes depend on gauge and QNM normalization conventions] → attach gauge and normalization status to every field dataset and keep raw/recoverable data separate from normalized profiles.
- [Continuation's local results can be mistaken for exhaustive branch discovery] → make partial/completion state and local-scope limits visible in API documentation and every report.

## Migration Plan

No migration is required. Add the continuation API and report alongside existing point-search APIs and examples. Existing point and plasmon report schemas remain readable; new branch reports use a separately versioned schema. Rollback consists of removing the new module/example/report assets without changing existing solver behavior.

## Open Questions

None that block implementation. Exact thesis fixture parameter transcription is a required implementation investigation and verification task, not a design choice to defer.
