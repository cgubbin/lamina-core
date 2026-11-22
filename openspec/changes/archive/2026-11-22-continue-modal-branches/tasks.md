# Tasks

## 1. Source and benchmark fixture

- [x] 1.1 Transcribe the stopped-light planar-waveguide parameters and reference points from Adam Freddie Page's Imperial thesis and cited source material; record citations, units, conversions, material assumptions, and a checked fixture provenance file.
- [x] 1.2 Define a reproducible planar avoided-crossing fixture with independently seeded roots and reference expectations; document its stack, units, branch choices, and tolerances alongside the fixture.

## 2. Continuation request and trace model

- [x] 2.1 Add the source descriptor, validated ordered real-coordinate schedule, backend declaration, and bounded continuation settings; verify validation rejects empty, non-finite, repeated, reversed, and invalid-budget inputs.
- [x] 2.2 Add owned branch trace, accepted station, per-step outcome, and typed stop diagnostics; verify a failed trace retains its seed and every accepted prefix point.
- [x] 2.3 Document propagation and spectral source-coordinate conventions and the local-only scope; verify public API examples compile and state that branch traces do not certify global identity or ZGVs.

## 3. Local predictor-corrector

- [x] 3.1 Factor or add the checked seeded-correction path for a newly constructed physical problem, retaining fresh numerical evidence and source provenance; verify accepted corrections satisfy residual and correction limits.
- [x] 3.2 Implement slope prediction and bounded interval subdivision with explicit backend and exterior-branch continuity; verify inserted stations remain ordered and no step invokes global discovery or silently changes backend/strategy.
- [x] 3.3 Add regular-step, derivative-degenerate, nearby-root, branch-cutoff, correction-failure, and retry-exhaustion coverage; verify each yields the specified accepted point or typed partial-trace outcome.

## 4. Physical branch reports

- [x] 4.1 Extend the maintained plasmon validation to emit both continued physical slices and compare each numerical branch with its independent reference; verify sign, decay, residual, backend, and tolerance fields are present in machine-readable output.
- [x] 4.2 Add the planar avoided-crossing report with two seeded traces; verify both traces retain their accepted evidence through the gap and the plot does not infer identity from coordinate ordering.
- [x] 4.3 Add the thesis stopped-light branch report and compare measured local slopes near both reference points; verify report metadata identifies source values, units, model conversions, and limitations.

## 5. Field datasets and visual report

- [x] 5.1 Recover transverse electric and magnetic profiles from exact accepted points before, near, and after the avoided crossing and around both stopped-light points; verify each dataset records physical coordinates, exterior branches, field gauge, and explicit QNM-normalization status.
- [x] 5.2 Emit versioned branch/field data and plots from a documented full-report command; verify the plotter rejects incomplete or stale inputs and recreates dispersion and field figures from a newly generated report directory.
- [x] 5.3 Document how to run the physical report, inspect its citations, machine-readable samples, dispersion plots, and field plots; verify all documented commands and output paths are reproducible outside a workstation-specific directory.

## 6. Integration

- [x] 6.1 Run the documented quick numerical report and focused continuation/physical benchmark checks, then run the full report and inspect all generated dispersion and field plots; verify failures identify the case and stop outcome rather than producing a success report.
- [x] 6.2 Review public API and report artifacts against both delta specs and existing point-search behavior; verify existing plasmon validation remains compatible and all supported benchmark cases are represented.
