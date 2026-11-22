## 1. Independent numerical fixtures

- [x] 1.1 Promote the vacuum/Drude interface model and independent spatial and temporal reference formulas into maintained lamina-mode benchmark support; document units and branch conventions, and verify representative reference values without invoking the solver.
- [x] 1.2 Add bounded public-API regression cases for both complex-coordinate formulations, checking root error, unsquared TM mismatch, exterior decay and attenuation signs; verify deliberate sign errors fail and unsupported classification remains explicit.
- [x] 1.3 Add TE negative controls and Scatter2/Transfer2 comparisons at roots and regular off-root points; verify reports distinguish checked controls from untested points and compare consistently scaled characteristics.

## 2. Reproducible validation artifacts

- [x] 2.1 Add a maintained public example with quick and full suites, explicit output directory and machine-readable settings, roots, reference errors and failures; document invocation and verify the full suite retains unsuccessful cases rather than silently dropping them.
- [x] 2.2 Promote the plotting/report helper with ordinary optional plotting dependencies, fresh-run completeness checks and PNG plus vector outputs; document setup and verify missing or mismatched inputs fail instead of reusing stale data.
- [x] 2.3 Document the benchmark's physical evidence and limits, including locally guided discovery and unsupported metallic labels; regenerate the report from a fresh full run and inspect axes, units, legends and recorded diagnostics.

## 3. Integration review

- [x] 3.1 Run the bounded regressions, locked lamina-mode checks and the documented full artifact workflow; review the diff for unrelated changes and record results and any existing warnings.

Follow-up: apply `recover-modal-fields` next, then `expose-modal-local-derivatives`. Sync and archive this change when requested; retain the longer-term sequence in the proposal and design.
