# Implementation validation

Implemented in sibling lamina-mode with no solver algorithm or core primitive additions:

- `examples/support/plasmon.rs`: shared model, independent references, public searches, physical checks and auditable output.
- `tests/plasmon_benchmarks.rs`: bounded physical regressions and wrong-sign controls.
- `examples/plasmon_validation.rs`: explicit quick/full command and fresh output directory.
- `scripts/plot_plasmon_validation.py` and optional requirements: portable plotting and complete-run checks.
- `docs/plasmon-validation.md` and README link: reproduction, conventions and interpretation limits.

Fresh full-run artifacts are retained at `artifacts/plasmon-regression-report` in lamina-core. The original `artifacts/plasmon-validation` is preserved.

## Results

- Locked offline tests: 261 library, 9 existing public integration and 2 new benchmark tests passed; doctests passed (one existing ignored example).
- Final quick regression rerun passed after adding machine-readable numerical settings.
- Full release harness: 33 accepted spatial and 33 temporal points, two accepted map roots, ten local TE controls with no accepted roots, and two 61×81 checked characteristic grids.
- Maximum relative root error: 1.31e-14 spatial, 1.98e-14 temporal. Maximum unsquared mismatch/kp: 5.22e-15 spatial, 1.17e-13 temporal. Acceptance limits are 1e-8, not snapshots of these observed errors.
- Scatter2/Transfer2 agree at roots and nearby regular points; all exterior decay and passive evolution signs pass.
- Forced-failure experiment used an isolated temporary copy with an intentionally impossible root tolerance. It returned nonzero, recorded failed completion status and retained all 66 identified point failures. Maintained tolerances were unchanged.
- Plotter rejected missing characteristic data, nonfinite coordinates, mismatched model metadata and failed completion status. It generated PNG/PDF/SVG outputs from the final full run; all three figure sets were visually inspected.
- `cargo clippy --locked --offline --all-targets --all-features` passed with the same six existing warnings (two partial-order comparisons, two complex types, two unit-struct defaults); no new benchmark warnings.
- All three existing public examples (`guided_mode`, `conformal_cavity`, `dispersive_substrate`) ran successfully.
- `cargo doc --locked --offline --no-deps`, package formatting and whitespace checks passed. The broader initial `cargo fmt --all` encountered existing trailing whitespace in core's `src/evaluate/real_axis/excitation.rs`; package-scoped formatting passes.

The plotting environment was supplied locally for verification; the maintained script contains no workstation package discovery or absolute dependency paths. Optional dependencies are documented normally. Neither numerical validation nor plotting changes unsupported metallic classification, radiation provenance, global coverage or QNM normalization claims.
