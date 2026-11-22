# Design

## Context

See proposal.md for motivation. The public solver change ended with 270 passing tests and six existing clippy warnings. `artifacts/plasmon-validation` in lamina-core now contains a standalone Rust harness, CSV/JSON data, Python plotting, three figure sets and a local HTML report. It uses public lamina-mode problems and core/material/unit primitives; no library solver was changed for the plots.

The measured reference case is vacuum / Drude metal with epsilon_infinity=1, kp=1e5 cm^-1 and gamma/omega_p=0.02. Both sets of 33 numerical roots agree with independent equations to approximately 2e-14 relative error in that run. This is observed agreement, not a portable acceptance guarantee. Analytic answers guide local chart placement. The original plotting helper discovers Nix store packages using workstation paths; that is not a maintained portable plotting interface. The original CSV also records zero TE accepted count at unchecked points, with checked indices only recoverable from source; the maintained format must distinguish those cases.

## Goals / Non-Goals

**Goals:** One maintained numerical benchmark implementation shared by routine integration tests and full public examples, plus portable plotting from fresh recorded data. Physical assertions supplement the characteristic residual.

**Non-Goals:** No benchmark-specific solver shortcut, blind global discovery, production branch tracking, general metal classification, material fitting, field/QNM implementation or plotting dependency in the Rust library.

## Decisions

### 1. Promote the public harness without changing numerics

Place shared fixture/reference helpers in lamina-mode's validation/example support area, with integration tests importing maintained support code through an explicit dev-only module arrangement. Ordinary callers and examples use public imports; no private solver test files are included. Keep the original core artifact directory as historical evidence rather than overwriting its measured values.

Provide documented quick/full modes on a runnable example, for example `cargo run --release --example plasmon_validation -- --suite quick|full --output <dir>`. Standard-library argument handling is sufficient; no runtime dependency or general-purpose solver CLI is needed. The quick automated suite uses three to five well-separated regular points per complex case; the full mode retains the 33+33 curves and two characteristic maps. Plotting is a separate optional Python command with ordinary NumPy/Matplotlib dependency instructions. Do not scan absolute Nix store paths in the portable script. The Rust suite must run without Python.

Alternative rejected: snapshotting PNGs as numerical tests. Image rendering is neither a root-accuracy check nor robust across plotting versions.

### 2. Keep the independent physical equations and branches explicit

Use u=k0/kp=omega/omega_p and v=beta/kp, with epsilon_m=1-1/[u(u+i eta)] and eta=0.02. The independent spatial root is v=u sqrt[epsilon_m/(1+epsilon_m)] with passive forward sign. Temporal roots satisfy

```text
u^4 + i eta u^3 - (1+2v^2)u^2 - 2i eta v^2 u + v^2 = 0.
```

Select the temporal branch from the positive lossless bound root; check the unsquared TM equation kappa_d/epsilon_d+kappa_m/epsilon_m=0 so squaring cannot hide a wrong sheet. Require Im kappa_j>0 on both sides; require Im beta>0 for spatial attenuation and Im omega<0 for temporal decay under the existing exp(+i beta z), exp(-i omega t) convention.

Reference formulas calculate benchmark windows and explicit branch anchors, not the values returned as numerical observations. Each marker must come from public discovery plus acceptance, with original termination and candidate evidence available. Local frequency rectangles stay away from Drude poles and both exterior q zeros; validate the fixed fixture geometry and retain the caller declarations. Do not claim finite samples certify arbitrary domains.

Use documented normalized root/boundary tolerances initially of order 1e-8 on regular quick fixtures, rather than asserting the previous machine-level observed agreement. Record the exact chosen limits and measured errors. Compare Scatter2/Transfer2 values at accepted and at selected nearby regular points so backend checks are not based solely on two nearly zero values. TE controls use the same local windows on a declared subset; record `te_checked` and optional outcome, never zero-as-untested.

### 3. Generate data and plots as auditable outputs

CSV rows identify complex case, prescribed coordinate, numerical root, independent root, residual, unsquared mismatch, ordered exterior wavevectors, backend checks and actual TE control status. JSON records model constants, units, chart/reference construction, numerical tolerances/budgets, per-case success/failure and summary errors. Preserve failed cases as diagnostics and return a failing validation exit status; do not silently shorten curves.

Plot real/imaginary dispersion separately, plus characteristic landscapes from public checked evaluation with real/imaginary zero contours and numerical/analytic root markers. A physical-check figure displays root error and exterior decay. State that real-frequency and real-wavevector curves are different slices at finite loss, that the metal is an illustrative model, and that the named classifier/radiation qualification limitations remain unchanged. Output PNG plus PDF or SVG with a local HTML index. The plot command consumes a named complete run, verifies required cases/fields and rejects missing/non-finite inputs rather than falling back to old CSVs. Do not automatically open a browser or require a display in tests.

Alternative rejected: silently regenerate missing numerical data from analytic curves inside the plotting script. Analytical lines and measured markers must remain distinguishable.

## Risks / Trade-offs

- [Analytic chart placement is mistaken for global discovery] -> Label local benchmark placement in reports and separate curves from future branch continuation.
- [Platform rounding makes a machine-precision snapshot brittle] -> Assert documented physical/relative tolerances and retain measured errors rather than exact CSV bytes.
- [Quick tests become slow plotting jobs] -> Fixed small quick fixtures; full grids and Python remain opt-in.
- [Stale artifacts conceal a failed run] -> Explicit run directory/metadata and failure records; plotting validates the requested run.

## Migration Plan

First change in the sequence. Promote the harness and add quick tests, portable plotting and documentation, then rerun the full report and compare with the historical artifact. Run current locked tests, doctests, clippy, rustdoc and formatting; execute quick/full generation and plotting explicitly. No core change is expected. `recover-modal-fields` then reuses the physical fixture conventions; `expose-modal-local-derivatives` reuses the independently calculated roots. Later milestones remain branches, QNM normalization, features/sensitivities and optimization in lamina-design.
