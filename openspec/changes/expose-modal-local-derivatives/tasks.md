## 1. Core thickness differentiation

- [ ] 1.1 Expose a narrow selected-thickness directional compiler in lamina-core using the existing parameter selector and internal derivative mapping; verify finite-layer bounds and slot-zero seeding while preserving the constant compiler's behavior.
- [ ] 1.2 Document caller-facing thickness units and verify equivalent centimetre and nanometre stacks obey the derivative conversion, without adding an independent parameter or stack representation.

## 2. Physical characteristic partials

- [ ] 2.1 Add explicit local analyticity assumptions and validated derivative settings for accepted-point evaluation; verify a fixed-frequency source assumption alone cannot authorize spectral differentiation and document held-fixed conventions.
- [ ] 2.2 Evaluate physical spectral and in-plane partials with existing first-order jets and stored exterior branches, including electric and magnetic constitutive dispersion; verify opposite-root signs, regular-point finite differences and typed cutoff failures without numerical differentiation in production.
- [ ] 2.3 Evaluate selected finite-layer thickness partials with the core directional compiler while holding physical coordinates and exterior roots fixed; verify independent finite differences and unit conversions on finite-stack fixtures.
- [ ] 2.4 Return owned per-selector results with source association, explicit evaluation backend, units, assumptions and derivative evidence; verify one unavailable partial preserves successful partials and the original accepted point, and distinguish globally invalid settings from selector failures.

## 3. Guarded implicit slopes

- [ ] 3.1 Add spectral and thickness implicit slopes from the physical partials with finite dimensioned denominator thresholds and recorded numerator, denominator and threshold; verify regular ratios and unavailable results near a degenerate denominator preserve the gradients.
- [ ] 3.2 Document k0 units, slope meanings and limits, including the absence of energy-velocity, exceptional-point or ZGV certification; verify examples expose units and unavailable slopes without silently substituting values.

## 4. Independent derivative validation

- [ ] 4.1 Add native-to-physical chain-rule regressions for propagation, conformal spectral and direct frequency formulations at roots and regular off-root points; verify coordinate factors and constitutive derivatives using decreasing finite-difference steps within a consistent local branch neighborhood.
- [ ] 4.2 Extend guided TE/TM, both interface-plasmon cases, the two finite-film modes and dispersive fixtures with independent derivative and perturbed-root slope checks; verify each perturbed solve has its own valid problem/domain association and compare physical results only.
- [ ] 4.3 Add a runnable derivative-report example and public API documentation covering assumptions, backend provenance, partial failures and stored thickness units; verify the example runs against maintained fixtures and reports successful and unavailable quantities clearly.

## 5. Integration review

- [ ] 5.1 Run focused and full core checks for the compiler addition and locked lamina-mode checks, doctests and examples; review cross-crate changes and tidy temporary artifacts, recording validation results and existing warnings.

Follow-up: propose branch continuation, followed by QNM normalization and observables, feature and sensitivity analysis, and optimization integration in lamina-design. Sync and archive when requested.
