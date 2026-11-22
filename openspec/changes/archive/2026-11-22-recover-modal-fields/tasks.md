## 1. Core unscaled field recovery

- [x] 1.1 Add a normalization-independent retained modal solution operation in lamina-core, sharing existing boundary-wave reconstruction and exposing a usable public capability for supported backends; verify recovery succeeds when QNM normalization is unavailable and document its arbitrary amplitude.
- [x] 1.2 Preserve existing normalized mode construction and seed-normalization behavior; add focused compatibility checks comparing existing normalized results and the new unscaled solution under one common complex scale.

## 2. Accepted-point field adapter

- [x] 2.1 Add explicit-backend field recovery from an owned accepted ModePoint using its exact physical coordinates, ordered exterior roots and stack association; verify no discovery or refinement is rerun and record source and recovery backend provenance.
- [x] 2.2 Return owned fields through existing core sampling, position and field containers; verify results outlive the temporary retained evaluator and preserve sampling order and opposite interface sides without introducing new stack primitives.
- [x] 2.3 Add the core seed gauge and finite nonzero complex rescaling settings, applying the same scale to electric and magnetic fields; verify component ratios and residual evidence remain consistent and document the gauge contract.
- [x] 2.4 Add explicit recovery residual settings and evidence for approximate accepted roots, checking detectable inconsistent or nonfinite retained data; verify excessive boundary defects fail independently of the original point's acceptance.
- [x] 2.5 Separate unsupported recovery, invalid sampling and overflow failures, preserving the accepted point; verify finite leaky exterior growth is retained and unsupported named classification does not block field recovery.

## 3. Physical validation and examples

- [x] 3.1 Add guided TE/TM and analytic interface-plasmon field comparisons in both complex-coordinate cases; verify core component conventions, tangential E/H and normal D/B continuity, and agreement between supported backends after a single gauge alignment.
- [x] 3.2 Add a symmetric finite metal-film fixture with both coupled TM modes, checking independent tanh/coth dispersion conditions and even/odd magnetic-field parity; verify interface continuity and decay for each mode and document the parity convention.
- [x] 3.3 Add a leaky-mode fixture and runnable field-profile examples using the maintained benchmark output conventions; verify finite exterior growth and inspect plots of interface and film modes, clearly identifying arbitrary amplitude and absent QNM normalization.

## 4. Integration review

- [x] 4.1 Run targeted and full core checks for the additive recovery capability and locked lamina-mode checks, doctests and examples; review cross-crate compatibility and tidy temporary changes, recording remaining limitations and existing warnings.

Follow-up: apply `expose-modal-local-derivatives`, then plan branch continuation. QNM normalization and observables remain a later change. Sync and archive when requested.
