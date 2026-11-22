# Proposal

## Why

The new plasmon plots demonstrate physically correct local roots in both complex cases, but the checks currently live in an ad hoc artifact crate. A maintained public validation harness is needed to catch branch, attenuation and constitutive-dispersion regressions before field and derivative work builds on these results.

## What Changes

- Promote the vacuum/Drude-metal interface benchmark into lamina-mode's maintained public examples and integration tests, preserving the existing source/data report as a reference artifact.
- Add a quick numerical regression suite and a separately requested full plotting suite for real-frequency/complex-wavevector and real-wavevector/complex-frequency modes.
- Check independent unsquared TM boundary matching, complex roots, exterior decay, passive spatial/temporal attenuation, sampled TE negative controls and supported backend agreement.
- Produce portable CSV/JSON data, characteristic maps and publication-friendly figures with model parameters, tolerances, reference assumptions and local-search scope recorded.
- Document that analytic solutions locate benchmark charts and that metallic named classification and built-in radiation qualification remain unchanged.

No root algorithm, global discovery, production branch tracking, fitted-metal database, fields, normalization or classification extension is included.

## Capabilities

### New Capabilities

- `modal-physical-benchmarks`: Reproducible public physical validation runs and diagnostic plots with independent reference equations and explicit scope.

### Modified Capabilities

None. Existing solver and classification contracts are exercised without changing their behavior.

## Impact

Primarily sibling lamina-mode: public example/harness, focused tests, optional Python plotting and guide links. OpenSpec remains in lamina-core; artifacts/plasmon-validation supplies the source benchmark. Numerical tests need no plotting packages; plotting dependencies stay development-only. No core primitive or runtime dependency changes are expected.

## Sequence

First of three planned changes. `recover-modal-fields` reuses these fixtures; `expose-modal-local-derivatives` later validates derivative predictions against them. Branch continuation, QNM normalization, feature sensitivities and lamina-design optimization remain later milestones.
