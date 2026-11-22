# Validation and implementation notes

Core adds `ModalFieldBackend::unscaled_fields` for Scatter2 and Transfer2. It uses existing retained modal seeds and shared electromagnetic reconstruction, validates physical exterior data, and returns owned core fields plus fresh characteristic and seed-defect evidence. It never evaluates QNM normalization. A zero-layer interface demonstrates recovery when normalized mode creation is unavailable; a cavity verifies compatibility up to one common normalization scale.

Mode adds `ModePoint::recover_fields_with`, owned `ModeFields`, explicit gauge and residual settings, and a failure carrying the original accepted point and structured cause. No stack, field or sampling representation is duplicated. Both source and evaluation backend names are retained. The backend seed gauge is arbitrary; comparisons align all components with one complex scale.

## Necessary core correction

The physical regressions exposed cancellation in Scatter2's existing modal layer reconstruction, which used the driven incoming-wave cut formula at a resonance. It could give discontinuous guided/film profiles even with a tiny accepted characteristic. Homogeneous reconstruction now uses the right-outgoing seed and projective suffix networks, cancelling the global resonant denominator algebraically. Driven reconstruction and normalization semantics remain unchanged. All 39 existing scattering workspace tests pass, including off-root compatibility.

## Physical evidence

- Interface plasmon in both complex cases: component ratios, unsquared reference roots, tangential E/H, normal D/B and independent exterior exponentials.
- Guided TE and TM slabs: component zeros, interface continuity, shared-gauge agreement of both backends.
- Symmetric 40 nm Drude film: independent tanh/coth references, two distinct accepted modes, even/odd Hy parity and interface matching. Both boundaries are explicitly sampled; default core uniform sampling excludes the right boundary.
- Leaky spectral cavity: finite exterior growth is retained and overflow fails diagnostically.
- Invalid gauge, invalid sampling, strict residual rejection and post-acceptance nonfinite material injection preserve the accepted source point.

The runnable `modal_fields` example and portable `scripts/plot_modal_fields.py` generate five profiles in CSV and PNG/PDF/SVG, with full source/gauge evidence. The generated report is retained in `artifacts/modal-field-profiles` and was visually inspected. No QNM normalization, mode volumes, global discovery or branch continuation are claimed.

Core locked tests pass: 1475 library tests, 8 integration tests and one doctest (two existing ignored). Core clippy passes without warnings; rustdoc succeeds with two existing unresolved links (`TransferStabilityCheck` and `Parameter::InPlane`). Changed core files pass formatting and whitespace checks. Existing unrelated work was preserved; snapshots were used to review the actual core edits independently of earlier development.

The completed benchmark was synced and archived separately. Draft PRs: https://github.com/cgubbin/tmm-core/pull/1 and https://github.com/cgubbin/lamina-mode/pull/1. Current field work is on `codex/recover-modal-fields` in each repository and is excluded from those benchmark PRs.


Final mode verification: 261 library tests, 15 integration tests (276 total) and two doctests pass, with one existing ignored doctest. After boxing the source-preserving recovery failure, affected field and material-corruption tests pass again. Mode clippy retains the same six existing warnings and no new warnings. Mode rustdoc, package formatting and whitespace checks pass. The public release field example ran successfully for all five profiles. Film exterior decay was additionally checked independently on both sides.
