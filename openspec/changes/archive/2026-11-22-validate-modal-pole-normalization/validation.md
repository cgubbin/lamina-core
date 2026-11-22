# Implementation and validation

## One-dimensional normalization clarification

The user's clarification was correct for the normal-incidence, nondispersive-exterior Maxwell bilinear form: outgoing electric and magnetic exterior densities cancel, so no exterior integral needs to be evaluated. It was incorrect to infer that core's finite-layer Maxwell normalization was necessarily incomplete solely because it had no explicit exterior term.

Independent TE/TM analytic cavity regressions now establish the precise convention mapping. For the symmetric vacuum / n=2, d=.3 cm / vacuum cavity with unit right-outgoing primary amplitude, the finite-layer Maxwell overlap is M=.9 cm. The scalar unit-delta Green pole factors are N_k0=+k0 M for TE and N_k0=-k0 M for TM. Adding scalar endpoint contributions to M again fails the analytic reference. Existing normalized core numerical behavior is preserved.

The new complete capability is explicitly the scalar pole-residue convention. It uses analytic outgoing endpoint derivatives, never infinite exterior quadrature, and does not assume the same cancellation for dispersive exteriors or beta poles. Normalized fields therefore have the declared scalar-pole scaling, not an implicitly relabelled legacy Maxwell M=1 scaling.

## Implementation

- Core: `ModalNormalizationBackend`, physical `PoleVariable`, validated settings, analytic finite-layer state-product integration, exact outgoing endpoint terms and finite exterior partition shifts. Reports retain contributions, units, conditioning, thresholds, fresh characteristic, seed defect and optional legacy comparison evidence. No new stack or sampling primitives or runtime dependencies.
- Mode: separately fallible `ModePoint::normalize_with`, owned normalization reports and normalized E/H sampling. Both source and evaluation backends, exact branches and explicit local declarations are retained. Incompatible variables, singular/ill-conditioned factors and invalid settings preserve the accepted point. Sampling detects tested post-evaluation material corruption.
- At the user's request normalization proceeds before the general local-gradient proposal. It uses existing analytic core material derivatives; the broader derivative API remains queued rather than being duplicated here.
- Public guide: sibling `lamina-mode/docs/pole-normalization.md`. Runnable `pole_normalization` example supports quick/full runs, fresh directories and failed-case records. Portable plotting reads the recorded data and rejects incomplete or failed runs.

## Independent physical evidence

Reference Green functions use value-only homogeneous-layer sine/cosine propagation and weighted Wronskians with a tested positive unit flux jump. They do not use production normalization, its characteristic derivatives or normalized fields. Checks cover source/observation pairs, decreasing offsets in two complex directions and independently refined contour residues.

Twelve fixtures cover guided TE/TM propagation, direct and conformal spectral cavities in both polarizations, both complex interface-plasmon cases with zero finite layers, both coupled metal-film modes and electric/magnetic dispersion in finite layers and exteriors.

Full-run measured maxima:

| Check | Maximum relative discrepancy |
| --- | ---: |
| Normalized prediction versus independent contour residue | 3.53e-13 |
| Contour refinement | 3.06e-13 |
| Scatter2 versus Transfer2 residue | 9.20e-15 |
| Complex seed-gauge invariance | 2.52e-15 |
| Exterior partition invariance | 5.32e-16 |

Acceptance is 1e-7 for independent residues and 1e-8 for contour refinement and invariance checks, not snapshots of observed floating-point agreement. Seven decreasing offsets show convergence followed by roundoff/root-location error. Deliberate controls omitting scalar boundary terms or dispersion, conjugating a bilinear residue factor or using the wrong variable are detected. Unnormalizable controls are labelled unavailable, without fabricated error values.

The full numerical report and inspected PNG/PDF/SVG plots are retained at `artifacts/pole-normalization-validation`. Plotting rejects missing contributions, failed completion status, nonfinite data and mismatched settings. Both quick and full public commands ran successfully; plotting is optional and absent from Rust tests.

## Repository checks

- Core locked offline suite: 1478 library and 8 integration tests pass (1486 total), plus one doctest; two existing doctests remain ignored.
- Mode locked offline suite: 261 library and 18 integration tests pass (279 total), plus three doctests; one existing doctest remains ignored.
- Core clippy passes without warnings. Mode clippy retains the six existing warnings, with none introduced by normalization.
- Rustdoc succeeds. Core retains its two existing unresolved links (`TransferStabilityCheck` and `Parameter::InPlane`); mode rustdoc has no new warnings.
- Changed core files and the mode package pass formatting and whitespace checks. Independent reference tests pass again after consolidating shared dev modules.

No mode volume, Purcell factor, exceptional-point certificate, global discovery or causal radiation promotion is claimed. Existing arbitrary-amplitude field recovery and legacy normalization remain available independently.
