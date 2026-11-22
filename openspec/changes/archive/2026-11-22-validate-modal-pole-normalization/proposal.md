# Proposal

## Why

Recovered modal fields currently have arbitrary amplitude. The older core normalization path aggregates finite-layer bilinear contributions, and its existing compatibility tests do not establish complete open-system normalization. We need a normalization whose physical convention and exterior treatment are explicit and whose predicted pole residues agree with independent solutions before using it for observables.

## What Changes

- Prove the normal-incidence, nondispersive-exterior cancellation in the existing Maxwell bilinear form; no explicit exterior integral or second counting of the equivalent scalar endpoint term is introduced.

- Add separately fallible complete normalization for simple frequency poles at fixed beta and simple propagation poles at fixed k0; distinguish the two in the public result rather than applying one spectral formula to both.
- Include finite-layer contributions and analytically regularized exterior contributions using the accepted branches and actual constitutive dispersion. State reciprocal pairing, core impedance conventions, canonical units and the planar Fourier/per-unit-area convention.
- Return owned normalization evidence and normalized field sampling, preserving the original accepted point and arbitrary-amplitude recovery. Record the remaining square-root sign freedom.
- Validate against an analytic planar cavity, independent scalar Green-function pole residues, complex-gauge invariance, exterior-boundary displacement and backend agreement. Include the zero-layer interface plasmon and both coupled film modes in the relevant complex cases.
- Document the existing finite-layer normalization path as a legacy convention, retain its numerical behavior, and expose the complete normalization additively. Do not infer outgoing radiation qualification from a successful normalization on a caller-selected sheet.

Mode volumes, Purcell factors, driven-response expansion APIs, branch continuation, exceptional-point normalization and arbitrary nonreciprocal/anisotropic systems are outside this change. The Green-function implementation is an independent validation fixture, not a new general public response solver.

## Capabilities

### New Capabilities

- `modal-pole-normalization`: Complete, convention-specific simple-pole normalization and normalized fields, with exterior regularization, source association and independently validated residue predictions.

### Modified Capabilities

None. Numerical acceptance, radiation provenance and arbitrary-amplitude field recovery keep their contracts; the legacy normalized core entry point keeps its numerical semantics.

## Impact

Lamina-core owns reusable bilinear integration, exterior terms and the high-level backend capability. Sibling lamina-mode owns accepted-point association, explicit pole-variable selection, normalized field results and diagnostics. Reuse core stack, material, field, sampling and unit primitives. Public documentation must distinguish legacy finite-layer scaling from complete normalization. No new runtime dependencies or parallel device representation are planned.

## Sequence

Apply after `recover-modal-fields`. At the user's request this normalization work proceeds before the general local-derivative API; use existing core analytic material derivatives and retain explicit local declarations, leaving that broader API queued. This provides a validation milestone before mode volumes or excitation strengths; branch continuation can proceed independently afterwards. Longer-term work remains branches, derived observables, feature sensitivities and optimization in lamina-design.
