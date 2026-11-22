# Tasks

## 1. Explicit radiation references

- [x] 1.1 Add validated per-side radiation choices, real-anchor inputs, tolerances and opaque verification evidence in lamina-mode; implement the opt-in problem-owned domain constructor using the existing coordinate-reference path. Verify real-coordinate consistency, open/closed outgoing and reversed branches, material restrictions, cutoff and invalid-reference errors for both map families in focused tests; document the convention and error conditions alongside the API.
- [x] 1.2 Preserve verified or unspecified radiation provenance through domain identities, clones, evaluation halos, search results, refinement and chart plans. Verify all existing construction paths retain unspecified provenance, foreign inputs fail before evaluation and coordinate-sign changes do not select exterior signs; document provenance scope and the requirement that the anchor lies within the requested rectangle.
- [x] 1.3 Add a TE/TM normal-incidence cavity regression using a valid explicit domain containing both a real all-outgoing anchor and the independent pole k0=(pi-i*ln(3))/(2*d). Verify the half-period incoming copy cannot satisfy the outgoing request, symmetric mixed requests report mismatch, and no chart translation, geometry expansion or pointwise square-root reselection occurs; add the example to the mathematical guide.

## 2. Traceable physical point descriptions

- [x] 2.1 Add finite physical point data and adapters for single-chart accepted roots and batch accepted observations, retaining exact problem/domain ownership, acceptance evidence, coordinates, polarization, ordered kappa/alpha and radiation provenance. Verify consistent descriptions across adapters, foreign-domain rejection before evaluation, and retained typed errors when point data or batch signatures are unavailable; document snapshot and identity semantics.
- [x] 2.2 Implement independently accessible exterior decay/growth, solve-relevant spatial/temporal behaviour and supported real-frequency outward flux coefficients with validated, recorded tolerance bands. Verify all sign directions, dimensional scaling, zero/boundary cases, invalid settings and non-finite inputs; document why complex-frequency flux is unsupported and why the coefficient is not integrated mode power.

## 3. Dielectric mode classification and interpretation limits

- [x] 3.1 Implement all eight strict named dielectric sectors from the design's physical sign table, with left/right labels, guided-mode handling and explicit boundary, inconsistent-evidence and unsupported reasons. Verify physically consistent alpha²/beta² fixtures against the paper's independent sign tuples, including reversed exterior ordering and symmetric limits; document the table and exact exterior/propagation assumptions.
- [x] 3.2 Exercise classification on the existing independently engineered TE/TM guided slab and selected-sheet propagation fixtures where applicable. Verify cutoff and improper real-axis cases remain explicit, and complex/magnetic exteriors or nonpositive propagation direction withhold named types while preserving available components; publish the support matrix in docs/modal-transform.md.
- [x] 3.3 Add spectral interpretation with temporal behaviour and radiation qualification from verified anchor provenance, withholding fixed-frequency paper labels. Verify the anchored outgoing cavity pole, arbitrary-sheet decaying points and verified mixed/incoming alternatives remain distinguishable; document the limited outgoing-continuation qualification and deferred QNM normalisation.

## 4. Preserve result evidence and document the handoff

- [x] 4.1 Add optional batch interpretation without changing discovery, refinement or duplicate grouping. Verify one outcome per accepted observation, independent continuation after interpretation failure, unchanged source observations/groups/coverage, and no promotion of clusters or other non-accepted candidates; document the separation of physical interpretation and numerical acceptance.
- [x] 4.2 Update README, included mathematical guide and roadmap with real-reference/private usage examples, interpretation tolerances, unsupported cases and the next public point-solver API milestone. Verify links and rustdoc render, and that fields/derivatives, branches, QNM normalisation, features and lamina-design optimisation retain their agreed order and ownership.

## 5. Integration verification

- [x] 5.1 Run locked lamina-mode all-target tests, doctests, clippy and no-dependency rustdoc plus changed-file formatting. Record passing results against the 240-test baseline and identify existing warnings separately; run relevant core tests if core code changes.
- [x] 5.2 Review the completed implementation against every new spec scenario and the support matrix, verify no unrelated dependency/upstream edits remain, and run strict OpenSpec validation. Record the explicit limits on radiation qualification, classification and coverage in the verification summary.
