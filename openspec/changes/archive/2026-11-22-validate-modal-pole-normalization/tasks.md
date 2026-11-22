# Tasks

## 1. Physical conventions and independent references

- [x] 1.1 Implement dev-only analytic TE/TM planar-cavity and matched fundamental-solution Green references with the specified positive unit flux jump; verify boundary matching, independent pole locations and the jump without using production normalization or characteristic derivatives.
- [x] 1.2 Publish the scalar operator, reciprocal vector pairing, k0/beta held-fixed conventions, units and planar Fourier convention; verify the analytic cavity residue, the normal-incidence finite-layer Maxwell equivalence without double counting and explicit angular-frequency Jacobian agree with those definitions.

## 2. Complete core normalization

- [x] 2.1 Add finite-layer operator-derivative bilinear contributions using core analytic integrators and actual epsilon/mu dispersion; verify independent analytic layer integrals for TE/TM and electric/magnetic dispersion, and document conversion from legacy electromagnetic overlaps.
- [x] 2.2 Add analytic left/right exterior boundary contributions on stored roots, with empty finite-layer sums allowed and optional finite exterior partition shifts; verify zero-layer interface normalization, finite leaky regularization and partition cancellation without inserting stack layers.
- [x] 2.3 Expose a high-level complete-normalization capability for supported core backends, retaining contributions, total, units and conditioning; verify external use without private workspace bounds and focused legacy normalized-mode compatibility, and document the additive complete versus legacy paths.

## 3. Owned accepted-point normalization and sampling

- [x] 3.1 Add explicit pole-variable normalization on owned accepted points, using explicit local analytic declarations, existing core material derivatives and exact branches; verify all three solve strategies, dropped problem/report ownership, incompatible-convention rejection and unchanged acceptance/radiation provenance.
- [x] 3.2 Add owned normalized field sampling with core positions and field containers, recorded seed gauge, factor and square-root scale; verify one common E/H scale, interface-side order and invariance up to a common sign under multiple nonzero complex seed multipliers.
- [x] 3.3 Add validated dimensioned factor/residual thresholds, conditioning evidence and structured source-preserving failures; verify singular/near-zero/nonfinite factors, invalid settings and missing analytic assumptions leave the accepted point and arbitrary-amplitude recovery intact, without legacy fallback or exceptional-point claims.

## 4. Independent physical correctness and reproducibility

- [x] 4.1 Add residue-convergence tests using independent Green values at several source/observation locations and at least six decreasing complex offsets in multiple directions; verify convergence before roundoff and cross-check with independent small-contour residues on isolated regular fixtures.
- [x] 4.2 Extend validation to guided propagation, both complex interface-plasmon cases, both coupled film modes and dispersive finite-layer/magnetic fixtures; verify complete residues, backend agreement and independently shifted exterior partitions under the same documented conventions.
- [x] 4.3 Add deliberate incorrect-term controls for missing exterior terms, missing dispersion, conjugation and wrong pole variable; verify each control fails against a fixture chosen to expose that error, and document actual reference-error and acceptance thresholds.
- [x] 4.4 Add bounded quick and full public validation runs with fresh machine-readable contribution/residue/failure records and optional convergence plots; verify commands reproduce the report, failed/incomplete data are rejected, and public docs distinguish arbitrary gauge, legacy scaling, complete normalization and deferred observables.

## 5. Integration review

- [x] 5.1 Run targeted and full locked core/mode tests, doctests, lint, rustdoc, formatting and existing/new public examples; review cross-crate compatibility and the full validation report, tidy temporary changes and record results and pre-existing warnings.

## Workflow follow-up

- Apply after `recover-modal-fields`; the general local-derivative API remains queued independently; sync and archive when requested after review.
- Branch continuation remains independently implementable. Propose mode volumes and excitation observables separately once this normalization milestone passes; retain feature sensitivities and lamina-design optimization on the roadmap.
