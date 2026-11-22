# Design

## Context

See proposal.md for motivation. `recover-modal-fields` now provides owned arbitrary-amplitude E/H, exact source branches and reconstruction-defect evidence. Its tests establish physical profiles, including a necessary correction to homogeneous scattering reconstruction. `expose-modal-local-derivatives` is planned but not implemented. The user has selected normalization now: use existing analytic core material derivatives for the required internal physical variations, with an explicit local declaration. Do not implement or duplicate the broader public gradient API in this change.

Core's `ComplexPlaneMode::new` obtains a seed, calls `raw_qnm_normalisation_unchecked` and scales by its inverse square root. That routine integrates finite layers and aggregates `BilinearLayerNormalization`: electric minus magnetic, with constitutive spectral derivatives and no conjugation. The aggregate fails on zero layers. No explicit exterior contribution is added along that call path. This is an observed finite-layer convention, not evidence of complete normalization. Existing normalized/unscaled compatibility only checks a common scale.

## Goals / Non-Goals

**Goals:** A complete, explicitly defined simple-pole normalization on existing isotropic reciprocal planar stacks; owned normalized fields; independent evidence that the normalization predicts residues. Handle frequency and propagation poles without conflating their derivatives or units.

**Non-Goals:** No generic Green-function public solver, PML mesh, driven-field expansion API, 3D mode volume, Purcell factor, exceptional-point scheme, new stack primitives or silent replacement of legacy numerical behavior. A normalization on a chosen sheet does not certify causality or outgoing radiation provenance.

## Decisions

### 1. Define normalization through a physical operator and its residue

Use the existing transverse scalar primary field: psi=Ey for TE and psi=Hy for TM. Define the canonical scalar operator in cm units by

```text
L(lambda) psi = d_z [p(lambda,z) d_z psi] + q(lambda,z) psi
p = 1/f; q = k0^2 g - beta^2/f
TE: f=mu,      g=epsilon
TM: f=epsilon, g=mu
```

The Green-function convention is `L G(z,z';lambda) = delta(z-z')`, with a positive unit flux jump, the stored exterior sheets, and a planar Fourier component without an extra transverse area or Fourier prefactor. This is a per-transverse-component convention, not a finite 3D volume normalization. Record its identifier and dimensions in every result. Reconstruct all E/H components with the existing core impedance-scaled conventions.

Select `lambda=k0` at fixed beta for either spectral strategy, or `lambda=beta` at fixed k0 for a propagation point. Native conformal z is never the normalization variable. Any omega presentation follows k0=omega/c with explicit derivative and residue Jacobians; do not relabel k0 as Hz. A raw characteristic derivative depends on characteristic gauge and is not by itself the normalization.

For reciprocal scalar problems the residue at a regular simple pole is represented by the primary-field product divided by a complete operator-derivative bilinear factor. The scalar equation is symmetric under beta reversal; document the corresponding reciprocal full-vector pairing and the component sign changes explicitly. Do not substitute a Hermitian absolute-square product or assume the same full E/H pairing at +beta suffices. The operator convention fixes the otherwise ambiguous scalar factors when translating the legacy electromagnetic overlap to a pole normalization.

Alternative rejected: define correctness by making the existing finite-layer total equal one. That is a scale identity for one implementation and says nothing about the exterior or pole residue.

### 2. Combine finite-region integrals with exact outgoing boundary terms

For finite truncation planes a and b in homogeneous left and right exteriors, the outgoing flux conditions are

```text
p psi'(a) = B_left psi(a),  B_left  = -i kappa_left/f_left
p psi'(b) = B_right psi(b), B_right = +i kappa_right/f_right.
```

Use the weak scalar-operator derivative convention

```text
N_lambda = integral_a^b [-p_lambda (psi')^2 + q_lambda psi^2] dz
           + psi(b)^2 (B_right)_lambda
           - psi(a)^2 (B_left)_lambda.
```

There is no complex conjugation. Derivatives hold the opposite physical coordinate fixed. Include constitutive epsilon/mu derivatives and `kappa_lambda=q_exterior,lambda/(2*kappa)` on the stored branches, using existing analytic core constitutive derivatives. The formula specifies the complete convention; map existing layer overlaps into it with documented core units and factors rather than treating them as equal by name. Validate the signs and Jacobians against the independent Green function before exposing successful complete results.

Integrate finite-layer products using existing analytic wave/state-product integrators. At default truncation planes use the physical stack boundaries: the boundary terms encode analytically continued infinite exterior tails, including for growing leaky modes. Empty finite-layer sums are zero, not errors. Retain each layer, both exterior terms, total, a contribution-magnitude scale and a cancellation/conditioning indicator. Zero/near-zero f, kappa or a singular local material neighborhood is diagnostically unsupported in this increment; no cutoff limiting prescription is invented.

For validation and audit, support moving either truncation plane a finite outward distance while keeping the original device unchanged. Integrate the added homogeneous exterior segments and evaluate the adjusted boundary terms with exact exponentials. This is a partition of one device, not insertion of fictitious stack layers. Reject overflow rather than claiming invariance beyond representable precision. Production normalization must not depend on a user field-sampling grid or numerical quadrature of a divergent real-axis tail.

Alternative rejected: truncate growing exterior fields at a large distance. The result would depend on an arbitrary box and generally fail partition invariance.

### 3. Add a core capability and an owned accepted-point adapter

Core owns complete operator factors, analytic integration and exterior terms, with a high-level capability for the supported built-in backends hiding private retained-workspace bounds. Reuse unscaled modal seeds and shared reconstruction. Use existing analytic core physical variations and explicit local analyticity declarations; a fixed-frequency search assumption does not authorize a k0 variation.

Mode owns a separately fallible normalization request on `ModePoint`, explicit backend and pole-variable convention, validated tolerances and source/evaluation provenance. Reject incompatible held-fixed conventions. Do not rerun a search or refiner. Return an owned report with the original accepted point, exact branches, source numerical evidence, fresh characteristic/field-defect evidence, complete contributions, units, conditioning, assumptions and any normalization failure.

Provide normalized sampling via existing core `FieldSampling` and `SpatialResponse`. The normalization report retains enough source context and its backend capability to reconstruct with the same evaluated convention, not a borrowed temporary workspace. Preserve the stable-material lifetime contract and check detectable source inconsistencies. Original arbitrary-amplitude fields remain available if normalization fails.

### 4. Make the gauge and failure contracts explicit

For a common seed multiplier a, the scalar bilinear total obeys `N -> a^2 N`; choose `scale=1/sqrt(N)` and apply it to E and H together. Record the square-root convention. A change of initial phase can change normalized fields by a common sign, so compare residue products exactly within tolerance and fields only up to that sign. Arbitrary gauge changes must not survive in predicted physical residues.

Use finite validated absolute/relative thresholds with a documented reference/contribution scale in the units of N, and a separate dimensionless conditioning limit. Do not apply an unlabelled numerical floor. Validate fresh characteristic and field defects with explicit dimensioned settings inherited from the field stage and explicit pole convention. Zero, near-zero or nonfinite factors, strongly cancelling unresolved totals, branch/material singularities and incompatible assumptions produce structured failures carrying all available evidence. A small factor does not establish an exceptional point; no automatic legacy fallback is allowed.

### 5. Verify against an independently evaluated Green function

Build dev-only reference solutions using analytic homogeneous-layer fundamental solutions and matching/weighted Wronskians, with `G=u_left(z_<)u_right(z_>)/W` and the specified unit flux jump. They must not call the production normalization integrator, use its characteristic derivative as their denominator, or reuse normalized fields as input. Use independent dispersion roots to cross-check the accepted source coordinates.

At several non-node source/observation locations, compare `(lambda-lambda_n) G(lambda)` with the normalized primary-field product, using at least six decreasing complex offsets in more than one direction, away from other poles and cuts. Track a convergence window before cancellation/roundoff; do not require arbitrarily tiny offsets to improve. Also compare an independent small-contour residue integral where a regular isolated circle is available, using reference Green values and retained convergence diagnostics. This second estimator helps separate root-position error from normalization error. No finite step sequence certifies a general neighborhood.

Quick numerical fixtures cover an analytic nondispersive dielectric cavity (TE/TM and direct/conformal spectral strategies), regular guided propagation modes, interface Drude plasmons in both complex cases (zero finite layers) and the even/odd film branches. Add dispersive finite layers and an analytic magnetic-dispersion fixture to test both constitutive derivatives. Full validation expands source positions, offsets, contours, gauges and partition shifts with configurable fixed budgets.

Require independent reference errors around 1e-7 or better on documented well-conditioned fixtures, with tighter analytic cavity/gauge checks where justified. Derive acceptance bounds from measured root-position and reference convergence errors; never loosen tolerances to accept an incorrect exterior or coordinate factor. Record actual thresholds. Negative controls omit an exterior contribution, omit material dispersion, introduce conjugation and mix the k0/beta convention on fixtures selected to distinguish those errors. Backend agreement is supplementary since both backends share normalization algebra.

### 6. Preserve legacy behavior and document exactly what is validated

Keep `ComplexPlaneState::mode()` and `seed_normalisation()` numerically compatible and identify their finite-layer convention in public docs. The new complete path uses a separately named API and report; do not silently upgrade old return values or claim that a passing common-scale regression validates them physically. An audit that finds a shared primitive mathematically incorrect must be resolved without masking the independent validation; any incompatible legacy correction requires an explicit follow-up migration rather than silently changing this compatibility decision.

Publish a runnable normalization example using core/material primitives and fresh output directories, with machine-readable model and convention metadata, per-case success/failure, contribution tables, residue errors, partition and gauge invariance checks. Optional plotting produces convergence figures in PNG and vector format and rejects incomplete/failed data. Clearly distinguish the basis-field residue convention from future EM-source observables and finite 3D mode volume definitions.

## Risks / Trade-offs

- [Matching two formulas that share an error] -> Independent unit-jump Green function, analytic cavity and deliberate negative controls.
- [Frequency and propagation poles receive the same label] -> Explicit lambda/held-fixed metadata and different operator derivatives, with wrong-variable controls.
- [Planar translational extent is mistaken for finite volume] -> Fixed Fourier-component convention and explicit absence of a 3D mode-volume claim.
- [Exterior terms or reciprocal signs are missing] -> Zero-layer fixtures, partition shifts and component/source-sign checks.
- [A cancelling factor is presented as a reliable scale] -> Contribution conditioning and dimensioned unavailable-result thresholds.
- [Normalization changes radiation interpretation] -> Retain original provenance; validate algebra on the selected branches without inventing causal certification.

## Migration Plan

Apply after field recovery; the general local-derivative API remains queued independently. Add core complete factors and analytic references first, then the owned mode adapter and normalized sampling, then independent residue/invariance validation and reproducible docs. Run targeted and full core checks for the numerical boundary, all mode checks and public examples; preserve legacy numerical behavior and arbitrary-amplitude access. Rollback removes the additive capability without changing existing solver acceptance. Branch continuation remains independent; mode volumes, excitation observables, sensitivities and lamina-design integration follow only after their own convention-specific proposals.

## References

- Sauvan et al., [Normalization, orthogonality and completeness of quasinormal modes of open systems: the case of electromagnetism](https://arxiv.org/abs/2112.08103), background on open-system normalization and regularization.
- Lalanne et al., [Quasinormal mode solvers for resonators with dispersive materials](https://arxiv.org/abs/1811.11751), background on dispersive nonlinear eigenproblems and benchmarking.

The scalar operator, boundary signs, unit-jump reference and canonical-unit mapping above are this project's explicit design convention, to be checked by independent fixtures; they are not quoted formulas attributed to those reviews.

## User clarification: no explicit exterior integral in 1D

Do not infer incompleteness merely from a finite-layer Maxwell bilinear expression. At normal incidence with nondispersive exteriors, the outgoing E/H bilinear exterior density cancels. The finite-layer Maxwell form already includes, through its derivative-field identity, the boundary contribution appearing in the scalar operator form. Add a TE/TM analytic cavity regression demonstrating `N_k0 = +k0 M` (TE) and `N_k0 = -k0 M` (TM) in core's conventions, with M the legacy E/H overlap. No exterior integral is evaluated and no scalar boundary term is added to M a second time. General dispersive/exterior and beta-pole cases use the explicitly defined scalar boundary-derivative convention and independent residue checks. Public docs must state this equivalence and its assumptions, rather than describing legacy finite-layer normalization as categorically incomplete.
