# Design

## Context

See proposal.md for motivation. The benchmark and field proposals establish reusable physical fixtures. Current lamina-mode exposes a checked characteristic and one native derivative: dF/dz for conformal charts or dF/dk0 for direct frequency. Its accepted points already retain physical coordinates, ordered exterior wavevectors, source analytical stack and numerical evidence. They retain a backend name rather than an executable backend instance.

Core exposes `ModeJet1`, constitutive lifting, checked determinant evaluation, `Parameter::{Spectral,InPlane,LayerThickness}` and `FiniteLayerIndex`. Internal stack compilation can seed selected thickness derivatives before caller-unit conversion. The public `compile_canonical_constant_stack` intentionally seeds no thickness derivative, while canonical layer construction/access is private. This is a demonstrated small missing public capability; mode must not recreate canonical stack internals or a parameter framework to work around it. First-order coordinate gradients need only repeated single-direction evaluations, not public second-order/bivariate jet exports.

## Goals / Non-Goals

**Goals:** Owned local physical first partials, actual constitutive/exterior contributions, explicit parameter and unit semantics, independently checked simple-root slopes and failures separate from root acceptance.

**Non-Goals:** No field derivatives, arbitrary material-coefficient seeding, Hessians, mixed second partials, feature/ZGV/exceptional-point solver, branch tracker, optimization or new normalization. Material dispersion with frequency is included; a material's oscillator parameters are held fixed. Source analytic search coverage is not expanded by local derivatives.

## Decisions

### 1. Seed physical variables directly at the accepted point

Use a public operation such as `point.local_derivatives_with(backend, request, local_assumption)` with an explicit supported backend value and core parameter selectors. The local declaration has no default and covers a holomorphic two-coordinate neighbourhood and the selected real thickness variations, with nonzero exterior q. An in-plane search's fixed-k0 material assumption is insufficient to establish a frequency derivative; state this additional obligation explicitly.

At the source point, obtain k0, beta, ordered kappa and actual core material handles. Recheck finite data and the squared exterior relations; reject missing/inconsistent branch evidence before derivative work. For each requested axis use core's first-order holomorphic jet, holding the other physical variables and material parameters constant. Evaluate the same characteristic gauge through `checked_determinant`. No raw chart inversion, numerical continuation differentiation, production finite differences or root correction is required.

Reuse backend/core evaluation helpers introduced by field work where suitable, but do not depend on recovered fields or normalization. Field work precedes this change as an implementation sequence and fixture dependency, not a mathematical requirement.

Alternative rejected: differentiate only the search chart and call that a full physical gradient. That cannot separate F_k0 from F_beta or provide geometry partials.

### 2. Lift each exterior derivative on its stored branch

For q_j=epsilon_j mu_j k0²-beta², use actual core constitutive jets and

```text
partial q_j / partial k0 = (epsilon'_j mu_j + epsilon_j mu'_j) k0²
                           + 2 epsilon_j mu_j k0
partial q_j / partial beta = -2 beta
partial kappa_j / partial a = (partial q_j / partial a)/(2 kappa_j).
```

For a thickness-only partial, exterior kappa and both coordinates are constant. Choose the local algebraic jet sign to match the stored kappa and verify primal/derivative consistency; endpoint principal roots cannot replace accepted branch data. Material dispersion in finite layers remains part of the core seeded optical computation. Validate finite derivative components and per-side cutoff/nondegeneracy bands; exact or near q zeros produce a typed unavailable derivative, even if the source root was accepted on a strategy allowing that boundary. Do not perturb a cutoff or silently change sheets/backends.

Native consistency tests check F_z=F_beta*beta_z for in-plane charts, F_z=F_k0*k0_z for constant-exterior spectral charts and F_native=F_k0 for direct frequency. Physical partials use k0=omega/c in cm^-1, not Hz or angular frequency in s^-1. If a caller later converts to omega, F_omega=F_k0/c with c in the corresponding units; no unlabeled substitution is made here.

### 3. Expose only the missing core thickness-compilation operation

Add a narrow public core helper to compile a canonical stack with a first derivative for one selected finite-layer thickness in slot zero. Reuse the existing internal derivative mapping and compiler, validate the public `FiniteLayerIndex` against the source stack and keep all unselected layers constant. Do not expose unchecked canonical-layer constructors, a parallel stack or a generic material-parameter system. Existing constant compilation remains unchanged.

The derivative is with respect to the thickness value in its stored `LengthUnit`, matching the existing compiler's seed-before-centimetre-conversion rule. Carry that unit and selected layer in the public result. For physically identical stacks F_d[nm]=1e-7 F_d[cm]. Coordinate partial units remain F per cm^-1; implicit thickness slopes use the same labelled source length unit. Core compilation errors and invalid layer selections remain typed.

Alternative rejected: seed a new canonical centimetre thickness in mode and label it with the original unit. That silently loses the core conversion chain rule.

### 4. Own successful and failed partials without hiding the source

Return an owned local report containing source point/provenance, original and evaluation backend labels, physical point, local declaration, request/settings, checked characteristic value and ordered per-parameter outcomes. Share source context; no discovery rerun or references to temporary evaluators. Global invalid request/settings fail before evaluation; parameter-specific failures remain associated with their selectors and do not erase successful partials. Preserve original numerical/backend/continuation causes through the public diagnostic model. No field reconstruction or material normalization is a hidden prerequisite.

An explicit different supported backend can be requested for comparison at the same physical point and branch, with both backend identities recorded; no automatic substitution occurs. Non-finite source data, bad constitutive derivatives and unsupported selectors are not converted to zero derivatives. Source acceptance, classification, groups and coverage remain untouched.

### 5. Compute conditional first-order implicit slopes separately

Offer ratios as optional results over the gradient, never as a new root acceptance test:

```text
d beta / d k0 = -F_k0/F_beta
d k0 / d beta = -F_beta/F_k0
d beta / d d_j = -F_dj/F_beta
d k0 / d d_j = -F_dj/F_k0.
```

Require finite denominator derivatives outside an explicitly configured nondegeneracy band. The denominator's absolute threshold and reference scale carry the derivative units of that axis and characteristic gauge; relative settings are dimensionless, finite and validated. Retain the denominator, threshold, numerator and resulting slope units. Zero, near-zero or unrepresentable ratios are unavailable diagnostics while the valid partials remain inspectable.

These are local simple-root implicit slopes conditional on acceptance and analytic assumptions. A complex dk0/dbeta is not automatically energy velocity or a real group velocity, and a small denominator is not an exceptional-point or ZGV diagnosis. Second-order feature equations and branch identity remain later work.

### 6. Validate independently, with finite differences confined to tests

Use guided TE/TM slabs, the interface plasmon in both complex cases and the two symmetric-film roots from field validation. Include actual finite-layer Lorentz dispersion, nonzero exterior epsilon/mu dispersion in a focused analytic fixture, thickness derivatives in equivalent nm/cm stacks and opposite branch sign checks. Compare F values and first partials with independent scalar characteristic derivatives on regular points as well as roots; root-only tests can hide derivative-of-gauge terms.

For finite-difference validation use a documented decreasing step sweep, the same physical branch anchored to the accepted kappa and valid local neighbourhoods, and demonstrate a convergence interval before roundoff. Independently recompute perturbed roots on their own source domains to check implicit slopes; compare physical numbers without pretending different-device domain identities are the same. No numerical derivative is used as the production implementation. A local jet derivative is not a derivative of the analytic reference root provider.

## Risks / Trade-offs

- [Native and physical derivatives are conflated] -> Distinct parameter labels/units and all three strategy chain-rule regressions.
- [Changed units silently rescale geometry sensitivities] -> Reuse core's existing seeding semantics and equivalent nm/cm device tests.
- [Frequency partials exceed the original search assumption] -> Require an additional explicit local declaration; do not extend coverage claims.
- [Tiny slope denominators produce misleading finite predictions] -> Validated dimensioned nondegeneracy checks and independent perturbed-root validation.
- [Geometry support becomes a new material framework] -> Only existing core selectors; arbitrary material-coefficient sensitivities remain deferred.

## Migration Plan

Third working change: after benchmark preservation, and after field recovery for shared fixtures/API conventions. Add the small core thickness compiler and targeted core tests first, then physical point jets, public owned outcomes, guarded slopes and runnable derivative examples. Run relevant core compilation/derivative tests and required full core checks for the actual edits, plus current mode tests, doctests, clippy, rustdoc and formatting. Verify every spec scenario and update the guide with units, local assumptions and derivative limits. This completes the current fields/local-derivatives milestone and hands off to adaptive branch continuation. QNM normalization, features/sensitivities and lamina-design optimization remain separate changes.
