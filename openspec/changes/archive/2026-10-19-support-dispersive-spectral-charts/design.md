# Design

## Context

See proposal.md for motivation. OpenSpec planning lives in lamina-core; implementation is primarily in sibling lamina-mode. The existing 250-test baseline covers checked optical evaluation, dimensionless conformal charts, bound problems, refinement, finite chart plans and physical descriptions with radiation provenance.

Core already supplies AnalyticalMaterialStack, AnalyticalMaterialHandle, complex-frequency material evaluation/derivatives, ConstitutiveLift and checked determinant evaluation accepting explicit exterior wavevector jets. Its material contract distinguishes analytical complex continuation from real-axis tabulated data. The constant-exterior restriction comes from lamina-mode's ComplexWavenumberMap coefficients and spectral constructor, not from core's stack model or checked evaluator.

The current chart implementation owns one analytically continued square root for the signed physical coordinate; raw map evaluation itself is infallible. Direct frequency search instead has an identity physical coordinate and two independently continued exterior roots, whose numerical continuation can fail. A fallible chart-level interface is the integration point; placing an unreported quadrature failure inside the existing infallible map operation is not acceptable.

## Goals / Non-Goals

**Goals:** Add a local stack-based spectral strategy with actual exterior dispersion and branch-preserving derivatives, sharing the numerical evidence and physical-description pipeline with existing strategies. Make its native coordinates, domain obligations and failure causes explicit before the public API is designed.

**Non-Goals:** Replace the constant-exterior transform; invent stack, layer, material, polarization or physical-unit primitives; certify unknown material singularities; automatically locate thresholds, split regions or integrate contours with holes; uniformise a general dispersive Riemann surface; implement homotopy, disjoint-chart radiation transport, public exports or additional physical classification labels. Fields, branches, QNM normalisation and optimisation remain later roadmap work.

## Decisions

### 1. Add a direct-frequency strategy over an unchanged core stack

Add a private problem construction path taking an existing AnalyticalMaterialStack<C>, fixed real beta, polarization and backend. Compile the same stack for checked optics and retain its actual exterior handles for the chart. Validate finite fixed inputs and existing core compilation constraints; do not evaluate exteriors once and store those values as constant materials. Preserve the current constant-exterior constructor and its construction evidence as the specialised conformal path. Strategy choice is explicit; erased material handles and finite samples do not justify automatically selecting the constant strategy.

The native direct coordinate is k=k0 in cm^-1, with dk0/dk=1 and fixed beta. Accept an explicit finite frequency rectangle using a mode-specific bounds interface backed by the existing rectangle machinery and core frequency-unit conventions. No parallel device representation is introduced. Frequency bounds and dimensionless conformal bounds must be distinguishable at construction and in provenance. Keep search direction Spectral while separately recording chart coordinate kind and units.

Alternative rejected: substitute epsilon(k) and mu(k) into the existing map coefficients. That makes the transform implicit and no longer establishes its original exterior uniformisation. A fixed-coefficient map could later be used as a preconditioner, but it is not needed for this direct-coordinate implementation.

### 2. Generalise the fallible chart boundary, not the stack model

Extract a private analytic-chart contract from the operations consumers actually need: requested/evaluation rectangles, native-coordinate kind, fixed reference data, scalar physical-coordinate projection and fallible mapping of a seeded native-coordinate jet to canonical coordinates and ordered exterior jets. Preserve the existing conformal chart through an adapter; implement a frequency chart separately. Include an operation to validate a numerical halo without rebuilding the anchor or expanding the requested search. Do not force k0 through sqrt(k0²), which would introduce a coordinate branch point at zero into an identity map.

Parameterise the bound domain and dispersion/search pipeline by that chart contract rather than assuming every chart is the existing raw-map-plus-coordinate-square-root construction. Refactor only private adapters and keep existing constructor aliases usable internally. Share characteristic evaluation, argument-principle error tracking, refinement eligibility, ownership tokens, accepted roots, finite chart execution and physical interpretation. Do not build a parallel discovery/refinement subsystem.

Change terminology that incorrectly hardcodes every native point as dimensionless z: root evidence, enclosures, correction tolerances and matching settings must retain their chart kind. Existing conformal defaults and meanings remain unchanged. Direct-chart native absolute tolerances are in cm^-1 and are recorded as such. Explicit plans contain domains of one problem/strategy; no mixed conformal/frequency plans or cross-strategy duplicate reconciliation in this change. Physical results from separate strategies can be compared in tests without pretending they share a domain identity.

Alternative rejected: inject non-finite sentinels into an infallible raw map and let checked optics detect them later. That loses side-specific continuation failures and cannot meet the error contract.

### 3. Continue the exterior square roots analytically from fixed references

For each ordered exterior define

```text
q_j(k) = epsilon_j(k) mu_j(k) k² - beta²
q'_j(k) = [epsilon'_j mu_j + epsilon_j mu'_j] k²
          + 2 epsilon_j mu_j k.
```

Choose a reference k_ref strictly inside the requested frequency rectangle and finite nonzero kappa_j_ref satisfying kappa_j_ref²=q_j(k_ref) within validated settings. Under the explicit holomorphic/nonzero-q domain assumption, continue along the straight segment inside the rectangular chart:

```text
kappa_j(k) = kappa_j_ref * exp( 1/2 integral[k_ref -> k] q'_j(s)/q_j(s) ds )
kappa'_j(k) = q'_j(k)/(2 kappa_j(k)).
```

Evaluate every point from the same fixed reference. This avoids dependence on adaptive-contour evaluation order or on the previous sample. On a valid simply connected region the analytic branch is path independent mathematically; numerical continuation remains tolerance limited. Reference consistency, reconstructed kappa²=q, finite quantities and proximity to sampled q zeros are checked pointwise without calling those checks domain certification.

Use actual core material values and first derivatives for both the integral and local q jet. Use ConstitutiveLift to retain the full exterior material dependence. Select the sign of a local algebraic sqrt(q_jet) by agreement with the independently continued scalar kappa, then verify its value and derivative against the continuation formulas. This is branch-preserving jet construction, not endpoint principal-root radiation selection. The physical characteristic derivative is dF/dk0. Constant-path comparisons convert conformal derivatives using dk0/dz before comparing.

Extend the existing adaptive Simpson machinery to a fallible integrand or factor an equivalent reusable implementation: preserve errors from material evaluation/numerical checks without encoding them as successful values. Validate dimensionless integration/consistency tolerances, physical cutoff bands and explicit recursion/evaluation budgets. Report which side failed, failure kind and available work diagnostics. Do not silently fall back to a principal root, constant model, alternate backend or alternate search strategy. Control repeated continuation cost through budgets first; caching is optional only if keyed by the immutable chart and cannot alter branch decisions.

Alternative deferred: continue roots through a homotopy from frozen exteriors. It may supply useful seeds later but does not replace a discovery method or establish missing-mode coverage.

### 4. Keep analytic-domain responsibility explicit and record known exclusions

Keep the existing explicit caller assumption covering a holomorphic, pole-free characteristic and valid materials on a neighbourhood of the closed evaluation region. Add a separate required, non-default declaration that each q_j is holomorphic and nonzero there. These assumptions are stronger than finite contour checks. They are obligations on the caller/model information, not certificates manufactured by domain construction.

Domain input also carries a list of known excluded frequency rectangles with typed reasons (material singularity, exterior zero/threshold with side, or caller-declared exclusion). Exclusion lists may be empty explicitly; that means no known exclusions supplied, not none exist. Validate their finite bounds and reject intersection with the entire closed requested rectangle or the computed halo, including boundary contact. Preserve declarations and finite truncations in the coverage ledger. Only rectangular charts disjoint from exclusions are executable; the caller supplies separate charts to cover additional regions. No subtraction of holes from an APM rectangle and no automatic subdivision.

Every q zero is excluded conservatively in this first implementation, even if an even-order zero could admit a removable analytic square root. Material/backend poles and invalid constitutive factors remain governed by the characteristic assumption as well. The native identity coordinate does not itself require k0 nonzero; any material or backend failure at zero remains an actual physical/evaluation limitation rather than an invented coordinate branch.

The numerical halo retains the exact reference and branch values and is checked against the same declarations. A halo intersecting a known exclusion fails before discovery; callers must choose a different domain. A continuation path encountering a non-finite model, sampled zero or exhausted budget produces a diagnostic failure, never an assertion that every other interior point is safe.

### 5. Separate explicit branch data from verified physical radiation

Offer two reference modes. An explicit complex reference supplies ordered kappa values checked against the actual q values; it creates branch-consistency provenance with physical radiation Unspecified. This supports complex or magnetic exterior models where the current real-anchor dielectric interpretation is unavailable. Negating either reference must preserve a distinct exterior sheet.

A real radiation-reference mode reuses the existing per-side outgoing/decaying or incoming/growing choices. At real positive k_ref, require actual positive real exterior epsilon and mu=1, real fixed beta and separation from cutoffs. Verify those actual anchor values and continue from them. The material can vary away from the anchor. Thus a lossless dispersive model can have verified radiation at the real anchor and complex exterior values at its resonance. Lossy exteriors that fail those anchor assumptions remain supported through explicit branch references but do not receive a fabricated physical outgoing qualification.

Retain anchor native-coordinate kind, material values at the anchor, verification settings, explicit assumptions and ordered branches in the same problem/domain provenance system. Physical descriptions at accepted roots continue to report temporal behaviour, exterior growth/decay and radiation provenance. The paper's fixed-frequency named types and instantaneous complex-frequency flux remain unsupported. A referenced outgoing continuation is still conditional on the chart/material assumptions and is not a normalised QNM or proof of global causal-sheet coverage.

### 6. Prove non-frozen behaviour and preserve existing evidence

Use three complementary regression groups:

- **Analytic branch fixture:** for example q(k)=exp(k) on a rectangle crossing its principal square-root cut but excluding k=0 if represented by epsilon=exp(k)/k² at beta=0. The known continuation is exp(k/2). Check both sides, opposite reference signs, shuffled evaluation order, physical derivatives and typed quadrature failures. Material values and derivatives come through a core analytical model/handle, not a mock branch-label lookup.
- **Dispersive TE/TM device:** a finite constant-index slab with at least one Lorentz exterior from lamina-materials (or a core analytical test model with the same explicit formula). Choose a region away from its poles and q zeros. Check roots against an independent scalar slab characteristic with verified branches, and check full dF/dk0 against independent derivatives. Use nonzero oscillator strength that shifts the root beyond acceptance tolerance relative to freezing at k_ref. Include nonzero permeability derivatives in a focused jet test even if the real-radiation device fixture is nonmagnetic. A lossless oscillator permits supported real-anchor verification without pretending epsilon remains constant at complex k.
- **Constant limits and integration:** zero oscillator strength agrees with existing constant-exterior conformal values and roots on matching sheets for both polarizations. Preserve 250 existing tests. Exercise single-chart refinement/physical descriptions, multiple explicit frequency charts, overlap grouping, opposite exterior signatures, foreign domains, continuation failures, budget-limited clusters and coverage exclusions. Refinement never accepts a root outside its originating frequency enclosure or into the halo.

Do not require both strategies' raw native corrections or enclosures to match: their coordinates differ. Compare physical values and derivatives after the appropriate chain rule, preserving each strategy's own evidence.

### 7. Document the strategy distinction and revise the roadmap

Update README and docs/modal-transform.md, which is included in crate rustdoc. State that constant exteriors are a conformal-strategy restriction, not a core-stack or universal frequency-solver restriction. Explain local frequency charts, native units, reference modes, stronger domain assumptions, known exclusions, numerical controls, derivative semantics and unsupported physical classifications. Provide a private core-stack usage example and name the independent dispersive regressions. Clearly distinguish caller-declared exclusions from discovered or certified singularities.

The roadmap becomes: this dispersive spectral-chart change, then public point-solver API and executable examples, then fields/derivatives, branch continuation, QNM normalisation, feature sensitivities and optimisation in lamina-design. Core retains stack/material/field primitives; mode owns search/branch interpretation; design owns optimisation. No public API is frozen here.

## Risks / Trade-offs

- [Unknown poles or zeros invalidate analytic continuation] -> Require explicit full-domain assumptions, reject known exclusions and preserve runtime failures; do not infer certification from sampling.
- [Nested exterior integration increases search cost] -> Set per-continuation work bounds and retain diagnostics; profile before adding caches or different quadrature.
- [Refactoring charts weakens existing branch/error contracts] -> Keep conformal adapters and the 250-test baseline, add typed failure injection, and make native-coordinate kinds explicit.
- [A complex reference is mistaken for physically outgoing radiation] -> Preserve explicit versus real-verified provenance separately; keep existing physical support restrictions.
- [Local charts miss modes near excluded thresholds or material poles] -> Report excluded scope and retain the constant-exterior global transform as a distinct specialised strategy; no global completeness claim.

## Migration Plan

First factor the fallible chart adapter and coordinate metadata with unchanged conformal behaviour. Add the direct stack-bound problem, domain declarations and two exterior continuations. Integrate the existing search/refinement/plan/description consumers, then add independent dispersive and constant-limit fixtures and publish the guide. Keep existing internal constructors supported; change no core stack representations or dependencies. Run focused regressions followed by locked all-target tests, doctests, clippy, rustdoc, changed-file formatting and strict OpenSpec validation. Record existing warnings separately. Any unexpected need for a core interface change must be justified by a demonstrated missing primitive and kept scoped.

## Implementation review and verification

See [verification.md](verification.md) for the 2026-10-19 scenario review, runnable
regressions, integration checks, retained limitations and public-API handoff.
The longer-term roadmap is also maintained in lamina-mode's public modal guide.
