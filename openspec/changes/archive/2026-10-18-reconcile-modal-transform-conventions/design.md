# Design

## Context

See proposal.md for motivation. The target implementation is the sibling lamina-mode repository; planning lives in lamina-core. The map modules, problem construction, refinement and chart orchestration are private. The existing explicit-domain contract already retains selected coordinate branches and ordered exterior wavevectors; overlap reconciliation requires compatible exterior signatures.

Primary reference: A. Bakhtazad, H. Abiri and R. Ghayour, "A General Transform for Regularizing Planar Open Waveguide Dispersion Relation", Journal of Lightwave Technology 15(2), 383–390 (1997), DOI 10.1109/50.554392. The supplied local PDF is ~/Downloads/bakhtazad1997.pdf. Public documentation must use the bibliographic reference and DOI, not depend on that local path or redistribute the paper.

Equations 6–9 on printed pages 384–385 establish a full exterior-sheet period 2*pi/eta. Figure 2 and section III classify fixed-real-frequency propagation modes under real exterior dielectric-index assumptions. Current map documentation acknowledges exterior sign reversal at i*pi but incorrectly identifies width-pi boundaries. Both canonical builders select that narrower strip; orchestration records OutsideFundamentalStrip without explaining the missing exterior sheets. Existing spectral tests already distinguish canonical mixed exterior sheets from an explicit all-outgoing normal-incidence cavity pole.

## Goals / Non-Goals

**Goals:** Establish an independently checked mathematical contract, correct the coverage vocabulary and ledger, and provide a durable public guide and ordered roadmap before freezing an API.

**Non-Goals:** Implement a classification API, generate exhaustive full-surface charts, recover seams, merge disjoint periodic observations, generalise to dispersive spectral exteriors, or implement fields, continuation, normalisation or optimisation.

## Decisions

### 1. Preserve the map formulas unless independent checks expose a defect

For the paper's nonmagnetic case identify L with cladding and R with substrate. Its R=(n_sub^2-n_cl^2)/4 is minus our r; its S is our s. With z=eta*Z, the paper's w=beta^2/k0^2 equals our in-plane squared-coordinate expression. Algebraically its normalised decay constants g and h correspond to kappa_L/(i*k0) and kappa_R/(i*k0). This establishes formula correspondence; the implementation must additionally trace the core's left/right field exponent choices before applying physical sign labels. Swapping exterior order requires an explicit dictionary, not silently reusing the same region labels.

For either current map, exp(z+i*pi)=-exp(z). Thus the squared physical coordinate repeats, both exterior wavevectors reverse, and only the ordered exterior data repeat after 2*i*pi. The separate analytic square root selecting the signed physical coordinate can have its own continuation behaviour: do not promise that a closed continuation returns the same signed coordinate solely from this exterior-period identity.

For generic asymmetric exteriors, the quadratic inverse for exp(2z) has two roots; each yields two exterior-sign representatives within a 2*pi strip, giving four exterior-sheet preimages of a generic squared coordinate. Test exceptional light-line and symmetric limits separately rather than requiring four distinct points everywhere. The equal-exterior r=0 limit has reduced sheet structure and does not justify the generic width-pi identification.

Alternative rejected: replace the uniformisation wholesale. The algebra points to a representation/coverage error, not a demonstrated defect in the exponential formulas.

### 2. Preserve finite canonical geometry and correct its declared scope

Keep current canonical rectangles and margins, since enlarging them automatically would change cost, branch validation and user-selected scope. Replace misleading fundamental-strip terminology in the exclusion model with a selected canonical strip and an explicit omission of the complementary exterior-sheet sectors modulo 2*pi. Represent the selected interval and full exterior period in structured coverage data, or equivalent typed bounds, so consumers can distinguish a sheet omission from an ordinary contour margin. Carry this information unchanged into batch reports. Update private call sites and tests together.

Explicit domains can search other sectors now. Prove that path with a physical fixture outside the selected strip; do not imply that the existing canonical plan finds all outgoing QNMs. Keep duplicate grouping conservative: no disjoint periodic merging in this change. Review all period-dependent validators, inverse enumeration, boundary wording and tests, correcting only uses that actually conflate squared coordinates with exterior sheets. An i*pi period can remain correct for squared-coordinate branch-point locations.

Alternative deferred: automatically add translated charts spanning 2*pi. That needs explicit sheet-selection and seam-coverage policy and deserves a subsequent change.

### 3. Publish the correspondence before introducing classification labels

Add lamina-mode/docs/modal-transform.md, link it from README, and repair both map modules' rustdoc and related geometry comments. Explain axes and physical-coordinate naming, angular-wavenumber units, time convention, ordered field exponents, alpha-to-kappa conversion, r/R sign, eta scaling, squared versus signed physical coordinates, inverse multiplicity, full versus selected strips, symmetric limits and paper assumptions.

Document the spectral map as an extension with its own coefficients and fixed square-root choices for constant exterior constitutive products; it is not the paper's fixed-frequency chart. Complex or magnetic exterior support in the algebra does not establish Figure 2's physical classification. Dispersive interior materials retain the existing caller-assumed analytic-domain obligation; the paper's regular dielectric characteristic does not prove generic material pole-freedom or eliminate floating-point failures.

Explain decay/growth and incoming/outgoing exterior behaviour separately from temporal versus spatial attenuation. The paper includes source sectors as well as leaky, absorptive and gain sectors. Numerical root acceptance establishes neither these labels nor outgoing-frequency QNM status. Classification implementation follows this convention contract.

### 4. Keep a durable roadmap with clear crate boundaries

Publish the following ordered milestones in the guide, explicitly superseding older README handoffs that put the public API immediately next. These are future directions, not implementation tasks of this change.

1. **This change: transform and sheet contract.** Formula correspondence, public conventions, corrected coverage and independent regression evidence.
2. **Physical point results and public API.** Physical-input solver, executable examples, numerical evidence and provenance; classification from ordered exterior behaviour under stated assumptions, with boundary/unsupported cases visible. Support both current solve directions without claiming full-domain coverage. Decide explicit outgoing-sheet selection and broader chart planning here.
3. **Fields and local derivatives.** Recover arbitrary-amplitude homogeneous mode fields from an accepted point, with reconstruction failures separate from root validity; expose material/coordinate/parameter derivatives needed for continuation and normalisation. Core owns reusable reconstruction, differentiation and integration primitives; mode owns mode interpretation. This milestone enables field overlap and later QNM work.
4. **Branches.** Adaptive predictor-corrector continuation, branch identity using coordinate/exterior/field evidence, crossings and unresolved cases, rediscovery and explicit coverage limits. Add arclength parameterisation when the chosen continuation coordinate fails. Profile before undertaking performance changes.
5. **QNM normalisation and derived observables.** State the outgoing frequency-QNM formulation, material dispersion derivatives, exterior regularisation and amplitude convention; retain fallible normalisation separately from field recovery. Address left/right fields where the formulation requires them. Planar structures require per-unit-length/area conventions; do not present a planar effective quantity as an unspecified finite 3D mode volume. Never use a naive integral of diverging exterior fields.
6. **Dispersion features and sensitivities.** Locate features on branches and continue them under design variation. For D(omega,k,p)=0, a regular ZGV feature obeys D_k=0 with D_omega nonzero; fixed-frequency root multiplicity does not itself imply an exceptional point. Distinguish zero complex slope, stationary real frequency and zero energy velocity. Feature sensitivities need the appropriate second derivatives and a nonsingular feature Jacobian.
7. **Optimisation integration.** lamina-mode supplies points, branches, features, sensitivities and quality evidence; lamina-design owns objectives, constraints and optimisation algorithms, including moving a specified ZGV feature. Dependencies run design -> mode -> core, with no mode -> design dependency.

Numerical follow-ups remain visible across milestones: seam recovery, repeated roots and cluster splitting, periodic equivalence, material poles, broader exterior dispersion and scaling/dynamic range. They must be prioritised against the target device rather than silently advertised as supported.

## Risks / Trade-offs

- [Correct algebra mistaken for verified physical sign correspondence] -> Derive exponents from core and check an independent physical slab before publishing labels.
- [Retaining narrower canonical charts still misses desired modes] -> Expose omitted sectors in structured reports and demonstrate explicit charts; defer automatic coverage honestly.
- [Tests mirror formulas without detecting the original error] -> Check inverse sheet multiplicity, sign reversal and a known outgoing slab pole, not just self-consistency.
- [A 2*pi exterior identity mistaken for every signed-coordinate continuation period] -> Keep coordinate-sheet and exterior-sheet provenance separate.
- [Public guide becomes stale] -> Link from README and rustdoc, pair mathematical claims with named regression tests, and update the roadmap at each milestone.

## Migration Plan

Add independent regressions, repair period-dependent semantics and the private exclusion representation, then publish the guide and replace outdated README handoffs. Preserve explicit-domain behaviour and existing roots within their selected charts. Run focused and full lamina-mode checks and relevant core checks only if core changes. No dependency upgrade is expected. If the exclusion representation must be revised, retain the corrected scientific documentation and adapt private consumers together rather than reverting to false period claims.

## Implementation verification (2026-10-18)

- lamina-mode: 240 all-target tests pass, with no failures or ignored tests (238 before this change). The two new map tests establish inverse sheet multiplicity and full exterior periods; existing orchestration tests now also assert omitted-sector bounds, report preservation, simultaneous exterior-sign incompatibility, and the nonzero incoming half-period copy of the independently known outgoing cavity pole.
- All 180 tests matching `conformal` and all 10 matching `orchestration` pass. Existing independent TE/TM guided and canonical selected-sheet slab roots remain accepted.
- Locked offline doctests pass (zero executable doctests); the public guide is included in crate-level rustdoc, and the no-dependencies documentation build succeeds without warnings. Local Markdown links resolve.
- Changed Rust files pass rustfmt checks. Clippy succeeds with the same eight pre-existing warnings (type complexity, negated partial-order comparisons and unit-struct defaults); no new warnings.
- No map formulas, canonical rectangles, coordinate validators, dependencies, lockfiles or core implementation changed. The period-dependent execution audit found no additional sheet-conflation defect: squared-coordinate inverse and branch-point enumeration legitimately retain i*pi spacing. Existing unrelated working-tree changes were preserved.
- Strict OpenSpec validation passes. Automatic full-surface search and physical classification remain future milestones; the canonical ledger now explicitly reports omitted exterior-sheet sectors rather than describing its width-pi selection as the full surface.
