# Design

## Context

See proposal.md for motivation. Planning is in lamina-core and implementation is in sibling lamina-mode. Existing private AcceptedRoot values carry numerical evidence and exact DomainProvenance. Batch AcceptedObservation values add chart/candidate indices and a fallible BranchSignature containing k0, k_parallel and ordered exterior kappa values. Single-chart refinement lacks that convenience signature but can evaluate the originating chart. Duplicate grouping already requires compatible exterior signatures.

ReferenceSelection currently selects only the signed linear coordinate, either by CoordinateSheet or an explicit coordinate at a conformal reference. It does not select exterior radiation signs. Exterior data are analytic functions of z. The previous change established their full 2*pi period and documented that canonical width-pi charts can select mixed radiation sheets. An outgoing normal-incidence cavity pole is already tested on an explicit chart outside that selection.

Scientific reference: Bakhtazad, Abiri and Ghayour, JLT 15(2), 383–390 (1997), DOI 10.1109/50.554392, Figure 2 and section III. The supplied PDF and lamina-mode/docs/modal-transform.md establish alpha=-i*kappa, exp(-i*omega*t), and positive in-plane propagation exp(i*beta*u). This change uses those physical sign relations; conformal position alone is insufficient when exterior order, material assumptions or solve direction change.

## Goals / Non-Goals

**Goals:** Supply internal physical point descriptions suitable for the next public API change, with auditable radiation references and explicit interpretation limits. Preserve all existing numerical contracts.

**Non-Goals:** Public solver exports, automatic choice of full-surface search geometry, branch following, global mode identifiers, classification for arbitrary complex/magnetic exteriors, a general causal-sheet construction, fields, normalisation or mode volumes. Existing complex/magnetic problems remain searchable and receive partial descriptions with unsupported named classifications.

## Decisions

### 1. Verify radiation at a real anchor; do not infer it from a coordinate sign

Add a distinct problem-owned domain construction entry point, conceptually domain_with_radiation_reference. It takes the unchanged rectangle, a supplied conformal reference strictly inside it, the real physical value of the varying coordinate, explicit ordered per-side choices, finite comparison/cutoff tolerances, existing chart settings and the required analytic assumption. Keep the present domain constructor unchanged for unspecified radiation provenance. The supplied reference determines the coordinate branch through the existing explicit-coordinate constructor.

For this verification path require positive real k0, real beta, exactly real positive exterior epsilon and mu=1 at the anchor, and no exterior cutoff within tolerance. For spectral problems these exteriors remain constant by the existing construction contract. Let q_j=epsilon_j*k0²-beta². The requested outgoing/decaying value is +sqrt(q_j) for an open channel and +i*sqrt(-q_j) for a closed channel; incoming/growing selects its negative. Compare the existing map's ordered kappa values against those values and verify both physical coordinates. Report side-specific mismatch, unsupported material, invalid reference, cutoff or existing chart failures before discovery. Never move z by i*pi or adjust a kappa sign to make a request pass. A symmetric map may not represent every mixed request; report the mismatch.

Record an opaque verified anchor descriptor (physical and conformal coordinates, material values, requested branches, observed kappa values and verification settings) on DomainProvenance. Constructor-only creation prevents callers from manufacturing verified evidence. Clones, halos, plans, discovery, refinement and point descriptions retain it. Existing domains carry Unspecified. An all-outgoing anchor establishes a particular analytic radiation continuation on that chart, not a proof of causal radiation on every possible complex contour or material holomorphy. The reference must lie inside this chart; transporting anchors between disjoint charts is future continuation work.

Alternative rejected: use Re(kappa)>0 at every complex-frequency point. That can switch analytic sheets and cannot establish a frequency-QNM radiation condition. Automatic anchor inversion and chart generation are deferred to future API conveniences, not required for this explicit construction path.

### 2. Make physical descriptions a separate fallible layer

Add a private physical interpretation module with a common point-data representation and adapters for a problem/domain/AcceptedRoot and for batch accepted observations. Validate exact originating problem and domain before evaluation. Use the chart's analytic coordinate and exterior data at the accepted z, preserving polarization, solve direction, root evidence, source observation identity and any radiation descriptor. A local description is not a globally unique branch ID.

Expose raw finite k0, beta and ordered kappa/alpha alongside independent interpretation components. A shared evaluator supplies the single-chart path; batch adapters reuse valid signatures and preserve existing signature errors, with no hidden backend fallback. Material values needed for interpretation come from the same bound physical context. These are snapshots under the existing immutable-material-behaviour obligation, not a new material freeze.

Return explicit classification statuses: classified, boundary/ambiguous, unsupported with reason, or unavailable with the source evaluation error. Low-level invalid configuration, non-finite point data and foreign ownership are errors. Unsupported named classification does not suppress available decay and attenuation information. Do not modify AcceptedRoot, discovery estimates, duplicate membership, coverage or refinement eligibility. Batch interpretation produces one result per accepted observation and continues after an interpretation failure; non-accepted candidates remain solely in their existing outcome reports.

Alternative rejected: extend numerical acceptance with a requirement for an outgoing or named type. Searches intentionally include source and mixed-sheet modes; numerical validity and physical interpretation answer different questions.

### 3. Use separate physical components and explicit tolerance bands

For each side, Re(alpha)=Im(kappa) determines outward field decay (>0), growth (<0), or a boundary band. Im(beta)>0 indicates attenuation along positive u; Im(beta)<0 amplification. Im(k0)<0 indicates temporal decay; Im(k0)>0 temporal growth. Report the solve-relevant attenuation separately; a fixed-real coordinate is real by construction, not an independently computed lifetime or propagation loss.

For fixed positive real k0 with positive real nonmagnetic dielectric exteriors, use the dimensionless outward flux coefficient Re(kappa/(k0*f)), where f=mu for TE and epsilon for TM, with the positive unit-system prefactor suppressed. The left wave's stack-normal sign is reversed but its outward normal also reverses, so this coefficient uses the same sign on both sides. It is a local directional coefficient per squared primary-field amplitude, not integrated modal power or a field-normalised leakage rate. For complex frequency or unsupported exteriors, flux direction is Unsupported; keep the raw kappa rather than call Re(kappa) a universal energy-flow test.

Use an explicit validated InterpretationConfig: finite nonnegative absolute and relative bands for dimensional wavenumber quantities and dimensionless flux coefficients, with a positive cutoff band. For a wavenumber sign use abs_tol+rel_tol*max(|k0|,|beta|,|kappa_L|,|kappa_R|); use a separate coefficient scale for flux. Never apply an implicit dimensionless floor of one to dimensional values. Record settings and boundary reasons. These bands are numerical conventions, not error bounds on roots; signs near their bands remain ambiguous. Exact material assumptions are tested as material assumptions, not softened by sign tolerances into support for arbitrary lossy exteriors.

### 4. Derive named types from physical signs under the paper's assumptions

Named types require fixed real positive k0, Re(beta) above its positive-direction band, and positive real nonmagnetic exterior dielectric constants. Interior material may be dispersive, absorptive or active under the existing analytic assumption. Use actual left/right labels throughout. The following translates Figure 2 into signs of Re(alpha), with O/I denoting outward/inward flux coefficients:

| Exterior field behaviour (left, right) | Flux (left, right) | Spatial behaviour | Named type |
| --- | --- | --- | --- |
| growth, decay | O, I | attenuation | left leaky |
| decay, growth | O, I | amplification | right source |
| decay, growth | I, O | attenuation | right leaky |
| growth, decay | I, O | amplification | left source |
| growth, growth | I, I | amplification | two-sided source |
| decay, decay | I, I | attenuation | absorptive |
| decay, decay | O, O | amplification | gain |
| growth, growth | O, O | attenuation | two-sided leaky |

The sign products obey Im(alpha_j²)=Im(beta²) for real exterior epsilon and frequency. Treat inconsistent evidence as ambiguous rather than overriding a component. A strict sector label requires all relevant signs outside their bands. For a real beta above both light lines with both alphas decaying and flux in its zero band, report Guided. Other real-axis cases, improper growing solutions, cutoff neighbourhoods and unresolved signs get boundary reasons rather than a forced guided/leaky label. Negative or near-zero Re(beta) gets an explicit unsupported direction reason in this first change; reversing the physical propagation direction is not silently performed.

The table is independent of which exterior has higher index. Test reversed ordering with left/right label permutation; the symmetric limit has linked exterior signs and only its represented subset of the table. Do not require asymmetric canonical bounds to interpret an explicit symmetric-domain root. Test all eight strict sectors from physically consistent sign data independently of locating an actual device root in each sector; complement them with real slab and cavity root regressions.

Alternative rejected: classify by Im(z) quadrants or call every positive Im(beta) point lossy. Both discard the exterior sign information needed to distinguish leakage, absorption and source behaviour.

### 5. Keep frequency-QNM qualification narrower than numerical acceptance

For spectral points report temporal behaviour, exterior field growth/decay and anchor provenance. The paper's named type is UnsupportedSolveDirection. A decaying root with an all-outgoing verified anchor is an outgoing-continuation frequency-mode candidate under the chart and material assumptions; it is not a normalised QNM or a completeness certificate. An unspecified anchor remains unspecified even if both endpoint Re(kappa) signs look outgoing. Non-all-outgoing verified choices remain explicit alternatives rather than receiving the outgoing qualification.

Test the existing vacuum/index-two/vacuum pole k0=(pi-i*ln(3))/(2*d) using a finite explicit rectangle containing both its z and a real-frequency outgoing anchor at Im(z)=-pi/2. This demonstrates outgoing continuation into the omitted canonical sector. The same physical k0 at z+i*pi has incoming signs and must fail an outgoing-anchor request or remain differently qualified. Add single-chart and batch adapters without changing their source reports.

### 6. Update the guide and hand off to the public point API

Extend docs/modal-transform.md, included in crate-level rustdoc, with the support matrix, component sign definitions, tolerance interpretation, real-anchor requirements, all eight strict types, guided/boundary cases and private usage examples. Update README and roadmap progress without claiming public exports. The next change exposes the resulting physical-input API and executable examples. Fields/derivatives, branch continuation, QNM normalisation, feature sensitivities and lamina-design optimisation retain their documented order and crate ownership.

## Risks / Trade-offs

- [A verified anchor is mistaken for a universal radiation certificate] -> Record exactly what was verified; document analytic chart continuation and material assumptions separately.
- [The anchor is outside a small root-isolation chart] -> Require an explicit larger valid domain containing it; disjoint-chart radiation transport is deferred and never inferred.
- [Strict labels hide borderline modes] -> Expose raw quantities, thresholds and boundary reasons; keep acceptance separate.
- [A local flux coefficient is mistaken for modal power] -> Name and document its amplitude convention and withhold it for complex-frequency problems.
- [Interpretation breaks batch ownership or grouping] -> Use separate result adapters and regression comparisons against unchanged source reports.
- [All eight labels appear tested only by duplicating the implementation] -> Build independent alpha²/beta²-consistent fixtures and compare to the paper's sign tuples, plus independent TE/TM device roots.

## Migration Plan

Land the radiation-reference constructor and provenance without changing existing domain calls. Add point-data and interpretation adapters, then the supported classifier and physical regressions. Publish the support matrix and updated roadmap with each stage. Run the 240-test baseline plus new focused tests, all-target tests, doctests, clippy, rustdoc, changed-file formatting and strict OpenSpec validation. No dependency or core changes are expected. All new interfaces remain private and opt-in, so existing arbitrary-sheet searches retain their behaviour if interpretation is unused.

## Implementation verification (2026-10-19)

- Added the opt-in real-anchor constructor and opaque verification evidence in problem/radiation.rs. Both solve directions check all four requested per-side combinations for asymmetric open, closed and mixed channels. Tests also cover invalid settings/references, material restrictions, cutoff rejection, symmetric mixed-sign mismatch and coordinate-sign independence. Requested geometry and existing unspecified-provenance construction remain unchanged.
- Added problem/physical.rs with borrowed accepted-root descriptions, finite physical-data snapshots, recorded component bands, supported flux coefficients, eight strict named sectors, guided/boundary handling and explicit unsupported reasons. Single-chart and batch adapters agree on independently engineered TE/TM guided slabs and normal-incidence cavity roots. No numerical acceptance, candidate eligibility or grouping logic changed.
- Outgoing TE/TM and incoming TE cavity continuations use independently known poles and domains containing real anchors. Their frequency descriptions retain outgoing versus other verification and temporal decay versus growth. Unsupported spectral named labels, arbitrary radiation provenance, complex/magnetic exteriors and nonpositive propagation direction remain explicit.
- Failure-injection regressions verify ownership checks before chart evaluation, retained root and original chart-error cause after acceptance, borrowed batch signature-error identity, continued interpretation of other observations, and unchanged group/discovery/exclusion evidence. Budget-limited reports retain unresolved candidates rather than promoting them through interpretation.
- All 250 locked offline all-target tests pass (240 before this change), with no ignored tests. Doctests pass (zero executable doctests); no-dependency rustdoc builds without warnings. Local documentation links/anchors and changed-file rustfmt checks pass. Clippy succeeds with the same eight pre-existing warnings and no new ones. Strict OpenSpec validation passes.
- No core implementation, dependency, lockfile, map formula or canonical-coverage geometry changed. Existing unrelated working-tree edits were preserved. README and the included public guide describe the exact support matrix, real-anchor private usage and the public point-solver API as the next milestone; later roadmap ordering and crate ownership are retained.
