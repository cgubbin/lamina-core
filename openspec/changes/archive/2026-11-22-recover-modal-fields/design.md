# Design

## Context

See proposal.md for motivation. Apply after the benchmark change. Core already provides `ComplexPlaneEvaluator::retain`, retained homogeneous boundary seeds, layer/exterior reconstruction, `FieldSampling`, resolved interface-side positions, `SpatialResponse` and electromagnetic field containers. These are the primitives to reuse.

A material implementation detail prevents a direct adapter today: `ComplexPlaneState::mode()` calls `ComplexPlaneMode::new`, which obtains `modal_boundary_solution()`, computes `raw_qnm_normalisation_unchecked`, applies `normalisation_scale`, and only then exposes fields. Thus the existing public complex-mode path is normalized even though its lower-level seed and reconstruction machinery are unscaled. Several reconstruction traits and response types are crate-private. Mode must not bypass this by duplicating field formulas or relying on inaccessible backend workspaces.

Mode's owned `ModePoint` retains an accepted root, source stack/context and a fallible ordered physical branch signature. It records a backend diagnostic type name, not an executable backend instance or retained field workspace. Metal named classification and explicit-complex radiation provenance can remain unsupported/unspecified even for the verified benchmark SPP.

## Goals / Non-Goals

**Goals:** Recover owned arbitrary-amplitude E/H samples from accepted points for all three existing strategies, preserve physical branches and interface-side identity, and expose reconstruction quality/errors without changing numerical acceptance.

**Non-Goals:** No new QNM normalization formula, mode volume, energy normalization guarantee, driven excitation, field derivatives, mode classification extension or branch tracker. Preserve the existing normalized core mode API; do not reinterpret its output as arbitrary-amplitude seed data.

## Decisions

### 1. Add an unscaled core operation, sharing existing reconstruction

Separate creation of the homogeneous seed from normalization inside core. Add an explicit unscaled view/entry point that reuses `modal_boundary_solution`, modal boundary waves and `FieldSamplingContext::reconstruct_from_boundary_waves`, without invoking bilinear energy aggregation or a normalization scale. Keep the existing normalized construction and `seed_normalisation` semantics intact; both paths share reconstruction helpers rather than separate Maxwell implementations.

Provide a public high-level owned field operation/capability for the existing built-in backends, so lamina-mode does not need to name private workspace traits. Implement the capability within core over the existing retained reconstruction bounds; do not expose unrelated backend internals or create a parallel stack/field representation. It accepts source coordinates, ordered exterior data, source stack and core sampling, and returns owned core sample values plus a public seed/defect summary. A narrow public built-in-backend capability bound is acceptable if needed to hide private bounds; no custom-backend plugin framework is part of this change.

Test core's new path on a zero-finite-layer interface and on a deliberately unusable normalization case: unscaled reconstruction must not require normalization to succeed. Existing normalized field/constitutive tests must continue to pass. Only the justified additive core boundary and shared factoring are authorized by this proposal.

Alternative rejected: call `state.mode()` then undo its scale. That still depends on normalization existing and would fail precisely where the new operation is needed.

### 2. Recover from a stored point using an explicit backend value

Offer a public operation such as `point.recover_fields_with(backend, sampling, settings)`. Exact naming is not a numerical contract. It uses the point's source stack, polarization, physical k0/beta and ordered kappa values, with finite/squared-relation consistency checks against actual constitutive values. Do not rerun chart discovery, refine the point, freeze materials or compute exterior principal roots.

An explicit supported core backend value is required; a diagnostic type-name string cannot recreate backend state. The result records original and recovery backends separately. Explicitly requesting another supported backend is allowed for physical comparison at the same point; automatic fallback is not. Use core's constant scalar evaluation for primal field recovery; no dummy derivatives or material-parameter frameworks are exposed to callers. Reject absent signatures, foreign associations in any problem-bound overload and detectable inconsistent or non-finite material/exterior data before reconstruction. Arbitrary interior mutation cannot be certified absent and remains outside the existing stable-model lifetime contract.

The public operation obtains a retained seed locally and samples it, returning owned arrays/context rather than a structure borrowing a locally created evaluator or workspace. A point that has outlived its report remains usable. Persistent caching and retained reusable field workspaces can be added later only if justified; they are not required here.

### 3. Treat gauge and reconstruction defect as evidence

Use core's documented backend seed gauge and allow an explicit finite nonzero complex scalar multiplier. Record the seed gauge identifier and the applied multiplier, and scale E and H together. Invalid/vanishing reference components do not justify silently changing to another backend or normalization rule. Where comparisons need a common gauge, anchor both solutions to the same finite nonzero primary field component and record that comparison scale.

At a numerically accepted approximate root the outgoing boundary seed has a finite defect. Return the original numerical acceptance evidence, the checked characteristic at recovery and a core homogeneous boundary-defect measure with its scale/units. Validate finite reconstruction data and explicit residual settings; a defect exceeding the requested reconstruction limit fails this postprocessing stage without revoking the root. Do not claim an exact null solution, normalised power or SI amplitude unless established by the existing core units and the recorded arbitrary gauge.

Alternative rejected: accepting only an exactly zero determinant. Floating-point roots are approximate; quality must be reported and bounded explicitly.

### 4. Preserve core spatial conventions and physical components

Accept `FieldSampling`/`LayerSampling`/`ExteriorSampling` and return owned core electromagnetic values with resolved sampling metadata. Preserve request order and both interface sides even at the same geometric position. Validate tangential E/H continuity and normal D/B continuity under source-free interface conditions; normal E/H can be discontinuous and must not be forced equal. Reuse core unit/axis/impedance conventions and document them from the existing API.

Sample exterior fields only at requested finite distances. Guided/SPP profiles should decay on the selected bound branches; leaky frequency profiles may grow. Preserve finite growth and diagnose overflow rather than clamping it or applying a naive L2 normalization. Add public diagnostic stages/causes for invalid requests, seed degeneracy, reconstruction, sampling and non-finite outputs. Original report, groups, multiplicities, coverage and interpretation remain unchanged.

### 5. Validate profiles with independent fixtures

Use the preserved interface benchmark in both complex cases and independent outgoing/decay exponentials with the correct TM E/H relationships. Add guided dielectric TE/TM slabs to check nonzero and vanishing Cartesian components. Compare both supported core backends up to a single recorded complex gauge, not only scalar intensities.

Add a vacuum / finite Drude metal / vacuum film, initially with kp*d=0.4 (40 nm for kp=1e5 cm^-1), at regular frequencies away from thresholds and degeneracies. Let alpha_j=-i*kappa_j with positive real decay on the bound fixture. Independent magnetic-primary-field equations are

```text
alpha_d/epsilon_d + alpha_m/epsilon_m * tanh(alpha_m*d/2) = 0  (even H)
alpha_d/epsilon_d + alpha_m/epsilon_m * coth(alpha_m*d/2) = 0  (odd H).
```

Locate both roots independently, then use public discovery/acceptance and field recovery. Test magnetic parity about the film centre and all interface conditions after one common gauge alignment. The parity names refer to the primary magnetic component, not an ambiguous generic 'symmetric mode'. Keep these helpers reusable for the derivative stage.

Include an outgoing dielectric cavity frequency pole with finite exterior samples to demonstrate valid growth and failure on unrepresentable sampling. Add a runnable public field-profile example/plot command; optional plotting remains outside numerical tests. Test invalid gauge/sampling, missing/unusable seed and non-finite model failures with the accepted point retained. No field result silently changes radiation qualification.

## Risks / Trade-offs

- [Unscaled API accidentally changes legacy normalized results] -> Additive entry point with shared helpers and existing normalized core regressions.
- [Private backend traits leak through the adapter] -> Core owns the narrow high-level operation; compile external mode examples for both built-in backends.
- [Gauge differences look like physical backend disagreement] -> Compare all components after one documented complex gauge alignment and retain seed metadata.
- [Approximate roots look like exact interface solutions] -> Preserve source residual and reconstruction defect; enforce explicit postprocessing tolerances.
- [Leaky growth is treated as invalid physical behavior] -> Finite-distance sampling preserves growth; only unsupported/unrepresentable evaluation fails.

## Migration Plan

After `preserve-plasmon-benchmarks`, add and test the core unscaled boundary first, then the public mode adapter/ownership, then physical profiles and plots. Run targeted core retained-mode tests plus the core suite as required by the actual core edits, and the current mode locked suite, doctests, clippy, rustdoc, formatting and runnable field example. Preserve pre-existing work and existing normalized semantics. The derivative change follows in the working sequence. Branch continuation and explicit QNM normalization/derived observables remain later; core owns reusable reconstruction and mode owns association with accepted points.

## Implementation finding: homogeneous scattering recovery

Physical profile comparisons exposed cancellation in the existing scattering modal reconstruction: it passed a near-zero incident wave to the driven cut formula and divided by a near-zero resonance denominator. Guided and film fields could violate interface continuity despite an accepted root. The shared homogeneous reconstruction now uses the seed's right outgoing amplitude and projective suffix networks, cancelling the global resonance denominator algebraically. It preserves the backend, gauge and normalized-mode semantics; driven field reconstruction is unchanged. Existing off-root modal tests and new interface/parity checks validate compatibility. This is a necessary core correction within the authorized field recovery work.
