# Spec Delta

## Purpose

Recover homogeneous electromagnetic mode fields from accepted modal points using their source device and exterior branches, with explicit amplitude gauge, physical sampling and failure evidence independent of root acceptance and QNM normalization.

## ADDED Requirements

### Requirement: Normalization-independent homogeneous recovery
Homogeneous electric and magnetic field recovery SHALL be available without computing an energy integral, a QNM normalization factor or a driven incident-wave solution. It SHALL retain a documented nonzero arbitrary complex amplitude and phase gauge. The existing normalized core mode path SHALL preserve its prior behavior.

#### Scenario: Normalization unavailable
- **WHEN** homogeneous reconstruction is usable but an energy integral or normalization scale is zero, undefined or unavailable
- **THEN** unscaled field recovery remains possible independently of that normalization failure

#### Scenario: Existing normalized consumer
- **WHEN** a caller uses the pre-existing normalized core mode entry point
- **THEN** it retains its established normalization semantics and diagnostics

### Requirement: Accepted point and branch association
Public mode-field recovery SHALL use an accepted point's exact physical coordinates, polarization, source stack and ordered exterior wavevectors. An explicit supported recovery backend SHALL be recorded separately from the source backend. Missing branch evidence, inconsistent source data or foreign associations SHALL fail before reconstruction without changing the accepted point.

#### Scenario: Stored point
- **WHEN** an owned accepted point outlives its original report and problem
- **THEN** its fields can be requested with an explicit supported backend using its retained source context and branch data

#### Scenario: Unsupported named classification
- **WHEN** an accepted metallic or spectral point cannot receive a supported named dielectric type
- **THEN** that classification limit alone does not prevent otherwise valid homogeneous field recovery

### Requirement: Core spatial and field conventions
Field requests SHALL reuse core device, sampling, position, unit and field conventions. Returned electric and magnetic samples SHALL preserve requested order, resolved positions and interface-side identity. Distinct sides of one geometric interface SHALL NOT be collapsed into a single sample or forced to have equal discontinuous components.

#### Scenario: Opposite interface samples
- **WHEN** a request samples both sides of the same interface coordinate
- **THEN** both samples remain identifiable, with tangential continuity and appropriate constitutive normal-component behavior evaluated separately

#### Scenario: Exterior and finite-layer sampling
- **WHEN** valid finite-distance exterior samples and finite-layer samples are requested
- **THEN** the result uses the source geometry and units and retains those resolved sampling locations alongside the field components

### Requirement: Owned fields with gauge and residual evidence
A successful field result SHALL own its sampled values and retain source-point provenance, source and recovery backends, amplitude/phase gauge and homogeneous reconstruction residual evidence. Rescaling SHALL apply one consistent complex factor to both electric and magnetic fields. Recovery SHALL NOT silently move or re-refine the source root.

#### Scenario: Approximate accepted root
- **WHEN** fields are reconstructed at a numerically accepted point with a finite residual
- **THEN** the result retains that point and the reconstruction defect rather than claiming an exact null solution or replacing the root

#### Scenario: Complex rescaling
- **WHEN** a valid nonzero complex amplitude multiplier is applied
- **THEN** every electric and magnetic component changes by that same multiplier and the recorded gauge changes consistently

### Requirement: Separate recovery and sampling failures
Recovery SHALL report structured failures for invalid sampling or gauge, unavailable branch data, unusable or degenerate homogeneous seeds, backend failures and non-finite sampled values. Such failures SHALL retain their available causes and the originating accepted point without changing discovery, acceptance, groups or coverage.

#### Scenario: Invalid recovery request
- **WHEN** a sampling location, gauge multiplier or reconstruction setting is invalid
- **THEN** the request fails diagnostically without mutating the point or substituting another backend or driven field

#### Scenario: Growing leaky exterior
- **WHEN** a valid leaky mode is sampled at finite exterior distances
- **THEN** finite growth is preserved rather than clamped, while unrepresentable values produce a diagnostic without asserting an integrable or normalized field

### Requirement: Independent physical field validation
Public recovery SHALL be validated on guided dielectric modes, a metal-interface plasmon in both complex cases, both coupled modes of a symmetric finite metal film and a leaky spectral cavity. Validation SHALL compare physical profiles, interface conditions and exterior behavior under a common gauge, not only root residuals or total field norms.

#### Scenario: Interface plasmon profile
- **WHEN** a benchmark plasmon field is recovered
- **THEN** its relative electric and magnetic components, interface conditions and exterior exponential factors agree with the independent TM solution under the same gauge

#### Scenario: Coupled film modes
- **WHEN** the two supported finite-film plasmon modes are recovered
- **THEN** both remain distinct and their magnetic-field parity and interface matching agree with the independent symmetric-film equations
