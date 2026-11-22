# Spec Delta

## Purpose

Expose physical first-order characteristic derivatives and conditional implicit root slopes at accepted modal points while preserving exterior branches, actual constitutive dispersion, parameter units and separate failure evidence.

## ADDED Requirements

### Requirement: Physical characteristic partial derivatives
The public interface SHALL evaluate the checked characteristic and its requested first partial derivatives with respect to physical vacuum angular wavenumber, physical in-plane angular wavenumber and selected finite-layer thicknesses. Each partial SHALL identify its parameter, held-fixed quantities and units separately from a native-chart derivative.

#### Scenario: Frequency and propagation partials
- **WHEN** both physical coordinate partials are requested at one accepted point
- **THEN** one holds beta fixed while varying k0 and the other holds k0 fixed while varying beta, without conflating either with dF/dz

#### Scenario: Thickness partial
- **WHEN** a valid finite-layer thickness parameter is selected
- **THEN** the returned partial varies only that thickness, preserving the other physical coordinates, geometry and material models

### Requirement: Explicit local derivative assumptions
Physical derivative evaluation SHALL require an explicit local analyticity declaration covering the requested coordinate and parameter variations near the accepted physical point. It SHALL retain the source domain assumptions without claiming they already certify variations not covered by the original search.

#### Scenario: Fixed-frequency source
- **WHEN** a frequency partial is requested for a point found at fixed real frequency
- **THEN** the operation requires the additional local constitutive and characteristic assumption rather than inferring complex-frequency validity from that search alone

### Requirement: Branch-consistent constitutive derivatives
Coordinate partials SHALL include actual exterior and interior permittivity and permeability dispersion as applicable and preserve each accepted exterior wavevector branch. A cutoff or non-finite derivative SHALL remain a diagnostic. Numerical continuation quadrature SHALL NOT be differentiated as a substitute for the analytic physical derivative.

#### Scenario: Dispersive magnetic exterior
- **WHEN** exterior permittivity and permeability both vary with k0
- **THEN** the frequency partial includes both contributions to its exterior wavevector derivative on the accepted branch

#### Scenario: Opposite exterior branch
- **WHEN** the source point has a negated exterior branch
- **THEN** that exterior's wavevector derivative has the corresponding sign instead of being reconstructed from an endpoint principal root

### Requirement: Canonical and thickness unit provenance
Physical coordinate partials SHALL use the documented inverse-centimetre coordinates. A thickness partial SHALL identify the selected layer and differentiate with respect to its stored caller-facing thickness unit, preserving core's unit-conversion semantics. Equivalent devices expressed in different thickness units SHALL obey the corresponding chain rule.

#### Scenario: Nanometre versus centimetre thickness
- **WHEN** equivalent physical stacks express the selected thickness in nanometres and centimetres
- **THEN** their characteristic values agree and their labelled derivatives differ by the documented unit-conversion factor

### Requirement: Guarded simple-root implicit slopes
The interface SHALL provide implicit first-order root slopes only when the requested denominator derivative passes validated finite nondegeneracy checks. It SHALL retain the contributing partials, units and source acceptance evidence. An unavailable slope SHALL NOT erase a valid gradient or be labelled an exceptional point, ZGV point or energy velocity.

#### Scenario: Regular spectral slope
- **WHEN** F_k0 is usable at an accepted spectral point
- **THEN** the reported derivative dk0/dbeta equals -F_beta/F_k0 under the declared held-fixed convention

#### Scenario: Degenerate denominator
- **WHEN** a required denominator is non-finite, zero or inside the configured nondegeneracy band
- **THEN** the slope is explicitly unavailable with a diagnostic rather than infinity, a fabricated finite value or a branch identity claim

### Requirement: Independent owned derivative evidence
Derivative results SHALL be owned and retain the source point, selected parameters, units, local declaration, evaluation backend and checked value/derivative outcomes. Derivative and slope failures SHALL preserve available causes without changing acceptance, classification, original reports or coverage. No backend fallback or hidden root correction SHALL occur.

#### Scenario: Stored accepted point
- **WHEN** an accepted point outlives its original execution objects
- **THEN** its local derivatives remain requestable using its source context and an explicit supported evaluation backend

#### Scenario: Partial evaluation failure
- **WHEN** one requested partial or derived slope cannot be evaluated
- **THEN** the result identifies that failure and retains available valid partials and source acceptance evidence

### Requirement: Independent derivative verification
Physical partials SHALL be checked against independent analytic derivatives, finite-difference step sweeps on the same branches and native-coordinate chain rules on regular fixtures. Supported implicit slopes SHALL be checked against independently perturbed roots. Verification SHALL include both solve directions and a dispersive finite-layer device.

#### Scenario: Conformal chain rule
- **WHEN** a supported conformal source chart represents the accepted physical point
- **THEN** its native derivative agrees with the appropriate physical partial multiplied by the chart-coordinate derivative within documented tolerance

#### Scenario: Layer sensitivity
- **WHEN** a regular film or slab thickness is perturbed independently
- **THEN** the measured nearby root shift agrees with the first-order implicit thickness slope over a documented convergence range
