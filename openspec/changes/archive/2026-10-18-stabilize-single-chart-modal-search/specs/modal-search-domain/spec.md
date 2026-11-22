# Spec Delta

## ADDED Requirements

### Requirement: Removable interface continuation
For finite nonzero constitutive factors and regular exterior and interior directional bases, scattering characteristic evaluation SHALL continue through opposite-admittance interface sums when the complete characteristic and carried derivatives are finite and representable.

#### Scenario: Negative-index finite layer
- **WHEN** a finite layer has epsilon and mu equal to minus one between vacuum exteriors at normal incidence
- **THEN** the checked scattering characteristic and its carried derivatives agree with the independent slab characteristic rather than failing because an interface admittance sum vanishes

### Requirement: Removable intermediate-pole continuation
Scattering composition SHALL preserve a usable representation through an intermediate network pole when the composed response has a finite representable analytic continuation. The characteristic SHALL retain its existing outgoing-boundary gauge and derivative semantics.

#### Scenario: Intermediate pole with finite composed response
- **WHEN** an intermediate scattering denominator vanishes but the complete network has a finite analytic continuation
- **THEN** the computed characteristic agrees with that continuation without an artificial all-zero projective collapse

#### Scenario: A true outgoing mode
- **WHEN** the complete scattering denominator is zero and its transmission numerator is nonzero
- **THEN** checked evaluation preserves the finite zero characteristic instead of rejecting the physical mode

### Requirement: Regular scattering compatibility
The scattering representation repair SHALL preserve ordinary amplitudes, powers, retained observables, and characteristic values and derivatives on previously regular well-conditioned inputs. Non-finite unresolved results SHALL continue to produce checked numerical failures.

#### Scenario: Existing driven and retained evaluation
- **WHEN** a regular stack is evaluated before and after the repair
- **THEN** both incident-side responses and retained observables agree within their established tolerances

#### Scenario: Unrepresentable characteristic
- **WHEN** overflow, underflow-induced division, or a remaining directional-basis singularity makes the characteristic non-finite
- **THEN** checked evaluation returns its numerical failure without replacing the characteristic with a projective denominator or silently changing backend
