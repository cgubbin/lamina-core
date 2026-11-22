# Spec Delta

## MODIFIED Requirements

### Requirement: Constant exterior spectral compatibility
Constant-exterior conformal spectral construction SHALL use exterior materials established as constant by construction, with the same constitutive values in the optical stack and spectral map. It SHALL NOT infer constancy from samples or freeze dispersive exteriors. This restriction SHALL NOT apply to a direct-frequency strategy that evaluates the supplied analytical stack throughout its domain.

#### Scenario: Dispersive finite layers
- **WHEN** a constant-exterior conformal spectral problem has explicit constant exteriors and analytical dispersive finite layers
- **THEN** the exterior map uses those constants while optical evaluation retains the finite layers' spectral dependence and derivatives

#### Scenario: Unproven exterior constancy
- **WHEN** constant-exterior conformal spectral construction is requested with arbitrary exterior material handles lacking constant construction provenance
- **THEN** construction is rejected or unavailable without silently sampling and freezing their values

#### Scenario: Stack-based dispersive spectral construction
- **WHEN** a direct-frequency spectral problem is constructed from an existing analytical core stack with dispersive exteriors
- **THEN** both exterior branch evaluation and checked optical evaluation use that same stack without requiring separate replacement exterior descriptors

## ADDED Requirements

### Requirement: Native chart coordinate provenance
Domain and result provenance SHALL distinguish dimensionless conformal coordinates from physical frequency-chart coordinates in inverse centimetres. Requested bounds, evaluation margins, root enclosures and native-coordinate tolerances SHALL retain that distinction without changing their physical interpretation or implying periodic sheet coverage for direct frequency charts.

#### Scenario: Direct frequency enclosure
- **WHEN** a direct-frequency search reports an enclosure and a refined native-coordinate correction
- **THEN** both are expressed in the recorded vacuum-angular-wavenumber chart coordinates with frequency units rather than being described as dimensionless conformal quantities

#### Scenario: Existing conformal result
- **WHEN** the constant-exterior conformal or fixed-frequency in-plane strategy is used
- **THEN** its existing dimensionless coordinate conventions, radiation references and finite coverage limitations remain preserved
