# modal-search-problem Specification

## Purpose

Bind a physical layered-medium mode problem to compatible conformal search domains, explicit analytic assumptions, and traceable single-chart results.

## Requirements

### Requirement: Coherent physical problem
Problem construction SHALL bind one stack, polarization, search direction, and fixed physical coordinate. Its chart data and checked optical evaluation SHALL use that same physical problem without accepting independently substituted exterior parameters.

#### Scenario: Fixed spectral coordinate
- **WHEN** an in-plane problem is constructed at a supported real spectral coordinate
- **THEN** both exterior constitutive pairs are evaluated from its stack at that coordinate and used in its map and optical evaluation

#### Scenario: Incompatible replacement
- **WHEN** a caller attempts to attach exterior data or a chart from another physical problem through the bound construction path
- **THEN** the combination is rejected or cannot be expressed through that path

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

### Requirement: Bound domain and branch provenance
Each problem-owned domain SHALL retain its requested rectangle, validated evaluation region, and analytic reference selection. It SHALL preserve the corresponding physical branch during evaluation and reject domains belonging to a different problem.

#### Scenario: Reference continuation
- **WHEN** a domain is built with a supported explicit reference coordinate or coordinate-sheet selection
- **THEN** its evaluation halo preserves that reference selection without changing the searched rectangle or reselecting a principal exterior root

#### Scenario: Foreign domain
- **WHEN** execution is requested with a domain constructed for another problem
- **THEN** it fails before optical search or is prevented by the ownership interface

#### Scenario: Invalid construction
- **WHEN** fixed inputs, exterior constitutive values, chart geometry, or its surrounding evaluation region fail their existing validity checks
- **THEN** construction returns a diagnostic error without beginning adaptive search

### Requirement: Explicit analytic domain assumption
Domain construction SHALL require an explicit caller assumption that the characteristic is holomorphic and pole-free on a neighborhood of the closed search domain containing its evaluation region. This assumption SHALL remain distinguishable from coordinate validation and pointwise numerical checks.

#### Scenario: Missing declaration
- **WHEN** a caller supplies only a rectangle and finite optical samples
- **THEN** those inputs do not create an executable bound domain without the explicit analytic-domain declaration

#### Scenario: Direction-specific material validity
- **WHEN** the assumption is recorded for an in-plane or spectral problem
- **THEN** it applies respectively at the fixed spectral point or throughout the mapped spectral domain and includes the materials and characteristic representation

### Requirement: Traceable single-chart execution
Bound execution SHALL associate its results with the problem's polarization, search direction, fixed coordinate, domain, branch reference, and caller assumption. It SHALL preserve existing estimates, multiplicities, leaves, run summary, termination, and checked evaluation failures without claiming accepted modes or complete coverage.

#### Scenario: Successful incomplete execution
- **WHEN** the underlying search succeeds with unresolved clusters or budget termination
- **THEN** the bound result retains those diagnostics together with its problem and domain provenance

#### Scenario: Optical failure
- **WHEN** checked optical evaluation fails through the bound path
- **THEN** the originating error remains identifiable and no successful empty result is substituted

### Requirement: Explicit radiation reference
Domain construction SHALL offer an opt-in real-reference path that verifies requested per-exterior outgoing/decaying or incoming/growing branches for positive real nonmagnetic dielectric exteriors. This SHALL remain separate from selecting the sign of the linear physical coordinate. Invalid, cutoff, unsupported or incompatible references SHALL fail before search without altering the requested geometry.

#### Scenario: All-outgoing real reference
- **WHEN** a supplied interior conformal reference maps to consistent real physical coordinates away from exterior cutoffs and matches the requested outgoing or decaying branches
- **THEN** construction retains the verified reference and continues those exterior data analytically within the domain without pointwise principal-root reselection

#### Scenario: Coordinate sign does not select radiation
- **WHEN** a caller changes only the signed-coordinate branch at fixed conformal reference
- **THEN** that change is not treated as selecting different exterior radiation signs

#### Scenario: Wrong exterior signs
- **WHEN** the supplied reference has the correct squared physical coordinates but one or both exterior wavevectors disagree with the requested radiation branches
- **THEN** construction rejects the mismatch without translating the chart or flipping exterior signs silently

#### Scenario: Unsupported or cutoff reference
- **WHEN** the reference has unsupported exterior constitutive factors, is outside the requested domain, is not real in physical coordinates, or is within a cutoff boundary band
- **THEN** construction reports the specific limitation before adaptive search

### Requirement: Radiation-reference provenance
Domains and their derived results SHALL retain radiation-reference verification separately from coordinate-branch provenance and the caller's analytic-domain assumption. Existing domain construction SHALL remain available with unspecified radiation provenance. A verified reference SHALL NOT certify arbitrary material analyticity, global outgoing-sheet coverage or a pointwise flux sign throughout a complex-frequency domain.

#### Scenario: Halo and plan preservation
- **WHEN** a radiation-referenced domain is cloned, receives its evaluation halo or enters a chart plan
- **THEN** its reference, requested per-side branches and verification evidence remain attached to the same domain identity without reselecting a sheet

#### Scenario: Existing arbitrary-sheet domain
- **WHEN** a domain is built using the existing coordinate-reference path alone
- **THEN** it remains executable and its results explicitly retain unspecified radiation provenance

### Requirement: Native chart coordinate provenance
Domain and result provenance SHALL distinguish dimensionless conformal coordinates from physical frequency-chart coordinates in inverse centimetres. Requested bounds, evaluation margins, root enclosures and native-coordinate tolerances SHALL retain that distinction without changing their physical interpretation or implying periodic sheet coverage for direct frequency charts.

#### Scenario: Direct frequency enclosure
- **WHEN** a direct-frequency search reports an enclosure and a refined native-coordinate correction
- **THEN** both are expressed in the recorded vacuum-angular-wavenumber chart coordinates with frequency units rather than being described as dimensionless conformal quantities

#### Scenario: Existing conformal result
- **WHEN** the constant-exterior conformal or fixed-frequency in-plane strategy is used
- **THEN** its existing dimensionless coordinate conventions, radiation references and finite coverage limitations remain preserved
