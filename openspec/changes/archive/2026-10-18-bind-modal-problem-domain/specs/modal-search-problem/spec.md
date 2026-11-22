# Spec Delta

## Purpose

Bind a physical layered-medium mode problem to compatible conformal search domains, explicit analytic assumptions, and traceable single-chart results.

## ADDED Requirements

### Requirement: Coherent physical problem
Problem construction SHALL bind one stack, polarization, search direction, and fixed physical coordinate. Its chart data and checked optical evaluation SHALL use that same physical problem without accepting independently substituted exterior parameters.

#### Scenario: Fixed spectral coordinate
- **WHEN** an in-plane problem is constructed at a supported real spectral coordinate
- **THEN** both exterior constitutive pairs are evaluated from its stack at that coordinate and used in its map and optical evaluation

#### Scenario: Incompatible replacement
- **WHEN** a caller attempts to attach exterior data or a chart from another physical problem through the bound construction path
- **THEN** the combination is rejected or cannot be expressed through that path

### Requirement: Constant exterior spectral compatibility
Spectral problem construction SHALL use exterior materials established as constant by construction, with the same constitutive values in the optical stack and spectral map. It SHALL NOT infer constancy from finite samples or freeze an arbitrary dispersive exterior.

#### Scenario: Dispersive finite layers
- **WHEN** a spectral problem has explicit constant exteriors and analytic dispersive finite layers
- **THEN** the exterior map uses those constants while optical evaluation retains the finite layers' spectral dependence and derivatives

#### Scenario: Unproven exterior constancy
- **WHEN** spectral construction is requested with arbitrary exterior material handles lacking constant construction provenance
- **THEN** construction is rejected or unavailable without silently sampling and freezing their values

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
