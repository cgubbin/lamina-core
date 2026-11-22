# Spec Delta

## ADDED Requirements

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
