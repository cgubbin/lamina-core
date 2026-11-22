# Spec Delta

## Purpose

Provide a dependable branch-preserving determinant and argument-principle search pipeline over one explicitly bounded analytic conformal chart.

## ADDED Requirements

### Requirement: Checked branch-preserving dispersion
Optical dispersion evaluation SHALL use checked characteristic evaluation with the chart-supplied canonical coordinates and exterior wavevector jets unchanged. It SHALL return the characteristic and its conformal derivative from the same seeded evaluation for either search direction.

#### Scenario: Complex in-plane or spectral search
- **WHEN** a chart evaluates a conformal point and its jets
- **THEN** the optical evaluation preserves its exterior radiation sheet and returns the value and derivative on that sheet without recomputing principal exterior roots

### Requirement: Preserve evaluation failures
Dispersion and single-chart search SHALL propagate checked input, backend, and output failures with their source information. They SHALL NOT convert these failures into a successful empty search or a fabricated mode estimate.

#### Scenario: Non-finite optical derivative
- **WHEN** a finite characteristic value carries a non-finite derivative
- **THEN** dispersion evaluation fails and the search retains that checked failure in its error chain

#### Scenario: Evaluation failure during adaptive search
- **WHEN** the argument-principle engine returns an error caused by optical evaluation
- **THEN** the modal search reports that error with the optical source preserved and makes no completeness claim

### Requirement: Unchanged contour coverage with validated evaluation margin
The search SHALL preserve its requested conformal contour rectangle while evaluating boundary roundoff inside a separately validated surrounding chart region. It SHALL preserve the chart reference sheet and reject an unsafe evaluation margin rather than clamping points or crossing a coordinate branch point.

#### Scenario: Rounded boundary evaluation
- **WHEN** adaptive contour arithmetic evaluates a point within the validated numerical margin around a requested boundary
- **THEN** evaluation succeeds on the original analytic sheet without moving the contour or clamping the point

#### Scenario: Unsafe chart extension
- **WHEN** constructing the numerical margin would include a coordinate branch point or non-finite bound
- **THEN** search construction returns a domain error before starting the argument-principle search

#### Scenario: Genuinely exterior point
- **WHEN** evaluation is requested outside the validated surrounding region
- **THEN** the chart continues to reject the point

### Requirement: Preserve search diagnostics
Successful single-chart results SHALL retain the argument-principle estimates, multiplicities, unresolved leaves, run summary, and termination reason. Numerical success SHALL NOT imply accepted individual roots, pole-free material domains, or complete coverage of excluded regions.

#### Scenario: Budget-limited search
- **WHEN** the engine returns a successful result with unresolved cells and a budget termination
- **THEN** the modal result preserves those cells and termination rather than presenting the returned estimates as an exhaustive mode set

#### Scenario: Cluster estimate
- **WHEN** an engine estimate represents multiple enclosed targets
- **THEN** its multiplicity and estimate metadata remain available and it is documented as an unresolved estimate rather than a validated individual mode

### Requirement: Consistent coordinate projection
Single-chart coordinate projection SHALL agree between scalar and jet evaluation paths within documented numerical tolerance while preserving the selected coordinate sheet. Coordinate conversion failures SHALL remain errors.

#### Scenario: Equivalent chart coordinate evaluations
- **WHEN** equivalent scalar and jet evaluation paths differ only within their stated numerical tolerance
- **THEN** both evaluations represent the same physical coordinate within that tolerance

#### Scenario: Opposite coordinate sheets
- **WHEN** a nonzero coordinate is evaluated on the two opposite analytic sheets
- **THEN** the evaluations preserve their opposite signs rather than merging them through a loose coordinate tolerance
