# modal-candidate-validation Specification

## Purpose

Refine isolated modal estimates on their bound analytic branch and distinguish numerically accepted roots from unresolved candidates while preserving discovery evidence.

## Requirements

### Requirement: Bound candidate provenance
Refinement SHALL consume estimates only with their originating physical problem and domain. It SHALL reject foreign problem or domain combinations before optical evaluation and preserve the original estimate metadata and branch provenance.

#### Scenario: Different domain of the same problem
- **WHEN** a search result is submitted for refinement with another domain of its problem
- **THEN** refinement rejects the combination before evaluating candidates

#### Scenario: Same bound branch
- **WHEN** a candidate is refined using its originating domain
- **THEN** all evaluations and physical coordinate projections use that domain's existing analytic reference and exterior branch

### Requirement: Isolated zero eligibility
Only estimates identified as one zero with multiplicity one SHALL be eligible for individual-root refinement. Clusters, unsupported targets, inconsistent metadata, and invalid enclosures SHALL produce explicit non-accepted outcomes.

#### Scenario: Cluster centroid
- **WHEN** an estimate represents multiple enclosed roots or is marked as a cluster
- **THEN** it remains unresolved with its original multiplicity and is not converted into one accepted root

#### Scenario: Invalid enclosure
- **WHEN** an estimate's enclosure is invalid or extends outside the requested search domain
- **THEN** it is not refined or accepted as a root of that domain

### Requirement: Bounded branch-preserving refinement
Refinement SHALL evaluate candidates only within both their original enclosures and the requested search rectangle on the existing branch. The surrounding evaluation halo SHALL NOT extend the root-acceptance region. Inability to make admissible progress SHALL remain an explicit non-accepted outcome.

#### Scenario: Outward step
- **WHEN** a proposed refinement step exits the candidate enclosure or searched rectangle
- **THEN** it is reduced or rejected before optical evaluation without clamping the point or switching charts

#### Scenario: Boundary limitation
- **WHEN** a candidate cannot be refined to a strictly interior accepted point within its budgets
- **THEN** the outcome records the limitation instead of accepting a point in the evaluation halo

### Requirement: Explicit numerical acceptance
An accepted isolated root SHALL satisfy configured finite residual and local correction tolerances at a freshly checked point, with a finite nonzero derivative and a branch-consistent physical coordinate. Configuration and convergence evidence SHALL remain available. Acceptance SHALL be described as numerical and conditional on the caller's analytic-domain assumption.

#### Scenario: Small residual with large correction
- **WHEN** a finite candidate has a small characteristic residual but its local correction exceeds the configured tolerance
- **THEN** residual size alone does not permit acceptance

#### Scenario: Derivative degeneracy
- **WHEN** the characteristic derivative is zero or the local correction is unrepresentable
- **THEN** no isolated simple root is accepted on that evidence

### Requirement: Explicit refinement outcomes
Refinement SHALL record one outcome per input estimate, including nonconvergence, exhausted budgets, unresolved clusters, and checked evaluation failures with their original causes. It SHALL NOT replace failure with an accepted root, omit the failed estimate, or reinterpret it as a successful empty search.

#### Scenario: Optical evaluation failure
- **WHEN** checked evaluation fails for a candidate or trial point
- **THEN** that candidate records the source error and remains non-accepted while the original search result is retained

#### Scenario: Bounded work
- **WHEN** refinement exhausts its iteration or trial-evaluation budget without meeting acceptance criteria
- **THEN** it returns an explicit budget outcome with the available numerical evidence

### Requirement: Preserve discovery and coverage evidence
The refinement report SHALL retain the complete originating search result and all candidate outcomes. Root acceptance SHALL NOT change estimates, multiplicities, unresolved leaves, run summary, termination, or caller assumptions, nor imply complete domain coverage or mathematical isolation certification.

#### Scenario: Accepted root from incomplete discovery
- **WHEN** a candidate converges from a budget-limited search
- **THEN** it can be numerically accepted while the original budget termination and unresolved coverage remain visible
