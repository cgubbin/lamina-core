# Spec Delta

## Purpose

Provide reliable pointwise characteristic evaluation and an explicit domain contract for downstream complex-plane mode searches.

## ADDED Requirements

### Requirement: Checked determinant evaluation
The evaluator SHALL offer checked determinant evaluation that rejects non-finite coordinate or exterior-wavevector jet components and non-finite determinant jet components. It SHALL preserve backend failures separately from input and output numerical failures.

#### Scenario: Derivative failure with finite primal
- **WHEN** the determinant primal value is finite but a carried derivative is non-finite
- **THEN** checked evaluation returns an output numerical failure rather than a successful determinant

#### Scenario: Invalid exterior input
- **WHEN** an exterior wavevector carries a non-finite primal or derivative component
- **THEN** checked evaluation returns an input numerical failure

#### Scenario: Backend failure
- **WHEN** a backend rejects the evaluation
- **THEN** checked evaluation preserves the backend error as its source

### Requirement: Characteristic semantics
For finite regular inputs the checked determinant SHALL represent the outgoing homogeneous boundary characteristic in the same gauge as raw determinant evaluation. Successful checked evaluation SHALL preserve a zero determinant as a valid numerical result.

#### Scenario: Accepted zero value
- **WHEN** evaluation at a mode produces a finite zero determinant and finite carried derivatives
- **THEN** checked evaluation returns that determinant successfully

#### Scenario: Regular point agreement
- **WHEN** raw and checked evaluations succeed at the same regular point
- **THEN** their determinant values and carried derivatives agree within numerical precision

### Requirement: Analytic transfer continuation through interior cutoff
For finite nonzero polarization constitutive factors and finite representable results, transfer determinant evaluation SHALL continue through an interior longitudinal-wavevector cutoff with finite values and carried derivatives. Crossing an interior square-root branch cut SHALL NOT change the physical characteristic.

#### Scenario: Exact TE or TM interior cutoff
- **WHEN** a finite layer has squared longitudinal wavevector zero and its polarization constitutive factor is finite and nonzero
- **THEN** the transfer determinant and its carried derivatives agree with the analytic cutoff limits

#### Scenario: Interior branch-cut crossing
- **WHEN** the squared interior wavevector crosses the principal square-root cut while the characteristic remains regular
- **THEN** the transfer determinant and derivatives follow the same analytic continuation on both sides

### Requirement: Explicit domain responsibility
Public checked-determinant documentation SHALL require caller-selected branch-consistent exterior data and a holomorphic characteristic on a neighborhood of the closed region for pole-free zero searches. It SHALL distinguish pointwise numerical checks from domain certification, and explain material singularities and removable backend failures.

#### Scenario: Finite contour samples
- **WHEN** a downstream caller obtains finite checked values on a contour
- **THEN** the documented contract provides no claim that the region is free of interior material singularities or determinant poles

#### Scenario: Fixed or varying spectral coordinate
- **WHEN** a caller constructs a modal search
- **THEN** the contract requires material validity at the fixed spectral point or throughout the mapped varying-spectral domain, respectively

### Requirement: Unresolved representation failures remain failures
When a scattering representation failure produces non-finite determinant components, checked evaluation SHALL return a numerical failure. It SHALL NOT substitute a raw projective denominator, perturb the requested coordinate, or silently switch backends to manufacture a successful result.

#### Scenario: Singular transmission chart
- **WHEN** scattering projection produces a non-finite ratio at a transmission chart singularity
- **THEN** checked evaluation returns an output numerical failure without claiming a physical determinant pole or an accepted mode

### Requirement: Existing evaluation compatibility
The existing raw determinant entry point and backend error types SHALL remain available. At previously regular inputs, the new numerical continuation SHALL preserve existing determinant semantics and driven response within numerical tolerance.

#### Scenario: Existing raw consumer
- **WHEN** a downstream consumer continues using the raw determinant entry point
- **THEN** it is not required to handle the new checked-evaluation error type
