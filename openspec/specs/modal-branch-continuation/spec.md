# modal-branch-continuation Specification

## Purpose

Let callers trace accepted modal roots into useful physical dispersion branches while preserving local numerical evidence, solver strategy and exterior branch provenance at every step.

## Requirements

### Requirement: Seeded physical branch request
The interface SHALL continue from one accepted modal point over a finite, ordered schedule of real physical control coordinates in inverse centimetres. Propagation sources SHALL vary real k0 and solve complex beta; direct-frequency and supported conformal-spectral sources SHALL vary real beta and solve complex k0.

#### Scenario: Propagation branch
- **WHEN** a caller continues a propagation point over real k0 values
- **THEN** the trace reports the corresponding complex beta root at each accepted station

#### Scenario: Spectral branch
- **WHEN** a caller continues a direct-frequency or supported conformal-spectral point over real beta values
- **THEN** the trace reports the corresponding complex k0 root at each accepted station

#### Scenario: Invalid schedule
- **WHEN** the schedule is empty, non-finite, repeated, or reverses direction
- **THEN** the request fails validation before advancing the seed point

### Requirement: Predictor-corrector continuation
Each next modal root SHALL be predicted from the previous accepted point's physical implicit slope and corrected locally on the same solution branch. A continuation step SHALL NOT repeat an unrestricted global discovery search, and an unavailable or degenerate slope SHALL produce a typed step outcome.

#### Scenario: Regular predicted step
- **WHEN** the previous point has a usable slope and the predicted point is locally regular
- **THEN** the corrector returns a newly accepted physical root associated with the requested control coordinate

#### Scenario: Local correction failure
- **WHEN** local correction fails within its configured work bounds
- **THEN** the step is retried only within the configured continuation policy and does not silently fall back to a global search

### Requirement: Bounded adaptive stepping
Continuation settings SHALL validate finite step limits, retry counts and correction tolerances. When an interval cannot be corrected at its requested size, the trace MAY subdivide it down to the configured minimum step and SHALL retain any inserted accepted points.

#### Scenario: Interval subdivision
- **WHEN** a requested interval fails at its current step size and remains above the minimum step
- **THEN** the continuation controller tries smaller intervals within its retry budget and records accepted intermediate points in order

#### Scenario: Minimum step exhausted
- **WHEN** the correction still fails at the minimum step or retry budget
- **THEN** the trace stops with an explicit failure associated with the requested station and its last accepted point

### Requirement: Branch and solver provenance
Every accepted continuation point SHALL retain its own problem/domain association, physical coordinates, numerical evidence, backend, source strategy, material stack, polarization and exterior branches. Continuation SHALL NOT switch strategy or silently reselect principal exterior roots.

#### Scenario: Dispersive exterior branch
- **WHEN** a spectral branch uses dispersive exterior materials
- **THEN** each corrected point retains branch-consistent exterior wavevectors evaluated from the source materials at that point's physical frequency

#### Scenario: Backend declaration
- **WHEN** a caller supplies an evaluation backend for continuation
- **THEN** the trace records that backend and does not substitute another backend after a step failure

### Requirement: Owned partial branch results
The branch trace SHALL own the seed, all accepted points and per-step outcomes. A failed continuation SHALL preserve earlier accepted points and its typed diagnostic without marking the branch complete or discarding the valid prefix.

#### Scenario: Failure after successful steps
- **WHEN** a later step reaches a cutoff, derivative failure or corrector limit
- **THEN** the result retains the seed and successful prefix alongside the failed station and cause

### Requirement: Local branch scope
A continuation trace SHALL describe only the numerical path reached from its seed under its settings. It SHALL NOT claim exhaustive mode discovery, global branch identity, crossing resolution, exceptional-point detection or ZGV certification.

#### Scenario: Nearby roots
- **WHEN** another mode lies near the predicted continuation point
- **THEN** the trace reports its local corrector and branch evidence without claiming that proximity alone proves global mode identity
