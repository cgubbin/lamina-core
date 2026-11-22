# modal-point-solver-api Specification

## Purpose

Allow downstream Rust callers to construct finite modal point searches from physical device inputs and retain usable results, numerical evidence, radiation provenance and diagnostic limitations through a public library interface.

## Requirements

### Requirement: Public physical problem construction
The library SHALL expose fixed-frequency propagation, constant-exterior conformal spectral and dispersive direct-frequency problem construction. Inputs SHALL reuse existing core device primitives and retain the chosen strategy and backend. Callers SHALL NOT need private chart, jet or compiled-solver types to construct or run a supported solve.

#### Scenario: Existing analytical stack
- **WHEN** an external caller supplies an analytical core stack for a fixed-frequency or direct-frequency solve
- **THEN** the public interface uses that stack without requiring a parallel device representation or replacement exterior models

#### Scenario: Constant conformal exterior evidence
- **WHEN** an external caller selects the constant-exterior conformal spectral strategy
- **THEN** construction requires constant exterior evidence, retains dispersive finite layers and does not infer exterior constancy from samples

### Requirement: Public domain and reference controls
The public interface SHALL expose distinguishable conformal and frequency bounds, branch references, supported radiation-reference choices, analytic assumptions and known exclusions. Fixed coordinates, native bounds, derivatives, corrections and tolerances SHALL have documented units. Required analytic declarations SHALL NOT be silently supplied by a convenience solve.

#### Scenario: Direct-frequency reference
- **WHEN** a caller constructs a direct-frequency domain with an explicit complex reference or a supported real radiation anchor
- **THEN** the corresponding branch data and radiation qualification remain distinct and inspectable with physical frequency units

#### Scenario: Required declaration omitted
- **WHEN** a caller has not supplied the domain assumptions required by its strategy
- **THEN** the public interface cannot execute the solve by inventing those assumptions or treating finite samples as certification

### Requirement: Public finite chart execution
Callers SHALL be able to execute explicit finite chart plans and the existing supported canonical conformal plans, configure validated discovery and refinement controls, and inspect discovery before refinement. Domain and plan ownership checks SHALL reject foreign associations before optical execution. Unsupported geometry SHALL remain an error requiring explicit bounds.

#### Scenario: Explicit dispersive frequency plan
- **WHEN** a caller supplies valid finite direct-frequency domains of one problem
- **THEN** execution uses those exact domains and settings and retains each chart outcome in plan order

#### Scenario: Foreign association
- **WHEN** a caller supplies a domain or discovery result from another problem, domain or plan
- **THEN** the public interface rejects the association before evaluating it

### Requirement: Owned point-solve reports
The public interface SHALL provide an owned solve report that can be returned or stored after temporary execution objects are dropped. It SHALL retain discovery outcomes, all candidate outcomes, accepted observations, group membership, coverage and source provenance. Owning a report SHALL NOT require leaked allocations, unsafe lifetime extension or caller-maintained self-references.

#### Scenario: Return from a helper
- **WHEN** a helper constructs a problem and finite plan, performs a solve and returns its report
- **THEN** the caller can inspect the report's points, failures, enclosures, assumptions and coverage after the helper's local objects have been dropped

#### Scenario: Separate physical result access
- **WHEN** a caller accesses an accepted point independently of the report's temporary iteration views
- **THEN** a supported owned point value retains its numerical evidence and exact source provenance without claiming to retain computed fields or normalization

### Requirement: Public evidence-preserving result inspection
Public reports SHALL distinguish estimates, accepted roots, unresolved or failed candidates and unavailable interpretation. They SHALL expose original multiplicities, leaves, summaries, terminations, exclusions and overlap contributors. Accepted points or duplicate groups SHALL NOT imply exhaustive coverage or merge incompatible exterior branches.

#### Scenario: Successful and failed overlapping charts
- **WHEN** one listed chart succeeds and another fails, including in overlapping scope
- **THEN** the owned report retains both outcomes and the failure's coverage limitation alongside accepted observations and their groups

#### Scenario: Budget-limited cluster
- **WHEN** discovery stops with a cluster or unresolved leaves
- **THEN** public inspection retains its multiplicity and termination and never lists the cluster centroid as an accepted individual mode

### Requirement: Public structured diagnostics
Construction, domain, discovery, refinement and interpretation failures SHALL remain distinguishable through public diagnostics. Available original numerical, continuation and backend causes SHALL remain inspectable through structured fields or error sources. The interface SHALL NOT reduce failures to strings, successful empty results or silent backend or strategy substitutions.

#### Scenario: Exterior continuation failure
- **WHEN** direct-frequency evaluation exhausts an exterior continuation budget
- **THEN** the public diagnostic retains the exterior identity, failure kind and available work count through discovery or candidate reporting

#### Scenario: Checked optical failure
- **WHEN** checked optical evaluation fails on an input, backend or output component
- **THEN** an external caller can identify that cause through the public error surface without importing a private module

### Requirement: Public optional physical interpretation
Public accepted points SHALL expose physical coordinates, ordered exterior data and radiation provenance, with supported interpretation available as a separate fallible result. Unsupported classification or unavailable physical data SHALL retain acceptance and evidence. Temporal decay alone SHALL NOT label a spectral root as an outgoing or normalized QNM.

#### Scenario: Supported propagation type
- **WHEN** an accepted fixed-frequency point satisfies the supported dielectric interpretation assumptions
- **THEN** public interpretation reports its supported type and component evidence with the applied tolerances

#### Scenario: Spectral interpretation limitation
- **WHEN** a spectral point has temporal decay and an explicit complex reference without verified radiation
- **THEN** the public result preserves unspecified radiation, reports the fixed-frequency named classification as unsupported and makes no normalization claim

#### Scenario: Physical data unavailable
- **WHEN** physical interpretation fails after numerical acceptance
- **THEN** its diagnostic remains attached to the accepted point without removing the point or changing discovery and coverage

### Requirement: Executable public usage
The project SHALL provide runnable guided-mode, constant-exterior cavity and dispersive-substrate examples using public imports and existing core/material models. Examples SHALL check expected physical results and show numerical evidence, radiation qualification and relevant coverage limits. Public quickstart code SHALL compile as a doctest or runnable example rather than relying only on ignored private snippets.

#### Scenario: Downstream example execution
- **WHEN** documented example commands run in the repository development environment
- **THEN** all three examples execute through the public interface and check their expected results without private-module access or test-only inclusion tricks

#### Scenario: Dispersive substrate example
- **WHEN** the frequency example uses a dispersive exterior model
- **THEN** it evaluates that model throughout the chart, retains its assumptions and reference and does not freeze the exterior at the anchor
