# dispersive-spectral-search Specification

## Purpose

Search finite complex-frequency regions using actual dispersive exterior materials and analytically continued radiation branches while preserving local numerical evidence and explicit domain limitations.

## Requirements

### Requirement: Actual dispersive stack evaluation
Direct frequency search SHALL accept an existing analytical core stack at fixed real in-plane wavenumber. Exterior constitutive values and derivatives SHALL come from that same stack throughout evaluation. It SHALL NOT freeze exterior values or require callers to rebuild the device in a separate stack representation.

#### Scenario: Dispersive substrate
- **WHEN** a stack has an exterior whose analytical constitutive model varies over the requested frequency chart
- **THEN** characteristic values and frequency derivatives include that variation at each evaluated frequency

#### Scenario: Material without complex continuation
- **WHEN** a supplied material lacks the required complex-frequency values or derivatives
- **THEN** it cannot enter the analytical stack path or produces a construction error rather than an implicit fit, extrapolation or frozen replacement

### Requirement: Analytic exterior branch continuation
Each exterior wavevector SHALL be continued from a finite nonzero reference satisfying its exterior dispersion relation. Values and derivatives SHALL preserve that branch independently of evaluation order. Continuation SHALL NOT use endpoint principal square roots or the previously sampled point to choose a new sheet.

#### Scenario: Principal cut crossing
- **WHEN** the squared exterior wavevector crosses the principal square-root cut inside a valid analytic chart without vanishing
- **THEN** its wavevector and frequency derivative remain on the reference-defined analytic continuation

#### Scenario: Opposite exterior reference
- **WHEN** one reference wavevector is negated while its squared relation remains unchanged
- **THEN** only that exterior's continued value and derivative reverse sign and the two references remain distinct in provenance

### Requirement: Explicit frequency domain obligations
Construction SHALL require caller declarations of material and characteristic analyticity and nonvanishing holomorphic squared exterior wavevectors on a neighbourhood of the evaluation region. It SHALL retain known pole and threshold exclusions and reject charts or halos intersecting them. Samples SHALL NOT be presented as domain certification.

#### Scenario: Pole or threshold in the declared domain
- **WHEN** a requested chart or its numerical halo intersects a caller-declared material-pole or exterior-threshold exclusion
- **THEN** construction rejects the domain before adaptive search, preserving the specific exclusion reason

#### Scenario: Unknown interior singularity
- **WHEN** finite checked samples succeed without independent knowledge of all material poles or exterior zeros
- **THEN** results retain the caller assumptions and make no claim that those samples established a valid analytic domain

### Requirement: Controlled fallible frequency evaluation
Evaluation SHALL return the checked characteristic and its derivative with respect to physical vacuum angular wavenumber. Material, exterior-continuation, numerical and backend failures SHALL remain typed diagnostic outcomes with their available causes. Continuation work SHALL have validated finite tolerances and bounds; failures SHALL NOT produce fabricated roots or trigger silent strategy switches.

#### Scenario: Exterior derivative contribution
- **WHEN** exterior permittivity or permeability has a nonzero frequency derivative
- **THEN** the characteristic derivative includes the resulting exterior-wavevector derivative as well as the remaining optical derivative terms

#### Scenario: Continuation budget exhaustion
- **WHEN** exterior continuation cannot meet its tolerance within configured work limits
- **THEN** the evaluation and search report that failure with its exterior identity rather than returning a principal-root substitute or successful empty result

### Requirement: Reference and radiation provenance
Frequency domains SHALL retain explicit exterior reference values and their verification separately from physical radiation qualification. They SHALL support real-anchor outgoing/decaying verification under the existing positive-dielectric assumptions using the actual anchor material values. A complex explicit reference SHALL NOT acquire outgoing qualification merely by satisfying its squared relation.

#### Scenario: Dispersive real radiation anchor
- **WHEN** the actual dispersive exterior factors at a real anchor satisfy the supported dielectric conditions and match requested radiation choices
- **THEN** that verified choice is retained through analytic frequency continuation without requiring exterior constancy away from the anchor

#### Scenario: Lossy complex reference
- **WHEN** finite nonzero wavevectors are explicitly supplied at a complex reference and satisfy the actual exterior dispersion relations
- **THEN** the domain can preserve those branches under its assumptions while physical radiation qualification remains unspecified

### Requirement: Integrated local search evidence
Direct frequency charts SHALL participate in discovery, bounded refinement, physical descriptions and finite chart plans under the existing evidence-preserving contracts. Reports SHALL retain physical chart bounds, units, references, exclusions and unresolved outcomes. Opposite exterior sheets SHALL remain distinct; local successes SHALL NOT imply global coverage or recover excluded regions.

#### Scenario: Accepted dispersive mode
- **WHEN** a simple root converges inside its originating frequency domain
- **THEN** its description retains frequency coordinates, ordered exterior data, material assumptions, local acceptance evidence and radiation provenance without receiving the paper's fixed-frequency named type

#### Scenario: Independent charts with a failure
- **WHEN** one finite frequency chart fails while another succeeds
- **THEN** both outcomes and declared exclusions remain visible in the aggregate

#### Scenario: Constant exterior limit
- **WHEN** direct frequency and constant-exterior conformal searches represent the same regular physical region and exterior branches
- **THEN** their characteristic values, physical-coordinate derivatives and accepted physical roots agree within documented numerical tolerances
