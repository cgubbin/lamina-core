# Spec Delta

## Purpose

Describe the physical identity and supported mode classifications of accepted modal roots while preserving their numerical evidence and selected radiation branches.

## ADDED Requirements

### Requirement: Traceable physical mode points
Physical interpretation SHALL consume accepted roots only with their originating problem and domain. Each description SHALL retain polarization, solve direction, physical coordinates, ordered exterior wavevectors, numerical acceptance evidence, domain identity and radiation-reference provenance. Foreign associations SHALL be rejected before optical evaluation.

#### Scenario: Accepted root on its originating branch
- **WHEN** a single-chart or batch accepted root is interpreted with its originating problem and domain
- **THEN** its physical data are evaluated on that domain's existing analytic branch and remain traceable to the accepted observation

#### Scenario: Foreign domain of the same problem
- **WHEN** interpretation is requested using another domain of the same problem
- **THEN** it is rejected before optical evaluation rather than reinterpreting the root on that domain's branch

### Requirement: Separate physical behaviour components
Descriptions SHALL distinguish left and right exterior decay, growth or boundary behaviour from spatial attenuation or amplification and temporal decay or growth. They SHALL retain finite raw quantities and interpretation tolerances. Exterior flux direction SHALL be reported only under supported physical assumptions, otherwise with an explicit unsupported reason.

#### Scenario: Spatially attenuating propagation mode
- **WHEN** fixed real frequency gives positive imaginary propagation wavenumber outside the configured boundary band
- **THEN** the description reports spatial attenuation along the positive propagation coordinate independently of either exterior's field decay or radiation direction

#### Scenario: Temporally decaying spectral mode
- **WHEN** a complex-frequency root has negative imaginary vacuum wavenumber outside the configured boundary band
- **THEN** it reports temporal decay without assigning a fixed-frequency paper label or inferring instantaneous flux direction from the sign of an exterior square root

### Requirement: Supported dielectric mode types
For positive real frequency, positive real propagation direction and positive real nonmagnetic dielectric exteriors, descriptions SHALL support guided modes and the paper's one-sided leaky/source, two-sided leaky/source, absorptive and gain types. Labels SHALL follow ordered exterior decay and flux evidence under the documented conventions, including reversed exterior ordering and the symmetric limit.

#### Scenario: Guided slab
- **WHEN** a supported fixed-frequency mode has real propagation coordinate above both exterior light lines and both exterior fields decay away from the stack
- **THEN** its type is guided within the configured tolerance convention

#### Scenario: Leaky and source sectors
- **WHEN** a supported point lies strictly within one of the paper's leaky, source, absorptive or gain sectors
- **THEN** its type agrees with the translated exterior sign evidence and identifies the actual left or right exterior rather than relying on a hardcoded conformal quadrant

#### Scenario: Symmetric exteriors
- **WHEN** the exterior dielectric constants are equal
- **THEN** interpretation supports the reduced sign sectors represented by the map without inventing unequal-index sectors or requiring finite canonical bounds

### Requirement: Explicit interpretation limits
Interpretation SHALL distinguish supported labels, boundary or ambiguous cases, unsupported assumptions and unavailable physical data. Non-finite inputs and invalid tolerances SHALL produce diagnostic errors. Near-cutoff or conflicting evidence SHALL NOT be forced into a definite named type.

#### Scenario: Exterior cutoff
- **WHEN** an exterior wavevector lies within the configured cutoff boundary band
- **THEN** the description identifies the boundary condition without asserting a strict leaky or source sector

#### Scenario: Complex or magnetic exterior
- **WHEN** the exterior assumptions required for the paper's named labels do not hold
- **THEN** the named classification records an unsupported reason while available field-growth and attenuation components remain accessible

#### Scenario: Unsupported propagation direction
- **WHEN** the real propagation wavenumber is nonpositive or within its direction boundary band
- **THEN** a positive-direction paper label is withheld with an explicit reason

### Requirement: Interpretation preserves numerical results
Classification SHALL be an optional layer over accepted roots. It SHALL NOT promote unresolved candidates, revoke numerical acceptance, filter observations or duplicate groups, alter multiplicities, or erase discovery and coverage evidence. A physical-data evaluation failure SHALL retain its source and the originating accepted root.

#### Scenario: Failed physical-data evaluation
- **WHEN** interpretation encounters an evaluation failure after numerical acceptance
- **THEN** it reports unavailable interpretation with the originating error while retaining acceptance and its numerical evidence

#### Scenario: Unsupported batch observation
- **WHEN** one observation cannot receive a supported physical label
- **THEN** other observations can still be interpreted and the original batch, group membership and coverage ledger remain unchanged

### Requirement: Radiation-qualified frequency interpretation
A spectral description SHALL distinguish an outgoing radiation continuation established by a validated real reference from unspecified or other radiation provenance. Numerical acceptance and temporal decay alone SHALL NOT establish an outgoing QNM. Documentation SHALL explain the reference-based qualification and defer field normalisation and mode volumes.

#### Scenario: Outgoing cavity pole
- **WHEN** a decaying spectral cavity pole is accepted on a domain with a validated all-outgoing real reference
- **THEN** its description retains that outgoing-continuation qualification and temporal decay without claiming a normalised QNM or complete spectral coverage

#### Scenario: Arbitrary spectral sheet
- **WHEN** a spectral root is accepted on a domain with no validated radiation reference
- **THEN** radiation provenance remains unspecified even if the root's imaginary frequency is negative
