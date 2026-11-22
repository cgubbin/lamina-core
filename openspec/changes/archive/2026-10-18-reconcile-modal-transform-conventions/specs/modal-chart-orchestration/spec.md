# Spec Delta

## ADDED Requirements

### Requirement: Sheet-resolved conformal coverage
Chart plans SHALL distinguish squared-coordinate repetition from repetition of ordered exterior wavevectors. Canonical coverage reports SHALL identify omitted exterior-sheet sectors in addition to finite truncations and contour margins. A shift preserving only the squared coordinate SHALL NOT establish physical sheet equivalence.

#### Scenario: Half-period shift
- **WHEN** a conformal shift by i*pi preserves the squared physical coordinate and reverses both nonzero exterior wavevectors
- **THEN** the shifted point remains distinct in exterior-sheet identity and is not treated as recovered coverage by periodicity

#### Scenario: Canonical restricted strip
- **WHEN** canonical charts cover only an imaginary strip of width pi within a full exterior-sheet period of width 2*pi
- **THEN** the coverage report identifies the unsearched sheet sectors and does not describe the selected strip as the complete sheet-resolved surface

#### Scenario: Explicit domain in another sector
- **WHEN** a valid explicit domain lies outside the canonical selected strip
- **THEN** it remains searchable under the existing domain contract and retains its own exterior-sheet provenance

### Requirement: Public transform correspondence
Public solver documentation SHALL explain coordinate and exterior-sheet periods, the correspondence to the cited fixed-frequency transform, and restrictions of the spectral extension. It SHALL state the assumptions behind borrowed mode labels and distinguish numerical acceptance from physical classification.

#### Scenario: Comparing with the reference paper
- **WHEN** a reader follows the mathematical guide from the crate README
- **THEN** the guide provides the reference, notation and sign dictionary, supported material assumptions, sheet coverage limits and explicit differences from the paper

#### Scenario: Interpreting an accepted spectral root
- **WHEN** numerical refinement accepts a root on a selected exterior sheet
- **THEN** documentation does not imply that the root is an all-outgoing frequency QNM or assign the paper's fixed-frequency classification without establishing the relevant assumptions
