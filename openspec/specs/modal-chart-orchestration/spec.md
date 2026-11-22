# modal-chart-orchestration Specification

## Purpose

Coordinate finite chart searches for one physical mode problem while preserving branch provenance, local numerical evidence, and explicit coverage limitations.

## Requirements

### Requirement: Problem-owned chart plans
A chart plan SHALL contain a nonempty ordered set of domains belonging to one physical problem, preserving each domain's exact identity, branch reference, requested bounds, evaluation bounds, and caller assumption. Invalid or foreign domains SHALL be rejected before execution.

#### Scenario: Foreign domain in a plan
- **WHEN** a domain from another physical problem is supplied
- **THEN** plan construction fails before any chart is searched

#### Scenario: Stable plan identity
- **WHEN** a valid plan is cloned or executed
- **THEN** its domain identities and order remain stable without rebuilding charts or expanding halos

### Requirement: Explicit finite search geometry
Canonical plans SHALL preserve supported geometry builders' finite bounds and margins. Unsupported or degenerate canonical geometry SHALL require explicit domains rather than silently inventing a truncation. Explicit plans SHALL support valid symmetric-exterior and normal-incidence charts.

#### Scenario: Canonical spectral geometry
- **WHEN** supported spectral extent inputs are supplied
- **THEN** the plan contains both lower and upper chart domains and records their excluded strip and finite truncations

#### Scenario: Degenerate canonical geometry
- **WHEN** canonical construction cannot determine finite bounds for symmetric exteriors or zero fixed in-plane wavenumber
- **THEN** it reports the need for explicit bounds and allows a separately validated explicit-domain plan

### Requirement: Independent chart execution
Execution SHALL retain one discovery outcome per planned chart, including complete successful results or typed failure causes. A failed chart SHALL NOT suppress execution or evidence from other charts, nor become a successful empty chart result.

#### Scenario: One chart fails
- **WHEN** discovery fails on one chart and succeeds on another
- **THEN** both outcomes remain identifiable in plan order and the failed chart remains a coverage limitation

### Requirement: Preserved refinement evidence
Batch refinement SHALL use each successful chart's exact originating domain and retain every local refinement outcome. It SHALL preserve all original discovery results, multiplicities, leaves, summaries, terminations, and branch assumptions.

#### Scenario: Unresolved local candidate
- **WHEN** a chart contains an unresolved cluster, failed refinement, or budget-limited discovery
- **THEN** aggregation retains that evidence even when another chart produces accepted roots

### Requirement: Conservative overlap reconciliation
Duplicate groups SHALL contain only accepted observations from overlapping charts with compatible physical coordinates, ordered exterior branches, and enclosure evidence. All source observations SHALL remain available. Ambiguous or incompatible observations SHALL remain separate without increasing or erasing their source multiplicities.

#### Scenario: Same root in overlapping charts
- **WHEN** accepted observations have mutually compatible isolation enclosures and matching physical branch data under the configured tolerances
- **THEN** they can be represented as one duplicate group retaining every chart and candidate reference

#### Scenario: Different radiation sheets
- **WHEN** physical coordinates are close but exterior radiation branches differ or cannot be distinguished reliably
- **THEN** observations are not merged and any ambiguity remains visible

#### Scenario: Chained proximity
- **WHEN** pairwise closeness forms a chain whose endpoints fail the compatibility criteria
- **THEN** the chain is not collapsed into one duplicate group

### Requirement: Explicit coverage ledger
An aggregate SHALL report coverage relative to its finite chart plan, retaining declared exclusions, failed charts, unresolved discovery evidence, and non-accepted candidates. It SHALL NOT infer global completeness, pole-free materials, or recovered seam coverage from accepted roots or duplicate reconciliation.

#### Scenario: Canonical margins exclude a root
- **WHEN** a root lies in an excluded strip or outside finite truncation bounds
- **THEN** the report retains that exclusion rather than claiming the physical search was exhaustive

#### Scenario: Overlapping success and failure
- **WHEN** one chart fails in a region also searched by a successful chart
- **THEN** the failure remains recorded without silently declaring its coverage recovered

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
