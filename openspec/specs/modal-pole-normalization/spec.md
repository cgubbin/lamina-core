# modal-pole-normalization Specification

## Purpose

Provide complete, convention-specific normalization and normalized fields for simple planar modal poles, with exterior regularization and independent physical validation suitable for later observable calculations.

## Requirements

### Requirement: Explicit physical pole convention
Normalization SHALL identify the differentiated physical pole variable, held-fixed coordinate, reciprocal pairing, units and planar Fourier convention. Frequency normalization at fixed beta and propagation normalization at fixed k0 SHALL remain distinct. A native chart derivative or temporal decay alone SHALL NOT establish either normalization.

#### Scenario: Frequency pole
- **WHEN** a supported simple spectral point is normalized
- **THEN** the result identifies physical k0 at fixed beta and records any conversion to angular frequency explicitly

#### Scenario: Propagation pole
- **WHEN** a supported fixed-frequency point is normalized
- **THEN** the result uses the propagation-pole convention and does not label a spectral energy factor as a beta-pole normalization

### Requirement: Complete regularized exterior treatment
Complete normalization SHALL include the finite-layer and both exterior contributions for the selected branches, with actual material dispersion where required. It SHALL retain contribution values, units and the applied exterior regularization. Finite-distance truncation of growing fields SHALL NOT substitute for a regularized exterior calculation.

#### Scenario: Zero finite layers
- **WHEN** a regular interface-plasmon pole has no finite layers
- **THEN** complete normalization remains computable from the exterior contributions rather than failing because the layer aggregate is empty

#### Scenario: Leaky spectral exterior
- **WHEN** outgoing fields grow along the real exterior coordinate
- **THEN** the normalization uses the documented analytic exterior continuation and preserves sampling growth independently

#### Scenario: One-dimensional Maxwell cancellation
- **WHEN** a normal-incidence cavity has nondispersive homogeneous exteriors
- **THEN** validation demonstrates the equivalent finite-layer Maxwell and scalar boundary forms without adding a second exterior contribution to the Maxwell form

### Requirement: Source-associated owned normalization
Public normalization SHALL be separately fallible, own its evidence and retain the accepted point, ordered branches, source and evaluation backends, local analytic assumptions and residual settings. It SHALL reuse core device and field primitives, without rerunning discovery or refinement, replacing branches or changing acceptance and radiation provenance.

#### Scenario: Stored accepted point
- **WHEN** the original problem and solve report have been dropped
- **THEN** a retained accepted point supports normalization with an explicit supported backend and the required local assumptions

#### Scenario: Unspecified radiation provenance
- **WHEN** normalization succeeds on branches selected by an explicit complex anchor
- **THEN** the result preserves unspecified radiation qualification rather than promoting it to verified outgoing provenance

### Requirement: Gauge-consistent normalized fields
Normalization SHALL retain the original complex gauge, its normalization factor and the applied square-root scale. E and H SHALL receive one consistent scale. Physical pole residues SHALL be invariant under any tested nonzero complex seed rescaling; normalized fields SHALL be equivalent up to the documented common sign freedom.

#### Scenario: Complex seed rescaling
- **WHEN** the same mode is evaluated with different finite nonzero complex seed multipliers
- **THEN** its predicted residue is unchanged and its normalized fields differ by at most the documented common sign

#### Scenario: Sampling normalized fields
- **WHEN** a caller requests normalized fields at valid core sampling positions
- **THEN** owned electric and magnetic samples preserve position order, interface sides, source association and normalization evidence

### Requirement: Independent physical normalization checks
Maintained validation SHALL compare normalized residue predictions with an independently constructed Green function and an analytic planar cavity reference. It SHALL cover frequency and propagation poles, TE/TM, dispersive materials and zero-layer exterior contributions. Agreement between implementations sharing one normalization formula SHALL NOT be the sole correctness test.

#### Scenario: Independent pole residue
- **WHEN** a simple benchmark pole is approached through a decreasing sequence of complex offsets
- **THEN** the directly evaluated response residue converges to the normalized-mode prediction within documented tolerances before roundoff dominates

#### Scenario: Incorrect normalization term
- **WHEN** a validation control omits an exterior or dispersive term, conjugates a required bilinear product, or uses the wrong pole variable
- **THEN** the independent reference detects the incorrect residue or normalization

### Requirement: Exterior partition and backend invariance
For regular supported fixtures, moving finite integration boundaries within homogeneous exteriors SHALL leave the complete normalization and physical residues unchanged within documented tolerances after gauge alignment. Supported backends SHALL agree on normalized predictions under the same conventions. Individual finite and exterior contributions MAY change with the partition.

#### Scenario: Displaced exterior boundaries
- **WHEN** both finite integration boundaries are displaced independently within the unchanged homogeneous exteriors
- **THEN** the altered finite-region and exterior-tail terms cancel in the complete result

#### Scenario: Backend comparison
- **WHEN** the same accepted physical point is normalized with Scatter2 and Transfer2
- **THEN** their complete residue predictions and normalized field components agree under one allowed common sign alignment

### Requirement: Diagnostic normalization limits
Invalid settings, incompatible pole assumptions, cutoff or material singularities, nonfinite contributions and small or ill-conditioned normalization factors SHALL produce structured failures preserving available causes and the original point. Near-degeneracy SHALL NOT trigger an exceptional-point claim, a fabricated finite scale or a silent fallback to legacy normalization.

#### Scenario: Vanishing normalization factor
- **WHEN** a normalization factor falls within the configured dimensioned exclusion band
- **THEN** normalization is unavailable with its component evidence and thresholds retained, while the accepted point and arbitrary-amplitude fields remain usable

#### Scenario: Unsupported local analytic neighborhood
- **WHEN** the required branch or constitutive variation is singular or lacks its required local declaration
- **THEN** normalization fails diagnostically without inferring analyticity from a finite set of samples

### Requirement: Legacy and observable claim boundaries
The existing finite-layer normalization entry point SHALL retain its established numerical behavior and be documented separately from complete pole normalization. Public examples SHALL distinguish arbitrary gauge, legacy scaling and complete validated normalization. This capability SHALL NOT by itself claim finite three-dimensional mode volumes, Purcell factors or complete modal expansions.

#### Scenario: Legacy normalized consumer
- **WHEN** an existing caller uses the legacy normalized core mode operation
- **THEN** its numerical semantics remain unchanged and documentation identifies the finite-layer convention and its limits

#### Scenario: Reproducible validation report
- **WHEN** the normalization validation example runs
- **THEN** it records model, coordinate and source conventions, contributions, settings, failed cases and residue errors, and produces reproducible convergence evidence without claiming unsupported observables
