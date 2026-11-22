# modal-physical-benchmarks Specification

## Purpose

Provide repeatable public mode-solver benchmarks that demonstrate physical branch selection, complex eigenvalue accuracy and diagnostic limits through independent equations, automated checks and reproducible plots.

## Requirements

### Requirement: Independent interface plasmon references
The benchmark SHALL exercise numerical discovery and acceptance for a passive metal-dielectric interface in both real-frequency/complex-wavevector and real-wavevector/complex-frequency cases. Reference eigenvalues and the unsquared boundary condition SHALL be calculated independently of the solver characteristic. Analytic chart placement SHALL be disclosed.

#### Scenario: Spatial plasmon
- **WHEN** a real frequency is prescribed in the supported passive interface benchmark
- **THEN** the accepted TM wavevector agrees with the independent interface relation, attenuates along positive propagation and has decaying exterior factors on both sides

#### Scenario: Temporal plasmon
- **WHEN** a real wavevector is prescribed in the supported passive interface benchmark
- **THEN** the accepted complex frequency agrees with the independently selected physical dispersion root, decays in time and retains the bound exterior branches

### Requirement: Physical controls and diagnostic assertions
Benchmark validation SHALL check finite root error, unsquared TM boundary mismatch, exterior decay, passive attenuation and supported backend consistency against documented tolerances. Local TE controls SHALL distinguish checked windows from untested points and SHALL NOT imply exhaustive absence of modes.

#### Scenario: Incorrect sheet or attenuation sign
- **WHEN** a computed benchmark root has the wrong exterior decay or passive evolution sign despite a small characteristic value
- **THEN** validation fails with the offending case and physical check identified

#### Scenario: TE negative control
- **WHEN** the configured TE control window is searched
- **THEN** its actual outcome is recorded separately from points where no TE check was run

### Requirement: Maintained quick and full validation runs
The project SHALL provide a bounded quick suite suitable for routine regression tests and a separately invoked full public validation run. The numerical quick suite SHALL run without plotting dependencies. Neither suite SHALL substitute reference values for accepted numerical roots or hide failed search outcomes.

#### Scenario: Routine regression
- **WHEN** the quick suite runs with the repository's numerical dependencies
- **THEN** representative points in both complex cases are independently searched and validated within a documented bounded workload

#### Scenario: Full plotted run
- **WHEN** the full validation command is requested
- **THEN** it generates the configured dispersion samples and characteristic evaluations through public solver entry points and preserves any failure with its case identity

### Requirement: Reproducible benchmark artifacts
Validation outputs SHALL retain model parameters, units, numerical settings, reference branches, physical checks and measured errors in machine-readable form. Figures SHALL distinguish numerical roots from analytic references and show attenuation and characteristic-zero evidence. Plotting failures SHALL NOT silently reuse stale data.

#### Scenario: Recreate figures
- **WHEN** a user follows documented generation and plotting commands
- **THEN** the figures can be recreated from the recorded numerical data without workstation-specific filesystem paths or hidden notebook state

#### Scenario: Unsupported interpretation
- **WHEN** a metal mode lacks a supported built-in named classification or verified dielectric radiation anchor
- **THEN** the report preserves that limitation while identifying the separately checked physical branch, without claiming global discovery or QNM normalization
