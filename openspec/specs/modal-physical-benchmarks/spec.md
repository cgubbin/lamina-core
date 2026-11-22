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

### Requirement: Maintained physical branch report cases
The full public validation report SHALL include continued physical branches for the maintained plasmon fixture, a planar avoided-crossing fixture, and a stopped-light planar waveguide case from Adam Freddie Page's Imperial thesis. Each case SHALL retain its model parameters, units, source references, accepted points and physical comparison tolerances.

#### Scenario: Plasmon branches
- **WHEN** the full branch report runs the supported interface-plasmon cases
- **THEN** it plots real-frequency/complex-wavevector and real-wavevector/complex-frequency branches against their independent physical references

#### Scenario: Avoided crossing
- **WHEN** the report continues both roots through the maintained planar avoided-crossing fixture
- **THEN** it retains both seeded traces through the near-degeneracy and shows the gap without claiming persistent global mode identity from coordinate ordering alone

#### Scenario: Page stopped-light waveguide
- **WHEN** the report evaluates the thesis stopped-light fixture
- **THEN** it compares the continued complex-frequency branch and its measured local slope near the two reference points using traceable published parameters and units

### Requirement: Reproducible branch dispersion and field plots
The full report SHALL emit machine-readable branch samples and reproducible dispersion figures, and SHALL plot transverse electric and magnetic field profiles at accepted points around the avoided crossing and stopped-light points. Field gauge and normalization status SHALL be identified for each profile; unavailable normalization SHALL remain explicit.

#### Scenario: Recreate branch figures
- **WHEN** a user runs the documented report and plotting commands
- **THEN** the dispersion figures and field profiles are recreated from the recorded branch data, without stale output or workstation-specific paths

#### Scenario: Fields near branch events
- **WHEN** representative accepted points are selected before, near and after an avoided crossing or stopped-light point
- **THEN** their transverse field profiles are associated with the exact point, physical coordinates, exterior branches and recorded gauge or QNM normalization
