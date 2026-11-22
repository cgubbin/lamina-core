# Spec Delta

## ADDED Requirements

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
