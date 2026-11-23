# guided-mode-sensitivities Specification

## Purpose

Provide total derivatives of propagation constants and power-normalized eigenfields along declared physical perturbations for downstream guided-wave models.

## Requirements

### Requirement: Physical perturbation direction
The API SHALL accept an accepted simple lossless guided mode and a real parameter direction with explicit units, geometry frame, wavelength rate and regional constitutive rates. It SHALL distinguish material dispersion induced by wavelength from explicit material variation at fixed wavelength and reject invalid geometry or unsupported directions.

#### Scenario: Moving one interior interface
- **WHEN** an interior interface moves with its neighbours fixed
- **THEN** the two adjacent thickness rates have opposite signs and all other thickness rates vanish

#### Scenario: Index variation
- **WHEN** a nonmagnetic dielectric index changes at fixed wavelength
- **THEN** the direction uses the corresponding relative-permittivity rate at the specified region and retains the parameter units

### Requirement: Branch-consistent eigenvalue derivative
The solver SHALL evaluate the total propagation-constant derivative by implicit differentiation on the accepted exterior branches. It SHALL retain the point, derivative units and conditioning evidence and report typed unavailability at cutoff or an unresolved simple-root denominator.

#### Scenario: Dispersive wavelength sweep
- **WHEN** vacuum wavelength varies in a dispersive stack
- **THEN** the result includes wavelength-induced changes of materials and exterior wavevectors in the physical derivative

### Requirement: Moving observations and normalized fields
An owned result SHALL give complex E and H derivatives along the same root in the declared observation frame and the existing signed unit-power and phase convention. It SHALL include changes of finite-layer and exterior power, phase and observation positions, with a typed diagnostic for a fixed lab point swept by a discontinuity or a singular phase reference.

#### Scenario: Rigid translation
- **WHEN** every interface and a material-attached observation translate together
- **THEN** the propagation constant and co-moving field have zero derivative while a stationary lab observation inside a smooth region follows the translated field gradient

#### Scenario: One-sided sample
- **WHEN** the observation follows a specified side of a moving interface
- **THEN** the result retains that side and differentiates the corresponding limiting field

#### Scenario: Downstream use
- **WHEN** the caller releases temporary solver objects
- **THEN** it can use the base profile, propagation-constant derivative, field derivatives and diagnostics without accessing internal scattering matrices
