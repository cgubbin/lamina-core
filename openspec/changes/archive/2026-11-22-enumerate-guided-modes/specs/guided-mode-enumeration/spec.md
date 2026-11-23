# Spec Delta

## Purpose

Provide a counted basis of guided TE or TM modes at a real vacuum wavelength for downstream planar propagation models.

## ADDED Requirements

### Requirement: Physical enumeration domain
The solver SHALL accept an existing planar stack, positive real vacuum wavelength, polarization, propagation direction, finite positive interval of |β| in inverse centimetres and explicit numerical controls. It SHALL check for finite positive real permittivity and permeability at the operating wavelength and record the requested and eligible intervals and cutoff margins.

#### Scenario: Unsupported material
- **WHEN** the stack contains a lossy layer at the operating wavelength
- **THEN** the result identifies the unsupported material and carries no completeness claim

#### Scenario: Exterior light line
- **WHEN** the request reaches a light line
- **THEN** the report identifies the threshold and any excluded or numerically unresolved neighbourhood

### Requirement: Count reconciliation
The report SHALL give an expected count from an independent physical counting method when available, the number of distinct accepted roots, and a completeness status scoped to the eligible interval. Complete status SHALL require an unambiguous count and one accepted guided root in each singly counted bracket. Numerical budget exhaustion SHALL retain the partial basis and identify unresolved brackets.

#### Scenario: Empty guided interval
- **WHEN** an eligible interval has a reliable count of zero
- **THEN** the report returns an empty basis with complete coverage of that interval

#### Scenario: Close pair
- **WHEN** two expected states cannot be isolated within the numerical budget
- **THEN** the bracket remains unresolved and the report preserves its expected multiplicity where known

### Requirement: Owned guided basis
The result SHALL own modes ordered by decreasing |β| with the requested propagation sign. Each accepted point SHALL retain checked root evidence and decaying exterior branches and support the existing field, guided-power and local derivative operations independently of the enumeration workspace.

#### Scenario: Higher-order slab modes
- **WHEN** a downstream consumer requests a multimode dielectric slab basis
- **THEN** it can inspect ordered propagation constants and obtain power-normalized fields for each accepted order after enumeration returns
