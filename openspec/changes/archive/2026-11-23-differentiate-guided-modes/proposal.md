# Proposal

## Why

The enumerator supplies guided modes and their unit-power profiles. A differentiable slab propagator also needs the change in each propagation constant and in the complex field at its observation points as the device or wavelength changes.

## What Changes

- Describe a real perturbation direction on the existing planar stack in terms of interface motion, layer thickness, constitutive parameters and vacuum wavelength.
- Compute the total propagation-constant derivative for a simple guided root and the corresponding derivatives of all six unit-power field components.
- Expose the observation frame, side, units, phase reference and numerical failures in an owned result usable by downstream code.

## Capabilities

### New Capabilities

- `guided-mode-sensitivities` provides directional derivatives of physical guided-mode results.

### Modified Capabilities

None.

## Impact

Extends lamina-mode and the core jet paths required for homogeneous field reconstruction. Existing local characteristic partials, enumeration and profile APIs remain available. The first consumer is a downstream multimode slab propagator.
