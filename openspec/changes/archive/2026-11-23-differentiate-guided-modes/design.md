# Design

## Context

The current local derivatives API computes physical partials at fixed coordinates and guarded implicit slopes. The profile API already owns the accepted field, sided samples and interface coordinates. Guided-power recovery integrates layer flux, includes both analytic exterior tails and fixes the left-interface primary phase. These are the base quantities for total differentiation.

## Goals / Non-Goals

Support first derivatives along real perturbations of simple modes in the lossless dielectric guided basis. Device propagation, cross-structure overlaps, degenerate-subspace derivatives and leaky-mode normalization use separate contracts.

## Decisions

1. Use one physical tangent with a named parameter and units. Represent transverse geometry by interface velocities in a declared lab frame. Obtain layer thickness rates by differencing adjacent velocities. A thickness convenience direction anchors the left interface and moves subsequent interfaces. The base profile's first interface is its local origin; retain the origin velocity to distinguish translation from deformation.
2. Convert vacuum wavelength to k₀ with a consistent length unit. Combine the spectral material derivative with the user's constitutive tangent at fixed wavelength exactly once, in both finite layers and exteriors. Provide the nonmagnetic index direction through dε = 2n dn. Use existing material handles, with explicit directional inputs when a material model exposes no parameter derivative.
3. Evaluate the characteristic directional derivative at fixed β and form β′ = −Fₚ/Fᵦ. Reuse the accepted sheet for both exterior square roots. Propagate the resulting physical jets through homogeneous field reconstruction, exposing only field-level outputs in lamina-mode. Reuse core primitives and extend a core jet entry point where needed.
4. Transport either lab-fixed or material-attached coordinates, retaining the side at an interface. Differentiate the moving layer integration measures and exponential tail factors. For the guided scale, include −P′/(2P) and the derivative of the phase anchor, −i Im(a′/a), in the logarithmic derivative. Conjugated power is differentiated with respect to a real parameter. Use a common converged integration grid for the primal and tangent evaluations.
5. Validate by correcting perturbed roots on both sides and applying the same power and phase convention. Compare multiple step sizes to distinguish truncation error from solver noise. Return base data and the specific failure when a sample derivative is undefined. Finite differences serve as an independent validation path.

## Risks / Trade-offs

Near cutoff or a close degeneracy, conditioning can invalidate an individual derivative. Preserve diagnostics and avoid assigning a derivative across an ambiguous modal swap. The normal electric or magnetic field can jump at an interface, so a fixed lab sample and a sided moving sample require different treatment.
