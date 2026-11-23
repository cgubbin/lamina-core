# Tasks

## 1. Parameter and coordinate contract

- [x] 1.1 Add physical tangent and observation types; document units, origin, interface side and supported lossless directions. Verify geometry mapping for a thickness change, internal-interface shift and rigid translation.
- [x] 1.2 Implement explicit wavelength and regional constitutive tangents. Check unit conversion, dispersion accounting and the nonmagnetic index direction against independent analytic material derivatives.

## 2. Total propagation-constant derivative

- [x] 2.1 Compose characteristic jets with the physical direction and accepted exterior roots. Verify both signs of propagation and agreement with existing single-thickness slopes and corrected neighbouring slab eigenvalues.
- [x] 2.2 Return owned conditioning evidence and typed failures. Exercise near-cutoff, degenerate denominator, invalid selectors and unsupported loss-inducing directions.

## 3. Power-normalized eigenfield derivative

- [x] 3.1 Carry the total eigenvalue, geometry and material tangent through homogeneous field recovery. Check complex TE/TM components for stationary lab samples, layer fractions and both sides of moving interfaces.
- [x] 3.2 Differentiate signed longitudinal power, tails and phase. Verify zero derivative of normalized power, phase consistency, invariance to seed amplitude and explicit phase-anchor failure.
- [x] 3.3 Publish the owned sensitivity result and a downstream slab example. Compare field and eigenvalue derivatives against a sequence of independently corrected perturbations, including dispersive wavelength, index, internal-interface motion and translation. Run the relevant core and mode regressions and document comparison tolerances.
