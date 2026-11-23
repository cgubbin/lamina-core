# Design

## Context

The existing physical point solver and local corrector can accept a candidate on an explicitly selected sheet. A row count from its discovery report cannot distinguish an absent mode from a missed one. The completeness statement belongs to the real, bound, lossless problem.

## Goals / Non-Goals

Count the bound spectrum of a piecewise homogeneous positive dielectric for one polarization, wavelength and propagation sign. Complex resonances and open radiation channels remain in the existing point solver. The result is a local basis at one device configuration.

## Decisions

1. Work internally in |β|. Intersect the requested window with the interval above both exterior light lines and below the maximum material light line. Record numerical guard bands separately from physically ineligible regions. A threshold state with a nondecaying tail is outside the guided basis. Use open numerical endpoints and flag an endpoint root within tolerance, so adjacent windows can be reconciled explicitly.
2. Use the scalar Sturm oscillation count with decaying exterior boundary conditions. Propagate an unwrapped phase of the primary field and its weighted normal derivative from the left decaying solution. For TE the continuous pair is Ey and Ey′/μ, for TM it is Hy and Hy′/ε. Account for the β-dependent right decay boundary phase. Preserve complete rotations within thick oscillatory layers, use scaled expressions in evanescent layers and the regular linear limit at zero transverse wavenumber. Document the endpoint count convention before coding. A scan for characteristic sign changes alone can skip a close pair.
3. Use endpoint count differences to split a bracket until its count is one. Correct on the selected bound sheet with the existing checked characteristic. Verify bracket membership, branch decay, numerical acceptance and unique occupancy. Return expected counts as unknown when endpoint phase uncertainty crosses an integer. Completeness is numerical and conditional on the material class, stated interval and tolerances.
4. Keep the count ledger separate from the owned physical modes. Callers can use the modes for power profiles and derivatives while inspecting any unresolved part of the spectrum. Sort by decreasing |β| for both signs; a sorted index has meaning only within the stated interval.

## Risks / Trade-offs

Nearly degenerate supermodes, near-cutoff phase uncertainty and thick barriers can exhaust floating-point resolution. Preserve the bracket, count uncertainty and valid modes. Validate the counter independently against analytic TE/TM ladders, and exercise the failure path with a deliberately restricted budget.
