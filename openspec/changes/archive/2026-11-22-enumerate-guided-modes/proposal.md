# Proposal

## Why

A multimode slab propagator needs every guided propagation constant in its operating interval. Existing finite chart searches return accepted points and diagnostics but give the consumer no independent count of the guided states.

## What Changes

- Enumerate TE or TM bound modes of a lossless isotropic dielectric stack at a fixed positive vacuum wavelength, over a declared real propagation-constant interval and propagation direction.
- Compare an independent real-axis mode count with distinct checked physical roots. Report cutoff margins, unresolved intervals and any discrepancy alongside the modes already found.
- Deliver ordered, owned points suitable for guided-power profiles and downstream propagation.

## Capabilities

### New Capabilities

- `guided-mode-enumeration` covers the supported interval, count evidence, accepted modes and completeness diagnostics.

### Modified Capabilities

None.

## Impact

Add an entry point and report in lamina-mode, using lamina-core stack and material handles. Validate with independent slab spectra and a pair of separated guides. The downstream propagation model can consume the resulting basis without running the search itself.
