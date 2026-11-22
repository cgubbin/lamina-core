# Proposal

## Why

Accepted points now have convincing physical eigenvalue checks, but callers cannot recover their actual electromagnetic profiles through lamina-mode. Field recovery is needed for device characterization and later branch identity and QNM work, and must remain usable when normalization is unavailable.

## What Changes

- Add normalization-independent homogeneous field reconstruction in lamina-core, reusing its retained modal seed, layer reconstruction and spatial sampling machinery. Preserve the existing normalized complex-mode API and behavior.
- Expose separately fallible electric/magnetic field recovery from an owned accepted lamina-mode point, its source stack and ordered exterior branches, with an explicit supported backend and documented arbitrary amplitude/phase gauge.
- Reuse core sampling, positions, units and field containers; return owned samples with interface-side identity, source provenance, recovery backend, gauge and reconstruction residual evidence.
- Validate TE/TM dielectric slabs, the benchmark plasmon in both complex cases, and the two coupled modes of a finite symmetric metal film. Include a leaky frequency cavity sampled only at finite exterior distances.
- Publish runnable field profiles and continuity checks without changing numerical acceptance, coverage or unsupported metal classifications.

No new QNM normalization, mode volume, field derivatives, branch tracking, new device representation or driven-excitation substitute is included.

## Capabilities

### New Capabilities

- `modal-field-recovery`: Normalization-independent homogeneous fields associated with accepted modal points, including physical sampling, gauge and separate failure evidence.

### Modified Capabilities

None. Root acceptance, optional classification and the existing normalized core entry point retain their existing contracts.

## Impact

Work spans lamina-core's complex retained-mode/field API and sibling lamina-mode's public accepted-point adapter, diagnostics, tests, examples and guides. Core already reconstructs fields but its public complex-mode construction currently performs bilinear QNM normalization; a small additive unscaled entry point is necessary. No parallel stack, field, sampling or backend framework is introduced.

## Sequence

Apply after `preserve-plasmon-benchmarks`. Fields precede the local-derivative proposal in the working sequence, although characteristic derivatives do not mathematically depend on field recovery. Later work remains branches, explicit QNM normalization/observables, feature sensitivities and optimization in lamina-design.
