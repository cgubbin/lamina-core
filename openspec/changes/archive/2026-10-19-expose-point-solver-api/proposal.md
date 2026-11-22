# Proposal

## Why

The internal solver now handles fixed-frequency propagation modes and both constant-exterior conformal and dispersive direct-frequency spectral modes, but downstream crates cannot call it. A public point-solver interface and runnable examples are the next step toward device characterization, fields, branch continuation and design optimization.

## What Changes

- Add a deliberate public facade for the three existing strategies, with checked physical construction, problem-owned domains and finite chart plans. Keep raw maps, seeded jets and compiled solver types internal.
- Reuse core analytical stacks, materials, layers, polarization, backends and unit conventions. Preserve construction evidence for constant conformal exteriors; do not sample arbitrary exterior handles to infer constancy.
- Expose explicit branch/radiation selection, analytic assumptions, exclusions and numerical controls with distinguishable conformal/frequency inputs and coordinate names/units.
- Provide owned point-solve reports callers can return or store, retaining chart discovery outcomes, all candidate outcomes, accepted observations, duplicate groups, coverage and exact provenance. Keep discovery/refinement inspectable and physical interpretation separately fallible.
- Expose structured, inspectable error causes and the existing supported physical descriptions without converting estimates to accepted modes or promising global coverage.
- Add executable guided-mode, constant-exterior cavity and dispersive-substrate examples using only public imports, plus external integration tests and a compiled quickstart.
- Update the roadmap: after this public interface come fields/local derivatives, branch continuation, QNM normalization, features/sensitivities and optimization in lamina-design.

Fields, sweeping/continuation, normalization, mode volumes, new root algorithms, automatic singularity discovery, mixed-strategy reconciliation, a CLI and package publication are outside this change.

## Capabilities

### New Capabilities

- `modal-point-solver-api`: External construction, execution and inspection of finite point solves through public physical inputs, owned evidence-bearing reports, accessible diagnostics and runnable examples.

### Modified Capabilities

None. Existing numerical, domain, classification and coverage contracts remain authoritative; this change exposes them without changing their meaning.

## Impact

Implementation is primarily in sibling lamina-mode: a public facade and re-exports in `src/lib.rs`, small ownership/accessor adaptations around `problem`, `refinement` and `orchestration`, public integration tests, examples, README and the included mathematical guide. OpenSpec planning remains in lamina-core. No new core primitives, upstream numerical changes or runtime dependencies are expected. This introduces the first supported public solver surface; private type names are not a compatibility obligation. The baseline is 260 passing lamina-mode tests and six existing clippy warnings.
