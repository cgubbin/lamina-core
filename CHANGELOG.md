# Changelog

The lamina crates (lamina-units, lamina-core, lamina-materials, lamina-mode) are
released in step until 1.0. Each entry names any shared convention it changes.

## Unreleased

Added:

- `Termination`, `Wall` and `WallKind`: either side of a `Stack` may be closed
  by a perfect electric or perfect magnetic wall at an absolute position in the
  stack frame (`z = 0` at the first interface). Open remains the default.
- Validation of walls at compilation (`StackCompileError::Termination`), and a
  typed `StackCompileError::OutwardThicknessDerivative` when a thickness
  derivative would move the last interface past a wall on it.
- Characteristic, modal reconstruction and field sampling of walled stacks in
  both backends; `FieldSamplingError::BeyondTermination`.
- `SolveRequestError::ClosedIncidence` for a plane wave incident from a walled
  side.
- `Pml` and `StretchProfile`: a side may be closed by a wall-backed perfectly
  matched layer by complex coordinate stretching, `s(ζ) = 1 + c (ζ/d)^m`, with
  order and strength required. It enters the characteristic as a layer of
  complex thickness `τ = d (1 + c/(m+1))`, normalised by the
  `k∥`-independent `exp(iκ_ref τ)`; fields inside it follow the stretched
  coordinate. `TerminationError::PmlScaleOutOfRange` reports a scale that
  cannot be represented.

Changed (breaking):

- `ExteriorWavevectors::left`/`right` return `Option<&A>`: a wavevector is
  supplied for open sides only (`closed_left`, `closed_right`, `closed`,
  `from_options`). A pattern that disagrees with the terminations is a
  `TerminationError::ExteriorBranchMismatch`.
- `Scatter2`'s error type is `Scatter2Error` (previously uninhabited);
  `Transfer2Error` gains `Termination`.
- `ModalFieldBackend::unscaled_fields` takes `[Option<C>; 2]` exterior roots.
- Real-axis `amplitudes`, `power` and `excitation` return `SolveRequestError`,
  which now wraps `ProjectionConstraintError`.

Open stacks keep their numerical results bit for bit.

## 0.1.0 (unreleased)

First release on crates.io.

Conventions fixed by this release:

- the canonical spectral coordinate is vacuum angular wavenumber `k₀ = 2π/λ₀`
  in `cm⁻¹`; in-plane wavenumbers use the same unit.
