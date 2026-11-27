# Changelog

The lamina crates (lamina-units, lamina-core, lamina-materials, lamina-mode) are
released in step until 1.0. Each entry names any shared convention it changes.

## Unreleased

Added:

- Declared spectral singularities: `MeromorphicMaterial::spectral_singularities`
  and `MeromorphicTensorMaterial::tensor_spectral_singularities` list every
  pole and zero of ε and μ in the complex k0 plane (`SpectralSingularity`,
  `SingularityKind`, `ConstitutiveQuantity`), or `None` when a model does not
  declare them (the default). The analytical handles forward the declaration,
  and `AnalyticalMaterialRegion::spectral_singularities` gives it in the stack
  frame (a full-class region declares none). `Constant` declares none unless
  ε or μ is zero.
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
- Tensor materials: `TensorMaterial`, `DifferentiableTensorMaterial`,
  `MeromorphicTensorMaterial` and `DifferentiableMeromorphicTensorMaterial`
  return ε and μ as 3×3 tensors in the material's own frame; every scalar
  material implements them as `ε I`, `μ I`. A model declares
  `PRINCIPAL_AXES` when its tensors are diagonal in that frame. Erased form:
  `AnalyticalTensorMaterialHandle`.
- `Rotation`: a proper rotation `R` from the material frame to the stack frame
  (`z` normal, `x` along `k∥`, `y` invariant), with `T_stack = R T Rᵀ`. Built
  from a matrix validated to `1e-12` (`RotationError`), intrinsic z-y′-z″ Euler
  angles, an optic-axis direction `(θ, φ)` or one of the 24 proper signed axis
  permutations; the matrix is what is recorded. Angles within an ulp of a
  multiple of `π/2` give exact quarter turns.
- `Constitutive` regions (`AnalyticalMaterialRegion`): a scalar material, or a
  tensor material with a rotation. The class (isotropic, diagonal, full) and
  separability are decided structurally from the material's declaration and
  the rotation's exact entries. `AnalyticalTensorStack` mixes isotropic and
  tensor regions in layers and exteriors, including the gaps behind walls;
  `Stack::is_separable`. A scalar stack lifts to a region stack with `From`
  (`StackRegion` names a region).
- Real-axis tensor regions: `TensorMaterialHandle` erases any
  `TensorMaterial`, such as a tabulated or interpolated model, and records its
  `PRINCIPAL_AXES` declaration. `MaterialRegion` and `TensorStack` hold such
  regions beside real-axis scalar materials, with the same class, separability
  and validation as analytical regions. They evaluate values on the real axis
  only (amplitudes and power); derivative and complex-plane evaluation of them
  does not compile. `Tensor3` is re-exported at the crate root.
- Class-preserving stack-frame evaluation of ε and μ through the jet algebra,
  refusing non-finite entries or a vanishing `ε_zz` or `μ_zz`
  (`ConstitutiveError`). A scalar material evaluates as an isotropic region.
- `Rotation::rotate`: the stack-frame image `R T Rᵀ` of a material-frame
  tensor.
- Transfer2 and Scatter2 solve stacks whose ε and μ are diagonal in the stack
  frame, isotropic regions included: TE with
  `κ² = μ_xx (ε_yy k₀² − k∥²/μ_zz)` and admittance `κ/μ_xx`, TM with
  `κ² = ε_xx (μ_yy k₀² − k∥²/ε_zz)` and admittance `κ/ε_xx`. Exterior
  wavevectors are evaluated per polarisation, so the guided thresholds are
  `k∥² = ε_yy μ_zz k₀²` (TE) and `ε_zz μ_yy k₀²` (TM). Walls and PMLs in a
  diagonal exterior use its per-polarisation `κ` and admittance. A full region
  is refused with `ConstitutiveError::Unsupported`, naming the region and its
  class.
- Fields of diagonal stacks take each eliminated component from its own entry
  (`H_z = k∥ E_y/(k₀ μ_zz)`, `E_z = −k∥ H_y/(k₀ ε_zz)`); `D = εE`, `B = μH`,
  energy and dissipation densities, and layer energy, dissipation,
  participation, confinement and QNM normalisation weight each field
  component by its own entry, with per-entry dispersive weights
  `∂(k₀ ε_ii)/∂k₀`. Overlaps and norms use the normal entry in the normal
  component.
- Modal field recovery, pole normalisation and directional derivatives accept
  region stacks; a modal direction takes per-entry stack-frame rates
  `[ε', μ']` (`ConstitutiveRates`) as well as isotropic `(ε', μ')`.

- `Transfer4` and `Scatter4`: 4×4 transfer- and scattering-matrix backends
  for any region, including full ε and μ tensors that couple TE and TM. They
  propagate the tangential state `ψ = (E_x, Z₀H_y, E_y, −Z₀H_x)` with
  `dψ/dz = i k₀ Δ ψ`, `Δ` the Berreman matrix obtained by eliminating `E_z`
  and `H_z` through `ε_zz` and `μ_zz`, and solve the coupled channel, or TE
  and TM of a separable stack. `Transfer4` uses the entire propagator
  `exp(i k₀ d Δ)` (scaling and squaring), exact through layer cutoffs and
  degenerate or defective layer eigenvalues; `Scatter4` composes bounded
  forward and backward eigen-wave amplitudes and refuses a layer whose
  forward and backward eigenvalues meet (`Scatter4Error::DefectiveLayer`).
- `Channel` (`TransverseElectric`, `TransverseMagnetic`, `Coupled`;
  `From<Polarisation>`) and `ChannelError` (`Unsupported`, `NotInvariant`,
  `BranchKind`). The 2×2 backends refuse the coupled channel; the 4×4
  backends refuse TE and TM on a stack with a region that couples them.
- `ExteriorBranch`: an open side's outward normal wavenumbers, `Scalar(κ)`
  for TE or TM and `Pair([κ₁, κ₂])` for the coupled channel (an outgoing wave
  varies as `exp(iκ|z|)` away from the stack, so an isotropic exterior gives
  `(κ, κ)`). `ExteriorWavevectors::from_branches` builds them. Each
  exterior's outgoing subspace is the spectral projector of the supplied
  roots applied to `E_y` and `Z₀H_y`, analytic where two outgoing roots
  coincide and singular where an outgoing root meets an incoming one; its
  conditioning is reported by `CoupledContext::basis_conditioning`.
- Coupled characteristic `det[V_L | M V_R]` of the two boundary subspaces
  (`Transfer4`), and the product of the pivot determinants of the same
  interface system in bounded amplitudes (`Scatter4`). For a separable stack
  it is a nonvanishing gauge times the TE and TM characteristics, so mode
  counts add. Walls close a side by the states they admit
  (`span{Z₀H_y, Z₀H_x}` for a perfect electric wall, `span{E_x, E_y}` for a
  perfect magnetic one), carried through the gap and the PML by the
  exterior's propagator; the PML is normalised by its `k∥`-independent
  reference phase.
- `JonesAmplitudes`: reflection and transmission matrices between exterior
  eigen-waves (TE then TM for an isotropic exterior), through
  `ProjectJonesAmplitudes`; `ProjectPower` gives `PlaneWavePower` of one
  value per incident eigen-wave. Scalar `PlaneWaveAmplitudes`
  stay with the 2×2 backends.
- `CoupledWorkspace::fields`: all six components of `E` and `Z₀H` in any
  region for a driven 4×4 solution.
- Closed-form layer integrals of the 4×4 backends: Hermitian and bilinear
  forms of the tangential state over a layer of thickness `d` from the block
  exponential `exp([[−X, Q], [0, Y]] d)` (Van Loan), entire in the layer
  generator and so exact at degenerate and defective layer eigenvalues,
  including a layer cutoff; over an exterior holding only decaying outgoing
  waves, from the Sylvester equation `X I + I Y = −Q`.
- `CoupledSolution` (`CoupledWorkspace::excitation`): one 4×4 solution with
  its six-component fields, weighted Hermitian and bilinear overlaps, field
  norms, longitudinal power `½ Re ∫ (E_y Z₀H_z* − E_z Z₀H_y*) dz` and interface
  flux, per finite layer, wall gap or decaying exterior.
- `CoupledWorkspace::layer_power`, `layer_energy` (with the Hermitian parts of
  caller-supplied weights: `[ε, μ]`, or `[∂(k₀ε)/∂k₀, ∂(k₀μ)/∂k₀]` for the
  dispersive energy) and `layer_dissipation` (from `(ε − ε†)/2i`), per unit
  incident flux; participation and confinement follow from the energy. On
  separable stacks they equal the 2×2 results.
- `Layers`, `LayerEnergy`, `LayerParticipation`, `EnergyConfinement`,
  `HermitianLayerOverlap` and `BilinearLayerOverlap` are exported.
- Modes of the 4×4 backends: `CoupledWorkspace::mode` gives the null vector
  of the homogeneous interface system (`[V_L | M V_R]` for `Transfer4`, the
  last pivot back-substituted for `Scatter4`) as the state at every
  interface, with the residual `σ_min/σ_max`, and refuses a two-dimensional
  null space (`ModeReconstructionError::DegenerateNullSpace`);
  `CoupledWorkspace::mode_solution` gives its fields, which keep only the
  outgoing waves outside, and its integrals including the exterior tails.
  On separable stacks the modes equal the 2×2 ones with no cross terms.
- Modal directions take full stack-frame rates `[ε', μ']` for the 4×4
  backends (`ModalDerivativeBackend` for `Transfer4` and `Scatter4`, any
  `Into<Channel>` as the direction's channel); the 2×2 backends refuse a full
  rate.
- `Scatter4::continued_from(k₀, k∥)`: each layer's forward/backward split is
  classified at a fixed complex reference point and continued to every
  evaluation point by tracking the eigenvalues along the segment between
  them, so the characteristic is holomorphic over a whole chart region that
  holds no layer branch point; `Scatter4Error::Continuation` refuses an
  evaluation whose segment passes near one. The default keeps the
  classification at the real point.
- Real-axis evaluation of the coupled channel selects the principal
  outgoing pair (decaying, or carrying flux away from the stack) of each open
  exterior; `ExteriorSelectionError`, `RealAxisEvaluationError::Exterior`.

Changed (breaking):

- The determinant and field evaluators (`checked_determinant`, `determinant` and
  friends) require `J: TensorConstitutiveLift<E, M>` in place of
  `ConstitutiveLift`. It holds for every `ScalarMaterial` (every `Material` and
  the erased scalar handles) and for `Constitutive` regions; generic callers
  replace the old bound with it.
- New error variants: `Transfer2Error::Constitutive`,
  `Scatter2Error::Constitutive`, `RealAxisEvaluationError::Constitutive`,
  `LayerEnergyError::Constitutive`, and
  `CoordinateCompileError::IncidentMedium` for an incident angle in an
  exterior that is not isotropic.
- `ModalFieldBackend` and `ModalNormalizationBackend` take the stack's region
  material as a type parameter (defaulting to `AnalyticalMaterialHandle`);
  `ModalCharacteristicDirection` takes the material and rate types, and the
  `ModalDerivativeBackend` methods are generic over them.

- `ExteriorWavevectors::left`/`right` return `Option<&A>`: a wavevector is
  supplied for open sides only (`closed_left`, `closed_right`, `closed`,
  `from_options`). A pattern that disagrees with the terminations is a
  `TerminationError::ExteriorBranchMismatch`.
- `Scatter2`'s error type is `Scatter2Error` (previously uninhabited);
  `Transfer2Error` gains `Termination`.
- `ModalFieldBackend::unscaled_fields` takes `[Option<C>; 2]` exterior roots.
- A scalar material, or an isotropic region, whose ε or μ is not finite at an
  evaluated point fails with `ConstitutiveError::NonFinite` naming the region,
  as diagonal and full regions already did; `Scatter2` previously returned NaN
  amplitudes and `Transfer2` `NonFiniteFinalMatrix`. For every class the check
  now covers each derivative part of the lifted jet as well as its value.
  `Transfer2Error::NonFiniteFinalMatrix` now says the inputs were finite.
- `Backend::solve` and `Backend::retain` take a `Channel` instead of a
  `Polarisation`, and the real-axis and complex-plane evaluators take
  `impl Into<Channel>`; a `Polarisation` argument still converts.
  `RealAxisState::polarisation` and `ComplexPlaneState::polarisation` are
  replaced by `channel`, and `RealAxisPairError::PolarisationMismatch` by
  `ChannelMismatch`.
- `ExteriorWavevectors::left`/`right` return `Option<&ExteriorBranch<A>>`
  (`ExteriorBranch::scalar` gives the TE/TM wavevector).
- `PlaneWaveEntries::ExteriorContext` no longer requires
  `ExteriorContextProvider`.
- `TensorConstitutiveLift` gains `is_separable`.
- `Transfer2Error` and `Scatter2Error` gain `Channel`.
- Real-axis `amplitudes`, `power` and `excitation` return `SolveRequestError`,
  which now wraps `ProjectionConstraintError`.

Open stacks keep their numerical results bit for bit, and isotropic stacks
through every observable.

Changed:

- Scatter4 and Transfer4 evaluate the characteristic with about a fifth fewer
  instructions and the same results: minors are read in place rather than
  copied, `solve` takes the determinant from the cofactors of its adjugate,
  the characteristic polynomial sums the trace of each product without
  forming it, and a projector's idempotency is checked on its values.
- Scatter4 and Transfer4 compute one-sample (`Ix0`) characteristics on point
  jets, whose value and derivatives are stored inline, instead of jets of
  one-element arrays that allocate on every operation: the material tensors,
  coordinates, thicknesses and exterior waves are converted once on entry and
  the results once on exit. A first-order characteristic of a one-layer
  rotated uniaxial film takes 1.1 M instead of 7.8 M instructions with Scatter4
  and 0.6 M instead of 5.9 M with Transfer4, with 20 times fewer allocations;
  results are bit for bit the same. The new `#[doc(hidden)]` trait `InlineJet`,
  implemented for every jet family and dimension, bounds the 4×4 `Backend`
  impls.
- The 4×4 jet matrices take minors, cofactors and blocks through index lists
  on the stack, read a one-sample matrix's values without collecting them per
  entry, and reuse the root orderings of 2×2 and 4×4 generators: a further 10 %
  fewer instructions and 40 % fewer allocations per Scatter4 characteristic.

## 0.1.0 (unreleased)

First release on crates.io.

Conventions fixed by this release:

- the canonical spectral coordinate is vacuum angular wavenumber `k₀ = 2π/λ₀`
  in `cm⁻¹`; in-plane wavenumbers use the same unit.
