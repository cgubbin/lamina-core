# Design: checked determinant foundation and continuation record

## Context

Source audit and planning on 2026-10-18; see proposal.md for motivation and specs/modal-search-domain/spec.md for the behavior contract. No Rust source was changed and no Cargo test suite was run during planning. Both sibling repositories were read; this change implements only the core foundation.

Confirmed user direction: develop the rough mode implementation into a usable solver; neither fixed real frequency nor fixed real in-plane wavenumber needs priority; preserve the roadmap between sessions. The user requested proceeding from the exploration. The decisions below define the first bounded change; later solver work stays on the roadmap.

## Goals / Non-Goals

Establish the characteristic function and its admissibility conditions before making root-count or completeness claims. Separate mathematical analyticity, numerical evaluability, root acceptance, and search coverage.

Goals: add a checked point contract and repair the analytically understood transfer cutoff. Non-goals: automatic material-domain certification, scattering chart repair, full search orchestration, or edits to sibling repositories. Remaining scattering representation failures are explicitly reported when they produce non-finite outputs; finite outputs alone are not certified correct search domains.

## Decisions

### Additive checked API

Add `ComplexPlaneEvaluator::checked_determinant` returning `Result<PlaneWaveDeterminant<J>, DeterminantEvaluationError<B::Error>>` for the built-in determinant projection. Error variants are `NonFiniteInput`, `Backend { source }`, and `NonFiniteOutput`. Use `ScalarAlgebra::all_finite` on both coordinates, both exterior wavevectors, and the full output jet; exact finite zeros succeed. Keep `determinant`, raw backend APIs, and `Scatter2::Error = Infallible` unchanged. The wrapper separates numerical checks from backend failure without changing every existing consumer. Do not silently fall back to Transfer2, shift coordinates, or normalize away a singularity.

This API does not check material analyticity, infer physical radiation sheets, or certify that supplied exterior wavevectors obey the material dispersion relation. Those remain explicit caller obligations. Material numerical failures reflected in the determinant are output failures, not evidence identifying a particular physical pole.

### Continue transfer propagation as functions of squared wavevector

Compute and retain q=epsilon*mu*k0^2-k_parallel^2 before its square root. Introduce internal entire-function jet operations C(x)=cos(sqrt(x)) and S(x)=sin(sqrt(x))/sqrt(x), x=q*d^2, with series around zero and stable scalar evaluation away from zero. Evaluate their first and second scalar derivatives and compose across all existing jet families, following the exprel implementation pattern. All coefficients must use the same analytic continuation; branch decisions select numerical formulas only and must agree across their switch boundary.

Construct the transfer layer as m11=m22=C(x), m12=-f*d*S(x), m21=(q*d/f)*S(x). Obtain q and f directly without algebraically squaring the already differentiated sqrt(q). At q=0 the unused cached kappa jet may contain non-finite derivatives; the transfer determinant path must not consume it. Field reconstruction at cutoff remains outside this contract. Zero f and material poles remain outside the regular continuation domain.

This makes the reference determinant insensitive to finite-layer square-root sign flips and repairs cutoff jets. A pointwise special case applied after sqrt(q) is rejected because it cannot recover lost derivatives. No blanket transfer-backend substitution is introduced for scattering consumers.

### Bound scattering work to characterization and diagnostics

Preserve the existing projection ratio and its gauge. Add regressions for the algebraic projective collapse and a regular finite physical stack with opposite-admittance interface chart failure, comparing to analytic transfer values. Require checked failure where the raw result is non-finite; do not require this change to make those scattering charts evaluable. A separate representation change must repair these cases before the broader mode solver can promise coverage across them. Do not add failing tests whose expected future repair is outside this scope.

### Document admissibility without inventing certification

Document the full-stack holomorphic-neighborhood and zero-free-contour obligations, fixed versus varying spectral domains, branch consistency, transmission-chart exclusions, material factor zeros, and overflow/underflow limitations. A finite sample is not an analytic-domain certificate or root-acceptance test. No domain-token type or material metadata trait is added here: those would require downstream orchestration and a separate API decision.

### Verification boundaries

Verify first/second directional and mixed jet continuation at exact cutoff against power-series identities, plus off-cutoff independent slab characteristics for both TE/TM and both varying coordinates. Compare built-in backends only on regular well-conditioned stacks. Test numerical failures, backend error preservation, finite zero acceptance, and input jets with finite values but non-finite derivatives. Use existing test support; no argloc production dependency or sibling test requirement is introduced.

## Audit evidence

### 1. Use the outgoing boundary residual as the mathematical reference

Sources: `src/backend/transfer2/entries.rs::right_gauged_mode_candidate`, `src/backend/transfer2/state.rs`, `src/backend/scatter2/projection.rs::project_determinant`.

Let xi_L = -i Y_L, xi_R = -i Y_R and let M map right boundary state to left. With unit right outgoing field the reference is

```
D = [xi_L, -1] M [1, -xi_R]^T.
```

The scattering projection implements D = 2 xi_L d / n21 = 2 xi_L / t_left, wherever those expressions are defined. Physical scattering poles are intended zeros of this characteristic; their existence alone does not violate the zero-finder contract.

Do not replace D by the raw projective denominator d without proving the removed factors are analytic and nowhere zero throughout the domain. Common projective factors can create spurious zeros; cancelling them can also be required to avoid 0/0. Holomorphic nonvanishing scaling preserves zero counts. A pointwise absolute-value normalization is not automatically a valid analytic search function.

Proposed role for Transfer2: a reference on well-conditioned problems, not the default cure for scattering failures, because its documented exponential instability remains.

### 2. Chart validity is only one part of domain validity

The sibling `lamina-mode/src/conformal/chart.rs` validates/continues the linear physical coordinate and explicit exterior wavevectors. It does not establish full-stack determinant analyticity. The dispersion adapter forwards the chart's wavevectors correctly, but maps can be constructed independently from the stack.

Require the selected exterior data to satisfy kappa_side^2 = epsilon_side mu_side k0^2 - k_parallel^2 and preserve the chosen radiation sheet and jet derivatives. Derive map coefficients from the stack where possible. Complex-frequency maps currently require constant exterior constitutive products; automatic physical region builders additionally require positive real products and reject some common degenerate geometries. Explicit charts and automatic regions have different support boundaries.

### 3. Several numerical failures are removable, not physical poles

- `src/backend/transfer2/entries.rs::from_layer` computes sin(kappa d) times factor/kappa. For finite nonzero factor its limit at kappa=0 is factor*d, but direct evaluation divides by zero. Jet sqrt derivatives are also singular there before later cancellations.
- With q=kappa^2 and factor f=mu for TE or epsilon for TM, the analytic layer entries can be written cos(d sqrt(q)), -f*d*sinc(d sqrt(q)), and (q*d/f)*sinc(d sqrt(q)). These are entire in q for fixed finite nonzero f. This explains why interior square-root sign changes need not be physical branch boundaries. Value and derivative evaluation must implement the removable continuation, not merely rely on symbolic cancellation. Zeros of f need separate analysis.
- `src/backend/scatter2/backend.rs::interface` divides by Y_left+Y_right before conversion to the projective representation. Projective composition therefore does not remove interface-level failures already materialized as infinities.
- The scattering determinant divides by n21. Its zeros can indicate a removable chart failure, a genuine singularity, or numerical underflow, depending on accompanying factors; classify rather than blanket-exclude them as physical poles.
- `src/backend/scatter2/projection.rs::cascade_projection` includes d_left*d_right in the output denominator and transmission numerators. At a pole of an intermediate network it can collapse to the all-zero homogeneous representation. Example using its exact formula: left (d,n11,n12,n21,n22)=(0,1,2,3,6), right=(1,.2,.7,.8,.1), produces all zeros. The corresponding rational cascade has a finite limiting response because 1*6-2*3=0 cancels the leading reflected term. This is an algebraic counterexample, not yet a physical stack regression.
- `Scatter2::Error = Infallible`; successful Result does not establish finite determinant values or finite derivatives.

Algebraic spot checks were evaluated directly: the cascade example produces (0,0,0,0,0), and transfer m12 at kappa=1e-7, f=2, d=.3 is approximately -.6, agreeing with its cutoff limit. These checks reproduce formulas; they do not test Rust implementations or establish all-stack behavior.

Recommendation: repair representational continuation before broad exclusions. Excluding every interior cutoff, interface denominator zero, or intermediate pole would omit otherwise valid searches and confuse numerical charts with physical singularities.

### 4. Meromorphic materials do not certify a holomorphic determinant

Sources: `src/material/mod.rs::MeromorphicMaterial`, sibling material models, and `argloc/src/lib.rs`.

The material traits supply values/derivatives but no complex-domain singularity certificate. Real fit-interval metadata is not such a certificate. Material poles, zeros of the polarization factor, exterior degeneracies, and their preimages need explicit treatment. Composition through propagation can turn a material pole into a more complicated singularity; do not assume every full-stack singularity is a finite-order pole removable by multiplying a known polynomial.

Initial admissibility contract should require a neighborhood of the closed search region on which the complete characteristic is holomorphic, plus zero-free integration contours. Material poles and other unhandled singularities must be outside that neighborhood. At fixed frequency, material validity is checked at that frequency; at varying frequency, check the entire mapped domain. Sampling only the contour is insufficient to certify an interior free of singularities.

If generic materials cannot provide domain information, distinguish caller-asserted validity from automatically established validity. Decide that API during proposal work; do not pretend finite sampling proves analyticity.

### 5. Keep zero counting separate from candidate validation

The local `argloc` documentation explicitly gives find_zeros a no-poles-in-domain assumption. A negative signed count may expose a violation, but a nonnegative count cannot certify its absence: zeros and poles can cancel. Multi-target estimates can be centroids rather than individual roots.

`lamina-mode/src/search.rs` currently maps every estimate to a mode; `mode.rs` retains only z, the varying coordinate, and raw residual. Future acceptance needs multiplicity/cluster handling, physical-coordinate convergence, boundary residual evidence, branch provenance, and retained unresolved regions/termination status. These are downstream work, not a substitute for a valid determinant domain.

## Implementation sequence

1. Establish the current core baseline without altering unrelated working-tree edits.
2. Land checked point evaluation, error tests, and its point-versus-domain documentation together.
3. Land entire-function algebra and derivative tests, then transfer propagation and physical cutoff/branch regressions.
4. Characterize unresolved scattering charts and document explicit restrictions; verify compatibility at regular points.
5. Run core integration checks and preserve the downstream handoff in the roadmap.

## Migration Plan

The checked API is opt-in; existing raw calls continue compiling. Regular determinant gauge and values remain unchanged. Improved transfer cutoff behavior extends the evaluated domain without changing the meaning of existing finite results. No release, deploy, sibling migration, or archive is required to verify this change. Revert the additive API and cutoff edits independently if regressions occur; preserve the investigation record.

## Future roadmap and resume point

Preserved order from the conversation:

1. Coherent build/dependency baseline.
2. Determinant and domain contract (current investigation).
3. Bind stack, polarization, fixed coordinate, charts, and radiation sheets into one consistent problem.
4. Trustworthy single-chart candidates, refinement, multiplicities, residuals, partial results.
5. Orchestrate charts for both search directions; reconcile overlap without merging distinct radiation sheets; support symmetric exteriors and normal incidence.
6. Public physical-input API, sensible defaults, complete examples, real documentation.
7. Reconstruct fields through core while separating root acceptance, reconstruction, and normalization failures.
8. Continuation/mode tracking, then profiling-led caching and performance work.

Milestones: trustworthy search kernel; usable solver; practical dispersion tool. Independent physical validation spans every milestone. Neither search direction dictates the shared architecture.

Implementation handoff: core now provides the additive checked determinant API and transfer continuation through interior cutoff, with value/derivative regressions, a physical opposite-admittance chart-failure fixture, the parameterized projective-collapse example, and thick evanescent cases. Exact-cutoff modal field reconstruction and finite accuracy/underflow certification remain outside these guarantees. Existing raw consumers retain their API.

The next change should restore the lamina-mode dependency baseline and plan stable scattering characteristic propagation using these regressions. Scattering interface and intermediate projective chart failures are characterized, not repaired; checked numerical success is not a domain certificate. Repair those representations before promising coverage across them. Then connect the checked evaluator and explicit caller domain assumptions to argloc, followed by the remaining roadmap above. Both complex spectral and complex in-plane searches use the same foundation.

## Risks / Trade-offs

- Domain restrictions can hide representation defects: remove proven removable failures and report remaining exclusions explicitly.
- Transfer reference can overflow: limit comparisons to well-conditioned cases and use independent physical equations too.
- Arbitrary projective rescaling can corrupt analytic derivatives/counts: prove scale properties and cancellation, including jet behavior.
- Automatic certification is unavailable for generic material traits: expose evidence/assumptions honestly.
- A finite checked result is weaker than a valid search domain: make this distinction prominent in API documentation and the downstream handoff.

## Implementation verification

Verified on 2026-10-18 using the repository's existing cached Nix development environment. Baseline: 1,459 unit tests and one retained-fields integration test passed. Final: 1,468 unit tests, seven checked-determinant integration tests, and the retained-fields integration test passed. Doctests: one passed, two existing ignored examples. `cargo clippy --all-targets`, `cargo doc --no-deps`, `git diff --check`, and strict OpenSpec validation succeeded.

All touched Rust files pass formatting checks. Repository-wide `cargo fmt --check` still reports existing differences in untouched tracked files: src/evaluate/real_axis/excitation.rs, src/evaluate/tests/retain/fields/electromagnetic_modal.rs, and src/observable/mod.rs. Documentation retains two existing unresolved links: TransferStabilityCheck in src/lib.rs and Parameter::InPlane in src/input/coordinate.rs. These unrelated issues were preserved. The pre-existing README, real-axis evaluator/state, input module, and spatial module/resolved edits were checked against the saved starting diff and remain unchanged; existing changes in shared export files were retained alongside the additive checked API exports.
