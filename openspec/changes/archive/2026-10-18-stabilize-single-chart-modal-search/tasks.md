# Tasks

## 1. Execution scope and baselines

- [x] 1.1 Establish an execution context authorized to edit both lamina-core and sibling lamina-mode before changing sibling source; verify the selected roots and writable scope explicitly include lamina-mode, preserve both working trees, and do not create workflow configuration as a side effect.
- [x] 1.2 Reproduce both locked dependency baselines and record the results in design.md; verify core's established passing tests and mode's two known failures against the planning record, with no unrequested lockfile changes.

## 2. Restore enforced single-chart regression behavior

- [x] 2.1 Replace the exact residual/coordinate comparison with the existing fixture's explicit numerical tolerance and assert all three discarded derivative-helper Results; verify evaluation errors fail those regressions and an opposite-sheet sign mismatch remains detectable.
- [x] 2.2 Add a test recorder for the failing adaptive contour evaluation, then separate the requested rectangle from the validated evaluation halo while preserving reference point/coordinate; verify endpoint/corner and affine-roundtrip cases, unchanged search bounds, branchpoint/non-finite halo rejection, genuine outside-point rejection, and the known asymmetric slab search. If the recorded point is materially outside the halo, report the different defect before changing argloc or broadening tolerance.
- [x] 2.3 Remove search dbg output and unused imports, retire uncompiled prototype paths only after verifying no consumers, and document the baseline and contour-versus-chart distinction in lamina-mode; verify its all-target tests enforce every fallible derivative comparison and no known baseline failure is suppressed.

## 3. Carry the complementary scattering determinant

- [x] 3.1 Extend the internal projective representation with h and construct identity, ordinary-entry, and propagation complements without division by a potentially zero projective denominator; verify `b*c-a*e=d*h` and component values/first/second jet derivatives, and document the invariant beside the type.
- [x] 3.2 Replace the redundant-factor cascade with the six reduced recurrences in design.md; verify ordinary Redheffer agreement on regular complex samples, homogeneous-scale invariance including h, associativity, mixed sampled arrays, and directional/bivariate derivative checks.
- [x] 3.3 Replace the archived all-zero-collapse characterization with the analytic intermediate-pole continuation assertion; verify the parameterized limiting response and its derivative at the pole against independently cancelled formulas and preserve a test for a genuine zero characteristic.

## 4. Construct direct projective interfaces and integrate workspaces

- [x] 4.1 Add direct polynomial interface projective components with their analytic complement, and route total response accumulation through them; verify exact opposite-admittance sums do not divide by zero and regular interface values/derivatives agree with Fresnel expressions for both polarizations.
- [x] 4.2 Adapt workspace append and retained component plumbing so optional ordinary reconstruction data cannot contaminate total characteristic accumulation; verify response-only/retained determinant equivalence at repaired exceptional points and existing amplitude, power, field, and layer-observable tests on regular stacks for both incident sides. Document exceptional-point reconstruction limitations with this change.

## 5. Establish physical characteristic evidence

- [x] 5.1 Turn the negative-index layer fixture in tests/checked_determinant.rs into a continuation-success test for Scatter2; verify TE/TM values and first/second spectral and in-plane derivatives against independent slab formulas and Transfer2 on well-conditioned cases.
- [x] 5.2 Add a physical multilayer case whose intermediate network reaches a pole while the total characteristic is regular; verify finite values and derivatives at and around that pole against an independent characteristic, keeping the algebraic pole test as separate evidence.
- [x] 5.3 Recheck thick evanescent, raw/checked gauge, and transfer-cutoff regressions; verify remaining non-finite scattering results are explicit errors, true finite mode zeros still succeed, and documentation accurately distinguishes repaired charts from exact-cutoff directional degeneracy and unrepresentable results.

## 6. Integrate checked optical evaluation in lamina-mode

- [x] 6.1 Forward CompiledPlaneWaveModeSolver through checked_determinant and propagate DeterminantEvaluationError as its associated error; verify a single seeded evaluation returns F and F_prime with the supplied exterior jets unchanged for both conformal map families, including a non-principal sheet case, and document the new internal error chain.
- [x] 6.2 Add injected checked input/output/backend failures at dispersion and adaptive-search boundaries; verify the originating error remains identifiable through DispersionError, argloc evaluation errors, and ModalSearchError, and no panic or successful empty result replaces a failure.
- [x] 6.3 Verify successful budget-limited results preserve estimates, multiplicities, unresolved leaves, run summary, and termination; document estimates versus accepted roots, cluster centroids, and caller pole-free-domain assumptions without promising partial successes that argloc returns only as errors.
- [x] 6.4 Run single-chart physical slab searches for TE/TM in both complex-in-plane and complex-spectral map families using checked evaluation; verify independently known roots, conformal derivative consistency, and branch provenance, and replace the template README with runnable developer commands and current limitations.

## 7. Integration verification and handoff

- [x] 7.1 Run all-target tests, doctests, Clippy, documentation builds, and formatting checks in both authorized repositories; verify new failures are resolved, known unrelated baseline warnings are recorded separately, dependency lockfiles change only if justified, and no source in argloc/materials/units was edited.
- [x] 7.2 Run `openspec validate stabilize-single-chart-modal-search --strict` and review the diff against both delta specs; verify existing core requirements remain intact and update design.md with delivered guarantees, final test counts, remaining exclusions, and the next bounded roadmap phase before marking the coordinated change complete.

## Workflow follow-up

- Review the complete cross-repository result, then sync delta specs and archive after implementation is verified.
- Next address the coherent problem/domain contract and candidate refinement/validation, followed by multi-chart coverage, the public physical-input API, fields, continuation, and profiling-led performance work as recorded in design.md.
