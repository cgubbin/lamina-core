# Tasks

## 1. Baseline and checked point evaluation

- [x] 1.1 Record the current core test baseline and pre-existing working-tree changes using `cargo test --all-targets` and `git status --short`; verify unrelated edits are preserved and any baseline failures are identified separately.
- [x] 1.2 Add and export the additive checked determinant method and typed input/backend/output error as designed; verify tests cover raw/checked agreement, finite zero acceptance, non-finite coordinate and exterior jet components, and output derivatives that fail despite a finite primal value.
- [x] 1.3 Verify backend error preservation using the existing Transfer2 non-finite-matrix failure path; assert the checked error retains the backend source rather than relabeling it as a material pole.
- [x] 1.4 Document the checked method's pointwise guarantee, raw API compatibility, and caller branch/domain obligations in its rustdoc; verify `cargo test --doc` and a public API usage test exercise the exported method and error.

## 2. Entire-function algebra for cutoff continuation

- [x] 2.1 Implement internal scalar C(x)=cos(sqrt(x)) and S(x)=sin(sqrt(x))/sqrt(x), including first and second derivatives, with a convergent small-argument series and stable regular formulas; verify exact zero coefficients, complex sample agreement, and numerical agreement across the formula switch.
- [x] 2.2 Lift both entire functions through value-only, first/second directional, and first/second bivariate jet families using the exprel composition pattern; verify chain-rule tests at zero and mixed sampled arrays with both zero and nonzero entries.
- [x] 2.3 Preserve the squared wavevector before taking its square root and expose the internal q and polarization factor needed by transfer propagation; verify q values and derivatives remain finite at cutoff without reconstructing q from a singular sqrt jet, and document this invariant next to the internal representation.

## 3. Transfer characteristic continuation

- [x] 3.1 Rewrite transfer layer construction using C(q*d^2) and S(q*d^2) as specified; verify TE/TM layer values and directional/mixed derivatives at exact interior cutoff against independent series limits, including thickness differentiation.
- [x] 3.2 Add evaluator-level single-layer characteristic regressions at and near cutoff for both polarizations and both varying canonical coordinates; verify against independent slab expressions and finite-difference derivatives away from cutoff, with regular exteriors to isolate the interior behavior.
- [x] 3.3 Add interior square-root-cut crossing and regular cross-backend comparison regressions; verify characteristic values and derivatives follow analytic continuation and existing driven-response tests retain their tolerances. Document that modal field reconstruction at exact cutoff remains outside this guarantee.

## 4. Scattering failure characterization and domain documentation

- [x] 4.1 Capture the audit's intermediate projective-pole collapse in an algebraic regression with a parameterized finite limiting response; verify the test distinguishes the homogeneous representation failure from a physical determinant singularity.
- [x] 4.2 Add a physical finite-layer fixture with finite nonzero constitutive factors and an opposite-admittance interface chart failure; verify the independent/transfer characteristic is finite and the checked Scatter2 path reports non-finite output, rather than substituting the projective denominator or changing the coordinate.
- [x] 4.3 Add a thick evanescent case that verifies checked evaluation either agrees with the independent characteristic when representable or reports an explicit numerical/backend failure; verify it never treats a non-finite component as success or silently switches backends.
- [x] 4.4 Update complex-plane and material-domain documentation with the boundary-residual gauge, inverse-transmission relation, material pole/factor-zero restrictions, removable chart failures, finite-sampling limitations, and fixed/varying spectral obligations; verify `cargo doc --no-deps` and cross-check each statement against the new regressions.
- [x] 4.5 Update the roadmap resume note with delivered guarantees, remaining scattering restrictions, and the next lamina-mode baseline/representation change; verify it retains both search directions, candidate validation, chart orchestration, common geometries, public API, fields, continuation, and performance work.

## 5. Integration verification

- [x] 5.1 Run `cargo test --all-targets`, `cargo test --doc`, `cargo clippy --all-targets`, and `cargo fmt --check`; verify new failures are resolved and separately report any pre-existing baseline issues without modifying unrelated files.
- [x] 5.2 Run `openspec validate establish-modal-search-contract --strict` and review the diff against each delta requirement; verify all in-scope guarantees have evidence and no sibling repository source or unsupported global search guarantee was introduced.

## Workflow follow-up

- Review and archive only after implementation and verification are complete.
- Restore the sibling lamina-mode baseline and plan stable scattering characteristic propagation as the next change before promising complete contour searches.
- Continue with the full roadmap in design.md; these follow-up phases are not tasks required to complete this bounded core change.
