# Tasks

## 1. Public problems, domains and controls

- [x] 1.1 Add the public facade and selected re-exports for propagation, constant-exterior conformal spectral and direct-frequency problems, hiding raw maps, jets and compiled evaluator types. Verify external integration tests construct all three from core inputs, preserve constant-exterior construction evidence and compile with Scatter2 and supported Transfer2 selection; document strategy and backend selection in public rustdoc.
- [x] 1.2 Expose distinguishable conformal/frequency regions, native coordinate/evidence terminology, branch and radiation references, required analytic declarations, known exclusions and continuation controls. Verify public tests check units, direct dF/dk0 versus conformal dF/dz, explicit-complex versus verified-real provenance, invalid settings/exclusions and omitted declarations; document units and assumptions on the public inputs.
- [x] 1.3 Add public explicit finite plans and existing supported canonical conformal builders with validated discovery/refinement/matching settings and staged discovery inspection. Verify public tests cover stable identities through clones, invalid/empty/repeated/foreign associations, unsupported canonical geometry requiring explicit bounds and unchanged requested/halo scope; document the staged flow and finite coverage.

## 2. Owned reports and accepted points

- [x] 2.1 Factor the minimal owned report payload needed to retain discovery and all refinement outcomes without borrowing local execution objects, duplicating solver execution or extending lifetimes unsafely. Verify tests compare retained settings, candidate outcomes, work/evidence, enclosures, multiplicities, leaves, summaries and terminations with staged execution; document ownership and source-index relationships.
- [x] 2.2 Expose a convenience finite-plan solve and owned accepted point values with controlled acceptance construction and exact source provenance. Verify an external helper can return a report after dropping its problem/plan/discovery locals and a point can outlive report views while retaining evidence and domain identity; document point lifetime and the distinction from branch identity or retained fields.
- [x] 2.3 Preserve per-chart errors, accepted observations, duplicate groups, diagnostics and exclusions in owned reports. Verify public integration tests cover overlapping success/failure, all overlap contributors, incompatible exterior sheets and budget-limited clusters without accepted centroids or inferred completeness; document report inspection and coverage limits.

## 3. Public diagnostics and interpretation

- [x] 3.1 Add public structured diagnostic adaptation for construction, domain/configuration, discovery, candidate and checked optical failures while preserving available original sources. Verify external tests inspect continuation side/kind/work and checked input/backend/output causes without private imports, and confirm failures do not become strings, empty success or silent fallback; document the stage/source error model.
- [x] 3.2 Expose separately fallible owned physical descriptions for accepted points/observations, retaining component data, applied tolerances, supported classification and radiation provenance. Verify public tests cover a guided point, outgoing and unspecified spectral references, unsupported labels and unavailable physical data without losing acceptance, groups or discovery; document the support matrix and absence of QNM normalization.

## 4. Runnable examples and public guide

- [x] 4.1 Add `guided_mode` and `conformal_cavity` executables using public imports and core primitives, with independent expected-root checks and concise acceptance/classification/radiation/coverage output. Run both documented locked offline example commands successfully and add source links plus usage explanations to the guide.
- [x] 4.2 Add `dispersive_substrate` using an actual analytical lamina-materials exterior model, direct-frequency bounds, a real outgoing anchor, explicit assumptions and known exclusions. Run the documented locked offline example command, verify its expected dispersive root without freezing the exterior and document continuation controls and local coverage.
- [x] 4.3 Add a complete non-ignored public quickstart doctest, update README/package description and the included guide to the public names and retained mathematical conventions, and advance the roadmap to fields/local derivatives with later milestones and crate ownership intact. Verify doctests, local source/heading links and rendered public rustdoc; retain labels on any remaining illustrative private snippets.

## 5. Integration verification and review

- [x] 5.1 Run locked offline all-target tests, doctests, clippy, no-dependency rustdoc, changed-file formatting and all three example executables. Record results against the 260-test/six-warning baseline, verify external clients need no private type annotations and run relevant core checks only if a demonstrated primitive change required them.
- [x] 5.2 Review every public API spec scenario and ownership/error conversion, confirm no second solver/device framework, unrelated upstream changes, dependencies or publication step were introduced, and run strict OpenSpec validation. Record the final exported API, verification evidence, remaining limits and fields/local-derivatives handoff.

## Workflow follow-up

- Sync the new public API spec and archive the completed change when requested.
- Continue with mode fields and local derivatives, followed by branches, QNM normalization, features/sensitivities and optimization in lamina-design.
