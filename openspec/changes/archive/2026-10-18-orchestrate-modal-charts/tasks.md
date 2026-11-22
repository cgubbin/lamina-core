# Tasks

## 1. Execution scope and baseline

- [x] 1.1 Reconcile the core-local CLI scope with the authorized sibling mode executor and preserve working-tree snapshots; run locked/offline baselines against core 1481 and mode 228 passing tests, record differences and inherited warnings, and verify no workflow configuration or lockfile changes.

## 2. Owned plans and finite geometry

- [x] 2.1 Add immutable nonempty ordered chart plans with problem/plan ownership and stable domain identities; test foreign domains, empty plans, duplicate exact-domain entries, clone reuse, and unchanged bounds/halos, documenting deliberate overlapping-domain support.
- [x] 2.2 Add problem-owned canonical plan constructors wrapping existing in-plane and spectral region builders with explicit per-chart reference selections and assumptions; verify one versus lower/upper chart order, exact builder geometry, typed invalid-input errors, and inspectable truncation/margin exclusions.
- [x] 2.3 Add explicit finite plan construction and declared-exclusion metadata; verify symmetric-exterior and normal-incidence charts are usable, unsupported/degenerate canonical construction reports its original cause, and no automatic finite truncation or physical-window completeness is invented, with matching module documentation.

## 3. Staged execution and coverage evidence

- [x] 3.1 Execute plans sequentially into owned indexed discovery batches retaining every successful result or original typed error; test a failure between successful charts, deterministic order, source traversal, foreign batch/plan rejection, and globally invalid configuration before any evaluation.
- [x] 3.2 Add batch refinement borrowing the discovery batch and using exact plan-owned domains; verify accepted, failed, clustered, empty, and budget-limited chart outcomes preserve every source estimate, multiplicity, leaf, summary, elapsed time, termination, and assumption without self-referential ownership or fabricated failed-discovery successes.
- [x] 3.3 Add an inspectable coverage ledger with canonical/declared exclusions and per-chart discovery/refinement evidence; test that overlapping success never erases a failed chart, budget termination, unresolved candidate, or excluded seam and document coverage as relative to the finite plan without a global completeness flag.

## 4. Conservative duplicate groups

- [x] 4.1 Add validated match configuration and accepted-observation branch signatures from original charts; test tolerance units/validation, same-sheet agreement, ordered exterior channels, opposite sheets, near-zero sign ambiguity, and typed signature-evaluation failures without changing local acceptance evidence.
- [x] 4.2 Group supported overlapping observations only with mutual enclosure, coordinate, and branch compatibility; test retained contributor indices/evidence, distinct close roots with incompatible enclosures, no within-chart merges, deterministic representatives, and nontransitive distance chains, documenting unmerged disjoint periodic representations and ambiguity diagnostics.

## 5. Physical integration and handoff

- [x] 5.1 Add TE/TM physical regressions for both search directions through plan discovery/refinement and coverage reporting; verify independent slab roots, supported canonical spectral lower/upper participation, duplicate reconciliation on intentionally overlapping charts, and explicit symmetric/normal-incidence cases while preserving single-chart tests.
- [x] 5.2 Update mode's developer README with executable commands and the problem-to-plan-to-discovery-to-refinement flow; verify commands run and examples distinguish duplicate groups, partial chart failures, finite coverage, and remaining seam/cluster limitations.
- [x] 5.3 Run all-target tests, doctests, Clippy, documentation builds, changed-file formatting, and whitespace checks in affected repositories; verify no new failures/warnings, record established exceptions, and confirm no unrequested upstream source or lockfile changes.
- [x] 5.4 Run `openspec validate orchestrate-modal-charts --strict`, review against the new and existing main contracts, and record final counts, exclusions, and the public physical-input API as the next bounded roadmap phase in design.md before completing implementation.

## Workflow follow-up

- Review the complete result, sync the delta spec, and archive when requested.
- Continue with the public API and executable examples, then fields/normalization, continuation, and profiling-led performance. Keep seam recovery, periodic duplicate equivalence, repeated-root/cluster work, and numerical exclusions visible in the roadmap.
