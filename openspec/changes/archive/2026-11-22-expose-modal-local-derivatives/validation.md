# Validation

Run from the corresponding crate directory; all Cargo commands used the lockfile and offline dependencies.

## lamina-core

- `cargo test --locked --offline`: passed, 1,481 unit tests, 8 integration tests, and one documentation test; two documentation tests are intentionally ignored.
- `cargo clippy --all-targets --locked --offline`: passed without warnings.
- `rustfmt --edition 2024 --config skip_children=true --check` on every modified Rust source: passed. `skip_children` avoids traversing an unrelated pre-existing trailing-whitespace error in `src/evaluate/real_axis/excitation.rs`.

## lamina-mode

- `cargo test --locked --offline`: passed, 264 unit tests, 28 integration tests, and four compile-fail/documentation tests; one documentation test is intentionally ignored.
- `cargo run --locked --offline --example local_derivatives`: passed and printed successful partials, backend provenance, a guarded slope, and an out-of-bounds thickness partial as unavailable while preserving the other results.
- `cargo clippy --all-targets --locked --offline`: passed. It reports existing warnings in conformal-region comparisons/type complexity and `Scatter2::default()` in existing tests; no warning points to the new derivative implementation.
- `rustfmt --edition 2024 --check` on every modified Rust source: passed.

## OpenSpec

- `openspec validate expose-modal-local-derivatives --strict --json`: passed with no issues.

Independent checks include branch-anchored finite differences for guided, dispersive, plasmon, film, and thickness derivatives; native-to-physical chain rules at accepted and off-root points; and perturbed-root slope comparisons using fresh propagation or frequency problem domains.
