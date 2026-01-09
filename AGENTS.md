# Codex Agent Guide — libsyn_tools (sim module correctness)

## Mission
Harden correctness of `libsyn_tools/sim/**` with small, robust fixes and targeted regression tests in `tests_sim/**`.
Avoid overengineering. Prefer fail-fast validation and consistent semantics.

## Mandatory workflow for every run
1) Read `codex_state.md` (repo state + task log + decisions).
2) Run baseline:
   - `pytest -q tests_sim`
   If failing, fix baseline before new work.
3) Implement the next task (from the task list) with:
   - minimal code changes,
   - clear error messages,
   - at least one regression test.
4) Run:
   - `pytest -q tests_sim`
5) Update `codex_state.md`:
   - mark task DONE,
   - list files changed,
   - summarize behavior change and tests added.

## Constraints
- You can only modify this repo. No external repos.
- Keep changes narrowly scoped and backwards compatible unless task says otherwise.
- Don’t use Python `assert` for runtime validation. Use explicit exceptions.

## Simulator invariants to preserve
- Selection + locking should not strand locks or “steal” pool items after interrupt.
- Presence (`is_present`) is availability; avoid mutating or selecting non-present objects unless explicitly intended.
- Precedents should be validated early; cycles must fail fast.
- SHACL remediation spawners must never fabricate invalid IDs.

## Testing guidance
- Prefer unique pool_type strings in tests to avoid cross-test leakage.
- Use small env timeouts and deterministic sequences.
- When validating deadlocks/cycles: raise early errors, don’t rely on env stalling.

## Commands
- Full sim tests: `pytest -q tests_sim`
- Single test: `pytest -q tests_sim/path/to/test_file.py::test_name`
