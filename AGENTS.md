# Codex Agent Guide — libsyn_tools (sim correctness polish)

## Scope
Work only in this repository.
Focus areas:
- `libsyn_tools/sim/**`
- `tests_sim/**`

## Mandatory workflow every run
1) Read `codex_state.md` (state + task log + decisions).
2) Run baseline:
   - `pytest -q tests_sim`
   If baseline fails, fix baseline first.
3) Implement the next task (in the task list).
4) Add/adjust tests in `tests_sim/` for each change.
5) Run:
   - `pytest -q tests_sim`
6) Update `codex_state.md` with:
   - task status DONE
   - files changed
   - tests added/updated
   - behavior/semantics notes

## Constraints / priorities
- Prefer small, targeted fixes.
- Do not overengineer: no complex deadlock detectors or schedulers unless explicitly required.
- Prefer explicit runtime validation over `assert`.
- Keep error messages descriptive (include object id / operation id where possible).
- Ensure report generation and basic simulator utilities behave correctly on edge cases (empty logs, all aborted, etc.).

## Testing guidance
- Add regression tests for each bug.
- Use unique pool types per test to prevent cross-test contamination.
- Keep tests deterministic with short env timeouts.

## Commands
- Full sim tests: `pytest -q tests_sim`
- Single test: `pytest -q tests_sim/test_file.py::test_name`
