# Codex Agent Guide — libsyn_tools (sim correctness fixes)

## Scope
Work only in this repository.
Primary directories:
- `libsyn_tools/sim/**`
- `tests_sim/**`

## Required workflow each run
1) Read `codex_state.md` (state + task list + log).
2) Run baseline:
   - `pytest -q tests_sim`
   If baseline fails, fix baseline first.
3) Implement the next task in the task list (below).
4) Add/extend regression tests in `tests_sim/` for each fix.
5) Run:
   - `pytest -q tests_sim`
6) Update `codex_state.md`:
   - mark task DONE
   - list files changed
   - list tests added/updated
   - summarize behavior change and any decisions

## Constraints
- Prefer small, surgical patches.
- Avoid overengineering: no new frameworks, no large refactors.
- Use explicit runtime validation (no asserts for correctness).
- Ensure fixes are robust (clear errors, no hidden AttributeError/TypeError).
- Preserve backwards compatibility unless explicitly stated.

## Testing guidance
- Add at least one regression test per issue.
- Use unique pool_type strings per test to avoid cross-test contamination.
- Tests should be deterministic; keep sim timeouts short.

## Commands
- Full sim tests: `pytest -q tests_sim`
- Single test: `pytest -q tests_sim/test_file.py::test_name`
