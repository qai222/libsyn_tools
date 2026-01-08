# Codex Agent Guide — libsyn_tools (focus: sim module)

## Scope
You are working **only** in this repository. Do not assume access to other repos or internal packages.
Primary focus is `libsyn_tools/sim/**` and `tests_sim/**`.

## First step every run
1) Read `codex_state.md` to understand:
   - what has been done,
   - what remains,
   - any decisions/constraints.
2) Run the simulator test suite to establish a baseline:
   - `pytest -q tests_sim`

If baseline fails, stop and fix baseline before starting new tasks.

## Repo layout (expected)
- `libsyn_tools/sim/` — simulator implementation
- `tests_sim/` — simulator tests (all currently pass)

## Coding principles for this repo
- Prefer **small, isolated changes** with targeted tests.
- Add tests that reproduce the bug/crash first (or at least in the same change).
- Preserve backwards compatibility unless a task explicitly changes semantics.
- Avoid direct mutation of SimPy internals unless absolutely necessary (e.g. `store.items.append/insert`).
- Ensure interrupts/cancellation paths do not:
  - crash the SimPy environment,
  - strand locks,
  - “steal” pool objects via stale `store.get()` events.

## Test strategy
- For each correctness fix, add a regression test under `tests_sim/`.
- Keep tests deterministic: use unique pool types per test to avoid cross-test leakage.
- Use short timeouts and minimal simulation steps.
- After each task:
  - Run `pytest -q tests_sim`
  - Ensure any new tests are stable across multiple runs if possible.

## Documentation & logging
- After completing each task, update `codex_state.md`:
  - Mark the task as done
  - Summarize changes (files, key logic)
  - Note any new decisions or follow-up items

## Style & safety
- Keep error messages actionable (include operation id, object id, precedent id, etc.).
- Prefer raising a domain error (e.g. `EngineMechanicalError`, `ContractViolationError`, or `ValueError`) over letting `KeyError/AttributeError` leak.
- If you introduce new event types or callback behaviors, add tests and document them in `codex_state.md`.

## Commands
- Run all sim tests:
  - `pytest -q tests_sim`
- Run a single test:
  - `pytest -q tests_sim/path/to/test_file.py::test_name`
