# Codex Agent Guide — libsyn_tools (sim correctness)

## Scope
Work only in this repository. Focus on:
- `libsyn_tools/sim/**`
- `tests_sim/**`

Do not assume access to external repos.

## Mandatory steps every run
1) Read `codex_state.md` (current state + decisions + task log).
2) Run baseline tests:
   - `pytest -q tests_sim`
   If baseline fails, stop and fix baseline first.
3) Implement the next task in the task list (below).
4) Add/extend regression tests in `tests_sim/`.
5) Run:
   - `pytest -q tests_sim`
6) Update `codex_state.md` with:
   - task status (DONE),
   - files changed,
   - tests added/updated,
   - any behavior/semantics decisions.

## Design constraints / invariants
- Never strand locks or pool items after interrupt/abort/exception.
- Presence (`is_present`) is availability. Do not select or mutate non-present objects unless explicitly creating them.
- Precedent dependency graphs must fail fast if invalid (missing IDs or cycles).
- Spawners must not fabricate invalid focus IDs (no "None", no blank nodes, no missing KG objects).
- Avoid Python `assert` for runtime validation; use explicit exceptions.

## Error handling policy
- Prefer domain errors with descriptive messages (include operation id, object id, precedent id, shape id).
- Unexpected exceptions inside operations should not deadlock the sim; always run cleanup.

## Testing guidance
- Use unique pool types per test to avoid cross-test contamination.
- Prefer deterministic small simulations and short timeouts.
- For concurrency/interrupt tests, explicitly schedule interrupts and insertions.
- Add one regression test per bug class.

## Useful commands
- All sim tests: `pytest -q tests_sim`
- Single test: `pytest -q tests_sim/test_file.py::test_name`
