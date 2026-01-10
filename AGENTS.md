# Codex Agent Guide — libsyn_tools (sim correctness)

## Scope
Work only in this repository.
Focus on:
- libsyn_tools/sim/**
- tests_sim/**

## Every run must do
1) Read codex_state.md.
2) Run baseline tests:
   pytest -q tests_sim
   If baseline fails, fix baseline first.
3) Execute the next task prompt in order.
4) Add/adjust regression tests in tests_sim for the task.
5) Run:
   pytest -q tests_sim
6) Update codex_state.md:
   - mark task DONE
   - files changed
   - tests added/updated
   - behavior notes/decisions

## Constraints
- Prefer small, surgical fixes; avoid large redesigns.
- Correctness > convenience: avoid silent failures; emit clear errors.
- Do not use assert for runtime validation.
- Keep semantics stable unless codex_state.md records a deliberate change.

## Testing guidance
- Use unique pool_type strings per test to avoid global registry bleed.
- Tests should be deterministic; keep SimPy timings small.
- Add at least 1 regression test per bug class.

## Useful commands
pytest -q tests_sim
pytest -q tests_sim/test_file.py::test_name
