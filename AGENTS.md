# Codex Agent Guide — libsyn_tools v0.3 (sim module)

These tasks target **only two correctness issues** in the simulator:
- (A) `LiteralSelector` can hang indefinitely when the FilterStore is out-of-sync.
- (B) Interrupt bookkeeping can mutate the KG without actually holding locks.

## Every run workflow
1) Read `codex_state_v0.3_AB.md`.
2) Run baseline:
   - `pytest -q tests_sim`
3) Execute the next task in `codes_task_prompts_v0.3_AB.md`.
4) Add/adjust regression tests in `tests_sim/`.
5) Run:
   - `pytest -q tests_sim`
6) Update `codex_state_v0.3_AB.md` with:
   - task status DONE
   - files changed
   - tests added/updated
   - notes/decisions

## Constraints
- Prefer **minimal, robust** changes (no redesign of the locking model).
- Do not introduce new blocking waits inside interrupt handlers.
- Preserve semantics:
  - If an object is legitimately unavailable (not present / locked), selectors may wait.
  - The fix must remove only the *unintended* hang when an object is available but missing from the store.
- No `assert` for correctness validation.

## Testing guidance
- Use unique pool_type strings per test.
- Keep tests deterministic with short sim timeouts.
