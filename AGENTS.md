# Codex Agent Guide — libsyn_tools (A/B hotfix + pre-lock checklist tests)

You are working **only** in this repository. Focus on:
- `libsyn_tools/sim/**`
- `tests_sim/**`

This work has two goals:
1) Fix two correctness issues:
   - (A) Cleanup/rollback can strand locks if FilterStore reinsertion raises.
   - (B) Interrupt bookkeeping may mutate KG without actually holding locks.
2) Add tests that validate the pre-lock checklist (time validation, periodic spawners termination rules, selector atomicity, no double-start, overlay strictness policy).

## Every run workflow
1) Read `codex_state_202601112041_AB_prelock.md`.
2) Run baseline tests:
   - `pytest -q tests_sim`
   If baseline fails, fix baseline first.
3) Execute the next incomplete task from `codes_task_prompts_202601112041_AB_prelock.md`.
4) Add/adjust regression tests in `tests_sim/`.
5) Run:
   - `pytest -q tests_sim`
6) Update the codex_state log and checklist.

## Constraints / rules
- Prefer minimal, robust fixes (no redesign of the locking model).
- Cleanup/interrupt paths must **never** strand locks and must **never** crash the sim.
- Do not add blocking waits inside interrupt handlers.
- When suppressing exceptions for robustness, always record diagnostics (history_log or logger).
- No `assert` for correctness validation.
- Keep tests deterministic and fast.

## Testing guidance
- Use unique `pool_type` strings per test to avoid global lookup bleed.
- Use `env.run(until=...)` with a small upper bound to avoid hangs.
