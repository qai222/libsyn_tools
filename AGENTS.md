# Codex Agent Guide — libsyn_tools (sim correctness hardening)

You are working **only** in this repository. Focus on **`libsyn_tools/sim/**`** and **`tests_sim/**`**.  
The repo already has `tests_sim` and **all tests currently pass** at the starting point.

## Run protocol (do this every run)
1. Read **`AGENTS.md`** (this file) and **`codex_state.md`**.
2. Run baseline tests:
   - `pytest -q tests_sim`
   If baseline fails, stop and fix baseline before doing new work.
3. Execute the **next incomplete task** from `codes_task_prompts.md`.
4. Add/adjust **regression tests** for the change.
5. Run:
   - `pytest -q tests_sim`
6. Update `codex_state.md`:
   - Mark the task DONE
   - List files changed
   - List tests added/updated
   - Note any behavior/semantics changes

## Constraints & priorities
- Prefer **simple, robust, non-overengineered** fixes.
- Avoid breaking APIs unless explicitly required; if you must, document in `codex_state.md`.
- **No `assert` for runtime validation.** Use explicit exceptions.
- Prevent:
  - leaked locks
  - drained pool stores
  - non-terminating sims (unless explicitly configured)
  - silent masking of correctness errors (make strictness configurable when needed)

## Design invariants to preserve
- **Presence = availability** (`is_present == {True}` required for selection and mutation unless creating).
- Operations must not proceed RUNNING if `pre_act` failed/cancelled.
- Validation/spawners must not fabricate invalid focus IDs.
- Reports/provenance should reflect terminal outcomes (END/ABORT/INTERRUPT).

## Testing guidance
- Use **unique pool types** in tests to avoid global registry bleed.
- Keep tests deterministic; short timeouts.
- Add at least one regression test per bug class.
- When a behavior is semantics-sensitive, encode the decision in a test.

## Commands
- Full suite: `pytest -q tests_sim`
- Single test: `pytest -q tests_sim/path/to/test_file.py::test_name`
