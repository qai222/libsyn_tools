# libsyn_tools.sim — Codex Agent Guide (Path A Refactor)

## Mission (what “success” means)
We are improving `libsyn_tools.sim` along **Path A**:
- Keep the atomic `UnitaryEdit` primitives as the *only* low-level graph mutation mechanism.
- Make the simulator stand out by being:
  1) **fast enough** to run many steps (avoid “rebuild the world” costs),
  2) **debuggable** (SHACL failures map to actionable context),
  3) **pleasant to extend** at the *Operation* level (effects DSL/macros), not at the RDF-triple-edit level.

## Non‑negotiables
- ✅ All tests must pass: `PYTHONPATH=. pytest tests_sim`
- ✅ No new third‑party dependencies.
- ✅ Do not break public APIs unless the task explicitly says so (prefer additive changes).
- ✅ Keep changes per task scoped and reviewable (avoid “mega refactors”).

## Where to work
Primary target modules:
- `libsyn_tools/sim/effect_engine.py`
- `libsyn_tools/sim/simulation.py`
- `libsyn_tools/sim/overlay/sppt_overlay.py`
- `libsyn_tools/sim/operation/*`
- `libsyn_tools/sim/operation_preset/*`
- `libsyn_tools/sim/report.py`

## Working style
- Prefer small helper functions over large rewrites.
- Preserve existing behavior unless explicitly improving it (performance refactors must be behavior‑preserving).
- Add/extend tests only when they reduce risk or lock in a new behavior.

## Commands
Run after each task:
- `PYTHONPATH=. pytest tests_sim`

Optional sanity checks (only if cheap):
- `python -m compileall libsyn_tools`

## Logging / Memory (critical)
Codex does NOT reliably remember prior runs. Therefore:
1) At the start of each run, read `codex_state.md`.
2) At the end of each run:
   - Update `codex_state.md`:
     - Mark the task as completed (checkbox).
     - Add a new entry to the Run Log table with date/time, summary, and test result.
     - Note any follow-ups or surprises.

## Definition of Done for each task
- Code compiles.
- Tests pass.
- `codex_state.md` updated.
