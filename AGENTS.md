# agents.md — Codex operating instructions (Phase 4: semantic selection + schedule integration)

## Mission (Phase 4)
Deliver the next differentiators:
1) **KgQuerySelector (Level 4 selector)**: select participants using KG/overlay queries (SPARQL) while preserving atomic selection + locking.
2) **Planning/Scheduling → Simulation bridge**: compile `chem_schema` / `opt` artifacts (OperationNetwork + SchedulerOutput) into runnable `sim.Operation` processes with correct timing, precedents, and module/resource locking.

Phase 4 should make it possible to:
- express “pick any available container satisfying semantic constraints” via KG query
- simulate a scheduled plan (start times, durations, assignments) in the DES kernel

Do NOT implement Phase 5 items (big chemistry pack, rich report UI, full replay system).
Do NOT change Phase 1–3 kernel semantics unless required to support Phase 4.

---

## Current repo assumptions
- Phase 3 merged. tests exist under `tests_sim/`.
- Recommended test command: `PYTHONPATH=. pytest tests_sim`.

---

## Required workflow (EVERY run)
1) Read `codex_state.md` first.
2) Execute exactly ONE task prompt (PH4-0, PH4-1, ...).
3) Run tests:
   - Minimum: `PYTHONPATH=. pytest tests_sim`
4) Update `codex_state.md`:
   - what changed and why
   - files touched
   - tests run + results
   - next prompt

If blocked, log the blocker with enough detail to continue later, then stop.

---

## Guardrails
- KgQuerySelector must remain **atomic**: selection removes object from pool (FilterStore) and acquires the lock before returning.
- KgQuerySelector must be deterministic under equal conditions: prefer FIFO order of FilterStore.
- Scheduling bridge must not require other repos. Use only `libsyn_tools.chem_schema` and `libsyn_tools.opt`.
- Keep adapters modular: put schedule bridge under `libsyn_tools/sim/adapters/` (new package).
- Avoid name confusion between `chem_schema.operation.Operation` and `sim.operation.operation.Operation` by aliasing imports in adapter code.

---

## Likely files
- `libsyn_tools/sim/effect_engine.py` (query graph helper)
- `libsyn_tools/sim/simulation.py` (attach effect_engine to env for selectors)
- `libsyn_tools/sim/operation/selector.py` (implement KgQuerySelector)
- `libsyn_tools/sim/adapters/schedule_bridge.py` (new)
- `tests_sim/*` (new regression tests)
