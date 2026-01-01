# agents.md — Codex operating instructions (Phase 5: productization + chemistry pack)

## Mission (Phase 5)
Turn libsyn_tools.sim into a *usable* and *friendly* digital-twin toolkit:

- Canonical state is a KG
- Correctness is enforced with SHACL contracts
- Endogenous spawners create a closed-loop environment that can diagnose and repair failures
- Phase 5 adds: **reports + UX + standard library ops + optional Chemistry Pack**

Phase 5 deliverables:
1) A first-class RunReport API (human-friendly outputs): summary metrics + CSV exports + optional markdown.
2) “Standard library” operations users will actually use (beyond transfer/drain): at minimum Wait/Hold and Mix/Combine.
3) Optional Chemistry Pack (overlay(s) + helper SHACL shapes + optional endogenous drift model) that makes chemistry queryable (not just JSON bookkeeping).
4) Examples/docs (5-minute quickstart scripts) and lightweight CLI entrypoint for running + producing a report.

Do NOT change Phase 1–4 semantics unless necessary. Avoid refactoring core engine logic.
This phase is about usability and value, not architecture rewrites.

---

## Repo constraints
- Only this repo is available; no other repos.
- tests_sim is the primary regression suite and must remain green:
  `PYTHONPATH=. pytest tests_sim`

---

## Required workflow (EVERY run)
1) Read `codex_state.md` first.
2) Execute exactly ONE Phase 5 prompt (PH5-0, PH5-1, ...).
3) Run tests:
   - Minimum: `PYTHONPATH=. pytest tests_sim`
4) Update `codex_state.md`:
   - what changed and why
   - files touched
   - tests run + results
   - next prompt to execute

If blocked, write the blocker + next-step hints in `codex_state.md` and stop.

---

## Guardrails
- Keep diffs small and reviewable.
- Preserve backwards compatibility:
  - Existing Simulation.export_* methods should keep working.
  - New APIs should be additive (RunReport, new ops, optional packs).
- Any Chemistry Pack must be OPTIONAL (not enabled by default) unless explicitly requested.
- Examples must be deterministic, short, and not depend on internet access.

---

## Likely files to touch
- `libsyn_tools/sim/simulation.py` (new report entrypoints)
- `libsyn_tools/sim/effect_engine.py` (read-only accessors; report helpers)
- `libsyn_tools/sim/report.py` (new)
- `libsyn_tools/sim/operation_preset/*` (new standard ops)
- `libsyn_tools/sim/overlay/*` (chemistry overlay; optional)
- `libsyn_tools/sim/cli.py` + `libsyn_tools/sim/__main__.py` (optional CLI)
- `examples/*` or `docs/*`
- `tests_sim/*` (new tests)
