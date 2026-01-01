# agents.md — Codex instructions (pre-Phase-2 hardening)

## Mission
Finalize Phase-1 kernel correctness so we can safely start Phase-2 (“SHACL contracts as enforceable policy + transactional semantics”).

This pre-Phase-2 patch set MUST address:
1) scheduled_start_time time scaling correctness
2) duplicate literal participant bindings (no self-deadlock)
3) mechanical precheck: write coverage requires lock for runtime-tracked objects (not only object properties)

Do NOT implement Phase-2 features (transactions, rollback, policy bundles, remediation libraries).

---

## Repo constraints
- Only this repo is available.
- `tests_sim/` exists and currently passes.
- Changes must keep the entire test suite passing.
- Prefer small diffs; avoid refactors unrelated to the prompt.

---

## Required workflow (EVERY run)
1) Read `codex_state.md` first.
2) Execute exactly one task prompt.
3) Run tests `PYTHONPATH=. pytest tests_sim`.
4) Update `codex_state.md`:
   - what changed + why
   - files touched
   - tests run + result
   - next prompt to run

If you hit a blocker, record it in `codex_state.md` with enough info to resume, then stop.

---

## Guardrails
- Do not introduce infinite waits in tests. If validating a “hang fix”, use `Simulation.run(until=...)` and assert completion/non-completion deterministically.
- Keep semantics consistent with docs:
  `simulation_speed_factor` multiplies base-time durations to produce sim-time timeouts.
- When enforcing lock coverage, apply it only to runtime-tracked objects (LabObject) to avoid breaking POM bookkeeping edits.

---

## Likely target files
- `libsyn_tools/sim/simulation.py` (docstring for speed factor)
- `libsyn_tools/sim/operation_process.py` or wherever `OperationProcess._run_core` lives
- `libsyn_tools/sim/operation/operation.py` (`Operation._pre_act_implementation`)
- `libsyn_tools/sim/operation/selector.py` (`LiteralSelector.resolve`)
- `libsyn_tools/sim/effect_engine.py` (`_precheck_mechanical`)
- `tests_sim/...` (new regression tests)
