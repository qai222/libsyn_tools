# AGENTS.md — Codex operating instructions (Phase 3: closed-loop spawners)

## Mission (Phase 3)
Make `libsyn_tools.sim` a *closed-loop* discrete-event lab digital twin:
- Canonical state is the KG
- Correctness is enforced with SHACL contracts (Phase 2 done)
- **Endogenous spawners diagnose + repair protocol failures** (Phase 3 focus)

Phase 3 deliverables:
1) Default overlay wiring (SPPT + currentVolume) so SHACL contracts can use derived facts by default.
2) A violation-driven repair spawner that consumes `LifecycleCallbacks.on_violation`
   and spawns remediation operations (with correct precedents).
3) A small remediation “standard library” (starting with overfill remediation using DrainExcess).
4) Optional: a validation-polling spawner that emits violation records (for invariants not tied to an op apply).

Do NOT implement Phase 4+ features (KG query selector, planning/scheduling adapters, big chemistry pack).

---

## Repo constraints
- Only this repo is available.
- tests exist under `tests_sim/` and currently pass with: `PYTHONPATH=. pytest tests_sim`.
- Your changes MUST keep `tests_sim` passing.

---

## Required workflow (EVERY run)
1) Read `codex_state.md` first.
2) Execute exactly ONE prompt (PH3-0, PH3-1, ...).
3) Run tests:
   - Minimum: `PYTHONPATH=. pytest tests_sim`
4) Update `codex_state.md`:
   - what changed and why
   - files touched
   - tests run + results
   - next prompt to execute

If blocked, log the blocker with enough detail to continue later, then stop.

---

## Guardrails
- Keep diffs small and goal-driven.
- Never monkey-patch Simulation/OperationProcess methods (use LifecycleCallbacks).
- Remediation ops spawned due to a violation MUST wait on the originating operation’s done_event:
  spawn with `precedents=[record.operation_id]` to avoid lock contention and ordering issues.
- Avoid infinite loops: repair spawner must dedupe repeated violations (by violation_id and/or op_id+shape_iri).

---

## Likely files to edit
- `libsyn_tools/sim/simulation.py`           (default overlay registration)
- `libsyn_tools/sim/spawner.py`              (new spawner(s))
- `libsyn_tools/sim/operation_preset/*`      (remediation helpers)
- `tests_sim/*`                              (new regression tests)
