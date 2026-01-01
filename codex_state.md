# codex_state.md — pre-Phase-2 patch memory + log (libsyn_tools.sim)

## Scope
This file is Codex’s persistent memory across runs for the **pre-Phase-2 patch set**.

Goal: finish Phase-1 kernel correctness so Phase-2 work (contracts/transactions/remediation) won’t sit on shaky semantics.

---

## Current repo state (starting point)
- Phase-1 refactor is already applied; user reports all tests pass.
- Runtime context exists per simpy.Environment.
- Filter stores are per env and enforce "available objects only".
- Spawner timeouts are scaled via speed factor.

---

## Known remaining issues to fix (pre-Phase-2)

### Issue A — scheduled_start_time uses inconsistent time units
Current code does:
- `delay = operation.scheduled_start_time - env.now`
- then `timeout(sim_time(delay))`

This mixes base-time and sim-time, especially when `simulation_speed_factor != 1`.
Reference: `OperationProcess._run_core` scheduled start gate.

Expected rule (documented):
- speed factor multiplies **base time** durations to obtain sim-time delays.

### Issue B — duplicate literal bindings can self-deadlock
- `LiteralSelector.resolve()` removes the specific object from its FilterStore via `store.get(candidate is obj)` and then locks it.
- If an Operation requests the same literal object twice (common when user passes a full IRI string twice), the second request can block forever waiting for the object to reappear in the pool — but it cannot reappear until the Operation finishes, so this is a self-deadlock.
- `Operation._pre_act_implementation` only short-circuits duplicates when `spec (raw string)` is in the acquired map; this is insufficient when the raw spec differs from the returned identifier (e.g., full IRI vs local identifier), or when specs are LiteralSelector instances.

### Issue C — Mechanical precheck does not enforce lock coverage for data edits
- `_precheck_mechanical` enforces "write coverage requires lock" only for ADD/REMOVE_OBJECT_PROPERTY, not for ADD/CHANGE_DATA_PROPERTY or ANNIHILATE.
- For Phase-2 contract work, unlocked writes on runtime-tracked objects must be prevented to avoid non-determinism.

---

## Phase exit criteria (must be true before starting Phase 2)
1) scheduled_start_time semantics are consistent with `simulation_speed_factor` and covered by tests.
2) duplicate literal participant bindings do NOT hang; they deterministically dedup or raise a clear error, covered by tests.
3) `_precheck_mechanical` requires locks for any edit that mutates an existing runtime-tracked object (LabObject), not only object properties, covered by tests.

---

## Status / Next prompt
- Next prompt to run: PREP0
- Completed prompts: (none yet)
- Blockers: (none yet)

---

## Execution log (append newest entries at top)

### [UNSTARTED] Initial state
- Tests: user reports all tests pass.
- No pre-Phase-2 patch work has been done.

### YYYY-MM-DD HH:MM
**Prompt executed:** PREP{N} — <title>
**Summary:**
- What changed:
- Why:
**Files changed:**
- ...
**Tests run:**
- ...
**Result:**
- ...
**Notes / follow-ups:**
- ...
**Next prompt:** PREP{N+1}
