# codes_task_prompts_202601112041_AB_prelock.md — Sequential tasks

Run baseline tests before starting each task:
`pytest -q tests_sim`

---

## Task 1 — Fix A: safe cleanup & rollback reinsertion (never strand locks)
**Goal**
Guarantee lock release completes even if pool reinsertion fails, and rollback cannot strand locks/pool items.

**Robust minimal strategy (two-phase + diagnostics)**
- Phase 1: release/cancel all lock requests with **no-throw** guarantees.
- Phase 2: attempt pool reinsertion per object under `try/except`, recording diagnostics; do not re-raise from cleanup paths.
- If pool_type is invalid (missing/multi), ensure object is removed from all stores so it cannot remain “ghost selectable”.

**Implementation requirements**
1) Add a helper used in cleanup/rollback, e.g.:
   - `FilterStoreRegistry.safe_put_obj_into_filter_store(obj, env, *, context: str, diagnostics_sink=...)`
   - OR a local helper in operation/effect_engine that wraps `put_obj_into_filter_store` and logs.
   It must:
   - never raise
   - record error details (object id, context, exception)
   - when `has_pool_type` is invalid, remove object from all stores by identifier (best effort)

2) In `Operation.post_act` and `Operation._release_all_locks`:
   - Ensure lock release loop cannot be cut short by reinsertion errors.
   - Release all locks first; then reinsertion in a separate loop (or continue-on-error per object).
   - Always clear `self.locks` at the end.

3) In rollback restore (`EffectEngine._restore_objects`):
   - Use the same safe reinsertion logic.
   - Ensure rollback does not re-raise due to reinsertion (it may record diagnostics).

**Tests**
Add targeted regression tests:
- `test_cleanup_releases_all_locks_even_if_pool_reinsert_raises`
  - Create an op that acquires two locks, then triggers a reinsertion failure for one (e.g., make `has_pool_type` invalid via a controlled edit or a crafted object).
  - After run/abort, assert both locks released (`lock.count==0`) and simulation continues.
  - Assert diagnostics recorded (history log contains a cleanup warning event or logger capture if you use caplog).

- `test_rollback_does_not_strand_locks_when_reinsert_fails`
  - Create an op that applies an edit that will fail SHACL/policy -> rollback.
  - Ensure rollback completes and no locks are stranded even if reinsertion fails for an object.
  - (If hard to force reinsertion fail in rollback, monkeypatch the reinsertion helper to raise and assert safe wrapper catches.)

Run `pytest -q tests_sim` and update `codex_state`.

---

## Task 2 — Fix B: interrupt bookkeeping uses held locks only
**Goal**
Prevent lockless KG mutation during interrupt bookkeeping. The interrupt handler must be crash-proof and non-blocking.

**Robust minimal strategy**
- Compute `held_iris` from granted lock requests actually held by the operation.
- Only apply bookkeeping edits whose targets are in `held_iris`.
- If `held_iris` is empty, skip engine apply (event log is authoritative).
- Wrap bookkeeping apply in `except Exception` and record diagnostics (do not re-raise).

**Implementation requirements**
1) In `OperationProcess._handle_interrupt` (or equivalent):
   - Build `held_iris`:
     - iterate `operation.locks`
     - include only granted/held requests
     - map `req.resource` -> object via existing reverse mapping
     - normalize via `identifier_from_iri`
     - guard failures (skip mapping errors)
   - Filter bookkeeping edits to `held_iris`.
   - Call `effect_engine.apply(..., locked_iris=held_iris)` only if non-empty.

2) Ensure handler never crashes:
   - wrap entire bookkeeping block in `try/except Exception as e` -> record diagnostic -> continue.

**Tests**
- `test_interrupt_bookkeeping_skips_when_no_locks_held`
  - Create an op that is interrupted before any locks are acquired (e.g., long scheduled_start wait or waiting on a precedent).
  - Monkeypatch `effect_engine.apply` to raise if called.
  - Assert sim completes and apply was not called.

- `test_interrupt_bookkeeping_only_targets_held_locks`
  - Create an op that acquires a lock on object A and has another participant B (not locked at interrupt time).
  - Interrupt after A lock acquired.
  - Monkeypatch apply to capture `locked_iris` and edit targets.
  - Assert: `locked_iris == {A}` (or subset) and edits only target A.

Run `pytest -q tests_sim` and update `codex_state`.

---

## Task 3 — Pre-lock checklist tests (and small wiring if needed)
**Goal**
Add tests that validate pre-lock checklist. Make small code changes only if required to expose a configuration or guarantee (avoid refactors).

### 3.1 Finite/valid time validation tests
Add tests that ensure invalid values are rejected with clear errors:
- speed_factor: 0, -1, NaN, inf
- temporal_cost: -1, NaN, inf
- scheduled_start_time: NaN, inf (and negative depending on current semantics; assert your chosen behavior)
- spawner intervals (TimerSpawner.interval/start_offset; inspector inspect_interval): NaN/inf and invalid sign

### 3.2 Periodic spawner termination rule tests
- Attach a periodic spawner (TimerSpawner or inspector with inspect_interval>0).
- Calling `Simulation.run(until=None)` must raise ValueError telling user to pass `until`.

### 3.3 Selector atomicity exception-safety tests
- Use an AttributeSelector predicate that raises for some candidates.
- Assert that after the failure:
  - the pool is not drained (object still in store or reinserted)
  - no locks are stranded (`lock.count==0` and no stuck queue)
  - sim does not hang

### 3.4 No double-start test
- Create a simple op that increments a counter on apply.
- Call `spawn_operation(op, start_immediately=True)` then call `sim.run(until=...)`.
- Assert the op ran exactly once (counter == 1).

### 3.5 Overlay strictness policy test
- Ensure strictness is explicitly configurable (if already present, test it; if not, add a minimal parameter on Simulation -> EffectEngine).
- Create an overlay provider that raises.
- In strict mode, validation/apply should fail with a controlled error.
- In tolerant mode, it should log and proceed.

Run `pytest -q tests_sim` and update `codex_state`.
