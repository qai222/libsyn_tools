# codex_state_202601112041_AB_prelock.md — v202601112041

Baseline expectation: `pytest -q tests_sim` passes before edits.

## Issues
### A) Cleanup/rollback can strand locks if FilterStore reinsertion raises
- In cleanup paths (`post_act`, `_release_all_locks`, rollback restore), a single exception during pool reinsertion can abort the loop early.
- That can leave lock requests held and/or pool objects missing, causing deadlocks.
- Fix must guarantee:
  1) Lock release completes even if reinsertion fails.
  2) Reinsertion failures are observable (diagnostics), not silently swallowed.
  3) If pool_type is invalid, ensure object is removed from all pool stores (no ghost entries).

### B) Interrupt bookkeeping may apply edits without holding locks
- Interrupt handler may call `effect_engine.apply` with `locked_iris` not corresponding to actually-held locks.
- Fix must guarantee:
  1) Only apply bookkeeping edits for objects whose locks are actually held at interrupt time.
  2) No lock acquisition inside interrupt handler.
  3) Interrupt handler never crashes; failures recorded as diagnostics.

## Pre-lock checklist tests to add
- Finite/valid time validation: reject NaN/inf/invalid for speed factor, temporal_cost, scheduled_start_time, spawner intervals.
- Periodic spawner termination rule: require `until` when periodic spawners are attached.
- Selector atomicity exception safety: predicate throws must not drain pools or strand locks.
- No double-start: starting ops pre-run must not cause a second process in `Simulation.run`.
- Overlay strictness policy is explicit and testable (strict vs tolerant).

## Task checklist
- [ ] Task 1 — Fix A (safe cleanup + safe rollback reinsertion) + tests
- [ ] Task 2 — Fix B (interrupt bookkeeping uses held locks only) + tests
- [ ] Task 3 — Add/adjust tests to validate pre-lock checklist + any small wiring changes needed

## Log
(append entries)

### Log entry template
- Date:
- Task:
- Status: DONE
- Summary:
- Files changed:
- Tests added/updated:
- Notes/decisions:
