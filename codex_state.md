# Codex State — libsyn_tools sim correctness hardening

Date: 2026-01-08
Owner: codex agent (sequential runs; agent does not retain memory between runs)

## Current repo state (as of starting point)
- `tests_sim/` exists
- All existing tests currently pass
- Several correctness issues were identified that are not fully covered by tests yet

## Known correctness issues to address (from issue list)
A) Precedent validation / dependency correctness
1. Missing/invalid precedents can crash at runtime (KeyError in operation_registry lookup).
2. ValidationAuditSpawner emits records with operation_id="VALIDATION_AUDIT" which is not in operation_registry.
   PolicyEnforcerSpawner uses record.operation_id as precedence, causing potential KeyError.

3. Dependents proceed even if precedent aborts/interrupts (no success gating).
   NOTE: This is semantics-sensitive. Decide and document intended behavior.

B) Selector correctness
4. LiteralSelector can select/lock non-present objects (bypasses presence gate).
5. LiteralSelector can crash on unknown IRI (AttributeError/None deref).
6. Interrupt/cancel paths must ensure:
   - pending FilterStore.get is cancelled,
   - pool objects are not “stolen” by stale get events,
   - locks are not stranded,
   - no orphan process crashes env.run.

C) EffectEngine presence semantics
7. Mechanical precheck allows edits on non-present objects (non-CREATE edits should likely require present objects).

D) Unitary edit idempotency
8. RemoveObjectProperty uses .remove and raises if absent (should likely be idempotent or at least consistent).

E) Lifecycle callback / SPPT correctness
9. Aborted/Interrupted operations do not emit operation_end callback, possibly leaving SPPT intervals open/incomplete.

F) Validation of parameters
10. simulation_speed_factor not validated; speed_factor==0 causes division by zero.

G) Reporting correctness
11. build_report can crash when history_log is empty (DataFrame has no columns).

H) Pool restore correctness
12. Rollback restore uses direct store.items mutation (may not wake waiting getters). Prefer store.put.

## Decisions needed
- Semantics for `required_precedents`:
  - Option 1 (completion-only): dependents start after precedent completes (success or failure).
  - Option 2 (success-gated): dependents require precedent success; otherwise they abort/skip.
  - Must choose one and encode in code + tests.
  - If changing semantics, ensure existing tests still pass or update them accordingly.

## Execution log
(append entries as tasks complete)

### Task log format
- Task N: <title>
  - Status: DONE / IN PROGRESS
  - Files changed:
  - Tests added/updated:
  - Notes / decisions:

## Completion checklist
- [ ] No runtime KeyError from missing precedents; validated early with clear message
- [ ] ValidationAuditSpawner + PolicyEnforcerSpawner cannot spawn invalid precedents
- [ ] Precedent semantics decided and tested
- [ ] LiteralSelector: unknown IRI handled, non-present selection behavior defined and tested
- [ ] EffectEngine presence gating defined and tested
- [ ] RemoveObjectProperty idempotent or explicitly mechanical error (documented + tested)
- [ ] Aborted/Interrupted operations close SPPT intervals (callbacks)
- [ ] simulation_speed_factor validated
- [ ] build_report safe on empty history
- [ ] Rollback restore uses store.put (wakes waiters); liveness regression test exists
- [ ] `pytest -q tests_sim` passes
