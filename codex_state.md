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
  - Decision: completion-only. Dependents wait for precedent completion regardless of abort/interrupt.
  - Rationale: preserves existing behavior and avoids introducing new failure modes without
    a clear policy for abort propagation; correctness is enforced via explicit tests.

## Execution log
(append entries as tasks complete)

### Task log format
- Task N: <title>
  - Status: DONE / IN PROGRESS
  - Files changed:
  - Tests added/updated:
  - Notes / decisions:

- Task 1: Validate precedents early (prevent runtime KeyError)
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/simulation.py (validate precedents in Simulation init + spawn_operation)
    - tests_sim/sim/test_missing_precedent.py (regression updated to expect early error)
    - tests_sim/conftest.py (ensure repo root is on sys.path for pytest)
  - Tests added/updated:
    - tests_sim/sim/test_missing_precedent.py::test_missing_precedent_id_fails_early
  - Notes / decisions:
    - required_precedents are validated against known operation IDs during Simulation __init__
      and against operation_registry during spawn_operation, raising a descriptive ValueError.

- Task 2: Fix ValidationAuditSpawner → PolicyEnforcerSpawner precedent crash
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/spawner.py (skip missing precedent ids for remediation spawns)
    - tests_sim/sim/test_spawners_validation_audit.py (regression coverage)
  - Tests added/updated:
    - tests_sim/sim/test_spawners_validation_audit.py::test_validation_audit_policy_enforcer_skips_missing_precedent
  - Notes / decisions:
    - PolicyEnforcerSpawner now only sets precedents when the violation's operation_id
      is registered in the simulation; audit-sourced violations omit precedents.

- Task 3: Decide and test precedent semantics (completion-only)
  - Status: DONE
  - Files changed:
    - tests_sim/sim/test_precedents_semantics.py (new regression tests for abort/interrupt precedents)
  - Tests added/updated:
    - tests_sim/sim/test_precedents_semantics.py::test_precedent_abort_does_not_block_dependents
    - tests_sim/sim/test_precedents_semantics.py::test_precedent_interrupt_does_not_block_dependents
  - Notes / decisions:
    - required_precedents are completion-only: dependent operations proceed once the precedent
      finishes, even if it aborts or is interrupted.

- Task 4: LiteralSelector unknown IRI fails cleanly
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation/selector.py (raise ValueError on unresolved literal IRI)
    - tests_sim/sim/test_literal_selector_unknown_iri.py (regression coverage)
  - Tests added/updated:
    - tests_sim/sim/test_literal_selector_unknown_iri.py::test_literal_selector_unknown_iri_raises_value_error
  - Notes / decisions:
    - LiteralSelector now fails fast with a descriptive ValueError when an IRI
      cannot be resolved in the KnowledgeGraph.

- Task 5: LiteralSelector presence gating (fail fast)
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation/selector.py (reject non-present objects)
    - tests_sim/sim/test_literal_selector_presence.py (regression coverage)
  - Tests added/updated:
    - tests_sim/sim/test_literal_selector_presence.py::test_literal_selector_rejects_non_present_object
  - Notes / decisions:
    - LiteralSelector treats non-present objects as unavailable and raises ValueError
      immediately (no blocking behavior).
    - Tests that rely on LiteralSelector now set participant objects to is_present={True}
      before selection to align with the new presence gate.

- Task 6: EffectEngine presence gating for non-CREATE edits
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/effect_engine.py (reject non-present subjects/objects for non-CREATE edits)
    - tests_sim/unit/test_effect_engine.py (regression coverage)
  - Tests added/updated:
    - tests_sim/unit/test_effect_engine.py::test_mechanical_abort_change_on_non_present
  - Notes / decisions:
    - Non-CREATE edits require present subjects (and object targets where applicable).
    - Edits against freshly created IRIs in the same batch are allowed.

- Task 7: RemoveObjectProperty idempotency
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation/unitary_edit.py (RemoveObjectProperty now no-ops if relation absent)
    - tests_sim/sim/test_remove_object_property_idempotent.py (regression coverage)
  - Tests added/updated:
    - tests_sim/sim/test_remove_object_property_idempotent.py::test_remove_object_property_is_idempotent
  - Notes / decisions:
    - RemoveObjectProperty now mirrors RemoveDataProperty: removing a missing relation is a no-op.

- Task 8: Close lifecycle intervals on abort/interrupt
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/simulation.py (emit operation_end on abort/interrupt)
    - tests_sim/sim/test_overlay_sppt_provider.py (regression coverage)
  - Tests added/updated:
    - tests_sim/sim/test_overlay_sppt_provider.py::test_sppt_overlay_closes_interval_on_abort
  - Notes / decisions:
    - Lifecycle callbacks now emit operation_end for abort/interrupt so SPPT intervals close.

- Task 9: Validate simulation_speed_factor
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/simulation.py (reject non-positive speed factors)
    - tests_sim/sim/test_simulation_speed_factor_validation.py (regression coverage)
  - Tests added/updated:
    - tests_sim/sim/test_simulation_speed_factor_validation.py::test_simulation_speed_factor_zero_rejected
    - tests_sim/sim/test_simulation_speed_factor_validation.py::test_simulation_speed_factor_negative_rejected
  - Notes / decisions:
    - simulation_speed_factor must be > 0; zero/negative values raise ValueError.

- Task 10: build_report handles empty history
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/simulation.py (seed empty event log DataFrame columns)
    - tests_sim/sim/test_build_report_empty.py (regression coverage)
  - Tests added/updated:
    - tests_sim/sim/test_build_report_empty.py::test_build_report_empty_history
  - Notes / decisions:
    - build_report now returns a sensible report when history_log is empty.

- Task 11: Rollback restore wakes waiters (FilterStore.put)
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/effect_engine.py (use store.put during rollback restore)
    - libsyn_tools/sim/operation/selector.py (use store.put when reinserting candidates)
    - tests_sim/sim/test_filterstore_restore_wakes_waiter.py (liveness regression)
  - Tests added/updated:
    - tests_sim/sim/test_filterstore_restore_wakes_waiter.py::test_restore_wakes_waiter_after_rollback
  - Notes / decisions:
    - Rollback/restore paths now use FilterStore.put to wake pending getters
      while avoiding duplicates.

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
