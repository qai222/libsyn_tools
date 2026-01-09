# Codex State — sim correctness hardening (2026-01-09)

## Baseline
- Repo: libsyn_tools (zip provided by user)
- tests: `pytest -q tests_sim` currently passes (baseline before changes)

## Semantics decisions (must remain consistent across tasks)
1) Presence semantics:
   - `is_present == {False}` means object is unavailable.
   - Non-CREATE edits on non-present objects are invalid (mechanical error).
2) Precedents semantics:
   - For now, `required_precedents` means "must reach a terminal state" (END/ABORT/INTERRUPT) before dependent may start.
   - Dependents are NOT success-gated unless explicitly implemented later.
3) Spawners:
   - Must not fabricate focus IDs. Only spawn remediation if focus is a valid IRI and maps to an existing KG object.

## Issues to fix (grouped)
A) Cleanup & lifecycle robustness
- Unexpected exceptions during pre_act / effects / callbacks can bypass cleanup and strand locks/pool items.
- pre_act interrupt currently can clear locks without re-queuing pool items.
- Terminal events ordering inconsistent; report misclassifies aborted/interrupted as "in progress".
- Callback exceptions can crash simulation.

B) SHACL / remediation ID validation
- focusNode can be blank node/literal/missing; current code stringifies into invalid IDs.
- KGInspectorSpawner may spawn with "None" focus/shape.
- PolicyEnforcer dedupe suppresses distinct focus violations (key too coarse).

C) Deterministic deadlock avoidance and selector robustness
- Lock ordering uses unstable selector string (repr(predicate) may include memory address).
- Sorting by selector string rather than resolved resources can diverge across operations.
- Selector interrupt around lock acquisition may strand pool item or queued lock.
- FilterStore reinsertion logic bypassed in rollback restore.
- put_obj_into_filter_store ignores queued lock requests.
- remove_obj_from_filter_store uses direct items removal; ensure no stale-get hazards in current logic.

D) Numeric / validation correctness
- PortionOfMaterial uses asserts; can be stripped and can divide by zero.
- make_drain_to_capacity can exceed capacity for tiny capacities.
- MaterialContainer.capacity can raise StopIteration; should raise clear error.
- temporal_cost not validated (negative can crash).

E) Graph / schedule correctness
- Precedent cycle validation can miss cycles when spawn_operation validates only new op.
- schedule bridge module id normalization (canonical vs identifier) can create duplicates.
- DrainExcess destination type not validated.

F) Transfer containment semantics
- Transfer annihilates original POMs without unlinking containment from src (stale links).

## Task execution log
(append entries)

### Log template
- Task N: <title>
  - Status: DONE / IN PROGRESS
  - Files changed:
  - Tests added/updated:
  - Notes/decisions:

- Task 1: Always-cleanup execution + callback isolation + terminal event/report correctness
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/lifecycle.py
    - libsyn_tools/sim/simulation.py
    - libsyn_tools/sim/report.py
    - tests_sim/sim/test_terminal_cleanup_callbacks.py
  - Tests added/updated:
    - tests_sim/sim/test_terminal_cleanup_callbacks.py
  - Notes/decisions:
    - Callback exceptions are recorded as CALLBACK_ERROR history events and do not crash the sim.
    - SHACLValidationError still raises after cleanup to preserve existing raise_shacl behavior.

- Task 2: Fix pre_act abnormal exits (interrupt/exception) to requeue pool items and propagate failure
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation/operation.py
    - tests_sim/sim/test_pre_act_interrupt_cleanup.py
  - Tests added/updated:
    - tests_sim/sim/test_pre_act_interrupt_cleanup.py
  - Notes/decisions:
    - pre_act now releases/requeues locks via a single helper and fails fast on aborted selector resolution.

- Task 3: SHACL remediation robustness: validate focus nodes + fix PolicyEnforcer dedupe key
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/spawner.py
    - tests_sim/sim/test_spawners_inspector.py
    - tests_sim/sim/test_spawners_policy_enforcer.py
    - tests_sim/sim/test_spawners_validation_audit.py
  - Tests added/updated:
    - tests_sim/sim/test_spawners_inspector.py
    - tests_sim/sim/test_spawners_policy_enforcer.py
    - tests_sim/sim/test_spawners_validation_audit.py
  - Notes/decisions:
    - Remediation spawners skip non-IRI or unknown focus nodes, and PolicyEnforcer dedupe now keys on focus.

- Task 4: Stable deadlock avoidance: deterministic lock ordering key + selector interrupt safety
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation/operation.py
    - libsyn_tools/sim/operation/selector.py
    - libsyn_tools/sim/effect_engine.py
    - tests_sim/sim/test_selector_stability_and_rollback.py
  - Tests added/updated:
    - tests_sim/sim/test_selector_stability_and_rollback.py
  - Notes/decisions:
    - Lock ordering now uses stable selector metadata instead of predicate repr.
    - Selector interrupt rollback always reinserts via FilterStoreRegistry.

- Task 5: Pool insertion rules: respect queued locks; consistent remove/put helpers
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation/selector.py
    - tests_sim/sim/test_filterstore_registry.py
  - Tests added/updated:
    - tests_sim/sim/test_filterstore_registry.py
  - Notes/decisions:
    - FilterStore reinsertion now treats queued lock requests as unavailable.

- Task 6: Numeric validation fixes: PortionOfMaterial, drain-to-capacity, container capacity, temporal_cost
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/knowledge_graph/ontology.py
    - libsyn_tools/sim/remediation.py
    - libsyn_tools/sim/operation/operation.py
    - tests_sim/sim/test_numeric_validation.py
  - Tests added/updated:
    - tests_sim/sim/test_numeric_validation.py
  - Notes/decisions:
    - Drain-to-capacity rejects capacities <= EPS; temporal_cost validates at init and pre_act.

- Task 7: Precedent cycle validation on spawn + schedule module id normalization + drain destination validation
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/simulation.py
    - libsyn_tools/sim/adapters/schedule_bridge.py
    - libsyn_tools/sim/operation_preset/drain.py
    - tests_sim/sim/test_spawn_and_serialization.py
    - tests_sim/sim/test_schedule_bridge.py
    - tests_sim/test_drain_excess.py
  - Tests added/updated:
    - tests_sim/sim/test_spawn_and_serialization.py
    - tests_sim/sim/test_schedule_bridge.py
    - tests_sim/test_drain_excess.py
  - Notes/decisions:
    - Spawn validation checks the full precedent graph; schedule bridge normalizes module IDs; DrainExcess rejects non-container destinations.

- Task 8: Schedule bridge normalization + DrainExcess destination validation
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/adapters/schedule_bridge.py
    - libsyn_tools/sim/operation_preset/drain.py
    - tests_sim/sim/test_schedule_bridge.py
    - tests_sim/test_drain_excess.py
  - Tests added/updated:
    - tests_sim/sim/test_schedule_bridge.py
    - tests_sim/test_drain_excess.py
  - Notes/decisions:
    - Module ID normalization checks both identifier and canonical forms; DrainExcess errors on missing/non-container destinations and releases locks on abort.

- Task 9: Transfer semantics unlink containment before annihilation
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation_preset/transfer.py
    - tests_sim/test_transfer_full_portion_size.py
  - Tests added/updated:
    - tests_sim/test_transfer_full_portion_size.py
  - Notes/decisions:
    - Transfer removes containment links (identifier/canonical match) before annihilating original POMs.

- Task 10: PortionOfMaterial ingredient dedupe correctness
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/knowledge_graph/ontology.py
    - tests_sim/unit/test_kg_entities.py
  - Tests added/updated:
    - tests_sim/unit/test_kg_entities.py
  - Notes/decisions:
    - Ingredients are merged by chemical identity (JSON sans mass/identifier) with quantities summed; missing mass/density raises.

## Done checklist
- [x] A: cleanup always runs on any exception; callbacks isolated; terminal events consistent; report correct
- [x] B: spawners validate focus nodes; policy dedupe includes focus; tests
- [x] C: stable lock ordering; selector interrupt safe during get+lock; pool reinsertion respects lock.queue; restore uses registry; tests
- [x] D: explicit validation replaces asserts; drain/capacity edge cases; temporal_cost; tests
- [x] E: cycle detection uses full graph on spawn; schedule id normalization; drain destination checks; tests
- [x] F: transfer unlinks containment before annihilation; tests
- [x] `pytest -q tests_sim` passes
