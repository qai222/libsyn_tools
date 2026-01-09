# Codex State — sim correctness fixes

Date: 2026-01-08

## Current baseline
- `tests_sim/` exists
- All tests currently pass before applying these tasks.

## Issues selected as REAL + worth fixing (holistic & minimal)
A) Transfer preset: stale containment link on annihilated POMs
- Transfer annihilates original POM(s) without unlinking from src container.
- Annihilate flips presence only; relationship persists.

B) PortionOfMaterial splitting: unsafe `assert` validation and possible division-by-zero
- `get_portion_by_volume` divides by `self.volume` and relies on asserts.
- Asserts can be stripped under `python -O`.

C) Pool availability check ignores queued lock requests
- `FilterStoreRegistry.put_obj_into_filter_store` checks `rs.lock.count` only.
- Should also consider `rs.lock.queue` to preserve atomic selection intent.

D) pre_act interrupt can silently “succeed”
- `Operation._pre_act_implementation` swallows Interrupt and returns.
- `_run_core` marks operation RUNNING after pre_act returns, potentially without resolved participants/locks.

E) KGInspectorSpawner can fabricate bogus "None" IDs
- It stringifies missing report fields, producing "None" focus/shape IDs.

F) Precedent cycles not detected
- Existence validated, but cycles can deadlock simulation (AllOf waits never resolve).

G) Schedule bridge module IDs not normalized → duplicate module objects
- `_ensure_modules_present` uses raw IDs; canonical vs identifier mismatch can cause duplicates.

(Also watch: remediation precedents from audit records must be valid operation IDs; ensure no missing precedents are introduced by spawners.)

## Decisions / semantics
- Presence is treated as availability: selecting or mutating non-present objects is invalid unless explicitly CREATE.
- pre_act interrupt should propagate cancellation upward (do not allow op to proceed RUNNING).
- Precedent graph must be a DAG (no cycles).

## Task log
(append as tasks complete)

### Template
- Task N: <title>
  - Status: DONE / IN PROGRESS
  - Files changed:
  - Tests added/updated:
  - Notes:

- Task 1: Transfer preset unlink containment before annihilation
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation_preset/transfer.py
    - tests_sim/test_transfer_containment_cleanup.py
    - tests_sim/unit/test_preset_transfer.py
    - codex_state.md
  - Tests added/updated:
    - tests_sim/test_transfer_containment_cleanup.py::test_transfer_removes_containment_link_from_annihilated_pom
    - tests_sim/unit/test_preset_transfer.py::test_transfer_by_volume_edits_shape
  - Notes:
    - Remove stale Is_directly_contained_by links before annihilating original POMs.

- Task 2: PortionOfMaterial validation guards
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/knowledge_graph/ontology.py
    - tests_sim/unit/test_kg_entities.py
    - codex_state.md
  - Tests added/updated:
    - tests_sim/unit/test_kg_entities.py::test_pom_get_portion_by_volume_zero_volume_raises
    - tests_sim/unit/test_kg_entities.py::test_pom_get_portion_rejects_invalid_size
  - Notes:
    - Replace asserts with explicit validation and guard against zero-volume splits.

- Task 3: Pool insertion respects queued locks
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation/selector.py
    - tests_sim/unit/test_selector.py
    - codex_state.md
  - Tests added/updated:
    - tests_sim/unit/test_selector.py::test_filter_store_rejects_queueing_locks
  - Notes:
    - Avoid reinserting pooled objects while lock requests are queued.

- Task 4: pre_act interrupt propagates cancellation
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/operation/operation.py
    - libsyn_tools/sim/simulation.py
    - tests_sim/unit/test_pre_act_interrupt.py
    - codex_state.md
  - Tests added/updated:
    - tests_sim/unit/test_pre_act_interrupt.py::test_pre_act_interrupt_propagates_and_does_not_run
  - Notes:
    - Re-raise pre_act interrupts after cleanup to prevent RUNNING state transitions.

- Task 5: KGInspectorSpawner skips missing focus/shape
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/spawner.py
    - tests_sim/sim/test_spawners_inspector.py
    - codex_state.md
  - Tests added/updated:
    - tests_sim/sim/test_spawners_inspector.py::test_inspector_skips_missing_focus_or_shape
  - Notes:
    - Ignore SHACL report entries missing focusNode or sourceShape.

- Task 6: Precedent cycle detection
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/simulation.py
    - tests_sim/sim/test_precedent_cycles.py
    - codex_state.md
  - Tests added/updated:
    - tests_sim/sim/test_precedent_cycles.py::test_precedent_cycle_fails_fast
  - Notes:
    - Detect precedent cycles during validation and raise with the cycle path.

- Task 7: Schedule bridge module ID normalization
  - Status: DONE
  - Files changed:
    - libsyn_tools/sim/adapters/schedule_bridge.py
    - tests_sim/sim/test_schedule_bridge.py
    - codex_state.md
  - Tests added/updated:
    - tests_sim/sim/test_schedule_bridge.py::test_schedule_bridge_reuses_module_identifier_form
  - Notes:
    - Normalize schedule module IDs to reuse existing module objects.

## Completion checklist
- [x] Transfer unlinks containment on annihilated POMs; regression test added
- [x] PortionOfMaterial uses explicit validation; zero-volume guarded; tests added
- [x] put_obj_into_filter_store checks lock queue; concurrency test added
- [x] pre_act interrupt cannot “silently succeed”; test added
- [x] KGInspectorSpawner ignores missing focus/shape rather than "None"; test added
- [x] precedent cycles detected with clear error; test added
- [x] schedule bridge normalizes module ids; test added
- [x] `pytest -q tests_sim` passes
