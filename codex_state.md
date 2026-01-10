# Codex State — sim correctness hardening (2026-01-09)

## Baseline
- `pytest -q tests_sim` passes before starting these tasks.

## Issues selected as "worth fixing" (simple + robust)
A) Operation/resources tracking
1) Duplicate resources can double-count provenance/utilization:
   - Operation.resources not deduplicated when same resolved resource bound to multiple roles.

B) Transfer preset validation & tolerance correctness
2) TransferMaterialByPortionSize missing source can crash with AttributeError in error message.
3) TransferMaterialByVolume missing/invalid source can crash (no None/type guard).
4) PortionOfMaterial.get_portion_by_volume tolerance allows values that then fail in get_portion.

C) Overlay completeness for subclassed types
5) CurrentVolumeOverlayProvider omits subclasses of MaterialContainer.
6) ChemistryOverlayProvider omits subclasses of PortionOfMaterial.
7) directly_contained_pom_volume undercounts if PortionOfMaterial subclasses exist.
8) get_directly_contained_individuals omits subclass instances when querying base class.

D) Reporting/export correctness
9) Instance history timestamps ignore abort/interrupt (only OPERATION_END).

E) Interrupt classification (optional-but-worth fixing for correctness observability)
10) Selector cancellation returning None can be misclassified as abort rather than interrupt.

## Non-goals (avoid overengineering)
- No global dead confirmation / advanced deadlock prevention.
- No large redesign of graph model.
- Keep semantics stable; only improve correctness/error clarity/coverage.

## Task log
(append entries)

### Template
- Task N: <title>
  - Status: DONE / IN PROGRESS
  - Files changed:
  - Tests added/updated:
  - Notes:

- Task 1: Operation.resources dedupe (order-preserving)
  - Status: DONE
  - Files changed: libsyn_tools/sim/operation/operation.py
  - Tests added/updated: tests_sim/sim/test_duplicate_literal_bindings.py
  - Notes: resources now dedupe before provenance tracking to avoid double-counting.

- Task 2: Transfer preset source validation (PortionSize + Volume)
  - Status: DONE
  - Files changed: libsyn_tools/sim/operation_preset/transfer.py
  - Tests added/updated: tests_sim/test_transfer_validation.py
  - Notes: added explicit None/type guards with clear errors before dereferencing.

- Task 3: PortionOfMaterial.get_portion_by_volume tolerance consistency
  - Status: DONE
  - Files changed: libsyn_tools/sim/knowledge_graph/ontology.py
  - Tests added/updated: tests_sim/unit/test_portion_by_volume_tolerance.py
  - Notes: clamp volume within epsilon to avoid post-validation failures.

- Task 4: Overlay & volume computations include subclass instances
  - Status: DONE
  - Files changed: libsyn_tools/sim/knowledge_graph/ontology.py; libsyn_tools/sim/overlay/current_volume_overlay.py; libsyn_tools/sim/overlay/chemistry_overlay.py
  - Tests added/updated: tests_sim/unit/test_subclass_instances.py
  - Notes: overlays and containment queries now scan subclass instances via all_instances().

- Task 5: Instance history terminal timestamps include ABORT/INTERRUPT
  - Status: DONE
  - Files changed: libsyn_tools/sim/simulation.py
  - Tests added/updated: tests_sim/sim/test_instance_history_terminal_events.py
  - Notes: terminal timestamp index now includes abort/interrupt events.

- Task 6: Preserve interrupt semantics for selector cancellation
  - Status: DONE
  - Files changed: libsyn_tools/sim/operation/operation.py
  - Tests added/updated: tests_sim/sim/test_selector_stability_and_rollback.py
  - Notes: selector cancellation now raises interrupt to avoid misclassified aborts.

## Completion checklist
- [x] Operation.resources deduped (order-preserving); utilization history not double-counted
- [x] Transfer presets: clear, fail-fast errors on missing/invalid source
- [x] get_portion_by_volume tolerance consistent with get_portion
- [x] Overlays include subclass instances (containers + POMs); volume computations include subclasses
- [x] Instance history timestamps include ABORT/INTERRUPT terminal times
- [x] Selector cancel is logged/classified as interrupt (not abort) when cancelled during acquisition
- [x] `pytest -q tests_sim` passes
