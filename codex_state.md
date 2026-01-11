# Codex State — libsyn_tools v0.2 sim correctness backlog

This file exists because Codex may not remember prior runs.  
Codex must update this file after each completed task.

## Baseline
- At start: `pytest -q tests_sim` passes.

## Scope
Fix all **unique, valid, worth-fixing** correctness issues identified for sim v0.2, grouped into merged tasks in `codes_task_prompts.md`.

## High-level semantics decisions (do not change without recording)
1. **Presence is availability**
   - Selecting or mutating non-present objects is invalid (except CREATE).
2. **Precedents**
   - `required_precedents` means “must reach a terminal state (END/ABORT/INTERRUPT) before dependents can start”.
   - Not success-gated unless explicitly added later.
3. **Overlay strictness**
   - Provide a strict mode for validation; default may remain tolerant only if explicitly documented.

## Task checklist
- [x] Task 1 — Finite/valid time validation everywhere
- [x] Task 2 — Termination/run semantics + spawner process lifecycle safety
- [x] Task 3 — Selector/pool robustness (exception-safe atomicity, store healing, fairness, cancellation semantics)
- [x] Task 4 — Operation spawn/precedents robustness (no double-start, NEW-only, cycles iterative, deferred-start semantics)
- [x] Task 5 — EffectEngine transactional correctness (no in-place edit mutation, rollback of runtime history, safe snapshotting, IRI normalization, lock rules)
- [x] Task 6 — SHACL/policy/spawner correctness (focus validation, dedupe, throttling, factory exception handling, strictness)
- [x] Task 7 — Presets/remediation validation fixes (transfer/drain/mix, containment unlink, capacity clamping)
- [x] Task 8 — Ontology/material math invariants + serialization stability
- [ ] Task 9 — Overlays/RDF graph correctness (union graph API, stable ingredient IRIs, owl:sameAs bridging, subclass coverage)
- [ ] Task 10 — Reporting/provenance correctness (terminal timestamps, ID collisions, utilization timing, SPPT safety)

## Progress log
Append entries under “Log” as tasks complete.

### Log
- Date: 2025-09-15
  Task: 1 — Finite/valid time validation everywhere
  Status: DONE
  Summary:
    - Behavior changes:
      - Non-finite or negative time inputs now raise ValueError across simulation speed factors, operation timing fields, spawner intervals, and schedule bridge inputs.
      - Schedule bridge now normalizes precedent identifiers via identifier_from_iri.
    - Files changed:
      - libsyn_tools/sim/operation/operation.py
      - libsyn_tools/sim/simulation.py
      - libsyn_tools/sim/spawner.py
      - libsyn_tools/sim/adapters/schedule_bridge.py
      - tests_sim/sim/test_numeric_validation.py
      - tests_sim/sim/test_schedule_bridge.py
      - tests_sim/sim/test_scheduled_start_time_speed.py
      - tests_sim/sim/test_simulation_speed_factor_validation.py
      - tests_sim/sim/test_spawner_time_validation.py
      - codex_state.md
    - Tests added/updated:
      - tests_sim/sim/test_numeric_validation.py
      - tests_sim/sim/test_schedule_bridge.py
      - tests_sim/sim/test_scheduled_start_time_speed.py
      - tests_sim/sim/test_simulation_speed_factor_validation.py
      - tests_sim/sim/test_spawner_time_validation.py
    - Notes/decisions:
      - scheduled_start_time requires finite and >= 0 (negative values rejected).

- Date: 2025-09-15
  Task: 2 — Termination/run semantics + spawner process lifecycle safety
  Status: DONE
  Summary:
    - Behavior changes:
      - Simulation.run now avoids double-starting operations and requires an explicit until when periodic spawners are attached.
      - Periodic spawners track their processes and stop cleanly on detach; spawner callbacks/logics are exception-shielded.
    - Files changed:
      - libsyn_tools/sim/simulation.py
      - libsyn_tools/sim/spawner.py
      - tests_sim/sim/test_spawner_lifecycle.py
      - codex_state.md
    - Tests added/updated:
      - tests_sim/sim/test_spawner_lifecycle.py
    - Notes/decisions:
      - Periodic spawners advertise requires_until and are interrupted on detach.

- Date: 2025-09-15
  Task: 3 — Selector/pool robustness (exception-safe atomicity, store healing, fairness, cancellation semantics)
  Status: DONE
  Summary:
    - Behavior changes:
      - Selector cancellation now yields OPERATION_INTERRUPT via SelectorCancelled handling.
      - Filter store removal scans all pools, and store insert de-dupes by identifier while respecting queued locks.
    - Files changed:
      - libsyn_tools/sim/operation/selector.py
      - libsyn_tools/sim/operation/operation.py
      - libsyn_tools/sim/simulation.py
      - tests_sim/sim/test_filterstore_registry.py
      - tests_sim/sim/test_selector_stability_and_rollback.py
      - tests_sim/unit/test_selector.py
      - codex_state.md
    - Tests added/updated:
      - tests_sim/sim/test_filterstore_registry.py
      - tests_sim/sim/test_selector_stability_and_rollback.py
      - tests_sim/unit/test_selector.py
  - Notes/decisions:
      - Selector cancellation is represented by SelectorCancelled and treated as interrupt-equivalent.

- Date: 2025-09-15
  Task: 4 — Operation spawn/precedents robustness (no double-start, NEW-only, cycles iterative, deferred-start semantics)
  Status: DONE
  Summary:
    - Behavior changes:
      - Operation identifiers must use identifier-form IDs; canonical IRIs are rejected to avoid log collisions.
      - spawn_operation normalizes/de-dupes precedents, validates cycles iteratively, and blocks dependencies on deferred ops.
      - Interrupts can be recorded for operations before their SimPy process is started.
    - Files changed:
      - libsyn_tools/sim/simulation.py
      - tests_sim/sim/test_spawn_precedent_robustness.py
      - tests_sim/sim/test_instance_history_terminal_events.py
      - codex_state.md
    - Tests added/updated:
      - tests_sim/sim/test_spawn_precedent_robustness.py
      - tests_sim/sim/test_instance_history_terminal_events.py
  - Notes/decisions:
      - Deferred-start operations are not allowed as precedents unless already terminal.

- Date: 2025-09-15
  Task: 5 — EffectEngine transactional correctness (no in-place edit mutation, rollback of runtime history, safe snapshotting, IRI normalization, lock rules)
  Status: DONE
  Summary:
    - Behavior changes:
      - EffectEngine now normalizes edits/locked IRIs without mutating inputs, enforces hashable data values, numeric capacity/time values, and functional cardinality.
      - Transaction rollback restores runtime history lengths and snapshotting reports invalid pool types as mechanical errors.
      - CREATE requires exclusivity or explicit locks when a runtime-tracked object appears in pool stores; interrupt bookkeeping uses held locks.
    - Files changed:
      - libsyn_tools/sim/effect_engine.py
      - libsyn_tools/sim/simulation.py
      - tests_sim/unit/test_effect_engine_transactional_correctness.py
      - codex_state.md
    - Tests added/updated:
      - tests_sim/unit/test_effect_engine_transactional_correctness.py
  - Notes/decisions:
      - CREATE on runtime-tracked objects is allowed without locks when the object is not present in any pool store.

- Date: 2025-09-15
  Task: 6 — SHACL/policy/spawner correctness (focus validation, dedupe, throttling, factory exception handling, strictness)
  Status: DONE
  Summary:
    - Behavior changes:
      - SHACL violations now ignore non-URIRef focus nodes and normalize focus IDs against KG objects before remediation.
      - Policy spawners cap remediation attempts per (shape, focus), count factory failures, and warn on unknown shapes.
      - Strict overlay validation defaults to enabled when SHACL shapes are provided (configurable via Simulation).
    - Files changed:
      - libsyn_tools/sim/effect_engine.py
      - libsyn_tools/sim/policy.py
      - libsyn_tools/sim/simulation.py
      - libsyn_tools/sim/spawner.py
      - tests_sim/sim/test_overlay_strictness.py
      - tests_sim/sim/test_spawners_policy_enforcer.py
      - codex_state.md
    - Tests added/updated:
      - tests_sim/sim/test_overlay_strictness.py
      - tests_sim/sim/test_spawners_policy_enforcer.py
    - Notes/decisions:
      - Policy bundle defaults now warn on unknown shapes unless explicitly disabled.

- Date: 2025-09-15
  Task: 7 — Presets/remediation validation fixes (transfer/drain/mix, containment unlink, capacity clamping)
  Status: DONE
  Summary:
    - Behavior changes:
      - Transfer presets now validate source/destination/device containers for existence, type, and presence before computing edits.
      - DrainExcess and MixInContainer provide clearer required-container errors; remediation capacity checks now fail fast on invalid capacity sets.
    - Files changed:
      - libsyn_tools/sim/operation_preset/transfer.py
      - libsyn_tools/sim/operation_preset/drain.py
      - libsyn_tools/sim/operation_preset/mix.py
      - libsyn_tools/sim/remediation.py
      - tests_sim/test_transfer_validation.py
      - tests_sim/sim/test_mix_preset.py
      - codex_state.md
    - Tests added/updated:
      - tests_sim/test_transfer_validation.py
      - tests_sim/sim/test_mix_preset.py
    - Notes/decisions:
      - Container validation errors now report role-specific messages for transfers.

- Date: 2025-09-15
  Task: 8 — Ontology/material math invariants + serialization stability
  Status: DONE
  Summary:
    - Behavior changes:
      - Chemical math now requires finite mass/density and strictly positive density; container capacity must be finite and > 0.
      - Chemistry overlay skips malformed ingredient blobs with a warning and rounds mass values consistently.
    - Files changed:
      - libsyn_tools/chem_schema/chemical.py
      - libsyn_tools/sim/knowledge_graph/ontology.py
      - libsyn_tools/sim/overlay/chemistry_overlay.py
      - tests_sim/sim/test_material_math_invariants.py
      - codex_state.md
    - Tests added/updated:
      - tests_sim/sim/test_material_math_invariants.py
    - Notes/decisions:
      - Capacity values of 0 are treated as invalid (must be finite and > 0).

### Log entry template
- Date: YYYY-MM-DD
- Task: <number + name>
- Status: DONE
- Summary:
  - Behavior changes:
  - Files changed:
  - Tests added/updated:
  - Notes/decisions:
