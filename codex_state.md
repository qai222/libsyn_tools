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
- [ ] Task 4 — Operation spawn/precedents robustness (no double-start, NEW-only, cycles iterative, deferred-start semantics)
- [ ] Task 5 — EffectEngine transactional correctness (no in-place edit mutation, rollback of runtime history, safe snapshotting, IRI normalization, lock rules)
- [ ] Task 6 — SHACL/policy/spawner correctness (focus validation, dedupe, throttling, factory exception handling, strictness)
- [ ] Task 7 — Presets/remediation validation fixes (transfer/drain/mix, containment unlink, capacity clamping)
- [ ] Task 8 — Ontology/material math invariants + serialization stability
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

### Log entry template
- Date: YYYY-MM-DD
- Task: <number + name>
- Status: DONE
- Summary:
  - Behavior changes:
  - Files changed:
  - Tests added/updated:
  - Notes/decisions:
