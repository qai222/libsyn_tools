# Codex State — sim correctness backlog (2026-01-09)

Baseline: tests_sim currently pass before work begins.

## Goals
Fix the consolidated correctness issues in a non-overengineered way:
- stop resource leaks / deadlocks from exceptions
- validate inputs early (finite times, missing objects, malformed sets)
- make spawners/overlays robust (no fabricated IDs, no silent masking unless configured)
- improve reporting/provenance correctness

## Decisions (recorded)
- Presence is availability: non-present objects should not be selected or mutated (except CREATE).
- Precedents mean "must reach terminal state" (END/ABORT/INTERRUPT) unless explicitly changed.
- Overlay strictness: add a switch; default remains tolerant unless strict requested.
- Remediation throttling: implement a simple, deterministic cap to prevent infinite storms (no complex state machines).

## Work plan (merged tasks)
Task 1: “Always-cleanup” lifecycle + interrupt bookkeeping + Simulation.run finalizers + terminal event ordering + report safety
Task 2: Selector & pool robustness (cleanup on any exception, lock.queue checks, LiteralSelector heal, HistorySelector predicate signature, cancellation classification)
Task 3: Operation pre_act/resource tracking (stable ordering key, resource dedupe, propagate selector-cancel as interrupt, robust lock release)
Task 4: EffectEngine correctness (normalize canonical IRIs, snapshot failures -> EngineMechanicalError, rollback reinsertion via registry, optional strict overlays)
Task 5: SHACL/spawners correctness (focus node validation, spawner exception shielding, TimerSpawner termination guard, detach spawners, PolicyEnforcer dedupe+throttle)
Task 6: Presets/remediation validation (transfer src guards, drain dst guards, drain-to-capacity clamp + clear errors, transfer containment unlink, mix-in presence checks)
Task 7: Ontology/material math robustness (portion_by_volume clamp, chemical None handling, reduce JSON drift with rounding, functional-set validation)
Task 8: Overlay/subclass coverage + containment helpers (all_instances traversal)
Task 9: Provenance/report exports (instance history timestamps for terminal events, stable IDs, SPPT IRI escaping + safe removal)

## Task log
(append as tasks complete)

### Task 1: “Always-cleanup” lifecycle + interrupt bookkeeping + Simulation.run finalizers + terminal event ordering + report safety
- Status: DONE
- Files changed:
  - libsyn_tools/sim/simulation.py
  - tests_sim/test_interrupt_cleanup.py
- Tests added/updated:
  - tests_sim/test_interrupt_cleanup.py::test_interrupt_bookkeeping_exception_does_not_crash
- Notes:
  - Cleanup now runs on any exception in the operation core; interrupt bookkeeping errors are swallowed.
  - Simulation.run always releases progress callbacks via try/finally.

### Task 2: Selector & pool robustness (cleanup on any exception, lock.queue checks, LiteralSelector heal, HistorySelector predicate signature, cancellation classification)
- Status: DONE
- Files changed:
  - libsyn_tools/sim/operation/selector.py
  - tests_sim/unit/test_selector.py
- Tests added/updated:
  - tests_sim/unit/test_selector.py::test_literal_selector_heals_missing_store_entry
  - tests_sim/unit/test_selector.py::test_history_selector_accepts_env_predicate
  - tests_sim/unit/test_selector.py::test_predicate_exception_returns_object_to_pool
- Notes:
  - Selector cleanup now handles any exception, preventing lost pool items or stranded locks.
  - HistorySelector now supports predicates that accept (obj, env).

### Task 3: Operation pre_act/resource tracking (stable ordering key, resource dedupe, propagate selector-cancel as interrupt, robust lock release)
- Status: DONE
- Files changed:
  - libsyn_tools/sim/operation/operation.py
  - tests_sim/sim/test_selector_stability_and_rollback.py
  - tests_sim/sim/test_operation_runtime_tracking.py
  - tests_sim/unit/test_operation_core.py
- Tests added/updated:
  - tests_sim/sim/test_selector_stability_and_rollback.py::test_selector_ordering_with_lambdas_is_stable_across_ops
  - tests_sim/sim/test_operation_runtime_tracking.py::test_recent_operations_dedupes_duplicate_bindings
  - tests_sim/unit/test_operation_core.py::test_release_all_locks_handles_missing_reverse_mapping
- Notes:
  - Pre-act ordering key now uses stable pool/role/class ordering and literal identifiers.
  - Lock cleanup no longer fails when reverse mappings are missing.

### Task 4: EffectEngine correctness (normalize canonical IRIs, snapshot failures -> EngineMechanicalError, rollback reinsertion via registry, optional strict overlays)
- Status: DONE
- Files changed:
  - libsyn_tools/sim/effect_engine.py
  - tests_sim/unit/test_effect_engine.py
- Tests added/updated:
  - tests_sim/unit/test_effect_engine.py::test_apply_accepts_canonical_iri
  - tests_sim/unit/test_effect_engine.py::test_snapshot_pool_type_error_becomes_mechanical
  - tests_sim/unit/test_effect_engine.py::test_overlay_provider_strictness
- Notes:
  - EffectEngine now resolves canonical identifiers for lookups and optionally fails fast on overlay provider errors.
  - Snapshot pool type validation errors are converted to EngineMechanicalError.

### Task 5: SHACL/spawners correctness (focus node validation, spawner exception shielding, TimerSpawner termination guard, detach spawners, PolicyEnforcer dedupe+throttle)
- Status: DONE
- Files changed:
  - libsyn_tools/sim/spawner.py
  - libsyn_tools/sim/simulation.py
  - tests_sim/sim/test_spawners_detach.py
  - tests_sim/sim/test_spawners_policy_enforcer.py
  - tests_sim/sim/test_spawners_timer_more.py
  - tests_sim/sim/test_spawners_validation_audit.py
  - tests_sim/sim/test_timer_spawner.py
- Tests added/updated:
  - tests_sim/sim/test_timer_spawner.py::test_timer_requires_until_for_unbounded_run
  - tests_sim/sim/test_spawners_timer_more.py::test_timer_factory_error_does_not_crash
  - tests_sim/sim/test_spawners_detach.py::test_spawner_detached_between_runs
  - tests_sim/sim/test_spawners_policy_enforcer.py::test_policy_enforcer_throttle_limits_repeats
  - tests_sim/sim/test_spawners_validation_audit.py::test_validation_audit_blank_focus_does_not_spawn_remediation
- Notes:
  - TimerSpawner errors are contained, and simulations require explicit until when timers are attached.
  - PolicyEnforcer now dedupes by (shape, focus) and enforces per-focus remediation caps; Simulation.run detaches spawners.

### Template
- Task N: <title>
  - Status: DONE / IN PROGRESS
  - Files changed:
  - Tests added/updated:
  - Notes:

## Completion checklist
- [ ] No leaked locks/pool items on any exception paths
- [ ] Selector/pool logic robust to predicate errors and lock waits
- [ ] Canonical IRIs accepted in EffectEngine edits
- [ ] Spawners never fabricate invalid focus IDs; remediation storms capped
- [ ] TimerSpawner cannot cause non-termination unless explicitly allowed
- [ ] Presets fail fast with clear errors; transfer unlink containment
- [ ] Material math safe for edge cases; functional sets validated
- [ ] Overlays include subclass instances
- [ ] Reports/provenance correct for abort/interrupt; SPPT IRIs safe
- [ ] pytest -q tests_sim passes
