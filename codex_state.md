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
