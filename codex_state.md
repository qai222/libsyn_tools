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

## Done checklist
- [ ] A: cleanup always runs on any exception; callbacks isolated; terminal events consistent; report correct
- [ ] B: spawners validate focus nodes; policy dedupe includes focus; tests
- [ ] C: stable lock ordering; selector interrupt safe during get+lock; pool reinsertion respects lock.queue; restore uses registry; tests
- [ ] D: explicit validation replaces asserts; drain/capacity edge cases; temporal_cost; tests
- [ ] E: cycle detection uses full graph on spawn; schedule id normalization; drain destination checks; tests
- [ ] F: transfer unlinks containment before annihilation; tests
- [ ] `pytest -q tests_sim` passes
