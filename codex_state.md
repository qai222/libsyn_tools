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

## Completion checklist
- [ ] Transfer unlinks containment on annihilated POMs; regression test added
- [ ] PortionOfMaterial uses explicit validation; zero-volume guarded; tests added
- [ ] put_obj_into_filter_store checks lock queue; concurrency test added
- [ ] pre_act interrupt cannot “silently succeed”; test added
- [ ] KGInspectorSpawner ignores missing focus/shape rather than "None"; test added
- [ ] precedent cycles detected with clear error; test added
- [ ] schedule bridge normalizes module ids; test added
- [ ] `pytest -q tests_sim` passes
