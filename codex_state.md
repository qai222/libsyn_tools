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
- [ ] Task 1 — Finite/valid time validation everywhere
- [ ] Task 2 — Termination/run semantics + spawner process lifecycle safety
- [ ] Task 3 — Selector/pool robustness (exception-safe atomicity, store healing, fairness, cancellation semantics)
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
- (none yet)

### Log entry template
- Date: YYYY-MM-DD
- Task: <number + name>
- Status: DONE
- Summary:
  - Behavior changes:
  - Files changed:
  - Tests added/updated:
  - Notes/decisions:
