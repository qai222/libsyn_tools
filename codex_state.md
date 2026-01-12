# codex_state_v0.3_AB.md — v0.3 sim correctness hotfixes (A & B)

Baseline expectation: `pytest -q tests_sim` passes before edits.

## Issues to fix

### A) LiteralSelector can hang indefinitely when FilterStore is out-of-sync
Scenario:
- target object is present/unlocked
- but it is missing from the pool store (`FilterStore.items`)
Current behavior: `LiteralSelector.resolve()` blocks on `store.get(...)` forever.

Fix goal:
- avoid blocking on `store.get` when object is actually available but not in store
- keep normal waiting behavior when object is legitimately unavailable (e.g., not present)

### B) Interrupt bookkeeping may apply edits without holding locks
Scenario:
- operation is interrupted before acquiring locks
- interrupt handler applies bookkeeping edits using `locked_iris` derived from participants/resources
- engine believes locks are held but they are not => lockless KG mutation

Fix goal:
- apply bookkeeping edits only for objects whose locks are actually held at interrupt time
- do not acquire locks during interrupt handling
- if none held, skip KG bookkeeping; event log remains authoritative

## Task checklist
- [x] Task 1 — Fix LiteralSelector out-of-sync hang (with tests)
- [x] Task 2 — Fix interrupt bookkeeping to use held locks only (with tests)

## Log
(append entries)

### Log entry template
- Date:
- Task:
- Status: DONE
- Summary:
- Files changed:
- Tests added/updated:
- Notes:

### Log entry
- Date: 2025-09-27
- Task: Task 1 — Fix LiteralSelector out-of-sync hang (with tests)
- Status: DONE
- Summary: Updated LiteralSelector pool handling to bypass filter store only when present/unlocked and to match by identifier, with immediate get handling to preserve ordering; added regression tests for store desync and presence waiting.
- Files changed: libsyn_tools/sim/operation/selector.py; tests_sim/unit/test_selector.py
- Tests added/updated: tests_sim/unit/test_selector.py (new LiteralSelector regression tests)
- Notes: Verified selector ordering by avoiding yield on already-triggered store.get events.

### Log entry
- Date: 2025-09-27
- Task: Task 2 — Fix interrupt bookkeeping to use held locks only (with tests)
- Status: DONE
- Summary: Avoided interrupt cleanup before interrupt bookkeeping, and limited bookkeeping edits to objects whose locks are still held.
- Files changed: libsyn_tools/sim/simulation.py; tests_sim/sim/test_interrupt_bookkeeping.py
- Tests added/updated: tests_sim/sim/test_interrupt_bookkeeping.py (new interrupt bookkeeping regression tests)
- Notes: Interrupt bookkeeping now uses held locks derived from operation requests and skips apply when none are held.
