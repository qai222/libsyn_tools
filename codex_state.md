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
- [ ] Task 1 — Fix LiteralSelector out-of-sync hang (with tests)
- [ ] Task 2 — Fix interrupt bookkeeping to use held locks only (with tests)

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
