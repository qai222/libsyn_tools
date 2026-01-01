# agents.md — Codex operating instructions (Phase 2: Contracts + Transactions)

## Mission (Phase 2)
Implement **enforceable SHACL contracts** and **transactional effect application** in `libsyn_tools.sim`:

- The canonical state is the KG (TWA BaseClass + KnowledgeGraph).
- Correctness is enforced using SHACL contracts (configurable hard/soft).
- Effect application becomes transactional: stage → apply → validate → (commit | rollback).

Phase 2 focuses on the EffectEngine and its interaction with OperationProcess.

Do NOT start Phase 3 (closed-loop remediation libraries, policy-driven spawners, KG query selectors).

---

## Current baseline assumptions
- Phase 1 refactor is already merged.
- Per-env RuntimeContext exists (resource_map/runtime_cache/filter_stores).
- Duplicate literal bindings and time scaling semantics are fixed.
- tests_sim passes under: `PYTHONPATH=. pytest tests_sim`.

---

## Required workflow (EVERY run)
1) Read `codex_state.md` first.
2) Execute exactly ONE task prompt (PH2-0, PH2-1, ...).
3) Run tests:
   - Minimum: `PYTHONPATH=. pytest tests_sim`
   - If there are other tests, also run: `PYTHONPATH=. pytest -q` (optional but recommended)
4) Update `codex_state.md`:
   - what changed and why
   - files touched
   - tests run + results
   - next prompt to execute

If blocked, log the blocker with enough detail to continue next run, then stop.

---

## Implementation rules / guardrails
- Keep diffs small and reversible.
- Preserve backwards compatibility where reasonable:
  - existing `EffectEngine.apply(...)` can remain, but introduce `apply_tx(...)` returning a result.
  - default behavior should remain "audit" (commit even if SHACL violations) unless enforcement is explicitly enabled.
- When enforcement is enabled and a transaction aborts, ensure:
  - KG state is rolled back to pre-transaction
  - runtime artifacts remain consistent (resource_map/filter stores not corrupted)
  - Operation locks are released (avoid hangs).
- Always emit structured violation records for BOTH:
  - engine/mechanical failures (origin="ENGINE")
  - SHACL failures (origin="SHACL")

---

## Likely files to touch
- `libsyn_tools/sim/effect_engine.py`        (transaction API + enforcement)
- `libsyn_tools/sim/effect_shacl.py`         (record + helpers; may add policy fields)
- `libsyn_tools/sim/simulation.py`           (catch/record aborted ops; ensure post_act executes)
- `libsyn_tools/sim/operation/unitary_edit.py` (optional: inverse helpers, if you choose that route)
- `tests_sim/...`                            (new tests for enforcement + rollback)

---

## Phase 2 exit criteria
1) A new transaction API exists (apply_tx / TransactionResult) with batch_id + fingerprints.
2) SHACL contracts support hard/soft modes and configurable disposition.
3) Enforced mode can abort a violating transaction AND roll back KG state.
4) Engine mechanical failures generate structured violation records.
5) tests_sim remains green.
