# codex_state.md — Phase 2 (Contracts + Transactions) memory + log

## Scope
This file is Codex’s persistent memory across runs for Phase 2.

Phase 2 goal:
Turn SHACL into an enforceable contract system and make effect application transactional
(stage → apply → validate → commit/rollback).

Do NOT implement Phase 3 (remediation spawners, policy enforcers, KG query selector).

---

## Repo state at Phase 2 start (as of latest snapshot)
- Phase 1 refactor applied.
- Per-env RuntimeContext exists (`env._libsyn_runtime_ctx`).
- Filter stores are per-env and represent available objects.
- OperationProcess uses speed_factor consistently (scheduled_start_time is base-time).
- tests_sim passes using: `PYTHONPATH=. pytest tests_sim` (user reported).

EffectEngine current behavior:
- Applies edits directly (micro-commit).
- Runs SHACL after applying edits (audit).
- Records SHACL violations in SHACLViolationRecord (currently always severity="soft", disposition="committed").
- Mechanical precheck raises RuntimeError on failure (not recorded as violation).
- No rollback / transactional abort semantics.

---

## Phase 2 design decisions (to be implemented)
### Transaction API
- Introduce `TransactionResult` (batch_id, fingerprints, committed bool, violations list).
- Add `EffectEngine.apply_tx(...) -> TransactionResult`.
- Keep `EffectEngine.apply(...)` as backwards-compatible wrapper:
  - call apply_tx and either return None or raise on hard-abort (depending on policy).

### Contract policy
- Introduce `PolicyBundle` (or similar) that maps shape IRIs (and optionally engine errors) to:
  - severity: soft/hard
  - disposition: committed/aborted/ignored
- Default policy: audit-only (commit) to preserve prior behavior unless enabled.

### Rollback strategy
Preferred approach for Phase 2:
- Snapshot and restore object fields for affected instances (is_present + mutated fields).
- Restore runtime artifacts for runtime-tracked objects that were registered/unregistered during apply.

(Do NOT attempt full RDF diffing; operate at object-field level.)

### OperationProcess behavior on abort
- Ensure locks are released even if transaction aborts.
- Record OPERATION_INTERRUPT or OPERATION_ABORT in history_log.
- Decide whether aborted ops mark done_event succeed() (recommended) vs fail.

---

## Status / Next prompt
- Next prompt: PH2-0
- Completed prompts: (none yet)
- Blockers: (none yet)

---

## Execution log (append newest entries at top)

### [UNSTARTED] Phase 2 start
- No Phase 2 work done yet.
- tests_sim reported passing.

### YYYY-MM-DD HH:MM
**Prompt executed:** PH2-{N} — <title>
**Summary:**
- What changed:
- Why:
**Files changed:**
- ...
**Tests run:**
- ...
**Result:**
- ...
**Notes / follow-ups:**
- ...
**Next prompt:** PH2-{N+1}
