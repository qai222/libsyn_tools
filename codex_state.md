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

### 2026-01-01 05:47
**Prompt executed:** None (user requested test run + logging only)
**Summary:**
- What changed: No code changes; recorded test execution per user request.
- Why: User asked to run `PYTHONPATH=. pytest tests_sim` and log results.
**Files changed:**
- codex_state.md
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (60 passed; warnings noted in test output).
**Notes / follow-ups:**
- No standard runner (no Makefile/pyproject/tox) detected at repo root.
**Next prompt:** PH2-0

### 2026-01-01 06:43
**Prompt executed:** PH2-5 — Docstrings + enforcement example
**Summary:**
- What changed: Documented policy/apply_tx semantics and added enforcement example; clarified scheduled_start_time wording.
- Why: Provide clear guidance for audit vs enforce behavior and base-time scheduling.
**Files changed:**
- libsyn_tools/sim/policy.py
- libsyn_tools/sim/effect_engine.py
- libsyn_tools/sim/operation/operation.py
- libsyn_tools/sim/simulation.py
- codex_state.md
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (64 passed; warnings noted in test output).
**Notes / follow-ups:**
- Example uses Simulation.effect_engine.policy to enable hard-abort enforcement.
**Next prompt:** PH2-6

### 2026-01-01 06:39
**Prompt executed:** PH2-4 — Operation abort handling
**Summary:**
- What changed: OperationProcess now catches mechanical/contract aborts, logs OPERATION_ABORT, and always calls post_act; added abort-flow test and adjusted mechanical tests for non-raising behavior.
- Why: Ensure aborted operations release locks, complete dependents deterministically, and avoid hangs.
**Files changed:**
- libsyn_tools/sim/effect_engine.py
- libsyn_tools/sim/simulation.py
- tests_sim/sim/test_operation_abort_flow.py
- tests_sim/sim/test_mechanical_abort_no_commit.py
- tests_sim/sim/test_mechanical_lock_coverage_data.py
- tests_sim/sim/test_remove_coverage.py
- tests_sim/sim/test_toctou_and_coverage.py
- tests_sim/unit/test_effect_engine.py
- codex_state.md
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (64 passed; warnings noted in test output).
**Notes / follow-ups:**
- Operation aborts now complete without exceptions in Simulation.run while preserving structured violations.
**Next prompt:** PH2-5

### 2026-01-01 06:34
**Prompt executed:** PH2-3 — Transaction rollback on SHACL abort
**Summary:**
- What changed: Added object-field snapshots and rollback in apply_tx; added rollback test with hard-abort policy.
- Why: Enforce abort disposition by restoring KG + runtime artifacts after SHACL violations.
**Files changed:**
- libsyn_tools/sim/effect_engine.py
- tests_sim/unit/test_effect_engine.py
- codex_state.md
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (63 passed; warnings noted in test output).
**Notes / follow-ups:**
- Rollback restores is_present and mutated property sets for affected objects plus runtime resource/filter-store membership.
**Next prompt:** PH2-4

### 2026-01-01 06:07
**Prompt executed:** PH2-2 — EngineMechanicalError + structured violations
**Summary:**
- What changed: Added EngineMechanicalError with violation payload; mechanical precheck emits structured ENGINE violations; tests updated for new error.
- Why: Ensure mechanical failures produce structured violations for audit and later enforcement.
**Files changed:**
- libsyn_tools/sim/effect_engine.py
- tests_sim/unit/test_effect_engine.py
- tests_sim/sim/test_mechanical_abort_no_commit.py
- tests_sim/sim/test_mechanical_lock_coverage_data.py
- tests_sim/sim/test_remove_coverage.py
- tests_sim/sim/test_toctou_and_coverage.py
- codex_state.md
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (62 passed; warnings noted in test output).
**Notes / follow-ups:**
- EngineMechanicalError now carries violations; apply_tx records and emits them before re-raising.
**Next prompt:** PH2-3

### 2026-01-01 06:04
**Prompt executed:** PH2-1 — TransactionResult + apply_tx scaffold
**Summary:**
- What changed: Added TransactionResult and apply_tx; kept apply as wrapper; added tests for transaction result/violations.
- Why: Establish transaction return structure and ensure violations are reported consistently.
**Files changed:**
- libsyn_tools/sim/effect_engine.py
- tests_sim/unit/test_effect_engine.py
- codex_state.md
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (62 passed; warnings noted in test output).
**Notes / follow-ups:**
- TransactionResult currently reports committed=True; rollback/enforcement still pending.
**Next prompt:** PH2-2

### 2026-01-01 05:50
**Prompt executed:** PH2-0 — Policy scaffolding for SHACL violations
**Summary:**
- What changed: Added policy bundle module; applied policy to SHACL violations; added unit test for policy severity/disposition.
- Why: Enable configurable severity/disposition per shape while preserving audit-only default behavior.
**Files changed:**
- libsyn_tools/sim/policy.py
- libsyn_tools/sim/effect_engine.py
- tests_sim/unit/test_effect_engine.py
- codex_state.md
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (61 passed; warnings noted in test output).
**Notes / follow-ups:**
- No enforcement/rollback yet; policy currently influences violation records only.
**Next prompt:** PH2-1

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
