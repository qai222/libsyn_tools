# codex_state.md — Phase 3 (closed-loop spawners) memory + log

## Scope
This file is Codex’s persistent memory across runs for Phase 3.

Phase 3 goal:
Implement a closed-loop environment where spawners can diagnose and repair failures using:
- LifecycleCallbacks.on_violation (records emitted by EffectEngine)
- Remediation operations spawned with correct precedents

---

## Current repo state (Phase 2 complete baseline)
- EffectEngine supports transactional apply (`apply_tx`) and rollback when any new SHACL violation
  has disposition="aborted".
- PolicyBundle exists and is used to label SHACLViolationRecord severity/disposition per sourceShape.
- Engine mechanical failures emit SHACLViolationRecord(origin="ENGINE", severity="hard", disposition="aborted").
- Simulation/OperationProcess catches EngineMechanicalError and ContractViolationError, releases locks, logs OPERATION_ABORT.
- Spawners exist:
  - TimerSpawner
  - KGInspectorSpawner (validates and spawns directly, not policy-aware, not record-driven)
  - ProcessInterruptSpawner
- Overlay providers exist:
  - SPPTOverlayProvider (lifecycle→events)
  - CurrentVolumeOverlayProvider
  BUT Simulation currently imports these overlays and does not clearly register them by default.

---

## Phase 3 targets
1) Default overlays registered automatically:
   - SPPTOverlayProvider(callbacks).snapshot
   - CurrentVolumeOverlayProvider().snapshot
2) Violation-driven repair spawner:
   - subscribes to callbacks.on_violation
   - dispatch by (origin, shape_iri) or shape_iri
   - spawns remediation with `precedents=[record.operation_id]`
   - dedup repeated events
3) Minimal remediation library:
   - start with "overfill → DrainExcess(to capacity)" recipe
4) Optional: validator spawner that emits SHACLViolationRecord from validate_now()
   for invariants that arise from endogenous drift (not tied to a specific op apply).

---

## Test command
Minimum: `PYTHONPATH=. pytest tests_sim`

---

## Status / Next prompt
- Next prompt: PH3-5
- Completed prompts: PH3-0, PH3-1, PH3-2, PH3-3, PH3-4
- Blockers: (none yet)

---

## Execution log (append newest entries at top)

### 2025-02-14 00:05
**Prompt executed:** PH3-4 — validation audit spawner + invariant test
**Summary:**
- What changed: Added ValidationAuditSpawner to emit SHACLViolationRecord from validate_now; added a polling audit test; noted preferred spawner composition in KGInspectorSpawner docstring.
- Why: Required to surface non-op SHACL violations and support closed-loop remediation wiring.
**Files changed:**
- libsyn_tools/sim/spawner.py
- tests_sim/sim/test_spawners_validation_audit.py
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (67 passed, 732 warnings)
**Notes / follow-ups:**
- None.
**Next prompt:** PH3-5

### 2025-02-14 00:04
**Prompt executed:** PH3-3 — remediation helper for drain-to-capacity
**Summary:**
- What changed: Added make_drain_to_capacity helper and updated policy-enforcer test to use it.
- Why: Centralize remediation recipe and keep spawner factory concise.
**Files changed:**
- libsyn_tools/sim/remediation.py
- tests_sim/sim/test_spawners_policy_enforcer.py
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (66 passed, 728 warnings)
**Notes / follow-ups:**
- None.
**Next prompt:** PH3-4

### 2025-02-14 00:03
**Prompt executed:** PH3-2 — violation-driven repair spawner + closed-loop test
**Summary:**
- What changed: Added PolicyEnforcerSpawner to react to SHACL violations and spawn remediation ops; added a closed-loop test using DrainExcess on committed overfill violations.
- Why: Required Phase 3 violation-driven repair wiring and regression coverage.
**Files changed:**
- libsyn_tools/sim/spawner.py
- tests_sim/sim/test_spawners_policy_enforcer.py
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (66 passed, 728 warnings)
**Notes / follow-ups:**
- None.
**Next prompt:** PH3-3

### 2025-02-14 00:02
**Prompt executed:** PH3-1 — default overlay registration + SHACL coverage test
**Summary:**
- What changed: Simulation now registers SPPT and currentVolume overlays by default; added a regression test ensuring currentVolume is available during SHACL validation.
- Why: Required to make derived facts available by default and to guard against regressions.
**Files changed:**
- libsyn_tools/sim/simulation.py
- tests_sim/sim/test_default_overlays.py
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (65 passed, 712 warnings)
**Notes / follow-ups:**
- None.
**Next prompt:** PH3-2

### 2025-02-14 00:01
**Prompt executed:** PH3-0 — baseline test run
**Summary:**
- What changed: No code changes.
- Why: User requested test run only.
**Files changed:**
- None
**Tests run:**
- `PYTHONPATH=. pytest tests_sim`
**Result:**
- Passed (64 passed, 714 warnings)
**Notes / follow-ups:**
- None.
**Next prompt:** PH3-1

### [UNSTARTED] Phase 3 start
- Phase 2 code present; tests_sim reported passing.
- No Phase 3 work done yet.

### YYYY-MM-DD HH:MM
**Prompt executed:** PH3-{N} — <title>
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
**Next prompt:** PH3-{N+1}
