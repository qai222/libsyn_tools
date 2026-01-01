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
- Next prompt: PH3-0
- Completed prompts: (none yet)
- Blockers: (none yet)

---

## Execution log (append newest entries at top)

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
