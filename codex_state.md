# codex_state.md — Phase 5 (productization + chemistry pack) memory + log

## Scope
Codex persistent memory for Phase 5.
Phase 5 goal: make the simulator module *useful and user-friendly* without destabilizing the kernel.

---

## Current repo state (Phase 4 complete)
The sim stack now has:
- Re-entrant runtime context per SimPy env (no global runtime artifacts).
- Transactional effect application with rollback on contract aborts.
- SHACL contract policy via PolicyBundle.
- Lifecycle callbacks and spawners (PolicyEnforcerSpawner, ValidationAuditSpawner, TimerSpawner, etc.).
- Default overlays registered in Simulation (SPPT + currentVolume).
- KgQuerySelector (SPARQL over KG + overlays) with atomic select+lock.
- schedule_bridge.compile_schedule_to_simulation for planned ops + SchedulerOutput.

Tests:
- `PYTHONPATH=. pytest tests_sim` is reported passing by user.

Known minor doc drift:
- selector module header says KgQuerySelector not implemented (but it is).

---

## Phase 5 goals (ordered)
1) RunReport: one-call “run + export + summarize” outputs.
2) Standard library operations: Wait/Hold + Mix/Combine + (optional) Measure.
3) Optional Chemistry Pack:
   - chemistry overlay making ingredients queryable (SMILES/mass/etc.)
   - optional SHACL shapes bundle for chemistry constraints
   - optional endogenous drift model (evaporation) as spawner
4) Examples/docs + optional CLI.

---

## Test command
Minimum: `PYTHONPATH=. pytest tests_sim`

---

## Status / Next prompt
- Next prompt: PH5-0
- Completed prompts: (none yet)
- Blockers: (none yet)

---

## Execution log (append newest entries at top)

### [UNSTARTED] Phase 5 start
- Phase 4 code present; tests_sim reported passing.
- No Phase 5 work done yet.

### YYYY-MM-DD HH:MM
**Prompt executed:** PH5-{N} — <title>
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
**Next prompt:** PH5-{N+1}
