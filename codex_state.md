# codex_state.md — libsyn_tools.sim Phase 1 state + execution log

## Scope
This file is Codex’s persistent memory across runs for Phase 1 only.

Phase 1 objective:
**Refactor the sim kernel to be re-entrant and deterministic by removing simulator-global mutable runtime artifacts and fixing lock/pool semantics.**

Do NOT start Phase 2/3 work here.

---

## Current architecture (as-is, before Phase 1 refactor)
Known Phase 1 issues in current codebase:
- Simulator-global mutable structures exist (resource map, runtime cache, filter store registry).
- FilterStoreRegistry uses a class-global `_stores` keyed only by `pool_type`, not by simulation/environment.
- EffectEngine imports and mutates global runtime structures (register/unregister).
- Operation pre_act tries to deduplicate locks after selection/locking; literal selection can deadlock if the same object is targeted twice.
- Overlays use inconsistent IRI strategies (e.g., current volume overlay uses `URIRef(instance_iri)` vs namespaced IRIs).

---

## Phase 1 target invariants
After Phase 1 refactor:
1) Per-simulation runtime isolation:
   - resources, runtime cache, filter stores are per `simpy.Environment` (or per Simulation context), not module globals.
2) Safe locking:
   - never request the same object lock twice in a single operation (duplicate participant bindings must reuse the same lock).
   - never release a lock twice.
3) Pool correctness:
   - filter stores represent *available* objects only (present + unlocked).
4) Deterministic behavior:
   - ordering and time scaling are consistent and documented.
5) Overlay IRI consistency:
   - overlay nodes align with the KG’s canonical IRI scheme.

---

## Status / Next prompt
- Current prompt to run next: (fill in after each run)
- Completed prompts: (fill in)
- Blockers / open questions: (fill in)

---

## Execution log (append newest entries at top)

### [UNSTARTED] Initial state
- Tests: user reports all tests pass.
- No Phase 1 refactor work has been done yet.

(When you run tasks, add entries like below)

### YYYY-MM-DD HH:MM (local)
**Prompt executed:** P{N} - <short name>
**Summary:**
- What changed:
- Key decisions:
**Files changed:**
- path/to/file.py
**Tests run:**
- command + result
**Notes / follow-ups:**
- ...
**Next prompt:** P{N+1}
