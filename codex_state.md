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
- Current prompt to run next: P2
- Completed prompts: P0, P1, P1.1, P1.2, P1.3, P1.4, P1.5, P1.6
- Blockers / open questions: None.

---

## Execution log (append newest entries at top)

### 2025-12-31 22:35 (local)
**Prompt executed:** P1.6 — Canonical overlay IRIs
**Summary:**
- What changed: Added canonical IRI helper and identifier normalization; updated current-volume overlay to emit canonical IRIs; normalized SHACL focus IRIs in KGInspectorSpawner; adjusted LiteralSelector to accept canonical IRIs; added overlay test asserting canonical IRIs align with SPPT participants.
- Key decisions: Normalize focus IRIs to identifiers before spawner dispatch; canonicalize current-volume subjects via shared helper.
**Files changed:**
- libsyn_tools/sim/knowledge_graph/base.py
- libsyn_tools/sim/overlay/current_volume_overlay.py
- libsyn_tools/sim/operation/selector.py
- libsyn_tools/sim/spawner.py
- tests_sim/unit/test_overlay_current_volume_provider.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim (passed; warnings present)
**Notes / follow-ups:**
- Next: continue Phase 1 refactors per prompt sequence.
**Next prompt:** P2

### 2025-12-31 22:20 (local)
**Prompt executed:** P1.5 — Spawner time scaling
**Summary:**
- What changed: Applied simulation_speed_factor to TimerSpawner and KGInspectorSpawner timeouts; documented scaling; added a timer spawner test for scaled timing.
- Key decisions: Spawner timeouts use the same dt * speed_factor rule as operations.
**Files changed:**
- libsyn_tools/sim/spawner.py
- libsyn_tools/sim/simulation.py
- tests_sim/sim/test_timer_spawner.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim (passed; warnings present)
**Notes / follow-ups:**
- Next: continue Phase 1 refactors per prompt sequence.
**Next prompt:** P2

### 2025-12-31 22:08 (local)
**Prompt executed:** P1.4 — Filter store availability semantics
**Summary:**
- What changed: Enforced filter stores to contain only present + unlocked objects; prevented literal selection from blocking on non-present objects; added tests for lock-aware store membership and EffectEngine sync behavior.
- Key decisions: Use lock count and is_present checks inside FilterStoreRegistry.put_obj_into_filter_store; literal selector waits on pool only for present objects.
**Files changed:**
- libsyn_tools/sim/operation/selector.py
- tests_sim/unit/test_selector.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim (passed; warnings present)
**Notes / follow-ups:**
- Next: continue Phase 1 refactors per prompt sequence.
**Next prompt:** P2

### 2025-12-31 22:01 (local)
**Prompt executed:** P1.3 — Participant lock dedupe + literal selector pooling
**Summary:**
- What changed: Deduplicated participant lock acquisition in Operation.pre_act; ensured LiteralSelector removes pooled objects from filter stores; added unit test for duplicate participant bindings and lock release integrity.
- Key decisions: Short-circuit literal duplicate bindings before lock request; keep dedupe for selectors by releasing redundant requests and not double-adding locks.
**Files changed:**
- libsyn_tools/sim/operation/operation.py
- libsyn_tools/sim/operation/selector.py
- tests_sim/unit/test_operation_core.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim (passed; warnings present)
**Notes / follow-ups:**
- Next: continue Phase 1 refactors per prompt sequence.
**Next prompt:** P2

### 2025-12-31 20:12 (local)
**Prompt executed:** P1.2 — EffectEngine runtime context confirmation
**Summary:**
- What changed: Confirmed EffectEngine/simulation are already free of module-global runtime caches; no code changes required beyond documentation.
- Key decisions: No refactor needed; runtime access already uses per-environment context.
**Files changed:**
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim (passed; warnings present)
**Notes / follow-ups:**
- Next: continue Phase 1 refactors per prompt sequence.
**Next prompt:** P2

### 2025-12-31 20:00 (local)
**Prompt executed:** P1.1 — Per-environment filter stores
**Summary:**
- What changed: Moved FilterStoreRegistry storage into per-environment RuntimeContext; updated filter store access/removal to require env; added isolation test for filter stores; updated tests that cleared class-global stores.
- Key decisions: Remove class-global registry entirely; rely on context for per-env isolation.
**Files changed:**
- libsyn_tools/sim/operation/runtime.py
- libsyn_tools/sim/operation/selector.py
- libsyn_tools/sim/effect_engine.py
- tests_sim/conftest.py
- tests_sim/unit/test_selector.py
- tests_sim/sim/test_exports_and_determinism.py
- tests_sim/sim/test_spawn_and_serialization.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim (passed; warnings present)
**Notes / follow-ups:**
- Next: continue Phase 1 refactors per prompt sequence.
**Next prompt:** P2

### 2025-12-31 19:54 (local)
**Prompt executed:** P1 — Per-environment runtime context
**Summary:**
- What changed: Introduced a per-environment RuntimeContext for resource/runtime caches and attached it to Simulation environments; refactored runtime usage to rely on the context; added isolation test for sequential simulations.
- Key decisions: Kept deprecated module globals as unused compatibility placeholders; runtime context is auto-created when accessed and explicitly attached in Simulation.
**Files changed:**
- libsyn_tools/sim/operation/runtime.py
- libsyn_tools/sim/effect_engine.py
- libsyn_tools/sim/simulation.py
- libsyn_tools/sim/operation/operation.py
- tests_sim/unit/test_runtime.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim (passed; warnings present)
**Notes / follow-ups:**
- Next: continue Phase 1 refactors (filter stores, locking, time scaling) per prompts.
**Next prompt:** P2

### 2025-12-31 19:32 (local)
**Prompt executed:** P0 — Baseline + tooling discovery (no behavior changes)
**Summary:**
- What changed: Recorded baseline tooling info and test run in codex_state.md.
- Key decisions: No repo lint/test config files found (no pyproject.toml, Makefile, tox.ini, or CI config). Continued to use required test command from agents.md.
**Files changed:**
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim (passed; warnings present)
**Notes / follow-ups:**
- Standard lint/test commands: none discovered in repo; only documented command is `PYTHONPATH=. pytest tests_sim`.
**Next prompt:** P1

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
