# agents.md — Codex operating instructions (libsyn_tools)

## Mission (Phase 1)
Harden `libsyn_tools.sim` into a **re-entrant, deterministic DES+KG kernel** by removing simulator-global mutable state and fixing lock/pool semantics.

Phase 1 goals map to this product statement:
> `libsyn_tools.sim` is a discrete-event lab digital twin where the canonical state is a knowledge graph, correctness is enforced with SHACL contracts, and endogenous spawners create a closed-loop environment that can diagnose and repair protocol failures.

In Phase 1 you are ONLY implementing the “harden the kernel” part:
- eliminate simulator-global mutable state for runtime artifacts (resources, runtime cache, filter stores)
- make locking/pooling correct and deterministic across runs
- align time scaling semantics (operations vs spawners)
- canonicalize KG IRI usage in overlays (minimum: current-volume overlay)

Do NOT start Phase 2/3 features (transactional rollback, policy bundles, remediation libraries, KG query selector).

---

## Repo constraints
- You only have access to THIS repo.
- The repo already has a `tests_sim/` folder and all tests currently pass.
- Your changes MUST keep the full test suite passing.
- Prefer minimal surface-area API breakage. Internal refactors are fine; public imports should remain stable when possible.

---

## Required workflow (every Codex run)
1) Read `codex_state.md` first.
2) Implement exactly the requested task prompt (no extra refactors).
3) Run tests:
   - `pytest tests_sim`.
4) Update `codex_state.md`:
   - Add a new log entry with: what changed, files touched, tests run, results, and what’s next.
   - Update the “Status / Next prompt” section.

If you discover a blocker, record it in `codex_state.md` with enough detail to continue later, then stop.

---

## Coding standards
- Keep changes small and reviewable.
- Add focused tests for new invariants:
  - multi-simulation runs in same Python process do not share runtime state
  - no double-release of locks; no deadlock on duplicate participant bindings
  - filter stores represent “available (present + unlocked) objects only” after the refactor
- Maintain type hints and docstrings.
- Avoid hidden behavior unless documented (e.g., attaching context to `simpy.Environment` is acceptable but must be documented and tested).

---

## Phase 1 target areas (likely files)
You should expect to touch these modules:
- `libsyn_tools/sim/operation/runtime.py`  (resource map + runtime cache)
- `libsyn_tools/sim/operation/selector.py` (FilterStoreRegistry)
- `libsyn_tools/sim/effect_engine.py`      (register/unregister + sync)
- `libsyn_tools/sim/operation/operation.py` (pre_act lock dedup; post_act release)
- `libsyn_tools/sim/simulation.py`         (resource build; exports using runtime cache)
- `libsyn_tools/sim/spawner.py`            (TimerSpawner / KGInspectorSpawner time scaling)
- `libsyn_tools/sim/overlay/current_volume_overlay.py` (IRI canonicalization)

---

## Definition of done for Phase 1
- No simulator-global mutable runtime artifacts shared across simulations:
  - resource maps
  - runtime caches
  - filter stores
- Two simulations created sequentially in the same process do not reuse each other’s SimPy objects.
- Locking is deterministic and safe:
  - no deadlock when the same object is bound to multiple roles (including literal participants)
  - no double-release or leaking locks
- Filter store semantics are coherent:
  - stores contain only present + unlocked objects
- Time scaling is consistent between operations and spawner timeouts.
- All tests in `tests_sim` pass.
