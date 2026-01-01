# codex_state.md — Phase 4 (KG query selectors + scheduling bridge)

## Scope
Codex memory across runs for Phase 4.

Phase 4 goals:
A) Implement KgQuerySelector (SPARQL-driven selection on KG + overlays).
B) Implement a schedule bridge: compile chem_schema/opt schedule outputs into sim Operations + Simulation.

---

## Current repo state (Phase 3 baseline)
- Simulation registers overlays by default:
  - SPPTOverlayProvider(callbacks).snapshot
  - CurrentVolumeOverlayProvider().snapshot
- Spawners are lifecycle-callback-driven:
  - PolicyEnforcerSpawner reacts to SHACLViolationRecord and spawns remediation ops with precedents
  - ValidationAuditSpawner can emit SHACLViolationRecord from validate_now()

- Phase 2 enforcement exists:
  - EffectEngine.apply_tx does snapshot → apply → SHACL → rollback if any disposition="aborted"
  - ContractViolationError raised by EffectEngine.apply(...) on aborted transactions

- Selector levels 1–3 exist:
  - LiteralSelector, AttributeSelector, HistorySelector
  - Level 4 KgQuerySelector is not implemented yet (doc says so)

---

## Phase 4 design decisions (planned)
### Query graph access
- Add a standard way for selectors to query KG + overlays:
  - EffectEngine.build_query_graph() -> rdflib.Graph (data + overlay)
  - Store a reference to the EffectEngine on simpy.Environment (e.g., env._libsyn_effect_engine)
  - KgQuerySelector uses that union graph for SPARQL.

### KgQuerySelector contract
- Inputs:
  - pool_type (FilterStore pool)
  - SPARQL SELECT query returning candidate IRIs
  - variable name (default "?s")
- Behavior:
  - run query on union graph
  - convert returned URIs to identifiers (strip base URL when needed)
  - atomically select + lock from FilterStore using FIFO
  - if query yields no candidates: raise a clear RuntimeError (deterministic)

### Schedule bridge contract
- New module: libsyn_tools/sim/adapters/schedule_bridge.py
- Compile:
  - Planned ops (chem_schema.Operation or OperationNetwork) + SchedulerOutput
  -> list of sim.Operation processes with:
    - identifier = planned op identifier (stable trace)
    - required_precedents from planned op precedents
    - scheduled_start_time and temporal_cost from SchedulerOutput times (base time)
    - participant_module locked via literal id from SchedulerOutput.assignments

- Ensure module objects exist in KG:
  - create LabObject(identifier=module_id, is_present={True}, has_pool_type={"MODULE"}) if missing

---

## Test command
Minimum: `PYTHONPATH=. pytest tests_sim`

---

## Status / Next prompt
- Next prompt: PH4-0
- Completed prompts: (none yet)
- Blockers: (none yet)

---

## Execution log (append newest entries at top)

### 2026-01-01 10:26
**Prompt executed:** PH4-5 — Doc snippet for selectors and schedule bridge
**Summary:**
- What changed: Added module docstring examples for KgQuerySelector and compile_schedule_to_simulation usage.
- Why: Provide minimal runnable usage snippets for Phase 4 features.
**Files changed:**
- libsyn_tools/sim/adapters/schedule_bridge.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim
**Result:**
- Passed (73 passed, 734 warnings).
**Notes / follow-ups:**
- None.
**Next prompt:** PH4-6

### 2026-01-01 10:20
**Prompt executed:** PH4-4 — Schedule bridge translator
**Summary:**
- What changed: Added op_translator support in schedule bridge and a translator test using TransferMaterialByVolume to mutate the KG.
- Why: Allow custom sim.Operation mapping from planned operations while preserving schedule semantics.
**Files changed:**
- libsyn_tools/sim/adapters/schedule_bridge.py
- tests_sim/sim/test_schedule_bridge.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim
**Result:**
- Passed (73 passed, 734 warnings).
**Notes / follow-ups:**
- None.
**Next prompt:** PH4-5

### 2026-01-01 10:15
**Prompt executed:** PH4-3 — Schedule bridge constraints
**Summary:**
- What changed: Added schedule bridge tests for precedents and module contention to ensure start times obey dependencies and locks.
- Why: Validate schedule bridge semantics under inconsistent schedules and shared-module contention.
**Files changed:**
- tests_sim/sim/test_schedule_bridge.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim
**Result:**
- Passed (72 passed, 734 warnings).
**Notes / follow-ups:**
- None.
**Next prompt:** PH4-4

### 2026-01-01 09:55
**Prompt executed:** PH4-2 — Schedule bridge adapter
**Summary:**
- What changed: Added schedule bridge adapter to compile scheduler output into sim operations and new tests for scheduled timing.
- Why: Enable planned schedule execution in the DES kernel per Phase 4.
**Files changed:**
- libsyn_tools/sim/adapters/__init__.py
- libsyn_tools/sim/adapters/schedule_bridge.py
- tests_sim/sim/test_schedule_bridge.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim
**Result:**
- Passed (70 passed, 734 warnings).
**Notes / follow-ups:**
- None.
**Next prompt:** PH4-3

### 2026-01-01 09:34
**Prompt executed:** PH4-1 — KgQuerySelector + query-based test
**Summary:**
- What changed: Implemented KgQuerySelector with KG/overlay SPARQL selection, normalized identifiers, and added a simulation test for currentVolume-based selection.
- Why: Enable deterministic Level 4 selection over KG + overlays per Phase 4 requirements.
**Files changed:**
- libsyn_tools/sim/operation/selector.py
- libsyn_tools/sim/env_utils.py
- tests_sim/sim/test_kg_query_selector.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim
**Result:**
- Passed (69 passed, 732 warnings).
**Notes / follow-ups:**
- None.
**Next prompt:** PH4-2

### 2026-01-01 08:41
**Prompt executed:** PH4-0 — Build query graph + env access
**Summary:**
- What changed: Added EffectEngine.build_query_graph, attached engine to Simulation env, added env helper, and new test for currentVolume overlay in query graph.
- Why: Enable KG query selection and safe access to the effect engine on the SimPy environment.
**Files changed:**
- libsyn_tools/sim/effect_engine.py
- libsyn_tools/sim/simulation.py
- libsyn_tools/sim/env_utils.py
- tests_sim/sim/test_effect_engine_query_graph.py
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim
**Result:**
- Passed (68 passed, 732 warnings).
**Notes / follow-ups:**
- None.
**Next prompt:** PH4-1

### 2025-09-24 06:25
**Prompt executed:** PH4-0 — Test run only
**Summary:**
- What changed: Ran the required test suite; no code changes.
- Why: User requested running tests and recording results.
**Files changed:**
- codex_state.md
**Tests run:**
- PYTHONPATH=. pytest tests_sim
**Result:**
- Passed (67 passed, 732 warnings).
**Notes / follow-ups:**
- None.
**Next prompt:** PH4-0

### [UNSTARTED] Phase 4 start
- Phase 3 baseline present.
- No Phase 4 work done yet.

### YYYY-MM-DD HH:MM
**Prompt executed:** PH4-{N} — <title>
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
**Next prompt:** PH4-{N+1}
