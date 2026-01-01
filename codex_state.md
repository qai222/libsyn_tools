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
