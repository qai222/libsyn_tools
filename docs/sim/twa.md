# `twa` in `libsyn_tools.sim`

This document tracks the latest state of how `twa` is used in the simulator code.

## What `twa` provides here
- Pydantic-based ontology/OGM primitives:
  - `BaseOntology`
  - `BaseClass`
  - `ObjectProperty`
  - `DatatypeProperty`
  - `KnowledgeGraph`
- In-memory global registries for ontology classes/properties/instances.
- Graph export from in-memory objects via `KnowledgeGraph.graph()`.

## How `sim` layers on top
- `libsyn_tools/sim/knowledge_graph/base.py` defines:
  - `SimOntology` (extends `BaseOntology`)
  - simulator property bases (`SimDataProperty`, `SimObjectProperty`)
  - `Individual` (extends `BaseClass`) with `identifier` aliasing `instance_iri`
- `libsyn_tools/sim/knowledge_graph/ontology.py` defines simulator entities and properties as `twa` classes.
- Instance identity is handled in both identifier and canonical-IRI forms via:
  - `identifier_from_iri(...)`
  - `canonical_iri(...)`

## State and mutation model in `sim`
- `twa` property fields are used as set-like containers (including functional properties with max-cardinality 1).
- The simulator mutates KG state through unitary edit primitives in:
  - `libsyn_tools/sim/operation/unitary_edit.py`
- `EffectEngine` applies edits transactionally and runs optional SHACL checks:
  - `libsyn_tools/sim/effect_engine.py`
- `Simulation` coordinates operation lifecycle and runtime execution:
  - `libsyn_tools/sim/simulation.py`

## Runtime behavior relevant to `twa`
- Object lookup is string-keyed via `KnowledgeGraph.get_object_from_lookup(...)`.
- Objects are auto-registered in lookup tables during model initialization (from `twa` `BaseClass` behavior).
- Selector/runtime layers (`operation/selector.py`, `operation/runtime.py`) build SimPy resource behavior on top of KG objects.

## Test isolation pattern
- `tests_sim/conftest.py` resets KG triples and lookup state between tests to avoid cross-test contamination from global registries.

## Keep this doc current
When changing `twa` integration points in `sim`, update this file in the same commit/PR.
