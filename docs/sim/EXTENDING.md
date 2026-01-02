## Extending libsyn_tools.sim (Path A)

This guide shows how to add new Operations and hook in overlays/spawners without touching
engine internals. The simulator expects you to express behavior as sequences of `UnitaryEdit`
primitives; the DSL in `EffectsBuilder` makes this ergonomic.

### Write a custom Operation
1. Subclass `libsyn_tools.sim.operation.operation.Operation`.
2. Add `participant_*` fields for resources you need.
3. Implement `get_operation_effects()` and return a `list[UnitaryEdit]`.
4. Prefer the DSL in `libsyn_tools.sim.operation.effects_dsl.EffectsBuilder`.

Example:
```python
from libsyn_tools.sim.operation import Operation
from libsyn_tools.sim.operation.effects_dsl import EffectsBuilder
from libsyn_tools.sim.knowledge_graph import Is_directly_contained_by
from twa.data_model.base_ontology import KnowledgeGraph

class MarkPresent(Operation):
    participant_target: str

    def get_operation_effects(self):
        target = KnowledgeGraph.get_object_from_lookup(self.participant_target)
        builder = EffectsBuilder()
        builder.create(target.identifier)
        return builder.build()
```

### Register overlays or spawners
- Overlays: call `EffectEngine.register_overlay_provider(...)` with a snapshot callback.
  The `Simulation` constructor wires in defaults (see `libsyn_tools/sim/simulation.py`).
- Spawners: attach to a simulation via `.attach(sim)` (see `libsyn_tools/sim/spawner/*`).

### Run a minimal simulation
See `examples/sim_minimal.py` for a small runnable script using presets and the effects DSL.

### Tests
Run the simulator tests with:
```
PYTHONPATH=. pytest tests_sim
```
