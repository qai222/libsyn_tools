## Extending libsyn_tools.sim (Path A)

This guide shows how to add new Operations and hook in overlays/spawners without touching
engine internals. The simulator executes `UnitaryEdit` transactions via `Simulation`;
`EffectsBuilder` is the easiest way to compose those edits.

### Write a custom Operation
1. Subclass `libsyn_tools.sim.operation.Operation`.
2. Add `participant_*` fields for resources you need.
3. Implement `get_operation_effects()` and return a `list[UnitaryEdit]`.
4. Prefer `libsyn_tools.sim.operation.effects_dsl.EffectsBuilder`.
5. Run the operation through `Simulation`; direct `Operation.execute()` is not supported.

Example:
```python
from libsyn_tools.sim.operation import Operation, UnitaryEdit
from libsyn_tools.sim.operation.effects_dsl import EffectsBuilder
from libsyn_tools.sim.knowledge_graph import Has_interrupt_events
from twa.data_model.base_ontology import KnowledgeGraph

class TagInterruptEvent(Operation):
    participant_target: str
    reason: str = "manual-tag"

    def get_operation_effects(self) -> list[UnitaryEdit]:
        target = KnowledgeGraph.get_object_from_lookup(self.participant_target)
        if target is None:
            raise RuntimeError(f"Unknown target {self.participant_target!r}")
        return (
            EffectsBuilder()
            .add_data(
                target.identifier,
                Has_interrupt_events.predicate_iri,
                self.reason,
            )
            .build()
        )
```

### Register overlays or spawners
- Overlays: call `sim.effect_engine.register_overlay_provider(...)` with a snapshot callback.
  The `Simulation` constructor already wires in default overlays (see `libsyn_tools/sim/simulation.py`).
- Spawners: subclass `Spawner` in `libsyn_tools/sim/spawner.py`, then attach with `.attach(sim)`.
  Note: periodic spawners require calling `sim.run(until=...)`.

### Run a minimal simulation
See `examples/sim_minimal.py` for a small runnable script using presets and direct `UnitaryEdit` primitives.

### Tests
Run the simulator tests with:
```
PYTHONPATH=. pytest -q tests_sim
```
