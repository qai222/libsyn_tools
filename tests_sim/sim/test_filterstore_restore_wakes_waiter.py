# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_filterstore_restore_wakes_waiter.py ###
from __future__ import annotations

import simpy

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph import LabObject, Has_interrupt_events
from libsyn_tools.sim.operation.selector import AttributeSelector, FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import AddDataProperty, Create


def test_restore_wakes_waiter_after_rollback():
    env = simpy.Environment()
    engine = EffectEngine()

    vial = LabObject()
    vial.has_pool_type.add("VIAL")
    vial.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(vial.identifier)
    Create(instance_1_iri=vial.identifier).apply()
    engine._register_if_new(vial, env)

    store = FilterStoreRegistry.get_filter_store("VIAL", env)
    assert vial in store.items

    edits = [
        AddDataProperty(
            instance_1_iri=vial.identifier,
            property_iri=Has_interrupt_events.predicate_iri,
            data_value="x",
        )
    ]
    snapshots = engine._snapshot_objects(env=env, edits=edits)

    store.items.remove(vial)

    selector = AttributeSelector(pool_type="VIAL", predicate=lambda obj: True)
    result: dict[str, str] = {}

    def _waiter():
        iri, req = yield from selector.resolve(env)
        result["iri"] = iri
        req.resource.release(req)

    env.process(_waiter())

    engine._restore_objects(env=env, snapshots=snapshots)
    env.run()

    assert result["iri"] == vial.identifier
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_filterstore_restore_wakes_waiter.py ###
