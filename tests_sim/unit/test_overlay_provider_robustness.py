"""Overlay provider robustness tests.

Overlay providers are user-extensible hooks; a faulty provider should not
crash the simulation core (it should be logged and ignored).
"""

import simpy
from rdflib import Graph

import libsyn_tools.sim.effect_engine as ee
from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph import Has_capacity, MaterialContainer
from libsyn_tools.sim.operation.runtime import get_runtime_context
from libsyn_tools.sim.operation.unitary_edit import ChangeDataProperty, Create
from twa.data_model.base_ontology import KnowledgeGraph


def test_faulty_overlay_provider_is_ignored_and_tx_can_commit(monkeypatch):
    env = simpy.Environment()
    get_runtime_context(env)

    engine = EffectEngine(shapes_graph=Graph())

    c = MaterialContainer()
    c.has_capacity.add(10.0)
    KnowledgeGraph.get_object_from_lookup(c.identifier)
    Create(instance_1_iri=c.identifier).apply()
    engine._register_if_new(c, env)

    # Make SHACL validate trivially succeed so this test is about overlay error handling.
    def ok_validate(*args, **kwargs):
        return True, Graph(), ""

    monkeypatch.setattr(ee, "validate", ok_validate)

    def bad_provider():
        raise RuntimeError("overlay boom")

    engine.register_overlay_provider(bad_provider)

    edits = [
        ChangeDataProperty(
            instance_1_iri=c.identifier,
            property_iri=Has_capacity.predicate_iri,
            data_value=5.0,
        )
    ]
    res = engine.apply_tx(edits, env, operation_id="op", locked_iris=[c.identifier])

    assert res.committed is True
