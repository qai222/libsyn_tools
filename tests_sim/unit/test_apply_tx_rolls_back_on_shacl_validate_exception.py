"""Regression: transaction atomicity must hold even if SHACL validation crashes.

This is distinct from a SHACL *violation* (handled via PolicyBundle / raise_shacl).
Here we simulate the SHACL engine itself throwing an unexpected exception.
"""

import pytest
import simpy
from rdflib import Graph

import libsyn_tools.sim.effect_engine as ee
from libsyn_tools.sim.effect_engine import EffectEngine, EngineMechanicalError
from libsyn_tools.sim.knowledge_graph import Has_capacity, MaterialContainer
from libsyn_tools.sim.operation.runtime import get_runtime_context
from libsyn_tools.sim.operation.unitary_edit import ChangeDataProperty, Create
from twa.data_model.base_ontology import KnowledgeGraph


def test_apply_tx_rolls_back_on_shacl_validate_exception(monkeypatch):
    env = simpy.Environment()
    get_runtime_context(env)  # attach per-env runtime context

    engine = EffectEngine(shapes_graph=Graph())

    c = MaterialContainer()
    c.has_capacity.add(10.0)

    # Ensure object exists in TWA lookup before applying edits.
    KnowledgeGraph.get_object_from_lookup(c.identifier)
    Create(instance_1_iri=c.identifier).apply()
    engine._register_if_new(c, env)

    def boom_validate(*args, **kwargs):
        raise RuntimeError("boom")

    # EffectEngine imports validate at module load time.
    monkeypatch.setattr(ee, "validate", boom_validate)

    edits = [
        ChangeDataProperty(
            instance_1_iri=c.identifier,
            property_iri=Has_capacity.predicate_iri,
            data_value=5.0,
        )
    ]

    with pytest.raises(EngineMechanicalError):
        engine.apply_tx(edits, env, operation_id="op")

    # Capacity must be restored (no partial state).
    assert c.has_capacity == {10.0}
