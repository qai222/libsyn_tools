# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_runtime.py ###
from __future__ import annotations

import simpy
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph   import LabObject
from libsyn_tools.sim.operation.runtime import (
    get_runtime_state,
    get_object_for_resource,
)


def test_runtime_state_memoization(env: simpy.Environment):
    o = LabObject()
    KnowledgeGraph.get_object_from_lookup(o.identifier)
    eng = EffectEngine()
    eng._register_if_new(o, env)

    rs1 = get_runtime_state(o, env)
    rs2 = get_runtime_state(o, env)
    assert rs1 is rs2
    assert getattr(o, "_runtime") is rs1


def test_reverse_resource_lookup(env: simpy.Environment):
    o = LabObject()
    KnowledgeGraph.get_object_from_lookup(o.identifier)
    eng = EffectEngine()
    eng._register_if_new(o, env)

    rs = get_runtime_state(o, env)
    back = get_object_for_resource(rs.lock)
    assert back.identifier == o.identifier
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_runtime.py ###
