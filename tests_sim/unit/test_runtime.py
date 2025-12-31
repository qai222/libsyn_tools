# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_runtime.py ###
from __future__ import annotations

import simpy
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph   import LabObject
from libsyn_tools.sim.operation.runtime import (
    get_runtime_context,
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


def test_runtime_context_isolated_between_simulations():
    obj = LabObject()
    KnowledgeGraph.get_object_from_lookup(obj.identifier)

    sim_one = Simulation([])
    ctx_one = get_runtime_context(sim_one.env, create=False)
    rs_one = get_runtime_state(obj, sim_one.env)
    res_one = ctx_one.resource_map[obj.instance_iri]

    sim_two = Simulation([])
    ctx_two = get_runtime_context(sim_two.env, create=False)
    rs_two = get_runtime_state(obj, sim_two.env)
    res_two = ctx_two.resource_map[obj.instance_iri]

    assert ctx_one is not ctx_two
    assert rs_one is not rs_two
    assert res_one is not res_two
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_runtime.py ###
