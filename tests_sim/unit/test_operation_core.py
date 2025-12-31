# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_operation_core.py ###
from __future__ import annotations

import simpy
from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph   import LabObject
from libsyn_tools.sim.operation.operation import Operation, _OpState
from libsyn_tools.sim.operation.selector import AttributeSelector, FilterStoreRegistry
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.unitary_edit import Create, UnitaryEdit


class _OpCreateTwice(Operation):
    """
    Test helper: create the *same* subject twice → duplicate CREATE IRIs
    must be detected in Operation.pre_act().
    """
    participant_obj: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            Create(instance_1_iri=self.participant_obj),
            Create(instance_1_iri=self.participant_obj),
        ]


class _OpCreateOnce(Operation):
    participant_obj: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [Create(instance_1_iri=self.participant_obj)]


class _OpDualBind(Operation):
    participant_literal: str = Field(...)
    participant_select: AttributeSelector = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_pre_act_lifecycle_and_resources(env: simpy.Environment):
    obj = LabObject()
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    # Register resource like Simulation._build_resources()
    EffectEngine()._register_if_new(obj, env)

    op = _OpCreateOnce(participant_obj=obj.identifier)

    # pre_act advances to PREPARED; locks held; resources populated
    evt = op.pre_act(env)
    env.run(evt)
    assert op.sim_state == _OpState.PREPARED
    assert op.resources == [obj.identifier]
    assert len(op.locks) == 1


def test_pre_act_duplicate_create_detected(env: simpy.Environment):
    obj = LabObject()
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    EffectEngine()._register_if_new(obj, env)

    op = _OpCreateTwice(participant_obj=obj.identifier)

    # pre_act should raise on duplicate CREATE IRIs
    try:
        env.run(op.pre_act(env))
        assert False, "Expected duplicate CREATE detection to raise"
    except RuntimeError as e:
        assert "Duplicate CREATE IRIs detected" in str(e)


def test_pre_act_dedupes_duplicate_participants(env: simpy.Environment):
    obj = LabObject()
    obj.has_pool_type.add("VIAL")
    KnowledgeGraph.get_object_from_lookup(obj.identifier)
    Create(instance_1_iri=obj.identifier).apply()
    EffectEngine()._register_if_new(obj, env)

    def pred(candidate: LabObject) -> bool:
        return candidate.identifier == obj.identifier

    op = _OpDualBind(
        participant_literal=obj.identifier,
        participant_select=AttributeSelector(pool_type="VIAL", predicate=pred),
    )

    env.run(op.pre_act(env))
    assert op.sim_state == _OpState.PREPARED
    assert op.resolved_resources["literal"] == obj.identifier
    assert op.resolved_resources["select"] == obj.identifier
    assert len(op.locks) == 1

    op._mark_running()
    op.post_act(env)
    rs = get_runtime_state(obj, env)
    assert rs.lock.count == 0
    store = FilterStoreRegistry.get_filter_store("VIAL", env)
    assert obj in store.items
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_operation_core.py ###
