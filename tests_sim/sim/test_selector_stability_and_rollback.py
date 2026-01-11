from __future__ import annotations

import simpy

from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.runtime import get_runtime_state
from libsyn_tools.sim.operation.selector import (
    AttributeSelector,
    FilterStoreRegistry,
    Selector,
    SelectorCancelled,
)
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, Annihilate, Create
from libsyn_tools.sim.policy import PolicyBundle, PolicyRule


class TwoSelectorOp(Operation):
    participant_alpha: AttributeSelector
    participant_zeta: AttributeSelector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


class SelectorWaitOp(Operation):
    participant_target: AttributeSelector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


class CancelSelector(Selector):
    def resolve(self, env: simpy.Environment):
        yield env.timeout(0)
        raise SelectorCancelled("cancelled")


class SelectorCancelOp(Operation):
    participant_target: Selector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


class AnnihilateOp(Operation):
    participant_container: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [Annihilate(instance_1_iri=self.participant_container)]


def test_selector_ordering_is_stable() -> None:
    alpha = MaterialContainer(identifier="alpha")
    zeta = MaterialContainer(identifier="zeta")
    for obj in (alpha, zeta):
        obj.has_pool_type.add("POOL_ORDERING")
        obj.is_present = {True}
        KnowledgeGraph.get_object_from_lookup(obj.identifier)

    op = TwoSelectorOp(
        participant_alpha=AttributeSelector("POOL_ORDERING", lambda obj: obj.identifier == alpha.identifier),
        participant_zeta=AttributeSelector("POOL_ORDERING", lambda obj: obj.identifier == zeta.identifier),
    )

    sim = Simulation([op])
    sim.run()

    assert op.resources == [alpha.identifier, zeta.identifier]


def test_selector_ordering_with_lambdas_is_stable_across_ops() -> None:
    alpha = MaterialContainer(identifier="alpha-lambda")
    zeta = MaterialContainer(identifier="zeta-lambda")
    for obj in (alpha, zeta):
        obj.has_pool_type.add("POOL_ORDERING_LAMBDA")
        obj.is_present = {True}
        KnowledgeGraph.get_object_from_lookup(obj.identifier)

    op_one = TwoSelectorOp(
        identifier="op-one",
        participant_alpha=AttributeSelector("POOL_ORDERING_LAMBDA", lambda obj: obj.identifier == alpha.identifier),
        participant_zeta=AttributeSelector("POOL_ORDERING_LAMBDA", lambda obj: obj.identifier == zeta.identifier),
    )
    op_two = TwoSelectorOp(
        identifier="op-two",
        participant_alpha=AttributeSelector("POOL_ORDERING_LAMBDA", lambda obj: obj.identifier == alpha.identifier),
        participant_zeta=AttributeSelector("POOL_ORDERING_LAMBDA", lambda obj: obj.identifier == zeta.identifier),
    )
    op_two.required_precedents.append(op_one.identifier)

    sim = Simulation([op_one, op_two])
    sim.run()

    assert op_one.resources == [alpha.identifier, zeta.identifier]
    assert op_two.resources == [alpha.identifier, zeta.identifier]


def test_interrupt_while_waiting_on_lock_request() -> None:
    container = MaterialContainer(identifier="lock-wait")
    container.has_pool_type.add("POOL_LOCK_WAIT")
    container.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(container.identifier)

    op = SelectorWaitOp(
        participant_target=AttributeSelector(
            "POOL_LOCK_WAIT", lambda obj: obj.identifier == container.identifier
        )
    )

    sim = Simulation([op])
    proc = sim.operation_registry[op.identifier]
    rs = get_runtime_state(container, sim.env)

    def locker(env: simpy.Environment) -> simpy.events.Event:
        req = rs.lock.request()
        yield req
        yield env.timeout(1.0)
        rs.lock.release(req)
        FilterStoreRegistry.put_obj_into_filter_store(container, env)

    def interrupter(env: simpy.Environment) -> simpy.events.Event:
        while proc._pre_act_process is None:
            yield env.timeout(0)
        while not rs.lock.queue:
            yield env.timeout(0)
        proc._pre_act_process.interrupt("lock-wait")

    sim.env.process(locker(sim.env))
    proc.simpy_process = sim.env.process(proc.run())
    sim.env.process(interrupter(sim.env))
    sim.env.run()

    store = FilterStoreRegistry.get_filter_store("POOL_LOCK_WAIT", sim.env)
    assert container in store.items
    assert rs.lock.count == 0
    assert not rs.lock.queue
    assert any(event.event_type == "OPERATION_INTERRUPT" for event in sim.history_log)
    assert not any(event.event_type == "OPERATION_START" for event in sim.history_log)


def test_selector_cancel_reports_interrupt() -> None:
    op = SelectorCancelOp(participant_target=CancelSelector())
    sim = Simulation([op])
    sim.run()

    assert any(event.event_type == "OPERATION_INTERRUPT" for event in sim.history_log)
    assert not any(event.event_type == "OPERATION_ABORT" for event in sim.history_log)


def test_rollback_restores_filter_store() -> None:
    lib = Namespace("https://libsyn-sim/kg/")
    shape_iri = str(lib.RollbackShape)
    ttl = f"""
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    lib:RollbackShape a sh:NodeShape ;
       sh:targetSubjectsOf lib:has_capacity ;
       sh:sparql [
         a sh:SPARQLConstraint ;
         sh:select \"\"\"
    SELECT ?this WHERE {{
      ?this lib:has_capacity ?cap .
      FILTER(xsd:double(?cap) > 0)
    }}
    \"\"\" ;
       ] .
    """
    shape_graph = Graph().parse(data=ttl, format="turtle")

    container = MaterialContainer(identifier="rollback-container")
    container.has_pool_type.add("POOL_ROLLBACK")
    container.is_present = {True}
    container.has_capacity.add(1.0)
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()

    op = AnnihilateOp(participant_container=container.identifier)

    sim = Simulation([op], shacl_shapes=shape_graph)
    sim.effect_engine.policy = PolicyBundle(
        per_shape={shape_iri: PolicyRule(severity="hard", disposition="aborted")}
    )
    get_runtime_state(container, sim.env)

    sim.run()

    rs = get_runtime_state(container, sim.env)
    assert container.is_present == {True}
    assert container.has_pool_type == {"POOL_ROLLBACK"}
    assert rs.lock.count == 0
    assert not rs.lock.queue

    FilterStoreRegistry.put_obj_into_filter_store(container, sim.env)
    store = FilterStoreRegistry.get_filter_store("POOL_ROLLBACK", sim.env)
    assert container in store.items
    assert container.is_present == {True}
