# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_interrupt_and_overlay.py ###
from __future__ import annotations

import simpy
from rdflib import Namespace
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    LabObject,
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.selector import LiteralSelector
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, Create, AddObjectProperty


# ---------- Interrupt path ----------

class LongOp(Operation):
    """Long-running op with a participant so interrupt writes a reason onto it."""

    participant_dev: LiteralSelector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_interrupt_writes_reason_and_no_end(env: simpy.Environment):
    dev = LabObject()
    dev.has_pool_type.add("DEV")
    KnowledgeGraph.get_object_from_lookup(dev.identifier)

    op = LongOp(
        identifier="L",
        participant_dev=LiteralSelector(dev.identifier),
        temporal_cost=100.0,  # long
    )
    sim = Simulation([op])

    # Schedule an interrupt at t=1
    def interrupter():
        yield sim.env.timeout(1.0)
        proc = sim.operation_registry["L"].simpy_process
        proc.interrupt("USR")

    sim.env.process(interrupter())
    sim.run(until=5.0)

    # Logs: interrupt recorded, no END for this op
    types = [r.event_type for r in sim.history_log if r.operation_id == "L"]
    assert "OPERATION_INTERRUPT" in types
    assert "OPERATION_END" not in types

    # Participant has interrupt reason recorded (via EffectEngine.apply inside except block)
    # Reason format: "<opid>:<cause>"
    assert "L:USR" in getattr(dev, "has_interrupt_events", set())


# ---------- Overlay ephemerality ----------

def _overflow_world():
    """Container with a 12 mL POM linked via is_directly_contained_by."""
    c = MaterialContainer(identifier="C")
    p = PortionOfMaterial(identifier="P")
    p.add_chemical(Chemical(mass=12.0, density=1.0))
    KnowledgeGraph.get_object_from_lookup(c.identifier)
    KnowledgeGraph.get_object_from_lookup(p.identifier)
    Create(instance_1_iri=c.identifier).apply()
    Create(instance_1_iri=p.identifier).apply()
    AddObjectProperty(
        instance_1_iri=p.identifier,
        instance_2_iri=c.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()


def test_overlay_not_persisted_in_base_graph():
    _overflow_world()
    sim = Simulation([])  # no ops needed; just want to inspect KG graph
    g = KnowledgeGraph.graph()

    lib = Namespace("https://libsyn-sim/kg/")
    # Overlay adds lib:currentVolume only in the ephemeral union during SHACL,
    # so the base graph should have zero such triples.
    assert len(list(g.triples((None, lib.currentVolume, None)))) == 0
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_interrupt_and_overlay.py ###
