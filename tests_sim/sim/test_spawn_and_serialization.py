# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_spawn_and_serialization.py ###
from __future__ import annotations

from rdflib import Graph
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph.physical_entities import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize


class CalOp(Operation):
    """No participants, no edits — used for spawn test."""

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _multi_resource_world():
    src = MaterialContainer(identifier="SRC")
    d1 = MaterialContainer(identifier="D1")
    d2 = MaterialContainer(identifier="D2")
    dev = MaterialContainer(identifier="PIP")
    pom = PortionOfMaterial(identifier="POM")
    pom.add_chemical(Chemical(mass=10.0, density=1.0))
    for o in (src, d1, d2, dev, pom):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=src.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return src, d1, d2, dev


def test_spawn_operation_midrun(tmp_path):
    sim = Simulation.compile_actions()  # empty

    # Spawn a simple op at t=1.0
    def spawner():
        yield sim.env.timeout(1.0)
        sim.spawn_operation(CalOp(identifier="SPN", temporal_cost=0.5))

    sim.env.process(spawner())
    sim.run(until=3.0)

    types = [r.event_type for r in sim.history_log if r.operation_id == "SPN"]
    assert "OPERATION_START" in types and "OPERATION_END" in types


def test_multi_resource_serialization_is_deterministic():
    # Two ops need the same source and device; they should serialize in a stable order
    src, d1, d2, dev = _multi_resource_world()
    op1 = TransferMaterialByPortionSize(
        identifier="A",
        participant_source=src.identifier,
        participant_destination=d1.identifier,
        participant_device=dev.identifier,
        portion_size=0.5,
        temporal_cost=1.0,
    )
    op2 = TransferMaterialByPortionSize(
        identifier="B",
        participant_source=src.identifier,
        participant_destination=d2.identifier,
        participant_device=dev.identifier,
        portion_size=0.5,
        temporal_cost=1.0,
    )

    sim1 = Simulation([op1, op2])
    sim1.run()
    log1 = [(r.event_type, r.timestamp, r.operation_id) for r in sim1.history_log]

    # Reset world
    g: Graph = KnowledgeGraph.graph()
    g.remove((None, None, None))
    from libsyn_tools.sim.operation.runtime import _RESOURCE_MAP, _RUNTIME_CACHE
    from libsyn_tools.sim.operation.selector import FilterStoreRegistry
    from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject
    for cls in (PortionOfMaterial, MaterialContainer, LabObject):
        try:
            cls.object_lookup.clear()
        except Exception:
            setattr(cls, "object_lookup", {})
    _RESOURCE_MAP.clear()
    _RUNTIME_CACHE.clear()
    FilterStoreRegistry._stores.clear()

    # Recreate world identically and run again
    src, d1, d2, dev = _multi_resource_world()
    op1 = TransferMaterialByPortionSize(
        identifier="A",
        participant_source=src.identifier,
        participant_destination=d1.identifier,
        participant_device=dev.identifier,
        portion_size=0.5,
        temporal_cost=1.0,
    )
    op2 = TransferMaterialByPortionSize(
        identifier="B",
        participant_source=src.identifier,
        participant_destination=d2.identifier,
        participant_device=dev.identifier,
        portion_size=0.5,
        temporal_cost=1.0,
    )
    sim2 = Simulation([op1, op2])
    sim2.run()
    log2 = [(r.event_type, r.timestamp, r.operation_id) for r in sim2.history_log]

    assert log1 == log2
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_spawn_and_serialization.py ###
