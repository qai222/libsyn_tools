from __future__ import annotations

from uuid import uuid4

from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.effect_shacl import SHACLViolationRecord
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty, UnitaryEdit
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByVolume
from libsyn_tools.sim.remediation import make_drain_to_capacity
from libsyn_tools.sim.spawner import PolicyEnforcerSpawner

_EPS = 1e-6


def _capacity_overfill_shape() -> tuple[Graph, str]:
    lib = Namespace("https://libsyn-sim/kg/")
    shape_iri = str(lib.CapacityCheckShape)
    ttl = f"""
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    lib:CapacityCheckShape a sh:NodeShape ;
       sh:targetSubjectsOf lib:currentVolume ;
       sh:sparql [
         a sh:SPARQLConstraint ;
         sh:select \"\"\"
    SELECT ?this WHERE {{
      ?this lib:currentVolume ?v ;
            lib:has_capacity ?cap .
      FILTER(xsd:double(?v) > xsd:double(?cap))
    }}
    \"\"\" ;
       ] .
    """
    return Graph().parse(data=ttl, format="turtle"), shape_iri


def _make_container_with_volume(identifier: str, volume: float) -> MaterialContainer:
    container = MaterialContainer(identifier=identifier)
    pom = PortionOfMaterial(identifier=f"{identifier}/pom")
    pom.add_chemical(Chemical(mass=volume, density=1.0))
    for obj in (container, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()
    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=container.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return container


def test_policy_enforcer_spawns_drain_on_committed_violation():
    shape_graph, shape_iri = _capacity_overfill_shape()
    base = "https://libsyn-sim/kg/"

    source = _make_container_with_volume(f"{base}source", 0.5)
    destination = _make_container_with_volume(f"{base}dest", 0.8)
    destination.has_capacity.add(1.0)
    waste = MaterialContainer(identifier=f"{base}waste")
    device = MaterialContainer(identifier=f"{base}device")
    for obj in (waste, device):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    op = TransferMaterialByVolume(
        identifier="fill",
        participant_source=source.identifier,
        participant_destination=destination.identifier,
        participant_device=device.identifier,
        transfer_volume=0.5,
        temporal_cost=0.0,
    )

    def _drain_factory(record: SHACLViolationRecord):
        if record.focus_iri is None:
            return None
        op = make_drain_to_capacity(record.focus_iri, waste.identifier)
        op.identifier = f"drain-{uuid4().hex[:6]}"
        op.temporal_cost = 0.0
        return op

    sim = Simulation([op], shacl_shapes=shape_graph)
    PolicyEnforcerSpawner(shape_dispatch={shape_iri: _drain_factory}).attach(sim)

    sim.run(until=5.0)

    assert destination.directly_contained_pom_volume <= destination.capacity + _EPS
    spawned = [r for r in sim.history_log if r.event_type == "OPERATION_END" and r.operation_id.startswith("drain-")]
    assert len(spawned) >= 1


class _NoOp(Operation):
    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_policy_enforcer_dedupe_includes_focus():
    shape_iri = "urn:shape:dedupe"

    def _factory(record: SHACLViolationRecord):
        if record.focus_iri is None:
            return None
        return _NoOp(identifier=f"remediate-{record.focus_iri}")

    sim = Simulation([])
    spawner = PolicyEnforcerSpawner(shape_dispatch={shape_iri: _factory}, dedupe=True)
    spawner.attach(sim)

    rec_a = SHACLViolationRecord(
        sim_time=0.0,
        operation_id="op-1",
        origin="SHACL",
        severity="soft",
        disposition="committed",
        shape_iri=shape_iri,
        focus_iri="focus-a",
    )
    rec_b = SHACLViolationRecord(
        sim_time=0.0,
        operation_id="op-1",
        origin="SHACL",
        severity="soft",
        disposition="committed",
        shape_iri=shape_iri,
        focus_iri="focus-b",
    )

    sim.callbacks.emit_violation(rec_a)
    sim.callbacks.emit_violation(rec_b)

    spawned = [op_id for op_id in sim.operation_registry if op_id.startswith("remediate-")]
    assert sorted(spawned) == ["remediate-focus-a", "remediate-focus-b"]
