# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_core_shapes.py ###
from __future__ import annotations

from uuid import uuid4

from rdflib import Namespace
from rdflib.namespace import SH
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.effect_shacl import _iter_validation_results, _first
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.remediation import make_drain_to_capacity
from libsyn_tools.sim.operation_preset.wait import Wait
from libsyn_tools.sim.shapes import load_core_shapes
from libsyn_tools.sim.spawner import PolicyEnforcerSpawner, ValidationAuditSpawner

_EPS = 1e-6


def _make_overfilled_container(identifier: str, volume: float, capacity: float) -> MaterialContainer:
    container = MaterialContainer(identifier=identifier)
    container.has_capacity.add(capacity)
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


def test_core_shapes_detect_overfill() -> None:
    shapes = load_core_shapes()
    lib = Namespace("https://libsyn-sim/kg/")
    _make_overfilled_container(f"{lib}overfill", volume=1.5, capacity=1.0)

    audit_anchor = Wait(identifier="VALIDATION_AUDIT", temporal_cost=0.0)
    sim = Simulation([audit_anchor], shacl_shapes=shapes)
    conforms, report, _ = sim.effect_engine.validate_now()

    assert not conforms
    assert any(
        _first(report, vr, SH.sourceShape) == str(lib.CapacityShape)
        for vr in _iter_validation_results(report)
    )


def test_core_shapes_remediation_via_audit_and_enforcer() -> None:
    shapes = load_core_shapes()
    lib = Namespace("https://libsyn-sim/kg/")

    target = _make_overfilled_container(f"{lib}target", volume=1.5, capacity=1.0)
    waste = MaterialContainer(identifier=f"{lib}waste")
    for obj in (waste,):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    def _drain_factory(record):
        if record.focus_iri is None:
            return None
        op = make_drain_to_capacity(record.focus_iri, waste.identifier)
        op.identifier = f"drain-{uuid4().hex[:6]}"
        op.temporal_cost = 0.0
        return op

    audit_anchor = Wait(identifier="VALIDATION_AUDIT", temporal_cost=0.0)
    sim = Simulation([audit_anchor], shacl_shapes=shapes)
    ValidationAuditSpawner(inspect_interval=1.0).attach(sim)
    PolicyEnforcerSpawner(shape_dispatch={str(lib.CapacityShape): _drain_factory}).attach(sim)

    sim.run(until=2.1)

    assert target.directly_contained_pom_volume <= target.capacity + _EPS
    spawned = [r for r in sim.history_log if r.event_type == "OPERATION_END" and r.operation_id.startswith("drain-")]
    assert len(spawned) >= 1
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_core_shapes.py ###
