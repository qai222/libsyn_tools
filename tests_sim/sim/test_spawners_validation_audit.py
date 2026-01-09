from __future__ import annotations

from uuid import uuid4

from pydantic import Field
from rdflib import Graph, Namespace, BNode, URIRef
from rdflib.namespace import SH, XSD, RDF
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.selector import FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty, UnitaryEdit
from libsyn_tools.sim.spawner import ValidationAuditSpawner, PolicyEnforcerSpawner


def _capacity_overflow_shape() -> tuple[Graph, str]:
    lib = Namespace("https://libsyn-sim/kg/")
    shape_iri = str(lib.CapacityAuditShape)
    ttl = f"""
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    lib:CapacityAuditShape a sh:NodeShape ;
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


def test_validation_audit_emits_violation_on_poll():
    shape_graph, shape_iri = _capacity_overflow_shape()
    base = "https://libsyn-sim/kg/"

    container = MaterialContainer(identifier=f"{base}audit-dest")
    container.has_capacity.add(1.0)
    pom = PortionOfMaterial(identifier=f"{base}audit-pom")
    pom.add_chemical(Chemical(mass=1.2, density=1.0))
    for obj in (container, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()
    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=container.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    sim = Simulation([], shacl_shapes=shape_graph)
    audit = ValidationAuditSpawner(inspect_interval=1.0)
    audit.attach(sim)

    records = []

    def _capture(record):
        records.append(record)

    sim.callbacks.on_violation.append(_capture)

    sim.run(until=1.1)

    assert any(r.shape_iri == shape_iri for r in records)
    assert any(r.operation_id == "VALIDATION_AUDIT" for r in records)


class _AuditRemediation(Operation):
    participant_container: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_validation_audit_policy_enforcer_skips_missing_precedent():
    shape_graph, shape_iri = _capacity_overflow_shape()
    base = "https://libsyn-sim/kg/"

    container = MaterialContainer(identifier=f"{base}audit-policy-dest")
    container.has_capacity.add(1.0)
    pom = PortionOfMaterial(identifier=f"{base}audit-policy-pom")
    pom.add_chemical(Chemical(mass=1.2, density=1.0))
    for obj in (container, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()
    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=container.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    def _remediation_factory(record):
        if record.focus_iri is None or record.operation_id != "VALIDATION_AUDIT":
            return None
        return _AuditRemediation(
            identifier=f"remediate-{uuid4().hex[:6]}",
            participant_container=record.focus_iri,
        )

    sim = Simulation([], shacl_shapes=shape_graph)
    FilterStoreRegistry.put_obj_into_filter_store(container, sim.env)
    ValidationAuditSpawner(inspect_interval=1.0).attach(sim)
    PolicyEnforcerSpawner(shape_dispatch={shape_iri: _remediation_factory}).attach(sim)

    sim.run(until=1.1)

    spawned = [
        r for r in sim.history_log
        if r.event_type == "OPERATION_END" and r.operation_id.startswith("remediate-")
    ]
    assert spawned


def test_validation_audit_skips_invalid_focus_nodes():
    sim = Simulation([])
    audit = ValidationAuditSpawner(inspect_interval=1.0)

    records = []

    def _capture(record):
        records.append(record)

    sim.callbacks.on_violation.append(_capture)

    report = Graph()
    vr_blank = URIRef("urn:vr:blank")
    vr_unknown = URIRef("urn:vr:unknown")
    report.add((vr_blank, RDF.type, SH.ValidationResult))
    report.add((vr_blank, SH.sourceShape, URIRef("urn:shape:blank")))
    report.add((vr_blank, SH.focusNode, BNode()))
    report.add((vr_unknown, RDF.type, SH.ValidationResult))
    report.add((vr_unknown, SH.sourceShape, URIRef("urn:shape:unknown")))
    report.add((vr_unknown, SH.focusNode, URIRef("urn:missing:focus")))

    audit._emit_violation_records(sim, report)

    assert records == []
