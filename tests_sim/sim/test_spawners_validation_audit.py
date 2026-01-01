from __future__ import annotations

from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.spawner import ValidationAuditSpawner


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
