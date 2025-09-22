# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_shacl_soft_violation_commits.py ###
from __future__ import annotations

import simpy
from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph.physical_entities import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize


def _overflow_shape() -> Graph:
    """Flag any container whose overlay currentVolume > 9.9."""
    lib = Namespace("https://libsyn-sim/kg/")
    ttl = f'''
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    [] a sh:NodeShape ;
       sh:targetSubjectsOf lib:currentVolume ;
       sh:sparql [
         a sh:SPARQLConstraint ;
         sh:select """
    SELECT ?this WHERE {{ ?this lib:currentVolume ?v . FILTER(xsd:double(?v) > 9.9) }}
    """ ;
       ] .
    '''
    return Graph().parse(data=ttl, format="turtle")


def _world_reservoir(amount: float):
    """Reservoir with a single POM of `amount` mL; empty v1; one pipette."""
    v1 = MaterialContainer(identifier="v1")
    res = MaterialContainer(identifier="reservoir")
    pip = MaterialContainer(identifier="pipette")

    pom = PortionOfMaterial(identifier="water")
    pom.add_chemical(Chemical(mass=amount, density=1.0))
    pom.is_directly_contained_by.add(res)

    for o in (v1, res, pip, pom):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=res.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return v1, res, pip


def test_shacl_violation_is_soft_commit(env: simpy.Environment):
    v1, res, pip = _world_reservoir(12.0)  # move all → v1 volume=12 (>9.9)

    op = TransferMaterialByPortionSize(
        identifier="fill",
        participant_source=res.identifier,
        participant_destination=v1.identifier,
        participant_device=pip.identifier,
        portion_size=1.0,
        temporal_cost=0.5,
    )
    sim = Simulation([op], shacl_shapes=_overflow_shape())
    sim.run()

    # Edits committed
    assert v1.directly_contained_pom_volume > 9.9

    # SHACL recorded a soft, committed violation
    viols = sim.effect_engine._shacl_violations  # test access
    assert any(v.origin == "SHACL" and v.severity == "soft" and v.disposition == "committed" for v in viols)
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_shacl_soft_violation_commits.py ###
