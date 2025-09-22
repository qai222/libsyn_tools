# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_shacl_raise_policy.py ###
from __future__ import annotations

import pytest
from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.effect_engine import SHACLValidationError
from libsyn_tools.sim.knowledge_graph.physical_entities import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize


def _overflow_shape() -> Graph:
    """Flag any container whose overlay currentVolume > 9.9 (soft policy)."""
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


def _world():
    v1 = MaterialContainer(identifier="v1")
    res = MaterialContainer(identifier="res")
    pip = MaterialContainer(identifier="pip")

    p = PortionOfMaterial(identifier="water")
    p.add_chemical(Chemical(mass=12.0, density=1.0))

    # seed base graph: present & link p -> res
    for o in (v1, res, pip, p):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=p.identifier,
        instance_2_iri=res.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return v1, res, pip


def test_shacl_raise_policy_commits_then_raises():
    v1, res, pip = _world()

    # Move all 12 mL from reservoir → v1 → overlay will violate (>9.9)
    op = TransferMaterialByPortionSize(
        identifier="fill",
        participant_source=res.identifier,
        participant_destination=v1.identifier,
        participant_device=pip.identifier,
        portion_size=1.0,
        temporal_cost=0.0,
    )
    sim = Simulation([op], shacl_shapes=_overflow_shape())
    # Turn on raise behavior
    sim.effect_engine.raise_shacl = True

    with pytest.raises(SHACLValidationError):
        sim.run()

    # State was committed (policy is validated after edits apply)
    assert v1.directly_contained_pom_volume > 9.9

    # We started, but the END event likely isn't logged because the raise
    # happens inside EffectEngine.apply (before post_act / END).
    types = [r.event_type for r in sim.history_log if r.operation_id == "fill"]
    assert "OPERATION_START" in types
    assert "OPERATION_END" not in types
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_shacl_raise_policy.py ###
