# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_operation_abort_flow.py ###
from __future__ import annotations

from pydantic import Field
from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize
from libsyn_tools.sim.overlay.current_volume_overlay import CurrentVolumeOverlayProvider
from libsyn_tools.sim.policy import PolicyBundle, PolicyRule
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty, UnitaryEdit


def _capacity_overflow_shape() -> tuple[Graph, str]:
    lib = Namespace("https://libsyn-sim/kg/")
    shape_iri = str(lib.CapacityOverflowShape)
    ttl = f'''
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    lib:CapacityOverflowShape a sh:NodeShape ;
       sh:targetSubjectsOf lib:currentVolume ;
       sh:sparql [
         a sh:SPARQLConstraint ;
         sh:select """
    SELECT ?this WHERE {{
      ?this lib:currentVolume ?v ;
            lib:has_capacity ?cap .
      FILTER(xsd:double(?v) > xsd:double(?cap))
    }}
    """ ;
       ] .
    '''
    return Graph().parse(data=ttl, format="turtle"), shape_iri


def _world_reservoir(amount: float):
    base = "https://libsyn-sim/kg/"
    v1 = MaterialContainer(identifier=f"{base}v1")
    res = MaterialContainer(identifier=f"{base}reservoir")
    pip = MaterialContainer(identifier=f"{base}pipette")
    p = PortionOfMaterial(identifier=f"{base}water")
    p.add_chemical(Chemical(mass=amount, density=1.0))
    for o in (v1, res, pip, p):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=p.identifier,
        instance_2_iri=res.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return v1, res, pip


class NoopOp(Operation):
    name: str = Field("noop")

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def test_abort_releases_locks_and_completes():
    shape_graph, shape_iri = _capacity_overflow_shape()
    policy = PolicyBundle(per_shape={shape_iri: PolicyRule(severity="hard", disposition="aborted")})

    v1, res, pip = _world_reservoir(12.0)
    v1.has_capacity.add(5.0)

    op_a = TransferMaterialByPortionSize(
        identifier="fill",
        participant_source=res.identifier,
        participant_destination=v1.identifier,
        participant_device=pip.identifier,
        portion_size=1.0,
        temporal_cost=0.0,
    )
    op_b = NoopOp(
        identifier="after",
        required_precedents=[op_a.identifier],
        temporal_cost=0.0,
    )

    sim = Simulation([op_a, op_b], shacl_shapes=shape_graph)
    sim.effect_engine.policy = policy
    sim.effect_engine.register_overlay_provider(CurrentVolumeOverlayProvider().snapshot)

    sim.run()

    types = [r.event_type for r in sim.history_log if r.operation_id == op_a.identifier]
    assert "OPERATION_ABORT" in types

    types_b = [r.event_type for r in sim.history_log if r.operation_id == op_b.identifier]
    assert "OPERATION_START" in types_b
    assert "OPERATION_END" in types_b
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_operation_abort_flow.py ###
