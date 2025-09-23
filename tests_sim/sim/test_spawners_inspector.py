# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_spawners_inspector.py ###
from __future__ import annotations

from uuid import uuid4

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
from libsyn_tools.sim.overlay.current_volume_overlay import CurrentVolumeOverlayProvider
from libsyn_tools.sim.spawner import KGInspectorSpawner

SHAPE_IRI = "https://libsyn-sim/shapes/OverflowShape"


def _overflow_shape_named() -> Graph:
    """A named NodeShape so KGInspectorSpawner can dispatch by IRI reliably."""
    lib = Namespace("https://libsyn-sim/kg/")
    exs = Namespace("https://libsyn-sim/shapes/")
    ttl = f'''
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    PREFIX exs: <{exs}>
    exs:OverflowShape a sh:NodeShape ;
       sh:targetSubjectsOf lib:currentVolume ;
       sh:sparql [
         a sh:SPARQLConstraint ;
         sh:select """
    SELECT ?this WHERE {{ ?this lib:currentVolume ?v . FILTER(xsd:double(?v) > 9.9) }}
    """ ;
       ] .
    '''
    return Graph().parse(data=ttl, format="turtle")


def _overflow_world():
    v1 = MaterialContainer(identifier="v1")
    res = MaterialContainer(identifier="res")
    pip = MaterialContainer(identifier="pip")

    p = PortionOfMaterial(identifier="water")
    p.add_chemical(Chemical(mass=12.0, density=1.0))

    for o in (v1, res, pip, p):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=p.identifier,
        instance_2_iri=res.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return v1, res, pip


def _conforming_world():
    v1 = MaterialContainer(identifier="v1c")
    res = MaterialContainer(identifier="resc")
    pip = MaterialContainer(identifier="pipc")

    p = PortionOfMaterial(identifier="waterc")
    p.add_chemical(Chemical(mass=8.0, density=1.0))

    for o in (v1, res, pip, p):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    # fix: use 'res' (not an undefined 'resc') as the container IRI here
    AddObjectProperty(
        instance_1_iri=p.identifier,
        instance_2_iri=res.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return v1, res, pip


def _fix_factory_transfer(v1_id: str, res_id: str, pip_id: str):
    """
    Return a dispatch function that **resolves** overflow by transferring
    a safe fraction from the focus container to the other container.
    Using TransferMaterialByPortionSize keeps coverage intact (new POMs are created).
    """

    def _factory(focus_iri: str):
        src = focus_iri
        dst = res_id if focus_iri == v1_id else v1_id
        return TransferMaterialByPortionSize(
            identifier=f"fix-{focus_iri}-{uuid4().hex[:6]}",
            participant_source=src,
            participant_destination=dst,
            participant_device=pip_id,
            portion_size=0.5,  # move half; enough to get below 9.9 in our test worlds
        )

    return _factory


def test_inspector_per_op_spawns_and_resolves_once():
    v1, res, pip = _overflow_world()

    viol = TransferMaterialByPortionSize(
        identifier="fill",
        participant_source=res.identifier,
        participant_destination=v1.identifier,
        participant_device=pip.identifier,
        portion_size=1.0,
        temporal_cost=0.0,
    )

    sim = Simulation([viol], shacl_shapes=_overflow_shape_named())
    sim.effect_engine.register_overlay_provider(CurrentVolumeOverlayProvider().snapshot)
    KGInspectorSpawner(
        shape_dispatch={SHAPE_IRI: _fix_factory_transfer(v1.identifier, res.identifier, pip.identifier)},
        inspect_interval=0.0,  # hook on OPERATION_END
    ).attach(sim)

    sim.run(until=5.0)

    # After fix is spawned, container should be under threshold.
    assert v1.directly_contained_pom_volume <= 9.9

    # There should be at least one fix op that ran to completion.
    spawned = [r for r in sim.history_log if r.event_type == "OPERATION_END" and r.operation_id.startswith("fix-")]
    assert len(spawned) >= 1


def test_inspector_no_violation_no_spawn():
    v1, res, pip = _conforming_world()

    ok = TransferMaterialByPortionSize(
        identifier="ok",
        participant_source=res.identifier,
        participant_destination=v1.identifier,
        participant_device=pip.identifier,
        portion_size=0.5,
        temporal_cost=0.0,
    )

    sim = Simulation([ok], shacl_shapes=_overflow_shape_named())
    sim.effect_engine.register_overlay_provider(CurrentVolumeOverlayProvider().snapshot)
    KGInspectorSpawner(
        shape_dispatch={SHAPE_IRI: _fix_factory_transfer(v1.identifier, res.identifier, pip.identifier)},
        inspect_interval=0.0,
    ).attach(sim)

    sim.run(until=5.0)
    spawned = [r for r in sim.history_log if r.event_type == "OPERATION_END" and r.operation_id.startswith("fix-")]
    assert len(spawned) == 0


def test_inspector_polling_interval_spawns_after_dt():
    v1, res, pip = _overflow_world()

    viol = TransferMaterialByPortionSize(
        identifier="over",
        participant_source=res.identifier,
        participant_destination=v1.identifier,
        participant_device=pip.identifier,
        portion_size=1.0,
        temporal_cost=0.0,
    )

    sim = Simulation([viol], shacl_shapes=_overflow_shape_named())
    sim.effect_engine.register_overlay_provider(CurrentVolumeOverlayProvider().snapshot)
    KGInspectorSpawner(
        shape_dispatch={SHAPE_IRI: _fix_factory_transfer(v1.identifier, res.identifier, pip.identifier)},
        inspect_interval=3.0,  # poll at t=3
    ).attach(sim)

    # Single run; assert earliest fix starts at or after t=3.0
    sim.run(until=6.0)
    fix_starts = [r.timestamp for r in sim.history_log if
                  r.event_type == "OPERATION_START" and r.operation_id.startswith("fix-")]
    assert fix_starts and min(fix_starts) >= 3.0
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_spawners_inspector.py ###
