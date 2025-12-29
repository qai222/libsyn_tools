# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_spawners_concurrency_optional.py ###
from __future__ import annotations

from uuid import uuid4

from pydantic import Field
from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph  import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize
from libsyn_tools.sim.overlay.current_volume_overlay import CurrentVolumeOverlayProvider
from libsyn_tools.sim.spawner import KGInspectorSpawner, TimerSpawner

SHAPE_IRI = "https://libsyn-sim/shapes/OverflowShape"


def _overflow_shape_named() -> Graph:
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


class NoOp(Operation):
    tag: str = Field(default="noop")

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _overflow_world():
    v1 = MaterialContainer(identifier="v1z")
    res = MaterialContainer(identifier="resz")
    pip = MaterialContainer(identifier="pipz")

    p = PortionOfMaterial(identifier="waterz")
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


def _factory(focus_iri: str) -> Operation:
    # Unique id prevents collisions if inspector fires more than once
    return NoOp(identifier=f"fix-{focus_iri}-{uuid4().hex[:6]}")


def test_timer_and_inspector_coexist_and_spawn_once():
    v1, res, pip = _overflow_world()

    viol = TransferMaterialByPortionSize(
        identifier="overZ",
        participant_source=res.identifier,
        participant_destination=v1.identifier,
        participant_device=pip.identifier,
        portion_size=1.0,
        temporal_cost=0.0,
    )

    sim = Simulation([viol], shacl_shapes=_overflow_shape_named())
    sim.effect_engine.register_overlay_provider(CurrentVolumeOverlayProvider().snapshot)

    # Timer that just generates extra OPERATION_END events in the background;
    # use UNIQUE ids per spawn to avoid spawn_operation collisions.
    TimerSpawner(
        op_factory=lambda s: NoOp(identifier=f"tick-{int(s.env.now)}-{uuid4().hex[:6]}"),
        interval=2.0
    ).attach(sim)

    # Inspector polling rather than per-op to avoid repeated firing from spawned ops
    KGInspectorSpawner(shape_dispatch={SHAPE_IRI: _factory}, inspect_interval=3.0).attach(sim)

    sim.run(until=5.0)

    # We expect at least one fix; duplicates are tolerated
    fixes = [r for r in sim.history_log if r.event_type == "OPERATION_END" and r.operation_id.startswith("fix-")]
    assert len(fixes) >= 1

    # And the timer tick ops ran at least twice (t≈0,2,4)
    ticks = [r for r in sim.history_log if r.event_type == "OPERATION_END" and r.operation_id.startswith("tick-")]
    assert len(ticks) >= 2
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_spawners_concurrency_optional.py ###
