# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_overlay_sppt_provider.py ###
from __future__ import annotations

from urllib.parse import quote

from rdflib import Namespace, URIRef
from rdflib.namespace import RDF
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph  import (
    MaterialContainer, PortionOfMaterial, Is_directly_contained_by
)
from libsyn_tools.sim.knowledge_graph.ontology import Has_begin_time_base, Has_end_time_base
from libsyn_tools.sim.lifecycle import LifecycleCallbacks
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import RemoveObjectProperty, UnitaryEdit
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize
from libsyn_tools.sim.overlay.sppt_overlay import SPPTOverlayProvider, _OpSpan

LIB = Namespace("https://libsyn-sim/kg/")


def test_sppt_overlay_emits_process_interval_and_participants():
    # Build a minimal world and a single transfer operation
    src = MaterialContainer()
    dst = MaterialContainer()
    dev = MaterialContainer()  # device
    pom = PortionOfMaterial()
    pom.add_chemical(Chemical(mass=10.0, density=1.0))

    for o in (src, dst, dev, pom):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=src.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    op = TransferMaterialByPortionSize(
        identifier="X",
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=dev.identifier,
        portion_size=0.5,
        temporal_cost=0.01,
    )

    sim = Simulation([op])
    # Explicitly attach the SPPT provider (already registered by default, but we want a handle)
    provider = SPPTOverlayProvider(sim.callbacks)
    sim.effect_engine.register_overlay_provider(provider.snapshot)

    sim.run()

    g = provider.snapshot()

    # Process node for op "X"
    proc_iri = LIB[f"Process/{op.identifier}"]
    assert (proc_iri, RDF.type, LIB.Process) in g

    # Interval node exists and is linked
    int_iri = LIB[f"Interval/{op.identifier}"]
    assert (int_iri, RDF.type, LIB.TimeInterval) in g
    assert (proc_iri, LIB.has_interval, int_iri) in g

    # Participants include src, dst, and device
    assert (proc_iri, LIB.has_participant, LIB[src.identifier]) in g
    assert (proc_iri, LIB.has_participant, LIB[dst.identifier]) in g
    assert (proc_iri, LIB.has_participant, LIB[dev.identifier]) in g


class BadRemove(Operation):
    participant_src: str
    bad_dst: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            RemoveObjectProperty(
                instance_1_iri=self.participant_src,
                instance_2_iri=self.bad_dst,
                property_iri=Is_directly_contained_by.predicate_iri,
            )
        ]


def test_sppt_overlay_closes_interval_on_abort():
    src = MaterialContainer()
    KnowledgeGraph.get_object_from_lookup(src.identifier)
    Create(instance_1_iri=src.identifier).apply()

    op = BadRemove(
        identifier="abort-op",
        participant_src=src.identifier,
        bad_dst="missing-dst",
    )

    sim = Simulation([op])
    provider = SPPTOverlayProvider(sim.callbacks)
    sim.effect_engine.register_overlay_provider(provider.snapshot)

    sim.run()

    g = provider.snapshot()
    int_iri = LIB[f"Interval/{op.identifier}"]
    assert (int_iri, LIB.has_begin_time, None) in g
    assert (int_iri, LIB.has_end_time, None) in g


def test_sppt_overlay_escapes_operation_identifier_and_preserves_metadata() -> None:
    op_id = "unsafe op/with space"
    span = _OpSpan(t0=0.0, t1=1.0, participants=[])
    provider = SPPTOverlayProvider(LifecycleCallbacks())
    escaped = quote(op_id, safe="")
    proc_iri = LIB[f"Process/{escaped}"]
    custom_pred = LIB.customNote
    provider._g.add((proc_iri, custom_pred, URIRef("https://example.com/meta")))

    provider._materialize_span(op_id, span)

    assert " " not in str(proc_iri)
    assert (proc_iri, RDF.type, LIB.Process) in provider._g
    assert (proc_iri, custom_pred, URIRef("https://example.com/meta")) in provider._g


def test_sppt_overlay_emits_base_time_when_speed_factor_used() -> None:
    src = MaterialContainer()
    KnowledgeGraph.get_object_from_lookup(src.identifier)
    Create(instance_1_iri=src.identifier).apply()

    op = BadRemove(
        identifier="base-time-op",
        participant_src=src.identifier,
        bad_dst="missing-dst",
    )

    sim = Simulation([op], simulation_speed_factor=2.0)
    provider = SPPTOverlayProvider(sim.callbacks)
    sim.effect_engine.register_overlay_provider(provider.snapshot)
    sim.run()

    g = provider.snapshot()
    int_iri = LIB[f"Interval/{op.identifier}"]
    base_begin = URIRef(Has_begin_time_base.predicate_iri)
    base_end = URIRef(Has_end_time_base.predicate_iri)

    assert (int_iri, base_begin, None) in g
    assert (int_iri, base_end, None) in g
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_overlay_sppt_provider.py ###
