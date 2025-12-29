# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_overlay_sppt_provider.py ###
from __future__ import annotations

from rdflib import Namespace
from rdflib.namespace import RDF
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph  import (
    MaterialContainer, PortionOfMaterial, Is_directly_contained_by
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize
from libsyn_tools.sim.overlay.sppt_overlay import SPPTOverlayProvider

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
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_overlay_sppt_provider.py ###
