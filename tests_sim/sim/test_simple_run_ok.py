# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_simple_run_ok.py ###
from __future__ import annotations

import simpy
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph  import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize


def _world_one_pom(vol_ml: float):
    """Source container with exactly one POM (known volume), empty dest, one pipette."""
    src = MaterialContainer(identifier="src")
    dst = MaterialContainer(identifier="dst")
    pip = MaterialContainer(identifier="pip")  # treat as device LabObject for locking

    pom = PortionOfMaterial(identifier="pomA")
    pom.add_chemical(Chemical(mass=vol_ml, density=1.0))  # 1 g/mL → volume=mass
    pom.is_directly_contained_by.add(src)

    # Seed KG: present everything & link pom→src
    for o in (src, dst, pip, pom):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=src.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return src, dst, pip


def test_simple_run_splits_and_moves(env: simpy.Environment):
    src, dst, pip = _world_one_pom(10.0)

    op = TransferMaterialByPortionSize(
        identifier="opX",
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=pip.identifier,
        portion_size=0.3,  # 30% of the only POM: 3 mL
        temporal_cost=1.0,
    )
    sim = Simulation([op])
    sim.run()  # full run

    # Volumes updated via POM reconstitution
    assert abs(src.directly_contained_pom_volume - 7.0) < 1e-9
    assert abs(dst.directly_contained_pom_volume - 3.0) < 1e-9

    # Event log integrity
    types = [r.event_type for r in sim.history_log]
    assert "OPERATION_START" in types and "OPERATION_END" in types
    # ordering: start before end
    t_start = next(r.timestamp for r in sim.history_log if r.event_type == "OPERATION_START")
    t_end = next(r.timestamp for r in sim.history_log if r.event_type == "OPERATION_END")
    assert t_start <= t_end
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_simple_run_ok.py ###
