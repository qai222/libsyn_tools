# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_precedents_ordering.py ###
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


def _seed_two_targets():
    src = MaterialContainer(identifier="srcB")
    d1 = MaterialContainer(identifier="d1")
    d2 = MaterialContainer(identifier="d2")
    pip = MaterialContainer(identifier="pipB")
    pom = PortionOfMaterial(identifier="pomB")
    pom.add_chemical(Chemical(mass=10.0, density=1.0))
    pom.is_directly_contained_by.add(src)

    for o in (src, d1, d2, pip, pom):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=src.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return src, d1, d2, pip


def test_precedence_enforced_in_run(env: simpy.Environment):
    src, d1, d2, pip = _seed_two_targets()

    op1 = TransferMaterialByPortionSize(
        identifier="first",
        participant_source=src.identifier,
        participant_destination=d1.identifier,
        participant_device=pip.identifier,
        portion_size=0.5,
        temporal_cost=1.0,
    )
    op2 = TransferMaterialByPortionSize(
        identifier="second",
        participant_source=src.identifier,
        participant_destination=d2.identifier,
        participant_device=pip.identifier,
        portion_size=0.5,
        temporal_cost=1.0,
        required_precedents=[op1.identifier],
    )
    sim = Simulation([op1, op2])
    sim.run()

    # START(first) < END(first) ≤ START(second) < END(second)
    def t(ev, opid):  # pick first timestamp for that op & event
        return next(r.timestamp for r in sim.history_log if r.event_type == ev and r.operation_id == opid)

    assert t("OPERATION_START", "first") <= t("OPERATION_END", "first") <= t("OPERATION_START", "second") <= t(
        "OPERATION_END", "second")
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_precedents_ordering.py ###
