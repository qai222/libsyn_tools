# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_exports_and_determinism.py ###
from __future__ import annotations

import pandas as pd
from rdflib import Graph
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize


def _build_sim():
    # stable identifiers for determinism
    src = MaterialContainer(identifier="S")
    dst = MaterialContainer(identifier="D")
    pip = MaterialContainer(identifier="P")
    pom = PortionOfMaterial(identifier="Q")
    pom.add_chemical(Chemical(mass=8.0, density=1.0))
    pom.is_directly_contained_by.add(src)

    for o in (src, dst, pip, pom):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=src.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    op = TransferMaterialByPortionSize(
        identifier="T",
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=pip.identifier,
        portion_size=0.25,
        temporal_cost=2.0,
    )
    return Simulation([op])


def test_exports_and_determinism(tmp_path):
    # Run 1
    sim1 = _build_sim()
    sim1.run()
    log1 = [(r.event_type, r.timestamp, r.operation_id) for r in sim1.history_log]

    # Export files
    csv1 = tmp_path / "log1.csv"
    ih1 = tmp_path / "inst1.csv"
    sim1.export_event_log(csv1)
    sim1.export_instance_history(ih1)
    assert csv1.exists() and ih1.exists()

    # Sanity: CSVs load with expected columns
    df_log = pd.read_csv(csv1)
    df_inst = pd.read_csv(ih1)
    assert {"operation_id", "timestamp", "event_type"} <= set(df_log.columns)
    assert {"instance_iri", "operation_id", "sim_timestamp"} <= set(df_inst.columns)

    # Reset world for run 2 (fresh KG/runtime)
    g: Graph = KnowledgeGraph.graph()
    g.remove((None, None, None))
    from libsyn_tools.sim.operation.runtime import _RESOURCE_MAP, _RUNTIME_CACHE
    from libsyn_tools.sim.knowledge_graph  import LabObject, MaterialContainer, PortionOfMaterial
    for cls in (PortionOfMaterial, MaterialContainer, LabObject):
        try:
            cls.object_lookup.clear()
        except Exception:
            pass
    _RESOURCE_MAP.clear()
    _RUNTIME_CACHE.clear()

    # Run 2 (identical)
    sim2 = _build_sim()
    sim2.run()
    log2 = [(r.event_type, r.timestamp, r.operation_id) for r in sim2.history_log]

    # Deterministic: byte-for-byte in this simplified projection
    assert log1 == log2
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_exports_and_determinism.py ###
