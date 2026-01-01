# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_report.py ###
from __future__ import annotations

from pathlib import Path

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


def _build_sim() -> Simulation:
    src = MaterialContainer(identifier="SRC")
    dst = MaterialContainer(identifier="DST")
    pip = MaterialContainer(identifier="PIP")
    pom = PortionOfMaterial(identifier="POM")
    pom.add_chemical(Chemical(mass=4.0, density=1.0))
    pom.is_directly_contained_by.add(src)

    for obj in (src, dst, pip, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=src.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    op = TransferMaterialByPortionSize(
        identifier="transfer-one",
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=pip.identifier,
        portion_size=0.5,
        temporal_cost=1.0,
    )
    return Simulation([op])


def test_run_and_report_outputs(tmp_path: Path) -> None:
    sim = _build_sim()
    out_dir = tmp_path / "report"
    report = sim.run_and_report(out_dir)

    assert report.summary.get("makespan") is not None
    assert (out_dir / "event_log.csv").exists()
    assert (out_dir / "summary.json").exists()
    summary_md = out_dir / "summary.md"
    assert summary_md.exists()
    content = summary_md.read_text(encoding="utf-8")
    assert "Run Summary" in content
    assert "Violations by Shape/Disposition" in content
    assert "| shape_iri | disposition | count |" in content
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_report.py ###
