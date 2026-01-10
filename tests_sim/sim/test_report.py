# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_report.py ###
from __future__ import annotations

from pathlib import Path

import simpy
from pydantic import Field

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import AddObjectProperty, Create, UnitaryEdit
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


class AbortOp(Operation):
    participant_container: str
    missing_iri: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            AddObjectProperty(
                instance_1_iri=self.participant_container,
                instance_2_iri=self.missing_iri,
                property_iri=Is_directly_contained_by.predicate_iri,
            )
        ]


class LongOp(Operation):
    temporal_cost: float = Field(default=1.0)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


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
    assert "effect_descriptions" in report.event_log.columns
    start_rows = report.event_log[report.event_log["event_type"] == "OPERATION_START"]
    assert start_rows["effect_descriptions"].iloc[0]


def test_build_report_handles_empty_history() -> None:
    sim = Simulation([])

    report = sim.build_report()

    assert report.summary["makespan"] is None
    assert report.summary["operation_counts"]["in_progress"] == 0
    assert report.event_log.empty


def test_build_report_abort_terminal_makespan() -> None:
    base = "https://libsyn-sim/kg/report-abort/"
    container = MaterialContainer(identifier=f"{base}container")
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()

    op_abort = AbortOp(
        identifier="abort-only",
        participant_container=container.identifier,
        missing_iri=f"{base}missing-object",
    )
    sim = Simulation([op_abort])
    sim.run()

    report = sim.build_report()

    counts = report.summary["operation_counts"]
    assert counts["abort"] == 1
    assert counts["in_progress"] == 0
    assert report.summary["makespan"] is not None


def test_build_report_interrupt_terminal_makespan() -> None:
    op_long = LongOp(identifier="interrupt-only", temporal_cost=1.0)
    sim = Simulation([op_long])
    proc = sim.operation_registry["interrupt-only"]

    def _interrupt(env: simpy.Environment) -> simpy.events.Event:
        yield env.timeout(0.1)
        proc.simpy_process.interrupt("boom")

    sim.env.process(_interrupt(sim.env))
    sim.run()

    report = sim.build_report()

    counts = report.summary["operation_counts"]
    assert counts["interrupt"] == 1
    assert counts["in_progress"] == 0
    assert report.summary["makespan"] is not None
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_report.py ###
