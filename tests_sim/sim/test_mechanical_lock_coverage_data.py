# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_mechanical_lock_coverage_data.py ###
from __future__ import annotations

from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer, PortionOfMaterial, Has_interrupt_events, Has_ingredient
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import Create, AddDataProperty, UnitaryEdit


class DataWriteOp(Operation):
    participant_locked: str = Field(...)
    target_iri: str = Field(...)
    pom_iri: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            AddDataProperty(
                instance_1_iri=self.target_iri,
                property_iri=Has_interrupt_events.predicate_iri,
                data_value="TOUCHED",
            ),
            AddDataProperty(
                instance_1_iri=self.pom_iri,
                property_iri=Has_ingredient.predicate_iri,
                data_value="H2O",
            ),
        ]


def _make_container(identifier: str) -> MaterialContainer:
    container = MaterialContainer(identifier=identifier)
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()
    return container


def _make_pom(identifier: str) -> PortionOfMaterial:
    pom = PortionOfMaterial(identifier=identifier)
    KnowledgeGraph.get_object_from_lookup(pom.identifier)
    Create(instance_1_iri=pom.identifier).apply()
    return pom


def test_unlocked_labobject_data_write_rejected():
    locked = _make_container("locked-container")
    target = _make_container("target-container")
    pom = _make_pom("pom-1")

    op = DataWriteOp(
        identifier="data-write-bad",
        participant_locked=locked.identifier,
        target_iri=target.identifier,
        pom_iri=pom.identifier,
    )
    sim = Simulation([op])

    sim.run()

    violations = sim.effect_engine._shacl_violations
    assert any(v.origin == "ENGINE" and v.disposition == "aborted" for v in violations)


def test_locked_labobject_data_write_allowed():
    locked = _make_container("locked-container-2")
    target = _make_container("target-container-2")
    pom = _make_pom("pom-2")

    class DataWriteOpWithTarget(DataWriteOp):
        participant_target: str = Field(...)

    op = DataWriteOpWithTarget(
        identifier="data-write-ok",
        participant_locked=locked.identifier,
        participant_target=target.identifier,
        target_iri=target.identifier,
        pom_iri=pom.identifier,
    )
    sim = Simulation([op])
    sim.run()

    types = [r.event_type for r in sim.history_log]
    assert "OPERATION_END" in types
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_mechanical_lock_coverage_data.py ###
