from __future__ import annotations

import simpy

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit


class AbortAfterPostAct(Operation):
    participant_obj: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []

    def post_act(self, env: simpy.Environment) -> None:
        super().post_act(env)
        raise RuntimeError("post-act-failed")


class InterruptAfterPostAct(Operation):
    participant_obj: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []

    def post_act(self, env: simpy.Environment) -> None:
        super().post_act(env)
        raise simpy.Interrupt("forced-interrupt")


def _build_container() -> MaterialContainer:
    container = MaterialContainer()
    container.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    return container


def test_instance_history_uses_abort_timestamp() -> None:
    container = _build_container()
    op = AbortAfterPostAct(participant_obj=container.identifier)

    sim = Simulation([op])
    sim.run()

    df = sim._build_instance_history_dataframe()
    row = df[df["operation_id"] == op.identifier].iloc[0]
    assert row["sim_timestamp"] is not None


def test_instance_history_uses_interrupt_timestamp() -> None:
    container = _build_container()
    op = InterruptAfterPostAct(participant_obj=container.identifier)

    sim = Simulation([op])
    sim.run()

    df = sim._build_instance_history_dataframe()
    row = df[df["operation_id"] == op.identifier].iloc[0]
    assert row["sim_timestamp"] is not None


def test_instance_history_records_identifier_form() -> None:
    container = _build_container()
    op_id = "op/with/iri"
    op = InterruptAfterPostAct(identifier=op_id, participant_obj=container.identifier)

    sim = Simulation([op])
    sim.run()

    df = sim._build_instance_history_dataframe()
    row = df[df["operation_id"] == op_id].iloc[0]
    assert row["operation_id"] == op_id
