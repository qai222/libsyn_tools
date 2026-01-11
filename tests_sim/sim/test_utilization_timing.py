# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_utilization_timing.py ###
from __future__ import annotations

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.selector import AttributeSelector
from libsyn_tools.sim.operation.unitary_edit import Create, UnitaryEdit


class UseOne(Operation):
    participant_obj: str

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


class UseTwo(Operation):
    participant_a: AttributeSelector
    participant_b: AttributeSelector

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _register_container(container: MaterialContainer) -> None:
    container.is_present = {True}
    KnowledgeGraph.get_object_from_lookup(container.identifier)
    Create(instance_1_iri=container.identifier).apply()


def test_instance_history_records_pool_type_at_usage_time() -> None:
    container = MaterialContainer(identifier="usage-pool-container")
    container.has_pool_type.add("POOL_USAGE_SNAPSHOT")
    _register_container(container)

    op = UseOne(participant_obj=container.identifier)
    sim = Simulation([op])
    sim.run()

    container.has_pool_type = {"POOL_USAGE_NEW"}

    report = sim.build_report()
    row = report.instance_history[report.instance_history["operation_id"] == op.identifier].iloc[0]
    assert row["pool_type"] == "POOL_USAGE_SNAPSHOT"


def test_lock_acquired_time_includes_pre_act_hold() -> None:
    a = MaterialContainer(identifier="preact-a")
    b = MaterialContainer(identifier="preact-b")
    for obj in (a, b):
        obj.has_pool_type.add("POOL_PREACT")
        _register_container(obj)

    op_hold = UseOne(
        identifier="hold-b",
        participant_obj=b.identifier,
        temporal_cost=1.0,
    )
    op_wait = UseTwo(
        identifier="wait-a-b",
        participant_a=AttributeSelector(
            pool_type="POOL_PREACT",
            predicate=lambda obj: obj.identifier == a.identifier,
        ),
        participant_b=AttributeSelector(
            pool_type="POOL_PREACT",
            predicate=lambda obj: obj.identifier == b.identifier,
        ),
        scheduled_start_time=0.1,
    )

    sim = Simulation([op_hold, op_wait])
    sim.run()

    start_time = next(
        r.timestamp
        for r in sim.history_log
        if r.operation_id == op_wait.identifier and r.event_type == "OPERATION_START"
    )
    df = sim._build_instance_history_dataframe()
    row = df[
        (df["instance_iri"] == a.identifier) & (df["operation_id"] == op_wait.identifier)
    ].iloc[0]

    assert row["lock_acquired_sim_timestamp"] is not None
    assert row["lock_acquired_sim_timestamp"] < start_time
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_utilization_timing.py ###
