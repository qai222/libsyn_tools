from __future__ import annotations

"""
Examples
--------
KgQuerySelector (SPARQL over KG + overlays):
    from libsyn_tools.sim.operation.selector import KgQuerySelector
    selector = KgQuerySelector(
        pool_type="VIAL",
        sparql=\"\"\"
        PREFIX lib: <https://libsyn-sim/kg/>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
        SELECT ?s WHERE { ?s lib:currentVolume ?v . FILTER(xsd:double(?v) > 0) }
        \"\"\",
    )

compile_schedule_to_simulation (planned ops + SchedulerOutput):
    from libsyn_tools.chem_schema import Operation, OperationType
    from libsyn_tools.opt import SchedulerOutput
    from libsyn_tools.sim.adapters.schedule_bridge import compile_schedule_to_simulation

    planned = [Operation(identifier="op1", type=OperationType.TransferLiquid)]
    schedule = SchedulerOutput(
        start_times={"op1": 0.0},
        end_times={"op1": 1.0},
        assignments={"op1": "module_1"},
    )
    sim = compile_schedule_to_simulation(planned, schedule)
    sim.run()
"""

from collections.abc import Callable, Iterable

from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import FunctionalModule
from libsyn_tools.chem_schema import Operation as PlannedOperation
from libsyn_tools.chem_schema import OperationNetwork
from libsyn_tools.opt import SchedulerOutput
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import LabObject
from libsyn_tools.sim.operation.operation import Operation as SimOperation
from libsyn_tools.sim.operation.operation import StrOrSelector
from libsyn_tools.sim.operation.unitary_edit import Create, UnitaryEdit


class ExecuteScheduledTask(SimOperation):
    participant_module: StrOrSelector
    planned_type: str | None = None
    planned_annotations: dict = Field(default_factory=dict)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


def _iter_planned_operations(
    planned_ops: list[PlannedOperation] | OperationNetwork,
) -> list[PlannedOperation]:
    if isinstance(planned_ops, OperationNetwork):
        return planned_ops.operations
    return list(planned_ops)


def _ensure_modules_present(module_ids: Iterable[str]) -> None:
    for module_id in module_ids:
        existing = KnowledgeGraph.get_object_from_lookup(module_id)
        if existing is None:
            module = LabObject(identifier=module_id)
            module.has_pool_type.add("MODULE")
            Create(instance_1_iri=module.identifier).apply()
        else:
            if existing.is_present != {True}:
                Create(instance_1_iri=existing.identifier).apply()
            if "MODULE" not in existing.has_pool_type:
                existing.has_pool_type.add("MODULE")


def compile_schedule_to_simulation(
    planned_ops: list[PlannedOperation] | OperationNetwork,
    schedule: SchedulerOutput,
    *,
    functional_modules: list[FunctionalModule] | None = None,
    op_translator: Callable[[PlannedOperation, str], SimOperation] | None = None,
    simulation_speed_factor: float = 1.0,
    random_seed: int | None = None,
    shacl_shapes=None,
    shacl_inference: str = "owlrl",
) -> Simulation:
    planned_list = _iter_planned_operations(planned_ops)
    planned_lookup = {op.identifier: op for op in planned_list}

    module_ids = set(schedule.assignments.values())
    if functional_modules is not None:
        module_ids.update(m.identifier for m in functional_modules)
    _ensure_modules_present(module_ids)

    operations: list[ExecuteScheduledTask] = []
    for op_id in schedule.start_times:
        if op_id not in schedule.end_times:
            raise KeyError(f"Missing end time for scheduled operation {op_id}")
        if op_id not in schedule.assignments:
            raise KeyError(f"Missing assignment for scheduled operation {op_id}")
        planned = planned_lookup.get(op_id)
        if planned is None:
            raise KeyError(f"Planned operation {op_id} not found in planned_ops")
        start = schedule.start_times[op_id]
        end = schedule.end_times[op_id]
        duration = end - start
        if duration < 0:
            raise ValueError(f"Negative duration for scheduled operation {op_id}: {duration}")

        module_id = schedule.assignments[op_id]
        if op_translator is not None:
            sim_op = op_translator(planned, module_id)
        else:
            sim_op = ExecuteScheduledTask(
                identifier=op_id,
                participant_module=module_id,
                planned_type=str(planned.type) if planned.type is not None else None,
                planned_annotations=dict(planned.annotations or {}),
            )

        sim_op.identifier = op_id
        sim_op.scheduled_start_time = start
        sim_op.temporal_cost = duration
        sim_op.required_precedents = list(planned.precedents)
        operations.append(sim_op)

    return Simulation(
        operations,
        simulation_speed_factor=simulation_speed_factor,
        random_seed=random_seed,
        shacl_shapes=shacl_shapes,
        shacl_inference=shacl_inference,
    )
