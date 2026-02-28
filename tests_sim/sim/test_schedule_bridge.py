from __future__ import annotations

import pytest
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical, Operation, OperationType
from libsyn_tools.opt import SchedulerOutput
from libsyn_tools.sim.adapters.schedule_bridge import (
    _ensure_modules_present,
    compile_schedule_to_simulation,
)
from libsyn_tools.sim.knowledge_graph import (
    LabObject,
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByVolume


def test_compile_schedule_to_simulation_runs_on_time():
    op_a = Operation(identifier="op_a", type=OperationType.TransferLiquid)
    op_b = Operation(identifier="op_b", type=OperationType.TransferSolid)

    schedule = SchedulerOutput(
        start_times={"op_a": 0.0, "op_b": 0.0},
        end_times={"op_a": 1.0, "op_b": 1.0},
        assignments={"op_a": "module_1", "op_b": "module_2"},
    )

    sim = compile_schedule_to_simulation([op_a, op_b], schedule)
    sim.run()

    starts = {
        r.operation_id: r.timestamp
        for r in sim.history_log
        if r.event_type == "OPERATION_START"
    }
    ends = {
        r.operation_id: r.timestamp
        for r in sim.history_log
        if r.event_type == "OPERATION_END"
    }

    assert starts["op_a"] == pytest.approx(0.0)
    assert starts["op_b"] == pytest.approx(0.0)
    assert ends["op_a"] == pytest.approx(1.0)
    assert ends["op_b"] == pytest.approx(1.0)


def test_schedule_bridge_enforces_precedents():
    op_b = Operation(identifier="op_b", type=OperationType.TransferLiquid)
    op_c = Operation(
        identifier="op_c",
        type=OperationType.TransferSolid,
        precedents=["op_b"],
    )

    schedule = SchedulerOutput(
        start_times={"op_b": 0.0, "op_c": 0.0},
        end_times={"op_b": 1.0, "op_c": 1.0},
        assignments={"op_b": "module_1", "op_c": "module_2"},
    )

    sim = compile_schedule_to_simulation([op_b, op_c], schedule)
    sim.run()

    starts = {
        r.operation_id: r.timestamp
        for r in sim.history_log
        if r.event_type == "OPERATION_START"
    }
    ends = {
        r.operation_id: r.timestamp
        for r in sim.history_log
        if r.event_type == "OPERATION_END"
    }

    assert starts["op_b"] == pytest.approx(0.0)
    assert ends["op_b"] == pytest.approx(1.0)
    assert starts["op_c"] == pytest.approx(1.0)
    assert ends["op_c"] == pytest.approx(2.0)


def test_schedule_bridge_serializes_module_contention():
    op_1 = Operation(identifier="op_1", type=OperationType.TransferLiquid)
    op_2 = Operation(identifier="op_2", type=OperationType.TransferSolid)

    schedule = SchedulerOutput(
        start_times={"op_1": 0.0, "op_2": 0.0},
        end_times={"op_1": 1.0, "op_2": 1.0},
        assignments={"op_1": "module_shared", "op_2": "module_shared"},
    )

    sim = compile_schedule_to_simulation([op_1, op_2], schedule)
    sim.run()

    starts = {
        r.operation_id: r.timestamp
        for r in sim.history_log
        if r.event_type == "OPERATION_START"
    }
    ends = {
        r.operation_id: r.timestamp
        for r in sim.history_log
        if r.event_type == "OPERATION_END"
    }

    assert starts["op_1"] == pytest.approx(0.0)
    assert ends["op_1"] == pytest.approx(1.0)
    assert starts["op_2"] == pytest.approx(1.0)
    assert ends["op_2"] == pytest.approx(2.0)


def test_schedule_bridge_translator_custom_operation_mutates_kg():
    source = MaterialContainer(identifier="https://libsyn-sim/kg/source")
    destination = MaterialContainer(identifier="https://libsyn-sim/kg/destination")
    device = MaterialContainer(identifier="https://libsyn-sim/kg/device")
    pom = PortionOfMaterial(identifier="https://libsyn-sim/kg/pom")
    pom.add_chemical(Chemical(mass=2.0, density=1.0))

    for obj in (source, destination, device, pom):
        KnowledgeGraph.get_object_from_lookup(obj.identifier)
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=source.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    planned = Operation(identifier="op_transfer", type=OperationType.TransferLiquid)
    schedule = SchedulerOutput(
        start_times={"op_transfer": 0.0},
        end_times={"op_transfer": 1.0},
        assignments={"op_transfer": "module_1"},
    )

    def _translate(_planned: Operation, _module_id: str) -> TransferMaterialByVolume:
        return TransferMaterialByVolume(
            participant_source=source.identifier,
            participant_destination=destination.identifier,
            participant_device=device.identifier,
            transfer_volume=1.0,
            temporal_cost=0.0,
        )

    sim = compile_schedule_to_simulation([planned], schedule, op_translator=_translate)
    sim.run()

    assert source.directly_contained_pom_volume == pytest.approx(1.0)
    assert destination.directly_contained_pom_volume == pytest.approx(1.0)


def test_schedule_bridge_reuses_module_identifier_form():
    module = LabObject(identifier="module_ident")
    module.has_pool_type.add("MODULE")
    KnowledgeGraph.get_object_from_lookup(module.identifier)
    Create(instance_1_iri=module.identifier).apply()

    canonical_iri = f"https://libsyn-sim/kg/{module.identifier}"
    _ensure_modules_present([canonical_iri])

    existing = KnowledgeGraph.get_object_from_lookup(module.identifier)
    assert existing is module
    assert canonical_iri not in LabObject.object_lookup


def test_schedule_bridge_normalizes_module_ids_for_ops():
    module = LabObject(identifier="module_norm")
    module.has_pool_type.add("MODULE")
    KnowledgeGraph.get_object_from_lookup(module.identifier)
    Create(instance_1_iri=module.identifier).apply()

    canonical_iri = f"https://libsyn-sim/kg/{module.identifier}"
    planned = Operation(identifier="op_mod", type=OperationType.TransferLiquid)
    schedule = SchedulerOutput(
        start_times={"op_mod": 0.0},
        end_times={"op_mod": 1.0},
        assignments={"op_mod": canonical_iri},
    )

    sim = compile_schedule_to_simulation([planned], schedule)

    assert sim.operations[0].participant_module == module.identifier


def test_schedule_bridge_rejects_extra_end_time_keys() -> None:
    planned = Operation(identifier="op_only_start", type=OperationType.TransferLiquid)
    schedule = SchedulerOutput(
        start_times={"op_only_start": 0.0},
        end_times={"op_only_start": 1.0, "op_extra": 1.0},
        assignments={"op_only_start": "module_start", "op_extra": "module_extra"},
    )

    with pytest.raises(ValueError, match="identical operation IDs"):
        compile_schedule_to_simulation([planned], schedule)


def test_schedule_bridge_rejects_extra_start_time_keys() -> None:
    planned = Operation(identifier="op_only", type=OperationType.TransferLiquid)
    schedule = SchedulerOutput(
        start_times={"op_only": 0.0, "op_extra_start": 0.0},
        end_times={"op_only": 1.0},
        assignments={"op_only": "module_only"},
    )

    with pytest.raises(ValueError, match="identical operation IDs"):
        compile_schedule_to_simulation([planned], schedule)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf"), -1.0])
def test_schedule_bridge_rejects_invalid_start_time(value: float) -> None:
    planned = Operation(identifier="op_bad_start", type=OperationType.TransferLiquid)
    schedule = SchedulerOutput(
        start_times={"op_bad_start": value},
        end_times={"op_bad_start": 1.0},
        assignments={"op_bad_start": "module_start"},
    )

    with pytest.raises(ValueError, match="start_time"):
        compile_schedule_to_simulation([planned], schedule)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf"), -1.0])
def test_schedule_bridge_rejects_invalid_end_time(value: float) -> None:
    planned = Operation(identifier="op_bad_end", type=OperationType.TransferLiquid)
    schedule = SchedulerOutput(
        start_times={"op_bad_end": 0.0},
        end_times={"op_bad_end": value},
        assignments={"op_bad_end": "module_end"},
    )

    with pytest.raises(ValueError, match="end_time"):
        compile_schedule_to_simulation([planned], schedule)


def test_schedule_bridge_normalizes_precedents():
    op_a = Operation(identifier="op_prec_a", type=OperationType.TransferLiquid)
    op_b = Operation(
        identifier="op_prec_b",
        type=OperationType.TransferSolid,
        precedents=[f"https://libsyn-sim/kg/{op_a.identifier}"],
    )
    schedule = SchedulerOutput(
        start_times={"op_prec_a": 0.0, "op_prec_b": 0.0},
        end_times={"op_prec_a": 1.0, "op_prec_b": 1.0},
        assignments={"op_prec_a": "module_prec_a", "op_prec_b": "module_prec_b"},
    )

    sim = compile_schedule_to_simulation([op_a, op_b], schedule)
    op_b_sim = next(op for op in sim.operations if op.identifier == "op_prec_b")

    assert op_b_sim.required_precedents == [op_a.identifier]
