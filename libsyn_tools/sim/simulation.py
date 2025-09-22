from __future__ import annotations

import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import simpy
from loguru import logger
from pandas._typing import FilePath
from pydantic import BaseModel
from rdflib import Graph, ConjunctiveGraph
from tqdm import tqdm

from .effect_engine import EffectEngine, KnowledgeGraph
from .knowledge_graph import LabObject, Has_interrupt_events
from .lifecycle import LifecycleCallbacks
from .operation.operation import Operation, _OpState
from .operation.runtime import _needs_runtime_tracking, get_runtime_state, _RUNTIME_CACHE
from .operation.unitary_edit import AddDataProperty


class OperationEventRecord(BaseModel):
    operation_id: str
    timestamp: float
    event_type: str
    operation_data: dict

    def __repr__(self) -> str:  # pragma: no cover
        return f"OperationEventRecord(operation_id={self.operation_id}, timestamp={self.timestamp:.3f}, event_type={self.event_type})"


class OperationProcess:
    """
    Wraps an Operation so that it can be executed as a SimPy process.
    """

    simpy_process: Optional[simpy.events.Process] = None

    def __init__(
            self,
            *,
            env: simpy.Environment,
            operation: Operation,
            operation_registry: Dict[str, "OperationProcess"],
            dependents: Dict[str, List[str]],
            history_log: List[OperationEventRecord],
            effect_engine: EffectEngine,
            speed_factor: float,
            callbacks: LifecycleCallbacks,  # NEW
    ):
        self.env = env
        self.operation = operation
        self.operation_registry = operation_registry
        self.dependents = dependents

        self.history_log = history_log
        self.effect_engine = effect_engine
        self.speed_factor = speed_factor
        self.callbacks = callbacks  # NEW

        self.done_event = env.event()

    def sim_time(self, dt: float) -> float:
        return dt * self.speed_factor

    def add_event_log(self, event_type: str, data: dict | None = None) -> None:
        self.history_log.append(
            OperationEventRecord(
                operation_id=self.operation.identifier,
                timestamp=self.env.now,
                event_type=event_type,
                operation_data=data or self.operation.model_dump(),
            )
        )

    def run(self):
        try:
            yield from self._run_core()
        except simpy.Interrupt as interrupt:
            edits = []
            reason_txt = f"{self.operation.identifier}:{interrupt.cause}"
            # record interrupt on participants (only if literal IRIs)
            for participant_name in self.operation.model_fields:
                if not participant_name.startswith("participant_"):
                    continue
                iri = getattr(self.operation, participant_name)
                if not isinstance(iri, str):
                    continue
                edits.append(
                    AddDataProperty(
                        instance_1_iri=iri,
                        property_iri=Has_interrupt_events.predicate_iri,
                        data_value=reason_txt,
                    )
                )
            if edits:
                self.effect_engine.apply(edits, self.env, operation_id=self.operation.identifier,
                                         locked_iris=self.operation.resources)
            self.operation.post_act(self.env)

            # NEW: lifecycle callback
            self.callbacks.emit_operation_interrupt(self, str(interrupt.cause))
            self.add_event_log("OPERATION_INTERRUPT", {"reason": str(interrupt.cause)})
            self.done_event.succeed()

    def _run_core(self):
        # scheduled start gate
        if self.operation.scheduled_start_time is not None:
            delay = self.operation.scheduled_start_time - self.env.now
            if delay > 0:
                yield self.env.timeout(self.sim_time(delay))

        # wait for precedents
        precedent_events = [self.operation_registry[pid].done_event for pid in self.operation.required_precedents]
        if precedent_events:
            yield simpy.events.AllOf(self.env, precedent_events)

        # pre-act
        yield self.operation.pre_act(self.env)
        self.operation._mark_running()

        # START
        self.callbacks.emit_operation_start(self)
        self.add_event_log("OPERATION_START")
        logger.debug(f"[t={self.env.now:.2f}] Start {self.operation.__class__.__name__}={self.operation.identifier}")

        # intrinsic duration
        if self.operation.temporal_cost:
            yield self.env.timeout(self.sim_time(self.operation.temporal_cost))

        # apply edits (mechanical checks + SHACL audit)
        staged = self.effect_engine.prepare(self.operation)
        self.effect_engine.apply(
            staged,
            self.env,
            operation_id=self.operation.identifier,
            locked_iris=self.operation.resources,
        )

        # provenance: remember which ops touched each runtime-tracked object
        # Guard against ANNIHILATE ⇒ resource removed; skip non-present objects
        for iri in self.operation.resources:
            obj = KnowledgeGraph.get_object_from_lookup(iri)
            if _needs_runtime_tracking(obj) and getattr(obj, "is_present", {False}) == {True}:
                get_runtime_state(obj, self.env).recent_operations.append(self.operation)

        # FINISH
        self.operation.post_act(self.env)
        self.done_event.succeed()
        self.callbacks.emit_operation_end(self)
        self.add_event_log("OPERATION_END")
        logger.debug(f"[t={self.env.now:.2f}] Finished {self.operation.identifier}")


class Simulation:
    """
    Orchestrates a collection of `Operation` instances inside a SimPy environment.
    """

    def __init__(
            self,
            operations: List[Operation],
            *,
            simulation_speed_factor: float = 1.0,
            random_seed: int | None = None,
            shacl_shapes: str | Path | Graph | None = None,
            shacl_inference: str = "owlrl"
    ):
        self.env = simpy.Environment()
        self.rng = random.Random(random_seed)

        self.operations = operations
        self.speed_factor = simulation_speed_factor

        if shacl_shapes is None:
            shapes_graph = None
        elif isinstance(shacl_shapes, (Graph, ConjunctiveGraph)):
            shapes_graph = shacl_shapes
        else:
            shapes_graph = Graph().parse(str(shacl_shapes), format="turtle")

        # NEW: lifecycle callbacks (shared among engine + processes)
        self.callbacks = LifecycleCallbacks()

        self.effect_engine = EffectEngine(
            shapes_graph=shapes_graph,
            inference=shacl_inference,
            callbacks=self.callbacks,  # NEW
        )

        self.operation_registry: Dict[str, OperationProcess] = {}
        self.dependents: Dict[str, List[str]] = defaultdict(list)
        self.history_log: List[OperationEventRecord] = []

        self._build_dependency_map()
        self._build_resources()
        self._build_processes()

    def _build_dependency_map(self) -> None:
        for op in self.operations:
            for pred in op.required_precedents:
                self.dependents[pred].append(op.identifier)

    def _build_resources(self) -> None:
        for obj in LabObject.all_instances():
            self.effect_engine._register_if_new(obj, self.env)

    def _build_processes(self) -> None:
        for op in self.operations:
            if op.sim_state is not _OpState.NEW:
                raise ValueError(f"Operation {op.identifier} already used in another Simulation")
            if op.identifier in self.operation_registry:
                raise ValueError(f"Duplicate Operation identifier {op.identifier}")
            self.operation_registry[op.identifier] = OperationProcess(
                env=self.env,
                operation=op,
                operation_registry=self.operation_registry,
                dependents=self.dependents,
                history_log=self.history_log,
                effect_engine=self.effect_engine,
                speed_factor=self.speed_factor,
                callbacks=self.callbacks,  # NEW
            )

    def run(self, until: float | None = None) -> None:
        """Start all processes and block until done or until time limit."""
        for proc in self.operation_registry.values():
            proc.simpy_process = self.env.process(proc.run())

        # progress bar via lifecycle callback (no monkey-patch)
        bar = tqdm(total=len(self.operation_registry), desc="Sim", unit="op")

        def _bar_on_end(proc: OperationProcess):
            bar.update()

        self.callbacks.on_operation_end.append(_bar_on_end)

        logger.info("Simulation start")
        self.env.run(until=until)
        logger.info(f"Simulation end @ t = {self.env.now}")

        self.callbacks.on_operation_end.remove(_bar_on_end)
        bar.close()

    def export_event_log(self, filename: FilePath) -> None:
        df_log = pd.DataFrame.from_records([r.model_dump() for r in self.history_log])
        df_log.to_csv(filename, index=False)
        logger.info(f"Event log exported → {filename}")

    def spawn_operation(
            self,
            op: Operation,
            *,
            precedents: list[str] | None = None,
            start_immediately: bool = True,
    ) -> "OperationProcess":
        if op.identifier in self.operation_registry:
            raise ValueError(f"Operation id {op.identifier!r} already exists")

        if precedents:
            op.required_precedents.extend(precedents)

        proc = OperationProcess(
            env=self.env,
            operation=op,
            operation_registry=self.operation_registry,
            dependents=self.dependents,
            history_log=self.history_log,
            effect_engine=self.effect_engine,
            speed_factor=self.speed_factor,
            callbacks=self.callbacks,  # NEW
        )
        self.operation_registry[op.identifier] = proc
        for pred in op.required_precedents:
            self.dependents[pred].append(op.identifier)

        if start_immediately:
            proc.simpy_process = self.env.process(proc.run())
        return proc

    @classmethod
    def compile_actions(cls, *actions: Operation, **kwargs) -> "Simulation":
        return cls(list(actions), **kwargs)

    def export_instance_history(self, filename: FilePath) -> None:
        end_time_index = {r.operation_id: r.timestamp for r in self.history_log if r.event_type == "OPERATION_END"}
        rows: list[dict] = []
        for rs in _RUNTIME_CACHE.values():
            if not _needs_runtime_tracking(rs.obj):
                continue
            for operation in rs.recent_operations:
                rows.append(
                    {
                        "instance_iri": rs.obj.identifier,
                        "instance_type": rs.obj.__class__.__name__,
                        "operation_id": operation.identifier,
                        "operation_type": operation.__class__.__name__,
                        "sim_timestamp": end_time_index.get(operation.identifier, None),
                    }
                )
        pd.DataFrame(rows).to_csv(filename, index=False)
        logger.info(f"Instance history exported → {filename}")
