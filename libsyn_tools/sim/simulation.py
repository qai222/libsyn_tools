from __future__ import annotations

import random
from collections import defaultdict
from typing import Dict, List, Optional

import pandas as pd
import simpy
from loguru import logger
from pandas._typing import FilePath
from pydantic import BaseModel

from .operation import Operation
from .effect_engine import EffectEngine
from .knowledge_graph import LabObject
from libsyn_tools.sim.operation.runtime import register_object_as_resource


class OperationEventRecord(BaseModel):
    operation_id: str
    timestamp: float
    event_type: str
    operation_data: dict

    def __repr__(self):
        return (
            f"OperationEventRecord(operation_id={self.operation_id}, timestamp={self.timestamp}, event_type={self.event_type})"
        )


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
            operation_registry: Dict[str, OperationProcess],
            dependents: Dict[str, List[str]],  # pre‑computed dependency map
            resource_map: Dict[str, simpy.Resource],
            history_log: List[OperationEventRecord],
            effect_engine: EffectEngine,
            speed_factor: float,
    ):
        """
        :param env: A SimPy Environment.
        :param operation: The Operation we want to simulate.
        :param operation_registry: Maps operation.identifier -> OperationProcess.
        :param resource_map: Maps resource IRI -> simpy.Resource for concurrency control.
        :param history_log: A list that will store ActionEventRecord objects.
        """
        self.dependents = dependents
        self.env = env
        self.operation = operation
        self.operation_registry = operation_registry
        self.resource_map = resource_map

        # We store the reference to a shared or global event log
        self.history_log = history_log

        # Signaling event to indicate this process has finished.
        self.done_event = env.event()

        self.effect_engine = effect_engine

        self.speed_factor = speed_factor

    def sim_time(self, dt: float) -> float:
        """Scale `dt` by the global *simulation_speed_factor* (P0‑5)."""
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
        """
        The main generator function describing the logic for running an Operation in the simulation.
        """

        resource_reqs: dict[simpy.Resource, simpy.events.Request] = dict()
        try:

            # if scheduled, minimum delay to the scheduled time
            if self.operation.scheduled_start_time is not None:
                delay = self.operation.scheduled_start_time - self.env.now
                if delay > 0:
                    yield self.env.timeout(self.sim_time(delay))

            # wait for all precedents at once
            precedent_events = [self.operation_registry[pid].done_event for pid in self.operation.required_precedents]
            if precedent_events:
                yield self.env.all_of(precedent_events)

            # check presumptions
            for pres in self.operation.presumptions:
                if not pres.validate(...):  # supply KG or context
                    raise ValueError(f"Presumption failed for {self.operation.identifier}: {pres}")

            # acquire resource
            for iri in sorted(set(self.operation.get_resources())):
                try:
                    res = self.resource_map[iri]
                except KeyError as exc:
                    raise RuntimeError(f"Resource {iri} is invalid; build_resources incomplete") from exc

                if res not in resource_reqs:
                    # choose request type based on resource class
                    resource_reqs[res] = res.request()
            yield self.env.all_of(resource_reqs.values())

            # operation start
            self.add_event_log(event_type="OPERATION_START")
            logger.debug(f"[t={self.env.now:.2f}] Start {self.operation.identifier}")

            # simulated execution delay
            if self.operation.temporal_cost:
                yield self.env.timeout(self.sim_time(self.operation.temporal_cost))

            # Stage → apply → finalize
            staged = self.effect_engine.prepare(self.operation)
            self.effect_engine.apply(staged)
            self.effect_engine.finalize(self.operation)

            # Normal completion
            self.done_event.succeed()
            self.add_event_log("ACTION_END")
            logger.debug(f"[t={self.env.now:.2f}] Finished {self.operation.identifier}")

        except simpy.Interrupt as intr:  # downstream abort (P0‑2)
            logger.warning(f"{self.operation.identifier} interrupted: {intr.cause!r}")
            self.done_event.fail(intr)
            self.add_event_log("ACTION_ABORTED")

        except Exception as err:
            self.done_event.fail(err)
            self.add_event_log("ACTION_ERROR")
            logger.error(f"[t={self.env.now:.2f}] {self.operation.identifier} failed: {err!r}")
            self._cascade_interrupt(err)
            raise
        finally:
            for req in resource_reqs.values():
                if req.triggered:
                    req.resource.release(req)

    def _cascade_interrupt(self, cause: Exception):
        """Interrupt all dependent SimPy processes (P0‑2)."""
        for dep_id in self.dependents[self.operation.identifier]:
            dep_proc = self.operation_registry[dep_id]
            if dep_proc.simpy_process and dep_proc.simpy_process.is_alive:
                try:
                    dep_proc.simpy_process.interrupt(cause)
                except RuntimeError:  # process already terminated
                    pass


class Simulation:
    def __init__(
            self,
            operations: List[Operation],
            *,
            simulation_speed_factor: float = 1.0,
            random_seed: int | None = None,
    ):
        self.env = simpy.Environment()
        self.operations = operations

        self.rng = random.Random(random_seed)
        self.speed_factor = simulation_speed_factor

        self.operation_registry: Dict[str, OperationProcess] = {}
        self.resource_map: Dict[str, simpy.Resource] = {}
        # A single list to store all event records from the entire simulation
        self.history_log: List[OperationEventRecord] = []

        self.effect_engine = EffectEngine()
        self.dependents: Dict[str, List[str]] = defaultdict(list)

        self.build_dependency_map()
        self.build_resources()
        self.build_processes()

    def build_dependency_map(self):
        for operation in self.operations:
            for pred in operation.required_precedents:
                self.dependents[pred].append(operation.identifier)

    def build_resources(self):
        for lab_obj in LabObject.object_lookup.values():
            register_object_as_resource(lab_obj, self.env)

    def build_processes(self):
        for act in self.operations:
            if act.identifier in self.operation_registry:
                raise ValueError(f"Duplicate Action identifier {act.identifier}.")
            self.operation_registry[act.identifier] = OperationProcess(
                env=self.env,
                operation=act,
                operation_registry=self.operation_registry,
                dependents=self.dependents,
                resource_map=self.resource_map,
                history_log=self.history_log,  # pass the shared log
                effect_engine=self.effect_engine,
                speed_factor=self.speed_factor,
            )

    def run(self, until: float = None):
        for proc in self.operation_registry.values():
            proc.simpy_process = self.env.process(proc.run())

        logger.info("Starting the simulation environment.")
        self.env.run(until=until)
        logger.info(f"Simulation ended at time: {self.env.now}")

    def export_event_log(self, filename: FilePath):
        """
        Example function to export the event log to a CSV or JSON file.
        """
        df_log = pd.DataFrame.from_records([r.model_dump() for r in self.history_log])
        df_log.to_csv(filename, index=False)
        logger.info(f"Event log exported to: {filename}")

    @classmethod
    def compile_actions(cls, *actions: "Operation", **kwargs) -> "Simulation":
        """
        Factory wrapper that instantiates :class:`ActionSimulation` directly.
        """
        return cls(list(actions), **kwargs)

    def get_resource_pool(self, iris: list[str]) -> list[simpy.Resource]:
        """
        Return the *Resource* objects corresponding to *iris* (ordered).
        """
        return [self.resource_map[i] for i in iris]
